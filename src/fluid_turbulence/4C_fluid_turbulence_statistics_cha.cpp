/*----------------------------------------------------------------------*/
/*! \file

\brief calculate mean values and fluctuations for turbulent channel


\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_turbulence_statistics_cha.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_fluid_xwall.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

#define NODETOL 1e-9

/*----------------------------------------------------------------------

                  Standard Constructor (public)

  ---------------------------------------------------------------------*/
FLD::TurbulenceStatisticsCha::TurbulenceStatisticsCha(Teuchos::RCP<Core::FE::Discretization> actdis,
    bool alefluid, Teuchos::RCP<Core::LinAlg::Vector<double>> dispnp,
    Teuchos::ParameterList& params, const std::string& statistics_outfilename,
    bool subgrid_dissipation, Teuchos::RCP<FLD::XWall> xwallobj)
    : discret_(actdis),
      scatradiscret_(Teuchos::null),
      alefluid_(alefluid),
      dispnp_(dispnp),
      params_(params),
      statistics_outfilename_(statistics_outfilename),
      subgrid_dissipation_(subgrid_dissipation),
      inflowchannel_(params_.sublist("TURBULENT INFLOW").get<bool>("TURBULENTINFLOW")),
      inflowmax_(params_.sublist("TURBULENT INFLOW").get<double>("INFLOW_CHA_SIDE", 0.0)),
      dens_(1.0),
      visc_(1.0),
      shc_(1.0),
      scnum_(1.0),
      myxwall_(xwallobj),
      numsubdivisions_(params_.sublist("TURBULENCE MODEL").get<int>("CHA_NUMSUBDIVISIONS"))
{
  //----------------------------------------------------------------------
  // plausibility check
  int numdim = params_.get<int>("number of velocity degrees of freedom");
  if (numdim != 3) FOUR_C_THROW("Evaluation of turbulence statistics only for 3d channel flow!");

  //----------------------------------------------------------------------
  // inflow channel check
  if (inflowchannel_)
  {
    if (discret_->get_comm().MyPID() == 0)
    {
      std::cout << "\n---------------------------------------------------------------------------"
                << std::endl;
      std::cout << "This is an additional statistics manager for turbulent inflow channels."
                << std::endl;
      std::cout << "Make sure to provide the outflow coordinate (INFLOW_CHA_SIDE)." << std::endl;
      std::cout << "Current coordinate is: " << inflowmax_ << std::endl;
      std::cout << "---------------------------------------------------------------------------\n"
                << std::endl;
    }

    // do not write any dissipation rates for inflow channels
    subgrid_dissipation_ = false;
  }

  //----------------------------------------------------------------------
  // switches, control parameters, material parameters

  // type of fluid flow solver: incompressible, Boussinesq approximation, varying density, loma
  physicaltype_ = Teuchos::getIntegralValue<Inpar::FLUID::PhysicalType>(params, "Physical Type");

  // get the plane normal direction from the parameterlist
  {
    std::string plainstring;
    if (inflowchannel_)
      plainstring =
          params_.sublist("TURBULENT INFLOW").get<std::string>("INFLOW_HOMDIR", "not_specified");
    else
      plainstring = params_.sublist("TURBULENCE MODEL").get<std::string>("HOMDIR", "not_specified");

    if (plainstring == "xz")
    {
      dim_ = 1;
    }
    else if (plainstring == "yz")
    {
      dim_ = 0;
    }
    else if (plainstring == "xy")
    {
      dim_ = 2;
    }
    else
    {
      FOUR_C_THROW("homogeneuous plane for channel flow was specified incorrectly.");
    }
  }

  // get turbulence model
  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));
  smagorinsky_ = false;
  multifractal_ = false;

  if (modelparams->get<std::string>("TURBULENCE_APPROACH", "DNS_OR_RESVMM_LES") ==
          "CLASSICAL_LES" and
      (!inflowchannel_))  // write model-related output only for pure turbulent channel flow
  {
    // check if we want to compute averages of Smagorinsky
    // constants, effective viscosities etc
    if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") == "Dynamic_Smagorinsky" ||
        params_.sublist("TURBULENCE MODEL").get<std::string>("PHYSICAL_MODEL", "no_model") ==
            "Smagorinsky_with_van_Driest_damping" ||
        params_.sublist("TURBULENCE MODEL").get<std::string>("PHYSICAL_MODEL", "no_model") ==
            "Smagorinsky")
    {
      if (discret_->get_comm().MyPID() == 0)
      {
        std::cout
            << "                             Initialising output for Smagorinsky type models\n\n\n";
        fflush(stdout);
      }

      smagorinsky_ = true;
    }
    // check if we want to compute averages of multifractal
    // quantities (N, B)
    else if (modelparams->get<std::string>("TURBULENCE_APPROACH", "DNS_OR_RESVMM_LES") ==
             "CLASSICAL_LES")
    {
      if (modelparams->get<std::string>("PHYSICAL_MODEL", "no_model") ==
          "Multifractal_Subgrid_Scales")
      {
        if (discret_->get_comm().MyPID() == 0)
        {
          std::cout << "                             Initializing output for multifractal subgrid "
                       "scales type models\n\n\n";
          fflush(stdout);
        }

        multifractal_ = true;
      }
    }
  }

  // not supported yet
  if (myxwall_ != Teuchos::null) multifractal_ = false;

  if (physicaltype_ == Inpar::FLUID::incompressible)
  {
    // get fluid viscosity from material definition --- for computation
    // of ltau
    int id = Global::Problem::instance()->materials()->first_id_by_type(Core::Materials::m_fluid);
    if (id == -1)
      FOUR_C_THROW("Could not find Newtonian fluid material");
    else
    {
      const Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance()->materials()->parameter_by_id(id);
      const Mat::PAR::NewtonianFluid* actmat = static_cast<const Mat::PAR::NewtonianFluid*>(mat);
      // we need the kinematic viscosity here
      dens_ = actmat->density_;
      visc_ = actmat->viscosity_ / dens_;
    }
  }
  else if (physicaltype_ == Inpar::FLUID::loma)
  {
    // get specific heat capacity --- for computation
    // of Temp_tau
    int id =
        Global::Problem::instance()->materials()->first_id_by_type(Core::Materials::m_sutherland);
    if (id == -1)
      FOUR_C_THROW("Could not find sutherland material");
    else
    {
      const Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance()->materials()->parameter_by_id(id);
      const Mat::PAR::Sutherland* actmat = static_cast<const Mat::PAR::Sutherland*>(mat);
      // we need the kinematic viscosity here
      shc_ = actmat->shc_;
    }
  }

  // ---------------------------------------------------------------------
  // up to now, there are no records written
  countrecord_ = 0;

  //----------------------------------------------------------------------
  // allocate some (toggle) vectors
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  meanvelnp_ = Core::LinAlg::create_vector(*dofrowmap, true);
  // this vector is only necessary for low-Mach-number flow or
  // turbulent passive scalar transport
  meanscanp_ = Core::LinAlg::create_vector(*dofrowmap, true);

  toggleu_ = Core::LinAlg::create_vector(*dofrowmap, true);
  togglev_ = Core::LinAlg::create_vector(*dofrowmap, true);
  togglew_ = Core::LinAlg::create_vector(*dofrowmap, true);
  togglep_ = Core::LinAlg::create_vector(*dofrowmap, true);

  // ---------------------------------------------------------------------
  // compute all planes for sampling

  // available planes of element nodes (polynomial)/corners
  // (Nurbs) of elements
  nodeplanes_ = Teuchos::RCP(new std::vector<double>);


  // available homogeneous (sampling) planes --- there are
  // numsubdivisions_ layers per element layer between two
  // nodes (Polynomial)/per element layer (Nurbs)
  planecoordinates_ = Teuchos::RCP(new std::vector<double>);


  // try to cast discretisation to nurbs variant
  // this tells you what kind of computation of
  // samples is required
  Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
      dynamic_cast<Core::FE::Nurbs::NurbsDiscretization*>(&(*actdis));

  // allocate array for bounding box
  //
  //          |  x  |  y  |  z
  //    ------+-----+-----+-----
  //      min |     |     |
  //    ------+-----+-----+-----
  //      max |     |     |
  //
  //
  boundingbox_ = Teuchos::RCP(new Core::LinAlg::SerialDenseMatrix(2, 3));
  for (int row = 0; row < 3; ++row)
  {
    (*boundingbox_)(0, row) = +10e+19;
    (*boundingbox_)(1, row) = -10e+19;
  }

  if (nurbsdis == nullptr)
  {
    // create set of available homogeneous planes. The normal direction
    // is read from the parameter list
    planecoordinates_ = Teuchos::RCP(new std::vector<double>);

    // the criterion allows differences in coordinates by 1e-9
    std::set<double, PlaneSortCriterion> availablecoords;

    // loop nodes, build set of planes accessible on this proc and
    // calculate bounding box
    for (int i = 0; i < discret_->num_my_row_nodes(); ++i)
    {
      Core::Nodes::Node* node = discret_->l_row_node(i);

      if (inflowchannel_ and node->x()[0] > inflowmax_ + NODETOL) continue;

      availablecoords.insert(node->x()[dim_]);

      for (int row = 0; row < 3; ++row)
      {
        if ((*boundingbox_)(0, row) > node->x()[row])
        {
          (*boundingbox_)(0, row) = node->x()[row];
        }
        if ((*boundingbox_)(1, row) < node->x()[row])
        {
          (*boundingbox_)(1, row) = node->x()[row];
        }
      }
    }

    // communicate mins
    for (int row = 0; row < 3; ++row)
    {
      double min;

      discret_->get_comm().MinAll(&((*boundingbox_)(0, row)), &min, 1);
      (*boundingbox_)(0, row) = min;
    }

    // communicate maxs
    for (int row = 0; row < 3; ++row)
    {
      double max;

      discret_->get_comm().MaxAll(&((*boundingbox_)(1, row)), &max, 1);
      (*boundingbox_)(1, row) = max;
    }

    //--------------------------------------------------------------------
    // round robin loop to communicate coordinates to all procs

    {
      int myrank = discret_->get_comm().MyPID();
      int numprocs = discret_->get_comm().NumProc();

      std::vector<char> sblock;
      std::vector<char> rblock;

      // create an exporter for point to point comunication
      Core::Communication::Exporter exporter(discret_->get_comm());

      for (int np = 0; np < numprocs; ++np)
      {
        Core::Communication::PackBuffer data;

        for (std::set<double, PlaneSortCriterion>::iterator plane = availablecoords.begin();
             plane != availablecoords.end(); ++plane)
        {
          add_to_pack(data, *plane);
        }
        swap(sblock, data());

        MPI_Request request;
        int tag = myrank;

        int frompid = myrank;
        int topid = (myrank + 1) % numprocs;

        int length = sblock.size();

        exporter.i_send(frompid, topid, sblock.data(), sblock.size(), tag, request);

        rblock.clear();

        // receive from predecessor
        frompid = (myrank + numprocs - 1) % numprocs;
        exporter.receive_any(frompid, tag, rblock, length);

        if (tag != (myrank + numprocs - 1) % numprocs)
        {
          FOUR_C_THROW("received wrong message (ReceiveAny)");
        }

        exporter.wait(request);

        {
          // for safety
          exporter.get_comm().Barrier();
        }

        // Unpack received block into set of all planes.
        {
          std::vector<double> coordsvec;

          coordsvec.clear();

          Core::Communication::UnpackBuffer buffer(rblock);
          while (!buffer.at_end())
          {
            double onecoord;
            extract_from_pack(buffer, onecoord);
            availablecoords.insert(onecoord);
          }
        }
      }
    }

    //----------------------------------------------------------------------
    // push coordinates of planes in a vector

    {
      nodeplanes_ = Teuchos::RCP(new std::vector<double>);


      for (std::set<double, PlaneSortCriterion>::iterator coord = availablecoords.begin();
           coord != availablecoords.end(); ++coord)
      {
        nodeplanes_->push_back(*coord);
      }

      // insert additional sampling planes (to show influence of quadratic
      // shape functions)

      for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
      {
        double delta = ((*nodeplanes_)[rr + 1] - (*nodeplanes_)[rr]) / ((double)numsubdivisions_);

        for (int mm = 0; mm < numsubdivisions_; ++mm)
        {
          planecoordinates_->push_back((*nodeplanes_)[rr] + delta * mm);
        }
      }
      planecoordinates_->push_back((*nodeplanes_)[(*nodeplanes_).size() - 1]);
    }
  }
  else
  {
    // pointwise sampling does not make any sense for Nurbs
    // discretisations since shape functions are not interpolating


    // planecoordinates are determined by the element (cartesian) number
    // in y direction and the number of sampling planes in between
    // and nodeplanes are kept as the corners of elements
    // (to be able to visualise stuff on the element center later on)

    // for nurbs discretisations, all vector sizes are already determined
    // by the knotvector size
    if (dim_ != 1)
    {
      FOUR_C_THROW("For the nurbs stuff, we require that xz is the hom. plane\n");
    }

    // get nurbs dis' knotvector sizes
    std::vector<int> n_x_m_x_l(nurbsdis->return_n_x_m_x_l(0));

    // get nurbs dis' element numbers
    std::vector<int> nele_x_mele_x_lele(nurbsdis->return_nele_x_mele_x_lele(0));

    // get the knotvector itself
    Teuchos::RCP<Core::FE::Nurbs::Knotvector> knots = nurbsdis->get_knot_vector();

    // resize and initialise to 0
    {
      (*nodeplanes_).resize(nele_x_mele_x_lele[1] + 1);
      (*planecoordinates_).resize(nele_x_mele_x_lele[1] * (numsubdivisions_ - 1) + 1);

      std::vector<double>::iterator coord;

      for (coord = (*nodeplanes_).begin(); coord != (*nodeplanes_).end(); ++coord)
      {
        *coord = 0;
      }
      for (coord = planecoordinates_->begin(); coord != planecoordinates_->end(); ++coord)
      {
        *coord = 0;
      }
    }

    // get element map
    const Epetra_Map* elementmap = nurbsdis->element_row_map();

    // loop all available elements
    for (int iele = 0; iele < elementmap->NumMyElements(); ++iele)
    {
      Core::Elements::Element* const actele = nurbsdis->g_element(elementmap->GID(iele));
      Core::Nodes::Node** nodes = actele->nodes();

      // get gid, location in the patch
      int gid = actele->id();

      int patchid = 0;

      std::vector<int> ele_cart_id(3);
      knots->convert_ele_gid_to_knot_ids(gid, patchid, ele_cart_id);

      // want to loop all control points of the element,
      // so get the number of points
      const int numnp = actele->num_node();

      // access elements knot span
      std::vector<Core::LinAlg::SerialDenseVector> knots(3);
      (*((*nurbsdis).get_knot_vector())).get_ele_knots(knots, actele->id());

      // aquire weights from nodes
      Core::LinAlg::SerialDenseVector weights(numnp);

      for (int inode = 0; inode < numnp; ++inode)
      {
        Core::FE::Nurbs::ControlPoint* cp =
            dynamic_cast<Core::FE::Nurbs::ControlPoint*>(nodes[inode]);

        weights(inode) = cp->w();
      }

      // get shapefunctions, compute all visualisation point positions
      Core::LinAlg::SerialDenseVector nurbs_shape_funct(numnp);

      switch (actele->shape())
      {
        case Core::FE::CellType::nurbs8:
        case Core::FE::CellType::nurbs27:
        {
          // element local point position
          Core::LinAlg::SerialDenseVector uv(3);

          {
            // standard

            //               v
            //              /
            //  w  7       /   8
            //  ^   +---------+
            //  |  /         /|
            //  | /         / |
            // 5|/        6/  |
            //  +---------+   |
            //  |         |   |
            //  |         |   +
            //  |         |  / 4
            //  |         | /
            //  |         |/
            //  +---------+ ----->u
            // 1           2
            // use v-coordinate of point 1 and 8
            // temporary x vector
            std::vector<double> x(3);

            // point 1
            uv(0) = -1.0;
            uv(1) = -1.0;
            uv(2) = -1.0;
            Core::FE::Nurbs::nurbs_get_3d_funct(
                nurbs_shape_funct, uv, knots, weights, actele->shape());
            for (int isd = 0; isd < 3; ++isd)
            {
              double val = 0;
              for (int inode = 0; inode < numnp; ++inode)
              {
                val += (((nodes[inode])->x())[isd]) * nurbs_shape_funct(inode);
              }
              x[isd] = val;
            }

            (*nodeplanes_)[ele_cart_id[1]] += x[1];
            (*planecoordinates_)[ele_cart_id[1] * (numsubdivisions_ - 1)] += x[1];

            for (int isd = 0; isd < 3; ++isd)
            {
              if ((*boundingbox_)(0, isd) > x[isd])
              {
                (*boundingbox_)(0, isd) = x[isd];
              }
              if ((*boundingbox_)(1, isd) < x[isd])
              {
                (*boundingbox_)(1, isd) = x[isd];
              }
            }

            for (int rr = 1; rr < numsubdivisions_ - 1; ++rr)
            {
              uv(1) += 2.0 / (numsubdivisions_ - 1);

              Core::FE::Nurbs::nurbs_get_3d_funct(
                  nurbs_shape_funct, uv, knots, weights, actele->shape());
              for (int isd = 0; isd < 3; ++isd)
              {
                double val = 0;
                for (int inode = 0; inode < numnp; ++inode)
                {
                  val += (((nodes[inode])->x())[isd]) * nurbs_shape_funct(inode);
                }
                x[isd] = val;
              }
              (*planecoordinates_)[ele_cart_id[1] * (numsubdivisions_ - 1) + rr] += x[1];
            }


            // set upper point of element, too (only for last layer)
            if (ele_cart_id[1] + 1 == nele_x_mele_x_lele[1])
            {
              // point 8
              uv(0) = 1.0;
              uv(1) = 1.0;
              uv(2) = 1.0;
              Core::FE::Nurbs::nurbs_get_3d_funct(
                  nurbs_shape_funct, uv, knots, weights, actele->shape());
              for (int isd = 0; isd < 3; ++isd)
              {
                double val = 0;
                for (int inode = 0; inode < numnp; ++inode)
                {
                  val += (((nodes[inode])->x())[isd]) * nurbs_shape_funct(inode);
                }
                x[isd] = val;
              }

              (*nodeplanes_)[ele_cart_id[1] + 1] += x[1];
              (*planecoordinates_)[(ele_cart_id[1] + 1) * (numsubdivisions_ - 1)] += x[1];

              for (int isd = 0; isd < 3; ++isd)
              {
                if ((*boundingbox_)(0, isd) > x[isd])
                {
                  (*boundingbox_)(0, isd) = x[isd];
                }
                if ((*boundingbox_)(1, isd) < x[isd])
                {
                  (*boundingbox_)(1, isd) = x[isd];
                }
              }
            }
          }
          break;
        }
        default:
          FOUR_C_THROW(
              "Unknown element shape for a nurbs element or nurbs type not valid for turbulence "
              "calculation\n");
          break;
      }
    }

    //----------------------------------------------------------------------
    // add contributions from all processors, normalize

    std::vector<double> lnodeplanes(*nodeplanes_);
    std::vector<double> lplanecoordinates(*planecoordinates_);

    discret_->get_comm().SumAll(lnodeplanes.data(), nodeplanes_->data(), nodeplanes_->size());
    discret_->get_comm().SumAll(
        lplanecoordinates.data(), planecoordinates_->data(), planecoordinates_->size());

    {
      (*nodeplanes_).resize(nele_x_mele_x_lele[1] + 1);
      (*planecoordinates_).resize(nele_x_mele_x_lele[1] * (numsubdivisions_ - 1) + 1);

      std::vector<double>::iterator coord;

      int nelelayer = (nele_x_mele_x_lele[0]) * (nele_x_mele_x_lele[2]);

      for (coord = (*nodeplanes_).begin(); coord != (*nodeplanes_).end(); ++coord)
      {
        *coord /= (double)(nelelayer);
      }
      for (coord = planecoordinates_->begin(); coord != planecoordinates_->end(); ++coord)
      {
        *coord /= (double)(nelelayer);
      }
    }

    // communicate mins
    for (int row = 0; row < 3; ++row)
    {
      double min;

      discret_->get_comm().MinAll(&((*boundingbox_)(0, row)), &min, 1);
      (*boundingbox_)(0, row) = min;
    }

    // communicate maxs
    for (int row = 0; row < 3; ++row)
    {
      double max;

      discret_->get_comm().MaxAll(&((*boundingbox_)(1, row)), &max, 1);
      (*boundingbox_)(1, row) = max;
    }
  }

  //----------------------------------------------------------------------
  // allocate arrays for sums of in plane mean values

  int size = planecoordinates_->size();

  // arrays for integration based averaging
  // --------------------------------------

  // first order moments
  sumu_ = Teuchos::RCP(new std::vector<double>);
  sumu_->resize(size, 0.0);

  sumv_ = Teuchos::RCP(new std::vector<double>);
  sumv_->resize(size, 0.0);

  sumw_ = Teuchos::RCP(new std::vector<double>);
  sumw_->resize(size, 0.0);

  sump_ = Teuchos::RCP(new std::vector<double>);
  sump_->resize(size, 0.0);

  sumrho_ = Teuchos::RCP(new std::vector<double>);
  sumrho_->resize(size, 0.0);

  sum_t_ = Teuchos::RCP(new std::vector<double>);
  sum_t_->resize(size, 0.0);

  sumrhou_ = Teuchos::RCP(new std::vector<double>);
  sumrhou_->resize(size, 0.0);

  sumrhou_t_ = Teuchos::RCP(new std::vector<double>);
  sumrhou_t_->resize(size, 0.0);

  // now the second order moments
  sumsqu_ = Teuchos::RCP(new std::vector<double>);
  sumsqu_->resize(size, 0.0);

  sumsqv_ = Teuchos::RCP(new std::vector<double>);
  sumsqv_->resize(size, 0.0);

  sumsqw_ = Teuchos::RCP(new std::vector<double>);
  sumsqw_->resize(size, 0.0);

  sumsqp_ = Teuchos::RCP(new std::vector<double>);
  sumsqp_->resize(size, 0.0);

  sumsqrho_ = Teuchos::RCP(new std::vector<double>);
  sumsqrho_->resize(size, 0.0);

  sumsq_t_ = Teuchos::RCP(new std::vector<double>);
  sumsq_t_->resize(size, 0.0);

  sumuv_ = Teuchos::RCP(new std::vector<double>);
  sumuv_->resize(size, 0.0);

  sumuw_ = Teuchos::RCP(new std::vector<double>);
  sumuw_->resize(size, 0.0);

  sumvw_ = Teuchos::RCP(new std::vector<double>);
  sumvw_->resize(size, 0.0);

  sumu_t_ = Teuchos::RCP(new std::vector<double>);
  sumu_t_->resize(size, 0.0);

  sumv_t_ = Teuchos::RCP(new std::vector<double>);
  sumv_t_->resize(size, 0.0);

  sumw_t_ = Teuchos::RCP(new std::vector<double>);
  sumw_t_->resize(size, 0.0);

  // arrays for point based averaging
  // --------------------------------

  pointsquaredvelnp_ = Core::LinAlg::create_vector(*dofrowmap, true);

  // first order moments
  pointsumu_ = Teuchos::RCP(new std::vector<double>);
  pointsumu_->resize(size, 0.0);

  pointsumv_ = Teuchos::RCP(new std::vector<double>);
  pointsumv_->resize(size, 0.0);

  pointsumw_ = Teuchos::RCP(new std::vector<double>);
  pointsumw_->resize(size, 0.0);

  pointsump_ = Teuchos::RCP(new std::vector<double>);
  pointsump_->resize(size, 0.0);

  // now the second order moments
  pointsumsqu_ = Teuchos::RCP(new std::vector<double>);
  pointsumsqu_->resize(size, 0.0);

  pointsumsqv_ = Teuchos::RCP(new std::vector<double>);
  pointsumsqv_->resize(size, 0.0);

  pointsumsqw_ = Teuchos::RCP(new std::vector<double>);
  pointsumsqw_->resize(size, 0.0);

  pointsumsqp_ = Teuchos::RCP(new std::vector<double>);
  pointsumsqp_->resize(size, 0.0);

  //----------------------------------------------------------------------
  // arrays for averaging of Smagorinsky constant etc.
  //

  // check if we want to compute averages of Smagorinsky
  // constants, effective viscosities etc
  if (smagorinsky_)
  {
    // extended statistics (plane average of Cs, (Cs_delta)^2, visceff)
    // for dynamic Smagorinsky model

    // vectors for element -> statistics communication
    // -----------------------------------------------

    // vectors containing processor local values to sum element
    // contributions on this proc
    Teuchos::RCP<std::vector<double>> local_Cs_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    Teuchos::RCP<std::vector<double>> local_Cs_delta_sq_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    Teuchos::RCP<std::vector<double>> local_visceff_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    Teuchos::RCP<std::vector<double>> local_Prt_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    Teuchos::RCP<std::vector<double>> local_Cs_delta_sq_Prt_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    Teuchos::RCP<std::vector<double>> local_diffeff_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    Teuchos::RCP<std::vector<double>> local_Ci_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    Teuchos::RCP<std::vector<double>> local_Ci_delta_sq_sum =
        Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));

    // store them in parameterlist for access on the element
    Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));


    modelparams->set<Teuchos::RCP<std::vector<double>>>("planecoords_", nodeplanes_);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Cs_sum", local_Cs_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>(
        "local_Cs_delta_sq_sum", local_Cs_delta_sq_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_visceff_sum", local_visceff_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Prt_sum", local_Prt_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>(
        "local_Cs_delta_sq_Prt_sum", local_Cs_delta_sq_Prt_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_diffeff_sum", local_diffeff_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Ci_sum", local_Ci_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>(
        "local_Ci_delta_sq_sum", local_Ci_delta_sq_sum);

    // vectors for statistics computation (sums and increments)
    // ----------------------------------

    // means for the Smagorinsky constant
    sum_cs_ = Teuchos::RCP(new std::vector<double>);
    sum_cs_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_cs_ = Teuchos::RCP(new std::vector<double>);
    incrsum_cs_->resize(nodeplanes_->size() - 1, 0.0);

    // means for (Cs*delta)^2
    sum_cs_delta_sq_ = Teuchos::RCP(new std::vector<double>);
    sum_cs_delta_sq_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_cs_delta_sq_ = Teuchos::RCP(new std::vector<double>);
    incrsum_cs_delta_sq_->resize(nodeplanes_->size() - 1, 0.0);

    // means for the effective viscosity
    sumvisceff_ = Teuchos::RCP(new std::vector<double>);
    sumvisceff_->resize(nodeplanes_->size() - 1, 0.0);

    incrsumvisceff_ = Teuchos::RCP(new std::vector<double>);
    incrsumvisceff_->resize(nodeplanes_->size() - 1, 0.0);

    // means for Prt
    sum_prt_ = Teuchos::RCP(new std::vector<double>);
    sum_prt_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_prt_ = Teuchos::RCP(new std::vector<double>);
    incrsum_prt_->resize(nodeplanes_->size() - 1, 0.0);

    // means for (Cs*delta)^2/Prt
    sum_cs_delta_sq_prt_ = Teuchos::RCP(new std::vector<double>);
    sum_cs_delta_sq_prt_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_cs_delta_sq_prt_ = Teuchos::RCP(new std::vector<double>);
    incrsum_cs_delta_sq_prt_->resize(nodeplanes_->size() - 1, 0.0);

    // means for the effective diffusivity
    sumdiffeff_ = Teuchos::RCP(new std::vector<double>);
    sumdiffeff_->resize(nodeplanes_->size() - 1, 0.0);

    incrsumdiffeff_ = Teuchos::RCP(new std::vector<double>);
    incrsumdiffeff_->resize(nodeplanes_->size() - 1, 0.0);

    // means for Ci
    sum_ci_ = Teuchos::RCP(new std::vector<double>);
    sum_ci_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_ci_ = Teuchos::RCP(new std::vector<double>);
    incrsum_ci_->resize(nodeplanes_->size() - 1, 0.0);

    // means for (Ci*delta)^2
    sum_ci_delta_sq_ = Teuchos::RCP(new std::vector<double>);
    sum_ci_delta_sq_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_ci_delta_sq_ = Teuchos::RCP(new std::vector<double>);
    incrsum_ci_delta_sq_->resize(nodeplanes_->size() - 1, 0.0);
  }

  //----------------------------------------------------------------------
  // arrays for averaging of parameters of multifractal subgrid-scales
  //

  // check if we want to compute averages of multifractal subgrid-scales
  if (multifractal_)
  {
    // vectors for statistics computation (sums and increments)
    // means for N
    sum_n_stream_ = Teuchos::RCP(new std::vector<double>);
    sum_n_stream_->resize(nodeplanes_->size() - 1, 0.0);
    sum_n_normal_ = Teuchos::RCP(new std::vector<double>);
    sum_n_normal_->resize(nodeplanes_->size() - 1, 0.0);
    sum_n_span_ = Teuchos::RCP(new std::vector<double>);
    sum_n_span_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_n_stream_ = Teuchos::RCP(new std::vector<double>);
    incrsum_n_stream_->resize(nodeplanes_->size() - 1, 0.0);
    incrsum_n_normal_ = Teuchos::RCP(new std::vector<double>);
    incrsum_n_normal_->resize(nodeplanes_->size() - 1, 0.0);
    incrsum_n_span_ = Teuchos::RCP(new std::vector<double>);
    incrsum_n_span_->resize(nodeplanes_->size() - 1, 0.0);

    // means for B
    sum_b_stream_ = Teuchos::RCP(new std::vector<double>);
    sum_b_stream_->resize(nodeplanes_->size() - 1, 0.0);
    sum_b_normal_ = Teuchos::RCP(new std::vector<double>);
    sum_b_normal_->resize(nodeplanes_->size() - 1, 0.0);
    sum_b_span_ = Teuchos::RCP(new std::vector<double>);
    sum_b_span_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_b_stream_ = Teuchos::RCP(new std::vector<double>);
    incrsum_b_stream_->resize(nodeplanes_->size() - 1, 0.0);
    incrsum_b_normal_ = Teuchos::RCP(new std::vector<double>);
    incrsum_b_normal_->resize(nodeplanes_->size() - 1, 0.0);
    incrsum_b_span_ = Teuchos::RCP(new std::vector<double>);
    incrsum_b_span_->resize(nodeplanes_->size() - 1, 0.0);

    // means for C_sgs
    sum_csgs_ = Teuchos::RCP(new std::vector<double>);
    sum_csgs_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_csgs_ = Teuchos::RCP(new std::vector<double>);
    incrsum_csgs_->resize(nodeplanes_->size() - 1, 0.0);

    // means for subgrid viscosity
    // if used in combination with eddy viscosity model
    sumsgvisc_ = Teuchos::RCP(new std::vector<double>);
    sumsgvisc_->resize(nodeplanes_->size() - 1, 0.0);

    incrsumsgvisc_ = Teuchos::RCP(new std::vector<double>);
    incrsumsgvisc_->resize(nodeplanes_->size() - 1, 0.0);

    // vectors for statistics computation (sums and increments)
    // means for Nphi
    sum_nphi_ = Teuchos::RCP(new std::vector<double>);
    sum_nphi_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_nphi_ = Teuchos::RCP(new std::vector<double>);
    incrsum_nphi_->resize(nodeplanes_->size() - 1, 0.0);

    // means for D
    sum_dphi_ = Teuchos::RCP(new std::vector<double>);
    sum_dphi_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_dphi_ = Teuchos::RCP(new std::vector<double>);
    incrsum_dphi_->resize(nodeplanes_->size() - 1, 0.0);

    // means for C_sgs_phi
    sum_csgs_phi_ = Teuchos::RCP(new std::vector<double>);
    sum_csgs_phi_->resize(nodeplanes_->size() - 1, 0.0);

    incrsum_csgs_phi_ = Teuchos::RCP(new std::vector<double>);
    incrsum_csgs_phi_->resize(nodeplanes_->size() - 1, 0.0);
  }

  //----------------------------------------------------------------------
  // arrays for averaging of residual, subscales etc.

  // prepare time averaging for subscales and residual
  if (subgrid_dissipation_)
  {
    //--------------------------------------------------
    // local_incrtauC            (in plane) averaged values of stabilization parameter tauC
    // local_incrtauM            (in plane) averaged values of stabilization parameter tauM
    // local_incrres(_sq,abs)    (in plane) averaged values of resM (^2) (||.||)
    // local_incrsacc(_sq,abs)   (in plane) averaged values of sacc (^2) (||.||)
    // local_incrsvelaf(_sq,abs) (in plane) averaged values of svelaf (^2) (||.||)
    // local_incrresC(_sq)       (in plane) averaged values of resC (^2)
    // local_incrspressnp(_sq)   (in plane) averaged values of spressnp (^2)
    //--------------------------------------------------
    Teuchos::RCP<std::vector<double>> local_incrvol =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrhk =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrhbazilevs =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrstrle =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrgradle =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incrtauC =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrtauM =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incrmk =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incrres =
        Teuchos::RCP(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrres_sq =
        Teuchos::RCP(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrabsres =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incrtauinvsvel =
        Teuchos::RCP(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incrsvelaf =
        Teuchos::RCP(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrsvelaf_sq =
        Teuchos::RCP(new std::vector<double>(3 * (nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrabssvelaf =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incrresC =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrresC_sq =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrspressnp =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrspressnp_sq =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incr_eps_pspg =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_supg =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_cross =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_rey =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_graddiv =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_eddyvisc =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_visc =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_conv =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_mfs =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_mfscross =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_mfsrey =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incr_eps_avm3 =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_incrcrossstress =
        Teuchos::RCP(new std::vector<double>(6 * (nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_incrreystress =
        Teuchos::RCP(new std::vector<double>(6 * (nodeplanes_->size() - 1), 0.0));

    // pass pointers to local sum vectors to the element
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrvol", local_incrvol);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhk", local_incrhk);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhbazilevs", local_incrhbazilevs);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrstrle", local_incrstrle);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrgradle", local_incrgradle);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrmk", local_incrmk);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("planecoords_", nodeplanes_);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauC", local_incrtauC);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauM", local_incrtauM);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres", local_incrres);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres_sq", local_incrres_sq);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabsres", local_incrabsres);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauinvsvel", local_incrtauinvsvel);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf", local_incrsvelaf);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf_sq", local_incrsvelaf_sq);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabssvelaf", local_incrabssvelaf);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC", local_incrresC);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC_sq", local_incrresC_sq);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp", local_incrspressnp);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp_sq", local_incrspressnp_sq);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_pspg", local_incr_eps_pspg);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_supg", local_incr_eps_supg);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_cross", local_incr_eps_cross);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_rey", local_incr_eps_rey);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_graddiv", local_incr_eps_graddiv);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_eddyvisc", local_incr_eps_eddyvisc);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_visc", local_incr_eps_visc);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_conv", local_incr_eps_conv);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfs", local_incr_eps_mfs);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfscross", local_incr_eps_mfscross);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfsrey", local_incr_eps_mfsrey);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_avm3", local_incr_eps_avm3);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrcrossstress", local_incrcrossstress);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrreystress", local_incrreystress);

    // means for comparison of of residual and subscale acceleration

    sumres_ = Teuchos::RCP(new std::vector<double>);
    sumres_->resize(3 * (nodeplanes_->size() - 1), 0.0);
    sumres_sq_ = Teuchos::RCP(new std::vector<double>);
    sumres_sq_->resize(3 * (nodeplanes_->size() - 1), 0.0);
    sumabsres_ = Teuchos::RCP(new std::vector<double>);
    sumabsres_->resize((nodeplanes_->size() - 1), 0.0);
    sumtauinvsvel_ = Teuchos::RCP(new std::vector<double>);
    sumtauinvsvel_->resize(3 * (nodeplanes_->size() - 1), 0.0);

    sumsvelaf_ = Teuchos::RCP(new std::vector<double>);
    sumsvelaf_->resize(3 * (nodeplanes_->size() - 1), 0.0);
    sumsvelaf_sq_ = Teuchos::RCP(new std::vector<double>);
    sumsvelaf_sq_->resize(3 * (nodeplanes_->size() - 1), 0.0);
    sumabssvelaf_ = Teuchos::RCP(new std::vector<double>);
    sumabssvelaf_->resize((nodeplanes_->size() - 1), 0.0);

    sumres_c_ = Teuchos::RCP(new std::vector<double>);
    sumres_c_->resize(nodeplanes_->size() - 1, 0.0);
    sumres_c_sq_ = Teuchos::RCP(new std::vector<double>);
    sumres_c_sq_->resize(nodeplanes_->size() - 1, 0.0);

    sumspressnp_ = Teuchos::RCP(new std::vector<double>);
    sumspressnp_->resize(nodeplanes_->size() - 1, 0.0);
    sumspressnp_sq_ = Teuchos::RCP(new std::vector<double>);
    sumspressnp_sq_->resize(nodeplanes_->size() - 1, 0.0);

    sumhk_ = Teuchos::RCP(new std::vector<double>);
    sumhk_->resize(nodeplanes_->size() - 1, 0.0);
    sumhbazilevs_ = Teuchos::RCP(new std::vector<double>);
    sumhbazilevs_->resize(nodeplanes_->size() - 1, 0.0);
    sumstrle_ = Teuchos::RCP(new std::vector<double>);
    sumstrle_->resize(nodeplanes_->size() - 1, 0.0);
    sumgradle_ = Teuchos::RCP(new std::vector<double>);
    sumgradle_->resize(nodeplanes_->size() - 1, 0.0);
    sumtau_m_ = Teuchos::RCP(new std::vector<double>);
    sumtau_m_->resize(nodeplanes_->size() - 1, 0.0);
    sumtau_c_ = Teuchos::RCP(new std::vector<double>);
    sumtau_c_->resize(nodeplanes_->size() - 1, 0.0);

    summk_ = Teuchos::RCP(new std::vector<double>);
    summk_->resize(nodeplanes_->size() - 1, 0.0);

    sum_eps_pspg_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_pspg_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_supg_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_supg_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_cross_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_cross_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_rey_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_rey_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_graddiv_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_graddiv_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_eddyvisc_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_eddyvisc_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_visc_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_visc_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_conv_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_conv_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_mfs_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_mfs_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_mfscross_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_mfscross_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_mfsrey_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_mfsrey_->resize(nodeplanes_->size() - 1, 0.0);
    sum_eps_avm3_ = Teuchos::RCP(new std::vector<double>);
    sum_eps_avm3_->resize(nodeplanes_->size() - 1, 0.0);

    sum_crossstress_ = Teuchos::RCP(new std::vector<double>);
    sum_crossstress_->resize(6 * (nodeplanes_->size() - 1), 0.0);
    sum_reystress_ = Teuchos::RCP(new std::vector<double>);
    sum_reystress_->resize(6 * (nodeplanes_->size() - 1), 0.0);

    // add quantities for scatra here
    Teuchos::RCP<std::vector<double>> local_scatra_incrvol =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_scatra_incrtauS =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incrresS =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incrresS_sq =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));

    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_supg =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_cross =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_rey =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_eddyvisc =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_visc =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_conv =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_mfs =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_mfscross =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_mfsrey =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));
    Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_avm3 =
        Teuchos::RCP(new std::vector<double>((nodeplanes_->size() - 1), 0.0));


    // pass pointers to local sum vectors to the element
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>("planecoords_", nodeplanes_);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>("incrvol", local_scatra_incrvol);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauS", local_scatra_incrtauS);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresS", local_scatra_incrresS);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incrresS_sq", local_scatra_incrresS_sq);

    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_supg", local_scatra_incr_eps_supg);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_cross", local_scatra_incr_eps_cross);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_rey", local_scatra_incr_eps_rey);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_eddyvisc", local_scatra_incr_eps_eddyvisc);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_visc", local_scatra_incr_eps_visc);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_conv", local_scatra_incr_eps_conv);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_mfs", local_scatra_incr_eps_mfs);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_mfscross", local_scatra_incr_eps_mfscross);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_mfsrey", local_scatra_incr_eps_mfsrey);
    scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
        "incr_scatra_eps_avm3", local_scatra_incr_eps_avm3);


    // means for comparison of of residual and dissipation
    sumres_s_ = Teuchos::RCP(new std::vector<double>);
    sumres_s_->resize(nodeplanes_->size() - 1, 0.0);
    sumres_s_sq_ = Teuchos::RCP(new std::vector<double>);
    sumres_s_sq_->resize(nodeplanes_->size() - 1, 0.0);
    sumtau_s_ = Teuchos::RCP(new std::vector<double>);
    sumtau_s_->resize(nodeplanes_->size() - 1, 0.0);

    sum_scatra_eps_supg_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_supg_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_cross_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_cross_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_rey_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_rey_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_eddyvisc_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_eddyvisc_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_visc_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_visc_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_conv_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_conv_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_mfs_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_mfs_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_mfscross_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_mfscross_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_mfsrey_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_mfsrey_->resize(nodeplanes_->size() - 1, 0.0);
    sum_scatra_eps_avm3_ = Teuchos::RCP(new std::vector<double>);
    sum_scatra_eps_avm3_->resize(nodeplanes_->size() - 1, 0.0);
  }


  //----------------------------------------------------------------------
  // initialize output and initially open respective statistics output file(s)

  Teuchos::RCP<std::ofstream> log;
  Teuchos::RCP<std::ofstream> log_Cs;
  Teuchos::RCP<std::ofstream> log_SSM;
  Teuchos::RCP<std::ofstream> log_MF;
  Teuchos::RCP<std::ofstream> log_res;
  Teuchos::RCP<std::ofstream> log_res_scatra;

  if (discret_->get_comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);

    if (physicaltype_ == Inpar::FLUID::loma)
    {
      if (inflowchannel_)
        s.append(".inflow.loma_statistics");
      else
        s.append(".loma_statistics");

      log = Teuchos::RCP(new std::ofstream(s.c_str(), std::ios::out));
      (*log) << "# Statistics for turbulent variable-density channel flow at low Mach number "
                "(first- and second-order moments)\n\n";

      log->flush();

      // additional output for dynamic Smagorinsky model
      if (smagorinsky_)
      {
        std::string s_smag(statistics_outfilename_);
        s_smag.append(".Cs_statistics");

        log_Cs = Teuchos::RCP(new std::ofstream(s_smag.c_str(), std::ios::out));
        (*log_Cs)
            << "# Statistics for turbulent incompressible channel flow (Smagorinsky constant)\n\n";
      }
    }
    else
    {
      if (inflowchannel_)
        s.append(".inflow.flow_statistics");
      else
        s.append(".flow_statistics");

      log = Teuchos::RCP(new std::ofstream(s.c_str(), std::ios::out));
      (*log) << "# Statistics for turbulent incompressible channel flow (first- and second-order "
                "moments)\n\n";

      log->flush();

      // additional output for dynamic Smagorinsky model
      if (smagorinsky_)
      {
        std::string s_smag(statistics_outfilename_);
        s_smag.append(".Cs_statistics");

        log_Cs = Teuchos::RCP(new std::ofstream(s_smag.c_str(), std::ios::out));
        (*log_Cs)
            << "# Statistics for turbulent incompressible channel flow (Smagorinsky constant)\n\n";
      }

      // additional output for multifractal subgrid scales
      if (multifractal_)
      {
        std::string s_mf(statistics_outfilename_);
        s_mf.append(".MF_statistics");

        log_MF = Teuchos::RCP(new std::ofstream(s_mf.c_str(), std::ios::out));
        (*log_MF) << "# Statistics for turbulent incompressible channel flow (parameter "
                     "multifractal subgrid scales)\n\n";
      }
    }

    // additional output of residuals and subscale quantities
    if (subgrid_dissipation_)
    {
      std::string s_res(statistics_outfilename_);
      s_res.append(".res_statistics");

      log_res = Teuchos::RCP(new std::ofstream(s_res.c_str(), std::ios::out));
      (*log_res) << "# Statistics for turbulent incompressible channel flow (residuals and "
                    "subscale quantities)\n";
      (*log_res) << "# All values are first averaged over the integration points in an element \n";
      (*log_res)
          << "# and after that averaged over a whole element layer in the homogeneous plane\n\n";

      std::string s_res_scatra(statistics_outfilename_);
      s_res_scatra.append(".res_scatra_statistics");

      log_res_scatra = Teuchos::RCP(new std::ofstream(s_res_scatra.c_str(), std::ios::out));
      (*log_res_scatra) << "# Statistics for turbulent incompressible channel flow with scalar "
                           "transport (residuals and subscale quantities)\n";
      (*log_res_scatra)
          << "# All values are first averaged over the integration points in an element \n";
      (*log_res_scatra)
          << "# and after that averaged over a whole element layer in the homogeneous plane\n\n";
      (*log_res_scatra)
          << "#                           THIS IS THE SCATRA FILE                          \n\n";
    }
  }

  // clear statistics
  this->clear_statistics();

  return;
}  // TurbulenceStatisticsCha::TurbulenceStatisticsCha



/*----------------------------------------------------------------------*

       Compute the in-plane mean values of first and second order
       moments for velocities, pressure and Cs are added to global
                            'sum' vectors.

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::do_time_sample(
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> velnp,
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> force)
{
  // we have an additional sample
  numsamp_++;

  // meanvelnp is a refcount copy of velnp
  meanvelnp_->Update(1.0, *velnp, 0.0);

  //----------------------------------------------------------------------
  // loop planes and calculate integral means in each plane

  this->evaluate_integral_mean_values_in_planes();

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane

  // try to cast discretisation to nurbs variant
  // this tells you whether pointwise computation of
  // samples is allowed
  Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
      dynamic_cast<Core::FE::Nurbs::NurbsDiscretization*>(&(*discret_));

  if (nurbsdis == nullptr)
  {
    this->evaluate_pointwise_mean_values_in_planes();
  }

  //----------------------------------------------------------------------
  // compute forces on top and bottom plate for normalization purposes

  for (std::vector<double>::iterator plane = planecoordinates_->begin();
       plane != planecoordinates_->end(); ++plane)
  {
    // only true for top and bottom plane
    if ((*plane - 2e-9 < (*planecoordinates_)[0] && *plane + 2e-9 > (*planecoordinates_)[0]) ||
        (*plane - 2e-9 < (*planecoordinates_)[planecoordinates_->size() - 1] &&
            *plane + 2e-9 > (*planecoordinates_)[planecoordinates_->size() - 1]))
    {
      // toggle vectors are one in the position of a dof in this plane,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);

      // activate toggles for in plane dofs
      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // if we have an inflow channel problem, the nodes outside the inflow discretization are
        // not in the bounding box -> we don't consider them for averaging
        if (node->x()[0] < (*boundingbox_)(1, 0) + NODETOL and
            node->x()[1] < (*boundingbox_)(1, 1) + NODETOL and
            node->x()[2] < (*boundingbox_)(1, 2) + NODETOL and
            node->x()[0] > (*boundingbox_)(0, 0) - NODETOL and
            node->x()[1] > (*boundingbox_)(0, 1) - NODETOL and
            node->x()[2] > (*boundingbox_)(0, 2) - NODETOL)
        {
          // this node belongs to the plane under consideration
          if (node->x()[dim_] < *plane + 2e-9 && node->x()[dim_] > *plane - 2e-9)
          {
            std::vector<int> dof = discret_->dof(node);
            double one = 1.0;

            toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
            togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
            togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
          }
        }
      }

      // compute forces by dot product
      double inc = 0.0;
      {
        double local_inc = 0.0;
        for (int rr = 0; rr < (*toggleu_).MyLength(); ++rr)
        {
          local_inc += (*toggleu_)[rr] * (*toggleu_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);

        if (abs(inc) < 1e-9)
        {
          FOUR_C_THROW("there are no forced nodes on the boundary\n");
        }

        local_inc = 0.0;
        for (int rr = 0; rr < force->MyLength(); ++rr)
        {
          local_inc += (*force)[rr] * (*toggleu_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        sumforceu_ += inc;

        local_inc = 0.0;
        for (int rr = 0; rr < force->MyLength(); ++rr)
        {
          local_inc += (*force)[rr] * (*togglev_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        sumforcev_ += inc;


        local_inc = 0.0;
        for (int rr = 0; rr < force->MyLength(); ++rr)
        {
          local_inc += (*force)[rr] * (*togglew_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        sumforcew_ += inc;
      }
    }
  }

  //----------------------------------------------------------------------
  // add increment of last iteration to the sum of Cs values
  // (statistics for dynamic Smagorinsky model)

  if (smagorinsky_)
  {
    for (unsigned rr = 0; rr < (*incrsum_cs_).size(); ++rr)
    {
      (*sum_cs_)[rr] += (*incrsum_cs_)[rr];
      (*sum_cs_delta_sq_)[rr] += (*incrsum_cs_delta_sq_)[rr];
      (*sumvisceff_)[rr] += (*incrsumvisceff_)[rr];
      (*sum_prt_)[rr] += (*incrsum_prt_)[rr];
      (*sum_cs_delta_sq_prt_)[rr] += (*incrsum_cs_delta_sq_prt_)[rr];
      (*sumdiffeff_)[rr] += (*incrsumdiffeff_)[rr];
      (*sum_ci_)[rr] += (*incrsum_ci_)[rr];
      (*sum_ci_delta_sq_)[rr] += (*incrsum_ci_delta_sq_)[rr];
    }
  }

  return;
}  // TurbulenceStatisticsCha::DoTimeSample


/*----------------------------------------------------------------------*

       Compute the in-plane mean values of first- and second-order
                 moments for low-Mach-number flow.

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::do_loma_time_sample(
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> velnp,
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> scanp,
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> force, const double eosfac)
{
  // we have an additional sample

  numsamp_++;

  // meanvelnp and meanscanp are refcount copies
  meanvelnp_->Update(1.0, *velnp, 0.0);
  meanscanp_->Update(1.0, *scanp, 0.0);

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane

  this->evaluate_loma_integral_mean_values_in_planes(eosfac);

  //----------------------------------------------------------------------
  // compute forces on top and bottom plate for normalization purposes

  for (std::vector<double>::iterator plane = planecoordinates_->begin();
       plane != planecoordinates_->end(); ++plane)
  {
    // only true for bottom plane
    if (*plane - 2e-9 < (*planecoordinates_)[0] && *plane + 2e-9 > (*planecoordinates_)[0])
    {
      // toggle vectors are one in the position of a dof in this plane,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // activate toggles for in plane dofs
      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // if we have an inflow channel problem, the nodes outside the inflow discretization are
        // not in the bounding box -> we don't consider them for averaging
        if (node->x()[0] < (*boundingbox_)(1, 0) + NODETOL and
            node->x()[1] < (*boundingbox_)(1, 1) + NODETOL and
            node->x()[2] < (*boundingbox_)(1, 2) + NODETOL and
            node->x()[0] > (*boundingbox_)(0, 0) - NODETOL and
            node->x()[1] > (*boundingbox_)(0, 1) - NODETOL and
            node->x()[2] > (*boundingbox_)(0, 2) - NODETOL)
        {
          // this node belongs to the plane under consideration
          if (node->x()[dim_] < *plane + 2e-9 && node->x()[dim_] > *plane - 2e-9)
          {
            std::vector<int> dof = discret_->dof(node);
            double one = 1.0;

            toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
            togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
            togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
            togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));
          }
        }
      }

      double inc = 0.0;
      force->Dot(*toggleu_, &inc);
      sumforcebu_ += inc;
      inc = 0.0;
      force->Dot(*togglev_, &inc);
      sumforcebv_ += inc;
      inc = 0.0;
      force->Dot(*togglew_, &inc);
      sumforcebw_ += inc;
      inc = 0.0;
      force->Dot(*togglep_, &inc);
      sumqwb_ += inc;
    }

    // only true for top plane
    if (*plane - 2e-9 < (*planecoordinates_)[planecoordinates_->size() - 1] &&
        *plane + 2e-9 > (*planecoordinates_)[planecoordinates_->size() - 1])
    {
      // toggle vectors are one in the position of a dof in this plane,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // activate toggles for in plane dofs
      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // if we have an inflow channel problem, the nodes outside the inflow discretization are
        // not in the bounding box -> we don't consider them for averaging
        if (node->x()[0] < (*boundingbox_)(1, 0) + NODETOL and
            node->x()[1] < (*boundingbox_)(1, 1) + NODETOL and
            node->x()[2] < (*boundingbox_)(1, 2) + NODETOL and
            node->x()[0] > (*boundingbox_)(0, 0) - NODETOL and
            node->x()[1] > (*boundingbox_)(0, 1) - NODETOL and
            node->x()[2] > (*boundingbox_)(0, 2) - NODETOL)
        {
          // this node belongs to the plane under consideration
          if (node->x()[dim_] < *plane + 2e-9 && node->x()[dim_] > *plane - 2e-9)
          {
            std::vector<int> dof = discret_->dof(node);
            double one = 1.0;

            toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
            togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
            togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
            togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));
          }
        }
      }

      double inc = 0.0;
      force->Dot(*toggleu_, &inc);
      sumforcetu_ += inc;
      inc = 0.0;
      force->Dot(*togglev_, &inc);
      sumforcetv_ += inc;
      inc = 0.0;
      force->Dot(*togglew_, &inc);
      sumforcetw_ += inc;
      inc = 0.0;
      force->Dot(*togglep_, &inc);
      sumqwt_ += inc;
    }
  }

  //----------------------------------------------------------------------
  // add increment of last iteration to the sum of Cs values
  // (statistics for dynamic Smagorinsky model)

  if (smagorinsky_)
  {
    for (unsigned rr = 0; rr < (*incrsum_cs_).size(); ++rr)
    {
      (*sum_cs_)[rr] += (*incrsum_cs_)[rr];
      (*sum_cs_delta_sq_)[rr] += (*incrsum_cs_delta_sq_)[rr];
      (*sumvisceff_)[rr] += (*incrsumvisceff_)[rr];
      (*sum_prt_)[rr] += (*incrsum_prt_)[rr];
      (*sum_cs_delta_sq_prt_)[rr] += (*incrsum_cs_delta_sq_prt_)[rr];
      (*sumdiffeff_)[rr] += (*incrsumdiffeff_)[rr];
      (*sum_ci_)[rr] += (*incrsum_ci_)[rr];
      (*sum_ci_delta_sq_)[rr] += (*incrsum_ci_delta_sq_)[rr];
    }
  }

  return;
}  // TurbulenceStatisticsCha::DoLomaTimeSample


/*----------------------------------------------------------------------*

       Compute the in-plane mean values of first- and second-order
                 moments for turbulent flow with passive scalar
                            transport.

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::do_scatra_time_sample(
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> velnp,
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> scanp,
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> force)
{
  // we have an additional sample

  numsamp_++;

  // meanvelnp and meanscanp are refcount copies
  meanvelnp_->Update(1.0, *velnp, 0.0);
  meanscanp_->Update(1.0, *scanp, 0.0);

  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane

  this->evaluate_scatra_integral_mean_values_in_planes();

  //----------------------------------------------------------------------
  // compute forces on top and bottom plate for normalization purposes

  for (std::vector<double>::iterator plane = planecoordinates_->begin();
       plane != planecoordinates_->end(); ++plane)
  {
    // only true for bottom plane
    if (*plane - 2e-9 < (*planecoordinates_)[0] && *plane + 2e-9 > (*planecoordinates_)[0])
    {
      // toggle vectors are one in the position of a dof in this plane,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // activate toggles for in plane dofs
      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // if we have an inflow channel problem, the nodes outside the inflow discretization are
        // not in the bounding box -> we don't consider them for averaging
        if (node->x()[0] < (*boundingbox_)(1, 0) + NODETOL and
            node->x()[1] < (*boundingbox_)(1, 1) + NODETOL and
            node->x()[2] < (*boundingbox_)(1, 2) + NODETOL and
            node->x()[0] > (*boundingbox_)(0, 0) - NODETOL and
            node->x()[1] > (*boundingbox_)(0, 1) - NODETOL and
            node->x()[2] > (*boundingbox_)(0, 2) - NODETOL)
        {
          // this node belongs to the plane under consideration
          if (node->x()[dim_] < *plane + 2e-9 && node->x()[dim_] > *plane - 2e-9)
          {
            std::vector<int> dof = discret_->dof(node);
            double one = 1.0;

            toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
            togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
            togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
            togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));
          }
        }
      }

      double inc = 0.0;
      force->Dot(*toggleu_, &inc);
      sumforcebu_ += inc;
      inc = 0.0;
      force->Dot(*togglev_, &inc);
      sumforcebv_ += inc;
      inc = 0.0;
      force->Dot(*togglew_, &inc);
      sumforcebw_ += inc;
      inc = 0.0;
      force->Dot(*togglep_, &inc);
      sumqwb_ += inc;
    }

    // only true for top plane
    if (*plane - 2e-9 < (*planecoordinates_)[planecoordinates_->size() - 1] &&
        *plane + 2e-9 > (*planecoordinates_)[planecoordinates_->size() - 1])
    {
      // toggle vectors are one in the position of a dof in this plane,
      // else 0
      toggleu_->PutScalar(0.0);
      togglev_->PutScalar(0.0);
      togglew_->PutScalar(0.0);
      togglep_->PutScalar(0.0);

      // activate toggles for in plane dofs
      for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
      {
        Core::Nodes::Node* node = discret_->l_row_node(nn);

        // if we have an inflow channel problem, the nodes outside the inflow discretization are
        // not in the bounding box -> we don't consider them for averaging
        if (node->x()[0] < (*boundingbox_)(1, 0) + NODETOL and
            node->x()[1] < (*boundingbox_)(1, 1) + NODETOL and
            node->x()[2] < (*boundingbox_)(1, 2) + NODETOL and
            node->x()[0] > (*boundingbox_)(0, 0) - NODETOL and
            node->x()[1] > (*boundingbox_)(0, 1) - NODETOL and
            node->x()[2] > (*boundingbox_)(0, 2) - NODETOL)
        {
          // this node belongs to the plane under consideration
          if (node->x()[dim_] < *plane + 2e-9 && node->x()[dim_] > *plane - 2e-9)
          {
            std::vector<int> dof = discret_->dof(node);
            double one = 1.0;

            toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
            togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
            togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
            togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));
          }
        }
      }

      double inc = 0.0;
      force->Dot(*toggleu_, &inc);
      sumforcetu_ += inc;
      inc = 0.0;
      force->Dot(*togglev_, &inc);
      sumforcetv_ += inc;
      inc = 0.0;
      force->Dot(*togglew_, &inc);
      sumforcetw_ += inc;
      inc = 0.0;
      force->Dot(*togglep_, &inc);
      sumqwt_ += inc;
    }
  }

  //----------------------------------------------------------------------
  // add increment of last iteration to the sum of Cs values
  // (statistics for dynamic Smagorinsky model)

  if (smagorinsky_)
  {
    for (unsigned rr = 0; rr < (*incrsum_cs_).size(); ++rr)
    {
      (*sum_cs_)[rr] += (*incrsum_cs_)[rr];
      (*sum_cs_delta_sq_)[rr] += (*incrsum_cs_delta_sq_)[rr];
      (*sumvisceff_)[rr] += (*incrsumvisceff_)[rr];
      (*sum_prt_)[rr] += (*incrsum_prt_)[rr];
      (*sum_cs_delta_sq_prt_)[rr] += (*incrsum_cs_delta_sq_prt_)[rr];
      (*sumdiffeff_)[rr] += (*incrsumdiffeff_)[rr];
      (*sum_ci_)[rr] += (*incrsum_ci_)[rr];
      (*sum_ci_delta_sq_)[rr] += (*incrsum_ci_delta_sq_)[rr];
    }
  }

  return;
}  // TurbulenceStatisticsCha::DoScatraTimeSample


/*----------------------------------------------------------------------*

          Compute in plane means of u,u^2 etc. (integral version)

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::evaluate_integral_mean_values_in_planes()
{
  //----------------------------------------------------------------------
  // loop elements and perform integration over homogeneous plane

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<FLD::Action>("action", FLD::calc_turbulence_statistics);

  // choose what to assemble
  eleparams.set("assemble matrix 1", false);
  eleparams.set("assemble matrix 2", false);
  eleparams.set("assemble vector 1", false);
  eleparams.set("assemble vector 2", false);
  eleparams.set("assemble vector 3", false);

  // set parameter list
  eleparams.set("normal direction to homogeneous plane", dim_);
  eleparams.set("coordinate vector for hom. planes", planecoordinates_);

  // in case of simultaneous assembly with inflow channel we have to decide on element level,
  // if this element is taken into account or not
  const double dummy = 99999.0;
  if (inflowchannel_)
    eleparams.set("INFLOW_CHA_SIDE", inflowmax_);
  else
    eleparams.set("INFLOW_CHA_SIDE", dummy);
  // set size of vectors
  int size = sumu_->size();

  // generate processor local result vectors
  Teuchos::RCP<std::vector<double>> locarea = Teuchos::RCP(new std::vector<double>);
  locarea->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumu = Teuchos::RCP(new std::vector<double>);
  locsumu->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumv = Teuchos::RCP(new std::vector<double>);
  locsumv->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumw = Teuchos::RCP(new std::vector<double>);
  locsumw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsump = Teuchos::RCP(new std::vector<double>);
  locsump->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumsqu = Teuchos::RCP(new std::vector<double>);
  locsumsqu->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumsqv = Teuchos::RCP(new std::vector<double>);
  locsumsqv->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumsqw = Teuchos::RCP(new std::vector<double>);
  locsumsqw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumuv = Teuchos::RCP(new std::vector<double>);
  locsumuv->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumuw = Teuchos::RCP(new std::vector<double>);
  locsumuw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumvw = Teuchos::RCP(new std::vector<double>);
  locsumvw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumsqp = Teuchos::RCP(new std::vector<double>);
  locsumsqp->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globarea = Teuchos::RCP(new std::vector<double>);
  globarea->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumu = Teuchos::RCP(new std::vector<double>);
  globsumu->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumv = Teuchos::RCP(new std::vector<double>);
  globsumv->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumw = Teuchos::RCP(new std::vector<double>);
  globsumw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsump = Teuchos::RCP(new std::vector<double>);
  globsump->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumsqu = Teuchos::RCP(new std::vector<double>);
  globsumsqu->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumsqv = Teuchos::RCP(new std::vector<double>);
  globsumsqv->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumsqw = Teuchos::RCP(new std::vector<double>);
  globsumsqw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumuv = Teuchos::RCP(new std::vector<double>);
  globsumuv->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumuw = Teuchos::RCP(new std::vector<double>);
  globsumuw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumvw = Teuchos::RCP(new std::vector<double>);
  globsumvw->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumsqp = Teuchos::RCP(new std::vector<double>);
  globsumsqp->resize(size, 0.0);

  // communicate pointers to the result vectors to the element
  eleparams.set("element layer area", locarea);
  eleparams.set("mean velocity u", locsumu);
  eleparams.set("mean velocity v", locsumv);
  eleparams.set("mean velocity w", locsumw);
  eleparams.set("mean pressure p", locsump);

  eleparams.set("mean value u^2", locsumsqu);
  eleparams.set("mean value v^2", locsumsqv);
  eleparams.set("mean value w^2", locsumsqw);
  eleparams.set("mean value uv", locsumuv);
  eleparams.set("mean value uw", locsumuw);
  eleparams.set("mean value vw", locsumvw);
  eleparams.set("mean value p^2", locsumsqp);

  // counts the number of elements in the lowest homogeneous plane
  // (the number is the same for all planes, since we use a structured
  //  cartesian mesh)
  int locprocessedeles = 0;

  eleparams.set("count processed elements", &locprocessedeles);

  // also set the xwall params if necessary

  if (myxwall_ != Teuchos::null) myxwall_->set_x_wall_params(eleparams);


  // set vector values needed by elements
  discret_->clear_state();
  discret_->set_state("u and p (n+1,converged)", meanvelnp_);
  if (alefluid_)
  {
    discret_->set_state("dispnp", dispnp_);
  }

  // call loop over elements
  discret_->evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->clear_state();


  //----------------------------------------------------------------------
  // add contributions from all processors
  discret_->get_comm().SumAll(locarea->data(), globarea->data(), size);

  discret_->get_comm().SumAll(locsumu->data(), globsumu->data(), size);
  discret_->get_comm().SumAll(locsumv->data(), globsumv->data(), size);
  discret_->get_comm().SumAll(locsumw->data(), globsumw->data(), size);
  discret_->get_comm().SumAll(locsump->data(), globsump->data(), size);

  discret_->get_comm().SumAll(locsumsqu->data(), globsumsqu->data(), size);
  discret_->get_comm().SumAll(locsumsqv->data(), globsumsqv->data(), size);
  discret_->get_comm().SumAll(locsumsqw->data(), globsumsqw->data(), size);
  discret_->get_comm().SumAll(locsumuv->data(), globsumuv->data(), size);
  discret_->get_comm().SumAll(locsumuw->data(), globsumuw->data(), size);
  discret_->get_comm().SumAll(locsumvw->data(), globsumvw->data(), size);
  discret_->get_comm().SumAll(locsumsqp->data(), globsumsqp->data(), size);


  //----------------------------------------------------------------------
  // the sums are divided by the layers area to get the area average

  Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
      dynamic_cast<Core::FE::Nurbs::NurbsDiscretization*>(&(*discret_));

  if (nurbsdis == nullptr)
  {
    discret_->get_comm().SumAll(&locprocessedeles, &numele_, 1);
  }
  else
  {
    // get nurbs dis' element numbers
    std::vector<int> nele_x_mele_x_lele(nurbsdis->return_nele_x_mele_x_lele(0));

    numele_ = nele_x_mele_x_lele[0] * nele_x_mele_x_lele[2];
  }


  for (unsigned i = 0; i < planecoordinates_->size(); ++i)
  {
    // get average element size
    (*globarea)[i] /= numele_;

    (*sumu_)[i] += (*globsumu)[i] / (*globarea)[i];
    (*sumv_)[i] += (*globsumv)[i] / (*globarea)[i];
    (*sumw_)[i] += (*globsumw)[i] / (*globarea)[i];
    (*sump_)[i] += (*globsump)[i] / (*globarea)[i];

    (*sumsqu_)[i] += (*globsumsqu)[i] / (*globarea)[i];
    (*sumsqv_)[i] += (*globsumsqv)[i] / (*globarea)[i];
    (*sumsqw_)[i] += (*globsumsqw)[i] / (*globarea)[i];
    (*sumuv_)[i] += (*globsumuv)[i] / (*globarea)[i];
    (*sumuw_)[i] += (*globsumuw)[i] / (*globarea)[i];
    (*sumvw_)[i] += (*globsumvw)[i] / (*globarea)[i];
    (*sumsqp_)[i] += (*globsumsqp)[i] / (*globarea)[i];
  }

  return;

}  // TurbulenceStatisticsCha::evaluate_integral_mean_values_in_planes()

/*----------------------------------------------------------------------*

         Compute in-plane means of u,u^2 etc. (integral version)
                     for low-Mach-number flow

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::evaluate_loma_integral_mean_values_in_planes(const double eosfac)
{
  //----------------------------------------------------------------------
  // loop elements and perform integration over homogeneous plane

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<FLD::Action>("action", FLD::calc_loma_statistics);

  // choose what to assemble
  eleparams.set("assemble matrix 1", false);
  eleparams.set("assemble matrix 2", false);
  eleparams.set("assemble vector 1", false);
  eleparams.set("assemble vector 2", false);
  eleparams.set("assemble vector 3", false);

  // set parameter list
  eleparams.set("normal direction to homogeneous plane", dim_);
  eleparams.set("coordinate vector for hom. planes", planecoordinates_);

  // set size of vectors
  int size = sumu_->size();

  // generate processor local result vectors
  Teuchos::RCP<std::vector<double>> locarea = Teuchos::RCP(new std::vector<double>);
  locarea->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumu = Teuchos::RCP(new std::vector<double>);
  locsumu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumv = Teuchos::RCP(new std::vector<double>);
  locsumv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumw = Teuchos::RCP(new std::vector<double>);
  locsumw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsump = Teuchos::RCP(new std::vector<double>);
  locsump->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumrho = Teuchos::RCP(new std::vector<double>);
  locsumrho->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumT = Teuchos::RCP(new std::vector<double>);
  locsumT->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumrhou = Teuchos::RCP(new std::vector<double>);
  locsumrhou->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumrhouT = Teuchos::RCP(new std::vector<double>);
  locsumrhouT->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumsqu = Teuchos::RCP(new std::vector<double>);
  locsumsqu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqv = Teuchos::RCP(new std::vector<double>);
  locsumsqv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqw = Teuchos::RCP(new std::vector<double>);
  locsumsqw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqp = Teuchos::RCP(new std::vector<double>);
  locsumsqp->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqrho = Teuchos::RCP(new std::vector<double>);
  locsumsqrho->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqT = Teuchos::RCP(new std::vector<double>);
  locsumsqT->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumuv = Teuchos::RCP(new std::vector<double>);
  locsumuv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumuw = Teuchos::RCP(new std::vector<double>);
  locsumuw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumvw = Teuchos::RCP(new std::vector<double>);
  locsumvw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumuT = Teuchos::RCP(new std::vector<double>);
  locsumuT->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumvT = Teuchos::RCP(new std::vector<double>);
  locsumvT->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumwT = Teuchos::RCP(new std::vector<double>);
  locsumwT->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globarea = Teuchos::RCP(new std::vector<double>);
  globarea->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumu = Teuchos::RCP(new std::vector<double>);
  globsumu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumv = Teuchos::RCP(new std::vector<double>);
  globsumv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumw = Teuchos::RCP(new std::vector<double>);
  globsumw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsump = Teuchos::RCP(new std::vector<double>);
  globsump->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumrho = Teuchos::RCP(new std::vector<double>);
  globsumrho->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumT = Teuchos::RCP(new std::vector<double>);
  globsumT->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumrhou = Teuchos::RCP(new std::vector<double>);
  globsumrhou->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumrhouT = Teuchos::RCP(new std::vector<double>);
  globsumrhouT->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumsqu = Teuchos::RCP(new std::vector<double>);
  globsumsqu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqv = Teuchos::RCP(new std::vector<double>);
  globsumsqv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqw = Teuchos::RCP(new std::vector<double>);
  globsumsqw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqp = Teuchos::RCP(new std::vector<double>);
  globsumsqp->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqrho = Teuchos::RCP(new std::vector<double>);
  globsumsqrho->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqT = Teuchos::RCP(new std::vector<double>);
  globsumsqT->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumuv = Teuchos::RCP(new std::vector<double>);
  globsumuv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumuw = Teuchos::RCP(new std::vector<double>);
  globsumuw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumvw = Teuchos::RCP(new std::vector<double>);
  globsumvw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumuT = Teuchos::RCP(new std::vector<double>);
  globsumuT->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumvT = Teuchos::RCP(new std::vector<double>);
  globsumvT->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumwT = Teuchos::RCP(new std::vector<double>);
  globsumwT->resize(size, 0.0);

  // communicate pointers to the result vectors to the element
  eleparams.set("element layer area", locarea);

  eleparams.set("mean velocity u", locsumu);
  eleparams.set("mean velocity v", locsumv);
  eleparams.set("mean velocity w", locsumw);
  eleparams.set("mean pressure p", locsump);
  eleparams.set("mean density rho", locsumrho);
  eleparams.set("mean temperature T", locsumT);
  eleparams.set("mean momentum rho*u", locsumrhou);
  eleparams.set("mean rho*u*T", locsumrhouT);

  eleparams.set("mean value u^2", locsumsqu);
  eleparams.set("mean value v^2", locsumsqv);
  eleparams.set("mean value w^2", locsumsqw);
  eleparams.set("mean value p^2", locsumsqp);
  eleparams.set("mean value rho^2", locsumsqrho);
  eleparams.set("mean value T^2", locsumsqT);

  eleparams.set("mean value uv", locsumuv);
  eleparams.set("mean value uw", locsumuw);
  eleparams.set("mean value vw", locsumvw);
  eleparams.set("mean value uT", locsumuT);
  eleparams.set("mean value vT", locsumvT);
  eleparams.set("mean value wT", locsumwT);

  // counts the number of elements in the lowest homogeneous plane
  // (the number is the same for all planes, since we use a structured
  //  cartesian mesh)
  int locprocessedeles = 0;

  eleparams.set("count processed elements", &locprocessedeles);

  // factor for equation of state
  eleparams.set("eos factor", eosfac);

  // set vector values needed by elements
  discret_->clear_state();
  discret_->set_state("u and p (n+1,converged)", meanvelnp_);
  discret_->set_state("scalar (n+1,converged)", meanscanp_);

  // call loop over elements
  discret_->evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->clear_state();


  //----------------------------------------------------------------------
  // add contributions from all processors
  discret_->get_comm().SumAll(locarea->data(), globarea->data(), size);

  discret_->get_comm().SumAll(locsumu->data(), globsumu->data(), size);
  discret_->get_comm().SumAll(locsumv->data(), globsumv->data(), size);
  discret_->get_comm().SumAll(locsumw->data(), globsumw->data(), size);
  discret_->get_comm().SumAll(locsump->data(), globsump->data(), size);
  discret_->get_comm().SumAll(locsumrho->data(), globsumrho->data(), size);
  discret_->get_comm().SumAll(locsumT->data(), globsumT->data(), size);
  discret_->get_comm().SumAll(locsumrhou->data(), globsumrhou->data(), size);
  discret_->get_comm().SumAll(locsumrhouT->data(), globsumrhouT->data(), size);

  discret_->get_comm().SumAll(locsumsqu->data(), globsumsqu->data(), size);
  discret_->get_comm().SumAll(locsumsqv->data(), globsumsqv->data(), size);
  discret_->get_comm().SumAll(locsumsqw->data(), globsumsqw->data(), size);
  discret_->get_comm().SumAll(locsumsqp->data(), globsumsqp->data(), size);
  discret_->get_comm().SumAll(locsumsqrho->data(), globsumsqrho->data(), size);
  discret_->get_comm().SumAll(locsumsqT->data(), globsumsqT->data(), size);

  discret_->get_comm().SumAll(locsumuv->data(), globsumuv->data(), size);
  discret_->get_comm().SumAll(locsumuw->data(), globsumuw->data(), size);
  discret_->get_comm().SumAll(locsumvw->data(), globsumvw->data(), size);
  discret_->get_comm().SumAll(locsumuT->data(), globsumuT->data(), size);
  discret_->get_comm().SumAll(locsumvT->data(), globsumvT->data(), size);
  discret_->get_comm().SumAll(locsumwT->data(), globsumwT->data(), size);


  //----------------------------------------------------------------------
  // the sums are divided by the layers area to get the area average
  discret_->get_comm().SumAll(&locprocessedeles, &numele_, 1);


  for (unsigned i = 0; i < planecoordinates_->size(); ++i)
  {
    // get average element size
    (*globarea)[i] /= numele_;

    (*sumu_)[i] += (*globsumu)[i] / (*globarea)[i];
    (*sumv_)[i] += (*globsumv)[i] / (*globarea)[i];
    (*sumw_)[i] += (*globsumw)[i] / (*globarea)[i];
    (*sump_)[i] += (*globsump)[i] / (*globarea)[i];
    (*sumrho_)[i] += (*globsumrho)[i] / (*globarea)[i];
    (*sum_t_)[i] += (*globsumT)[i] / (*globarea)[i];
    (*sumrhou_)[i] += (*globsumrhou)[i] / (*globarea)[i];
    (*sumrhou_t_)[i] += (*globsumrhouT)[i] / (*globarea)[i];

    (*sumsqu_)[i] += (*globsumsqu)[i] / (*globarea)[i];
    (*sumsqv_)[i] += (*globsumsqv)[i] / (*globarea)[i];
    (*sumsqw_)[i] += (*globsumsqw)[i] / (*globarea)[i];
    (*sumsqp_)[i] += (*globsumsqp)[i] / (*globarea)[i];
    (*sumsqrho_)[i] += (*globsumsqrho)[i] / (*globarea)[i];
    (*sumsq_t_)[i] += (*globsumsqT)[i] / (*globarea)[i];

    (*sumuv_)[i] += (*globsumuv)[i] / (*globarea)[i];
    (*sumuw_)[i] += (*globsumuw)[i] / (*globarea)[i];
    (*sumvw_)[i] += (*globsumvw)[i] / (*globarea)[i];
    (*sumu_t_)[i] += (*globsumuT)[i] / (*globarea)[i];
    (*sumv_t_)[i] += (*globsumvT)[i] / (*globarea)[i];
    (*sumw_t_)[i] += (*globsumwT)[i] / (*globarea)[i];
  }

  return;

}  // TurbulenceStatisticsCha::evaluate_loma_integral_mean_values_in_planes()


/*----------------------------------------------------------------------*

         Compute in-plane means of u,u^2 etc. (integral version)
                     for turbulent passive scalar transport

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::evaluate_scatra_integral_mean_values_in_planes()
{
  //----------------------------------------------------------------------
  // loop elements and perform integration over homogeneous plane

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<FLD::Action>("action", FLD::calc_turbscatra_statistics);

  // choose what to assemble
  eleparams.set("assemble matrix 1", false);
  eleparams.set("assemble matrix 2", false);
  eleparams.set("assemble vector 1", false);
  eleparams.set("assemble vector 2", false);
  eleparams.set("assemble vector 3", false);

  // set parameter list
  eleparams.set("normal direction to homogeneous plane", dim_);
  eleparams.set("coordinate vector for hom. planes", planecoordinates_);

  // set size of vectors
  int size = sumu_->size();

  // generate processor local result vectors
  Teuchos::RCP<std::vector<double>> locarea = Teuchos::RCP(new std::vector<double>);
  locarea->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumu = Teuchos::RCP(new std::vector<double>);
  locsumu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumv = Teuchos::RCP(new std::vector<double>);
  locsumv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumw = Teuchos::RCP(new std::vector<double>);
  locsumw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsump = Teuchos::RCP(new std::vector<double>);
  locsump->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumphi = Teuchos::RCP(new std::vector<double>);
  locsumphi->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumsqu = Teuchos::RCP(new std::vector<double>);
  locsumsqu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqv = Teuchos::RCP(new std::vector<double>);
  locsumsqv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqw = Teuchos::RCP(new std::vector<double>);
  locsumsqw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqp = Teuchos::RCP(new std::vector<double>);
  locsumsqp->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumsqphi = Teuchos::RCP(new std::vector<double>);
  locsumsqphi->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> locsumuv = Teuchos::RCP(new std::vector<double>);
  locsumuv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumuw = Teuchos::RCP(new std::vector<double>);
  locsumuw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumvw = Teuchos::RCP(new std::vector<double>);
  locsumvw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumuphi = Teuchos::RCP(new std::vector<double>);
  locsumuphi->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumvphi = Teuchos::RCP(new std::vector<double>);
  locsumvphi->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> locsumwphi = Teuchos::RCP(new std::vector<double>);
  locsumwphi->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globarea = Teuchos::RCP(new std::vector<double>);
  globarea->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumu = Teuchos::RCP(new std::vector<double>);
  globsumu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumv = Teuchos::RCP(new std::vector<double>);
  globsumv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumw = Teuchos::RCP(new std::vector<double>);
  globsumw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsump = Teuchos::RCP(new std::vector<double>);
  globsump->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumphi = Teuchos::RCP(new std::vector<double>);
  globsumphi->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumsqu = Teuchos::RCP(new std::vector<double>);
  globsumsqu->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqv = Teuchos::RCP(new std::vector<double>);
  globsumsqv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqw = Teuchos::RCP(new std::vector<double>);
  globsumsqw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqp = Teuchos::RCP(new std::vector<double>);
  globsumsqp->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumsqphi = Teuchos::RCP(new std::vector<double>);
  globsumsqphi->resize(size, 0.0);

  Teuchos::RCP<std::vector<double>> globsumuv = Teuchos::RCP(new std::vector<double>);
  globsumuv->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumuw = Teuchos::RCP(new std::vector<double>);
  globsumuw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumvw = Teuchos::RCP(new std::vector<double>);
  globsumvw->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumuphi = Teuchos::RCP(new std::vector<double>);
  globsumuphi->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumvphi = Teuchos::RCP(new std::vector<double>);
  globsumvphi->resize(size, 0.0);
  Teuchos::RCP<std::vector<double>> globsumwphi = Teuchos::RCP(new std::vector<double>);
  globsumwphi->resize(size, 0.0);

  // communicate pointers to the result vectors to the element
  eleparams.set("element layer area", locarea);

  eleparams.set("mean velocity u", locsumu);
  eleparams.set("mean velocity v", locsumv);
  eleparams.set("mean velocity w", locsumw);
  eleparams.set("mean pressure p", locsump);
  eleparams.set("mean scalar phi", locsumphi);

  eleparams.set("mean value u^2", locsumsqu);
  eleparams.set("mean value v^2", locsumsqv);
  eleparams.set("mean value w^2", locsumsqw);
  eleparams.set("mean value p^2", locsumsqp);
  eleparams.set("mean value phi^2", locsumsqphi);

  eleparams.set("mean value uv", locsumuv);
  eleparams.set("mean value uw", locsumuw);
  eleparams.set("mean value vw", locsumvw);
  eleparams.set("mean value uphi", locsumuphi);
  eleparams.set("mean value vphi", locsumvphi);
  eleparams.set("mean value wphi", locsumwphi);

  // counts the number of elements in the lowest homogeneous plane
  // (the number is the same for all planes, since we use a structured
  //  cartesian mesh)
  int locprocessedeles = 0;

  eleparams.set("count processed elements", &locprocessedeles);

  // set vector values needed by elements
  discret_->clear_state();
  discret_->set_state("u and p (n+1,converged)", meanvelnp_);
  discret_->set_state("scalar (n+1,converged)", meanscanp_);

  // call loop over elements
  discret_->evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->clear_state();


  //----------------------------------------------------------------------
  // add contributions from all processors
  discret_->get_comm().SumAll(locarea->data(), globarea->data(), size);

  discret_->get_comm().SumAll(locsumu->data(), globsumu->data(), size);
  discret_->get_comm().SumAll(locsumv->data(), globsumv->data(), size);
  discret_->get_comm().SumAll(locsumw->data(), globsumw->data(), size);
  discret_->get_comm().SumAll(locsump->data(), globsump->data(), size);
  discret_->get_comm().SumAll(locsumphi->data(), globsumphi->data(), size);

  discret_->get_comm().SumAll(locsumsqu->data(), globsumsqu->data(), size);
  discret_->get_comm().SumAll(locsumsqv->data(), globsumsqv->data(), size);
  discret_->get_comm().SumAll(locsumsqw->data(), globsumsqw->data(), size);
  discret_->get_comm().SumAll(locsumsqp->data(), globsumsqp->data(), size);
  discret_->get_comm().SumAll(locsumsqphi->data(), globsumsqphi->data(), size);

  discret_->get_comm().SumAll(locsumuv->data(), globsumuv->data(), size);
  discret_->get_comm().SumAll(locsumuw->data(), globsumuw->data(), size);
  discret_->get_comm().SumAll(locsumvw->data(), globsumvw->data(), size);
  discret_->get_comm().SumAll(locsumuphi->data(), globsumuphi->data(), size);
  discret_->get_comm().SumAll(locsumvphi->data(), globsumvphi->data(), size);
  discret_->get_comm().SumAll(locsumwphi->data(), globsumwphi->data(), size);


  //----------------------------------------------------------------------
  // the sums are divided by the layers area to get the area average
  discret_->get_comm().SumAll(&locprocessedeles, &numele_, 1);


  for (unsigned i = 0; i < planecoordinates_->size(); ++i)
  {
    // get average element size
    (*globarea)[i] /= numele_;

    // renmark:
    // scalar values named phi are stored in values named T (which are
    // taken form the loma case)
    (*sumu_)[i] += (*globsumu)[i] / (*globarea)[i];
    (*sumv_)[i] += (*globsumv)[i] / (*globarea)[i];
    (*sumw_)[i] += (*globsumw)[i] / (*globarea)[i];
    (*sump_)[i] += (*globsump)[i] / (*globarea)[i];
    (*sum_t_)[i] += (*globsumphi)[i] / (*globarea)[i];

    (*sumsqu_)[i] += (*globsumsqu)[i] / (*globarea)[i];
    (*sumsqv_)[i] += (*globsumsqv)[i] / (*globarea)[i];
    (*sumsqw_)[i] += (*globsumsqw)[i] / (*globarea)[i];
    (*sumsqp_)[i] += (*globsumsqp)[i] / (*globarea)[i];
    (*sumsq_t_)[i] += (*globsumsqphi)[i] / (*globarea)[i];

    (*sumuv_)[i] += (*globsumuv)[i] / (*globarea)[i];
    (*sumuw_)[i] += (*globsumuw)[i] / (*globarea)[i];
    (*sumvw_)[i] += (*globsumvw)[i] / (*globarea)[i];
    (*sumu_t_)[i] += (*globsumuphi)[i] / (*globarea)[i];
    (*sumv_t_)[i] += (*globsumvphi)[i] / (*globarea)[i];
    (*sumw_t_)[i] += (*globsumwphi)[i] / (*globarea)[i];
  }

  return;

}  // TurbulenceStatisticsCha::evaluate_scatra_integral_mean_values_in_planes()


/*----------------------------------------------------------------------*

          Compute in plane means of u,u^2 etc. (nodal quantities)

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::evaluate_pointwise_mean_values_in_planes()
{
  int planenum = 0;

  //----------------------------------------------------------------------
  // pointwise multiplication to get squared values

  pointsquaredvelnp_->Multiply(1.0, *meanvelnp_, *meanvelnp_, 0.0);


  //----------------------------------------------------------------------
  // loop planes and calculate pointwise means in each plane

  for (std::vector<double>::iterator plane = planecoordinates_->begin();
       plane != planecoordinates_->end(); ++plane)
  {
    // toggle vectors are one in the position of a dof in this plane,
    // else 0
    toggleu_->PutScalar(0.0);
    togglev_->PutScalar(0.0);
    togglew_->PutScalar(0.0);
    togglep_->PutScalar(0.0);

    // count the number of nodes in plane (required to calc. in plane mean)
    int countnodesinplane = 0;

    //----------------------------------------------------------------------
    // activate toggles for in plane dofs

    for (int nn = 0; nn < discret_->num_my_row_nodes(); ++nn)
    {
      Core::Nodes::Node* node = discret_->l_row_node(nn);

      // if we have an inflow channel problem, the nodes outside the inflow discretization are
      // not in the bounding box -> we don't consider them for averaging
      if (node->x()[0] < (*boundingbox_)(1, 0) + NODETOL and
          node->x()[1] < (*boundingbox_)(1, 1) + NODETOL and
          node->x()[2] < (*boundingbox_)(1, 2) + NODETOL and
          node->x()[0] > (*boundingbox_)(0, 0) - NODETOL and
          node->x()[1] > (*boundingbox_)(0, 1) - NODETOL and
          node->x()[2] > (*boundingbox_)(0, 2) - NODETOL)
      {
        // this node belongs to the plane under consideration
        if (node->x()[dim_] < *plane + 2e-9 && node->x()[dim_] > *plane - 2e-9)
        {
          std::vector<int> dof = discret_->dof(node);
          double one = 1.0;

          toggleu_->ReplaceGlobalValues(1, &one, &(dof[0]));
          togglev_->ReplaceGlobalValues(1, &one, &(dof[1]));
          togglew_->ReplaceGlobalValues(1, &one, &(dof[2]));
          togglep_->ReplaceGlobalValues(1, &one, &(dof[3]));

          // now check whether we have a pbc condition on this node
          std::vector<Core::Conditions::Condition*> mypbc;

          node->get_condition("SurfacePeriodic", mypbc);

          // yes, we have a pbc
          if (mypbc.size() > 0)
          {
            // loop them and check, whether this is a pbc pure master node
            // for all previous conditions
            unsigned ntimesmaster = 0;
            for (auto& numcond : mypbc)
            {
              const auto mymasterslavetoggle =
                  numcond->parameters().get<std::string>("Is slave periodic boundary condition");

              if (mymasterslavetoggle == "Master")
              {
                ++ntimesmaster;
              }  // end is slave?
            }    // end loop this conditions

            if (ntimesmaster != mypbc.size())
            {
              continue;
            }
            // we have a master. Remember this cause we have to extend the patch
            // to the other side...
          }
          countnodesinplane++;
        }
      }
    }

    int countnodesinplaneonallprocs = 0;

    discret_->get_comm().SumAll(&countnodesinplane, &countnodesinplaneonallprocs, 1);

    if (countnodesinplaneonallprocs)
    {
      double inc = 0.0;
      {
        //----------------------------------------------------------------------
        // compute scalar products from velnp and toggle vec to sum up
        // values in this plane
        double local_inc = 0.0;
        for (int rr = 0; rr < meanvelnp_->MyLength(); ++rr)
        {
          local_inc += (*meanvelnp_)[rr] * (*toggleu_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsumu_)[planenum] += inc / countnodesinplaneonallprocs;

        local_inc = 0.0;
        for (int rr = 0; rr < meanvelnp_->MyLength(); ++rr)
        {
          local_inc += (*meanvelnp_)[rr] * (*togglev_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsumv_)[planenum] += inc / countnodesinplaneonallprocs;

        local_inc = 0.0;
        for (int rr = 0; rr < meanvelnp_->MyLength(); ++rr)
        {
          local_inc += (*meanvelnp_)[rr] * (*togglew_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsumw_)[planenum] += inc / countnodesinplaneonallprocs;

        local_inc = 0.0;
        for (int rr = 0; rr < meanvelnp_->MyLength(); ++rr)
        {
          local_inc += (*meanvelnp_)[rr] * (*togglep_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsump_)[planenum] += inc / countnodesinplaneonallprocs;

        //----------------------------------------------------------------------
        // compute scalar products from squaredvelnp and toggle vec to
        // sum up values for second order moments in this plane
        local_inc = 0.0;
        for (int rr = 0; rr < pointsquaredvelnp_->MyLength(); ++rr)
        {
          local_inc += (*pointsquaredvelnp_)[rr] * (*toggleu_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsumsqu_)[planenum] += inc / countnodesinplaneonallprocs;

        local_inc = 0.0;
        for (int rr = 0; rr < pointsquaredvelnp_->MyLength(); ++rr)
        {
          local_inc += (*pointsquaredvelnp_)[rr] * (*togglev_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsumsqv_)[planenum] += inc / countnodesinplaneonallprocs;

        local_inc = 0.0;
        for (int rr = 0; rr < pointsquaredvelnp_->MyLength(); ++rr)
        {
          local_inc += (*pointsquaredvelnp_)[rr] * (*togglew_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsumsqw_)[planenum] += inc / countnodesinplaneonallprocs;

        local_inc = 0.0;
        for (int rr = 0; rr < pointsquaredvelnp_->MyLength(); ++rr)
        {
          local_inc += (*pointsquaredvelnp_)[rr] * (*togglep_)[rr];
        }
        discret_->get_comm().SumAll(&local_inc, &inc, 1);
        (*pointsumsqp_)[planenum] += inc / countnodesinplaneonallprocs;
      }
    }
    planenum++;
  }

  return;
}  // TurbulenceStatisticsCha::evaluate_pointwise_mean_values_in_planes()


/*----------------------------------------------------------------------*

        Add computed dynamic Smagorinsky quantities (Smagorinsky
           constant, effective viscosity and (Cs_delta)^2 used
                      during the computation)

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::add_dynamic_smagorinsky_quantities()
{
  // get sublist of turbulence parameters from the fluid dynamic
  // parameter list --- it is used to transfer data between element
  // and statistics method
  Teuchos::ParameterList* modelparams = &(params_.sublist("TURBULENCE MODEL"));

  // extract values for Cs, Cs_delta_sq_ and visceff from parameterlist
  // the values are stored in vectors --- each component corresponds to
  // one element layer
  Teuchos::RCP<std::vector<double>> global_incr_Cs_sum;
  Teuchos::RCP<std::vector<double>> local_Cs_sum;
  global_incr_Cs_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Cs_sum = modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_sum", Teuchos::null);
  if (local_Cs_sum == Teuchos::null) FOUR_C_THROW("local_Cs_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_Cs_delta_sq_sum;
  Teuchos::RCP<std::vector<double>> local_Cs_delta_sq_sum;
  global_incr_Cs_delta_sq_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Cs_delta_sq_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Cs_delta_sq_sum", Teuchos::null);
  if (local_Cs_delta_sq_sum == Teuchos::null)
    FOUR_C_THROW("local_Cs_delta_sq_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_visceff_sum;
  Teuchos::RCP<std::vector<double>> local_visceff_sum;
  global_incr_visceff_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_visceff_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_visceff_sum", Teuchos::null);
  if (local_visceff_sum == Teuchos::null)
    FOUR_C_THROW("local_visceff_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_Prt_sum;
  Teuchos::RCP<std::vector<double>> local_Prt_sum;
  global_incr_Prt_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Prt_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Prt_sum", Teuchos::null);
  if (local_Prt_sum == Teuchos::null) FOUR_C_THROW("local_Prt_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_Cs_delta_sq_Prt_sum;
  Teuchos::RCP<std::vector<double>> local_Cs_delta_sq_Prt_sum;
  global_incr_Cs_delta_sq_Prt_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Cs_delta_sq_Prt_sum = modelparams->get<Teuchos::RCP<std::vector<double>>>(
      "local_Cs_delta_sq_Prt_sum", Teuchos::null);
  if (local_Cs_delta_sq_Prt_sum == Teuchos::null)
    FOUR_C_THROW("local_Cs_delta_sq_Prt_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_diffeff_sum;
  Teuchos::RCP<std::vector<double>> local_diffeff_sum;
  global_incr_diffeff_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_diffeff_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_diffeff_sum", Teuchos::null);
  if (local_diffeff_sum == Teuchos::null)
    FOUR_C_THROW("local_diffeff_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_Ci_sum;
  Teuchos::RCP<std::vector<double>> local_Ci_sum;
  global_incr_Ci_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Ci_sum = modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Ci_sum", Teuchos::null);
  if (local_Ci_sum == Teuchos::null) FOUR_C_THROW("local_Ci_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_Ci_delta_sq_sum;
  Teuchos::RCP<std::vector<double>> local_Ci_delta_sq_sum;
  global_incr_Ci_delta_sq_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Ci_delta_sq_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Ci_delta_sq_sum", Teuchos::null);
  if (local_Ci_delta_sq_sum == Teuchos::null)
    FOUR_C_THROW("local_Ci_delta_sq_sum==null from parameterlist");

  // now add all the stuff from the different processors
  discret_->get_comm().SumAll(
      local_Cs_sum->data(), global_incr_Cs_sum->data(), local_Cs_sum->size());
  discret_->get_comm().SumAll(local_Cs_delta_sq_sum->data(), global_incr_Cs_delta_sq_sum->data(),
      local_Cs_delta_sq_sum->size());
  discret_->get_comm().SumAll(
      local_visceff_sum->data(), global_incr_visceff_sum->data(), local_visceff_sum->size());
  discret_->get_comm().SumAll(
      local_Prt_sum->data(), global_incr_Prt_sum->data(), local_Prt_sum->size());
  discret_->get_comm().SumAll(local_Cs_delta_sq_Prt_sum->data(),
      global_incr_Cs_delta_sq_Prt_sum->data(), local_Cs_delta_sq_Prt_sum->size());
  discret_->get_comm().SumAll(
      local_diffeff_sum->data(), global_incr_diffeff_sum->data(), local_diffeff_sum->size());
  discret_->get_comm().SumAll(
      local_Ci_sum->data(), global_incr_Ci_sum->data(), local_Ci_sum->size());
  discret_->get_comm().SumAll(local_Ci_delta_sq_sum->data(), global_incr_Ci_delta_sq_sum->data(),
      local_Ci_delta_sq_sum->size());

  // Replace increment to compute average of Smagorinsky Constant, effective
  // viscosity and (Cs_delta)^2
  for (unsigned rr = 0; rr < global_incr_Cs_sum->size(); ++rr)
  {
    (*incrsum_cs_)[rr] = (*global_incr_Cs_sum)[rr];
    (*incrsum_cs_delta_sq_)[rr] = (*global_incr_Cs_delta_sq_sum)[rr];
    (*incrsumvisceff_)[rr] = (*global_incr_visceff_sum)[rr];
    (*incrsum_prt_)[rr] = (*global_incr_Prt_sum)[rr];
    (*incrsum_cs_delta_sq_prt_)[rr] = (*global_incr_Cs_delta_sq_Prt_sum)[rr];
    (*incrsumdiffeff_)[rr] = (*global_incr_diffeff_sum)[rr];
    (*incrsum_ci_)[rr] = (*global_incr_Ci_sum)[rr];
    (*incrsum_ci_delta_sq_)[rr] = (*global_incr_Ci_delta_sq_sum)[rr];
  }

  // reinitialise to zero for next element call
  local_Cs_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Cs_delta_sq_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_visceff_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Prt_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Cs_delta_sq_Prt_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_diffeff_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Ci_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Ci_delta_sq_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));

  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Cs_sum", local_Cs_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>(
      "local_Cs_delta_sq_sum", local_Cs_delta_sq_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_visceff_sum", local_visceff_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Prt_sum", local_Prt_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>(
      "local_Cs_delta_sq_Prt_sum", local_Cs_delta_sq_Prt_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_diffeff_sum", local_diffeff_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Ci_sum", local_Ci_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>(
      "local_Ci_delta_sq_sum", local_Ci_delta_sq_sum);

  return;
}  // FLD::TurbulenceStatisticsCha::add_dynamic_smagorinsky_quantities


/*----------------------------------------------------------------------*

            Add parameters of multifractal
                 subgrid-scales model

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::add_model_params_multifractal(
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> velnp,
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> fsvelnp, const bool withscatra)
{
  // action for elements
  Teuchos::ParameterList paramsele;
  paramsele.set<FLD::Action>("action", FLD::calc_model_params_mfsubgr_scales);
  paramsele.sublist("MULTIFRACTAL SUBGRID SCALES") = params_.sublist("MULTIFRACTAL SUBGRID SCALES");
  paramsele.set("scalar", withscatra);
  if (withscatra)
  {
    paramsele.set("scnum", scnum_);
  }

  // vectors containing processor local values to sum element
  // contributions on this proc
  Teuchos::RCP<std::vector<double>> local_N_stream_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_N_normal_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_N_span_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_B_stream_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_B_normal_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_B_span_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_Csgs_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_sgvisc_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_Nphi_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_Dphi_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  Teuchos::RCP<std::vector<double>> local_Csgs_phi_sum =
      Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));


  // store them in parameterlist for access on the element
  Teuchos::ParameterList* modelparams = &(paramsele.sublist("TURBULENCE MODEL"));

  modelparams->set<Teuchos::RCP<std::vector<double>>>("planecoords", nodeplanes_);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_N_stream_sum", local_N_stream_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_N_normal_sum", local_N_normal_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_N_span_sum", local_N_span_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_B_stream_sum", local_B_stream_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_B_normal_sum", local_B_normal_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_B_span_sum", local_B_span_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Csgs_sum", local_Csgs_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_sgvisc_sum", local_sgvisc_sum);
  if (withscatra)
  {
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Nphi_sum", local_Nphi_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Dphi_sum", local_Dphi_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Csgs_phi_sum", local_Csgs_phi_sum);
  }

  // set state vectors for element call
  discret_->clear_state();
  discret_->set_state("velnp", velnp);
  if (fsvelnp == Teuchos::null) FOUR_C_THROW("Haven't got fine-scale velocity!");
  discret_->set_state("fsvelnp", fsvelnp);

  // call loop over elements to compute means
  discret_->evaluate(
      paramsele, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  discret_->clear_state();

  // extract values for N, B, Csgs and sgvisc from parameter list
  // the values are stored in vectors --- each component corresponds to
  // one element layer
  Teuchos::RCP<std::vector<double>> global_incr_N_stream_sum;
  global_incr_N_stream_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_N_stream_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_N_stream_sum", Teuchos::null);
  if (local_N_stream_sum == Teuchos::null)
    FOUR_C_THROW("local_N_stream_sum==null from parameterlist");
  Teuchos::RCP<std::vector<double>> global_incr_N_normal_sum;
  global_incr_N_normal_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_N_normal_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_N_normal_sum", Teuchos::null);
  if (local_N_normal_sum == Teuchos::null)
    FOUR_C_THROW("local_N_normal_sum==null from parameterlist");
  Teuchos::RCP<std::vector<double>> global_incr_N_span_sum;
  global_incr_N_span_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_N_span_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_N_span_sum", Teuchos::null);
  if (local_N_span_sum == Teuchos::null) FOUR_C_THROW("local_N_span_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_B_stream_sum;
  global_incr_B_stream_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_B_stream_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_B_stream_sum", Teuchos::null);
  if (local_B_stream_sum == Teuchos::null)
    FOUR_C_THROW("local_B_stream_sum==null from parameterlist");
  Teuchos::RCP<std::vector<double>> global_incr_B_normal_sum;
  global_incr_B_normal_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_B_normal_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_B_normal_sum", Teuchos::null);
  if (local_B_normal_sum == Teuchos::null)
    FOUR_C_THROW("local_B_normal_sum==null from parameterlist");
  Teuchos::RCP<std::vector<double>> global_incr_B_span_sum;
  global_incr_B_span_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_B_span_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_B_span_sum", Teuchos::null);
  if (local_B_span_sum == Teuchos::null) FOUR_C_THROW("local_B_span_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_Csgs_sum;
  global_incr_Csgs_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Csgs_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Csgs_sum", Teuchos::null);
  if (local_Csgs_sum == Teuchos::null) FOUR_C_THROW("local_Csgs_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_sgvisc_sum;
  global_incr_sgvisc_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_sgvisc_sum =
      modelparams->get<Teuchos::RCP<std::vector<double>>>("local_sgvisc_sum", Teuchos::null);
  if (local_sgvisc_sum == Teuchos::null) FOUR_C_THROW("local_sgvsic_sum==null from parameterlist");

  Teuchos::RCP<std::vector<double>> global_incr_Nphi_sum;
  global_incr_Nphi_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));

  Teuchos::RCP<std::vector<double>> global_incr_Dphi_sum;
  global_incr_Dphi_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));

  Teuchos::RCP<std::vector<double>> global_incr_Csgs_phi_sum;
  global_incr_Csgs_phi_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));

  if (withscatra)
  {
    local_Nphi_sum =
        modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Nphi_sum", Teuchos::null);
    if (local_Nphi_sum == Teuchos::null) FOUR_C_THROW("local_Nphi_sum==null from parameterlist");

    local_Dphi_sum =
        modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Dphi_sum", Teuchos::null);
    if (local_Dphi_sum == Teuchos::null) FOUR_C_THROW("local_Dphi_sum==null from parameterlist");

    local_Csgs_phi_sum =
        modelparams->get<Teuchos::RCP<std::vector<double>>>("local_Csgs_phi_sum", Teuchos::null);
    if (local_Csgs_phi_sum == Teuchos::null)
      FOUR_C_THROW("local_Csgs_phi_sum==null from parameterlist");
  }

  // now add all the stuff from the different processors
  discret_->get_comm().SumAll(
      local_N_stream_sum->data(), global_incr_N_stream_sum->data(), local_N_stream_sum->size());
  discret_->get_comm().SumAll(
      local_N_normal_sum->data(), global_incr_N_normal_sum->data(), local_N_normal_sum->size());
  discret_->get_comm().SumAll(
      local_N_span_sum->data(), global_incr_N_span_sum->data(), local_N_span_sum->size());
  discret_->get_comm().SumAll(
      local_B_stream_sum->data(), global_incr_B_stream_sum->data(), local_B_stream_sum->size());
  discret_->get_comm().SumAll(
      local_B_normal_sum->data(), global_incr_B_normal_sum->data(), local_B_normal_sum->size());
  discret_->get_comm().SumAll(
      local_B_span_sum->data(), global_incr_B_span_sum->data(), local_B_span_sum->size());
  discret_->get_comm().SumAll(
      local_Csgs_sum->data(), global_incr_Csgs_sum->data(), local_Csgs_sum->size());
  discret_->get_comm().SumAll(
      local_sgvisc_sum->data(), global_incr_sgvisc_sum->data(), local_sgvisc_sum->size());
  if (withscatra)
  {
    discret_->get_comm().SumAll(
        local_Nphi_sum->data(), global_incr_Nphi_sum->data(), local_Nphi_sum->size());
    discret_->get_comm().SumAll(
        local_Dphi_sum->data(), global_incr_Dphi_sum->data(), local_Dphi_sum->size());
    discret_->get_comm().SumAll(
        local_Csgs_phi_sum->data(), global_incr_Csgs_phi_sum->data(), local_Csgs_phi_sum->size());
  }

  // Replace increment to compute average of parameters N and B as well as
  // subgrid viscosity
  for (unsigned rr = 0; rr < global_incr_N_stream_sum->size(); ++rr)
  {
    (*incrsum_n_stream_)[rr] = (*global_incr_N_stream_sum)[rr];
    (*incrsum_n_normal_)[rr] = (*global_incr_N_normal_sum)[rr];
    (*incrsum_n_span_)[rr] = (*global_incr_N_span_sum)[rr];
    (*incrsum_b_stream_)[rr] = (*global_incr_B_normal_sum)[rr];
    (*incrsum_b_normal_)[rr] = (*global_incr_B_normal_sum)[rr];
    (*incrsum_b_span_)[rr] = (*global_incr_B_span_sum)[rr];
    (*incrsum_csgs_)[rr] = (*global_incr_Csgs_sum)[rr];
    (*incrsumsgvisc_)[rr] = (*global_incr_sgvisc_sum)[rr];
    if (withscatra)
    {
      (*incrsum_nphi_)[rr] = (*global_incr_Nphi_sum)[rr];
      (*incrsum_dphi_)[rr] = (*global_incr_Dphi_sum)[rr];
      (*incrsum_csgs_phi_)[rr] = (*global_incr_Csgs_phi_sum)[rr];
    }
  }

  // reinitialize to zero for next element call
  local_N_stream_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_N_normal_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_N_span_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_B_stream_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_B_normal_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_B_span_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_Csgs_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  local_sgvisc_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));

  if (withscatra)
  {
    local_Nphi_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    local_Dphi_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
    local_Csgs_phi_sum = Teuchos::RCP(new std::vector<double>(nodeplanes_->size() - 1, 0.0));
  }

  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_N_stream_sum", local_N_stream_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_N_normal_sum", local_N_normal_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_N_span_sum", local_N_span_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_B_stream_sum", local_B_stream_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_B_normal_sum", local_B_normal_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_B_span_sum", local_B_span_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Csgs_sum", local_Csgs_sum);
  modelparams->set<Teuchos::RCP<std::vector<double>>>("local_sgvisc_sum", local_sgvisc_sum);
  if (withscatra)
  {
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Nphi_sum", local_Nphi_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Dphi_sum", local_Dphi_sum);
    modelparams->set<Teuchos::RCP<std::vector<double>>>("local_Csgs_phi_sum", local_Csgs_phi_sum);
  }

  // add increment of last iteration to the sum of all values
  for (unsigned rr = 0; rr < (*incrsum_n_stream_).size(); ++rr)
  {
    (*sum_n_stream_)[rr] += (*incrsum_n_stream_)[rr];
    (*sum_n_normal_)[rr] += (*incrsum_n_normal_)[rr];
    (*sum_n_span_)[rr] += (*incrsum_n_span_)[rr];
    (*sum_b_stream_)[rr] += (*incrsum_b_stream_)[rr];
    (*sum_b_normal_)[rr] += (*incrsum_b_normal_)[rr];
    (*sum_b_span_)[rr] += (*incrsum_b_span_)[rr];
    (*sum_csgs_)[rr] += (*incrsum_csgs_)[rr];
    (*sumsgvisc_)[rr] += (*incrsumsgvisc_)[rr];
    if (withscatra)
    {
      (*sum_nphi_)[rr] += (*incrsum_nphi_)[rr];
      (*sum_dphi_)[rr] += (*incrsum_dphi_)[rr];
      (*sum_csgs_phi_)[rr] += (*incrsum_csgs_phi_)[rr];
    }
  }

  return;
}  // FLD::TurbulenceStatisticsCha::add_model_params_multifractal();


void FLD::TurbulenceStatisticsCha::evaluate_residuals(
    std::map<std::string, Teuchos::RCP<Core::LinAlg::Vector<double>>> statevecs,
    std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> statetenss, const double thermpressaf,
    const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
    std::map<std::string, Teuchos::RCP<Core::LinAlg::Vector<double>>> scatrastatevecs,
    std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> scatrafieldvecs)
{
  if (subgrid_dissipation_)
  {
    //--------------------------------------------------------------------
    // set parameter list (time integration)

    // action for elements
    eleparams_.set<FLD::Action>("action", FLD::calc_dissipation);

    // add velafgrad
    Teuchos::ParameterList* stabparams = &(params_.sublist("RESIDUAL-BASED STABILIZATION"));
    if (stabparams->get<bool>("Reconstruct_Sec_Der"))
    {
      for (std::map<std::string, Teuchos::RCP<Core::LinAlg::Vector<double>>>::iterator state =
               statevecs.begin();
           state != statevecs.end(); ++state)
      {
        if (state->first == "velaf")
        {  // project_gradient_and_set_param decides, if we want to project something or not
          FLD::UTILS::project_gradient_and_set_param(
              discret_, eleparams_, state->second, "velafgrad", alefluid_);
          break;
        }
      }
    }
    eleparams_.sublist("TURBULENCE MODEL") = params_.sublist("TURBULENCE MODEL");
    eleparams_.set<double>("thermpress at n+alpha_F/n+1", thermpressaf);
    eleparams_.set<double>("thermpress at n+alpha_M/n", thermpressam);
    eleparams_.set<double>("thermpressderiv at n+alpha_F/n+1", thermpressdtaf);
    eleparams_.set<double>("thermpressderiv at n+alpha_M/n+1", thermpressdtam);

    // set state vectors for element call
    for (std::map<std::string, Teuchos::RCP<Core::LinAlg::Vector<double>>>::iterator state =
             statevecs.begin();
         state != statevecs.end(); ++state)
    {
      discret_->set_state(state->first, state->second);
    }

    if (myxwall_ != Teuchos::null) myxwall_->set_x_wall_params(eleparams_);

    // call loop over elements to compute means
    discret_->evaluate(
        eleparams_, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

    discret_->clear_state();

    // do we also have a scalar field
    if (scatradiscret_ != Teuchos::null)
    {
      // add dissipation and residuals of scalar field

      // set action for elements
      Core::UTILS::add_enum_class_to_parameter_list<ScaTra::Action>(
          "action", ScaTra::Action::calc_dissipation, scatraeleparams_);
      // set parameters required for evaluation of residuals, etc.
      scatraeleparams_.set<double>("time-step length", scatraparams_->get<double>("TIMESTEP"));
      scatraeleparams_.set<int>("fs subgrid diffusivity",
          Teuchos::getIntegralValue<Inpar::ScaTra::FSSUGRDIFF>(*scatraparams_, "FSSUGRDIFF"));
      scatraeleparams_.sublist("TURBULENCE MODEL") =
          scatraextraparams_->sublist("TURBULENCE MODEL");
      scatraeleparams_.sublist("SUBGRID VISCOSITY") =
          scatraextraparams_->sublist("SUBGRID VISCOSITY");
      scatraeleparams_.sublist("MULTIFRACTAL SUBGRID SCALES") =
          scatraextraparams_->sublist("MULTIFRACTAL SUBGRID SCALES");
      scatraeleparams_.set<bool>("update material",
          (&(scatraextraparams_->sublist("LOMA")))->get<bool>("update material", false));
      scatraeleparams_.sublist("STABILIZATION") = scatraparams_->sublist("STABILIZATION");
      scatraeleparams_.sublist("TIME INTEGRATION") = *scatratimeparams_;
      // remark: this is the thermodynamic pressure taken form the fluid field
      //         since the scatra field is solved after the fluid field, the thermodynamic pressure
      //         may slightly change
      //         however, the error is expected to be small
      scatraeleparams_.set<double>("time derivative of thermodynamic pressure", thermpressdtaf);
      scatraeleparams_.set<double>("thermodynamic pressure", thermpressaf);
      scatraeleparams_.set<double>("thermodynamic pressure at n+alpha_M", thermpressam);

      // set state vectors for element call
      for (std::map<std::string, Teuchos::RCP<Core::LinAlg::Vector<double>>>::iterator state =
               scatrastatevecs.begin();
           state != scatrastatevecs.end(); ++state)
      {
        scatradiscret_->set_state(state->first, state->second);
      }

      // call loop over elements to compute means
      scatradiscret_->evaluate(scatraeleparams_, Teuchos::null, Teuchos::null, Teuchos::null,
          Teuchos::null, Teuchos::null);

      scatradiscret_->clear_state();
    }

    // ------------------------------------------------
    // get results from element call via parameter list
    Teuchos::RCP<std::vector<double>> local_vol =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrvol");

    Teuchos::RCP<std::vector<double>> local_incrhk =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrhk");
    Teuchos::RCP<std::vector<double>> local_incrhbazilevs =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrhbazilevs");
    Teuchos::RCP<std::vector<double>> local_incrstrle =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrstrle");
    Teuchos::RCP<std::vector<double>> local_incrgradle =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrgradle");

    Teuchos::RCP<std::vector<double>> local_incrtauC =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrtauC");
    Teuchos::RCP<std::vector<double>> local_incrtauM =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrtauM");

    Teuchos::RCP<std::vector<double>> local_incrmk =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrmk");

    Teuchos::RCP<std::vector<double>> local_incrres =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrres");
    Teuchos::RCP<std::vector<double>> local_incrres_sq =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrres_sq");
    Teuchos::RCP<std::vector<double>> local_incrabsres =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrabsres");
    Teuchos::RCP<std::vector<double>> local_incrtauinvsvel =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrtauinvsvel");
    Teuchos::RCP<std::vector<double>> local_incrsvelaf =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrsvelaf");
    Teuchos::RCP<std::vector<double>> local_incrsvelaf_sq =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrsvelaf_sq");
    Teuchos::RCP<std::vector<double>> local_incrabssvelaf =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrabssvelaf");
    Teuchos::RCP<std::vector<double>> local_incrresC =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrresC");
    Teuchos::RCP<std::vector<double>> local_incrresC_sq =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrresC_sq");
    Teuchos::RCP<std::vector<double>> local_incrspressnp =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrspressnp");
    Teuchos::RCP<std::vector<double>> local_incrspressnp_sq =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrspressnp_sq");

    Teuchos::RCP<std::vector<double>> local_incr_eps_visc =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_visc");
    Teuchos::RCP<std::vector<double>> local_incr_eps_conv =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_conv");
    Teuchos::RCP<std::vector<double>> local_incr_eps_eddyvisc =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_eddyvisc");
    Teuchos::RCP<std::vector<double>> local_incr_eps_avm3 =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_avm3");
    Teuchos::RCP<std::vector<double>> local_incr_eps_mfs =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_mfs");
    Teuchos::RCP<std::vector<double>> local_incr_eps_mfscross =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_mfscross");
    Teuchos::RCP<std::vector<double>> local_incr_eps_mfsrey =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_mfsrey");
    Teuchos::RCP<std::vector<double>> local_incr_eps_supg =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_supg");
    Teuchos::RCP<std::vector<double>> local_incr_eps_cross =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_cross");
    Teuchos::RCP<std::vector<double>> local_incr_eps_rey =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_rey");
    Teuchos::RCP<std::vector<double>> local_incr_eps_graddiv =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_graddiv");
    Teuchos::RCP<std::vector<double>> local_incr_eps_pspg =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_eps_pspg");

    Teuchos::RCP<std::vector<double>> local_incrcrossstress =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrcrossstress");
    Teuchos::RCP<std::vector<double>> local_incrreystress =
        eleparams_.get<Teuchos::RCP<std::vector<double>>>("incrreystress");

    int presize = local_incrresC->size();
    int velsize = local_incrres->size();
    int stresssize = local_incrcrossstress->size();

    //--------------------------------------------------
    // vectors to sum over all procs

    // volume of element layers
    Teuchos::RCP<std::vector<double>> global_vol;
    global_vol = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // element sizes of element layers
    Teuchos::RCP<std::vector<double>> global_incrhk;
    global_incrhk = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // element sizes in Bazilevs parameter, viscous regime in element layers
    Teuchos::RCP<std::vector<double>> global_incrhbazilevs;
    global_incrhbazilevs = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // element sizes of element stream length
    Teuchos::RCP<std::vector<double>> global_incrstrle;
    global_incrstrle = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // element sizes based on gradient length
    Teuchos::RCP<std::vector<double>> global_incrgradle;
    global_incrgradle = Teuchos::RCP(new std::vector<double>(presize, 0.0));


    // (in plane) averaged values of tauM/tauC

    Teuchos::RCP<std::vector<double>> global_incrtauM;
    global_incrtauM = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrtauC;
    global_incrtauC = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // mk of element layers
    Teuchos::RCP<std::vector<double>> global_incrmk;
    global_incrmk = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of resM (^2) (abs)

    Teuchos::RCP<std::vector<double>> global_incrres;
    global_incrres = Teuchos::RCP(new std::vector<double>(velsize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrres_sq;
    global_incrres_sq = Teuchos::RCP(new std::vector<double>(velsize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrtauinvsvel;
    global_incrtauinvsvel = Teuchos::RCP(new std::vector<double>(velsize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrabsres;
    global_incrabsres = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of svelaf (^2) (abs)

    Teuchos::RCP<std::vector<double>> global_incrsvelaf;
    global_incrsvelaf = Teuchos::RCP(new std::vector<double>(velsize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrsvelaf_sq;
    global_incrsvelaf_sq = Teuchos::RCP(new std::vector<double>(velsize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrabssvelaf;
    global_incrabssvelaf = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of resC (^2)

    Teuchos::RCP<std::vector<double>> global_incrresC;
    global_incrresC = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrresC_sq;
    global_incrresC_sq = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of spressnp (^2)

    Teuchos::RCP<std::vector<double>> global_incrspressnp;
    global_incrspressnp = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incrspressnp_sq;
    global_incrspressnp_sq = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation by pspg term

    Teuchos::RCP<std::vector<double>> global_incr_eps_pspg;
    global_incr_eps_pspg = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation by supg term

    Teuchos::RCP<std::vector<double>> global_incr_eps_supg;
    global_incr_eps_supg = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation by cross term

    Teuchos::RCP<std::vector<double>> global_incr_eps_cross;
    global_incr_eps_cross = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation by reynolds term

    Teuchos::RCP<std::vector<double>> global_incr_eps_rey;
    global_incr_eps_rey = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation by continuity stabilisation

    Teuchos::RCP<std::vector<double>> global_incr_eps_graddiv;
    global_incr_eps_graddiv = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation by eddy viscosity

    Teuchos::RCP<std::vector<double>> global_incr_eps_eddyvisc;
    global_incr_eps_eddyvisc = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    Teuchos::RCP<std::vector<double>> global_incr_eps_avm3;
    global_incr_eps_avm3 = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    // (in plane) averaged values of dissipation by mfs
    Teuchos::RCP<std::vector<double>> global_incr_eps_mfs;
    global_incr_eps_mfs = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    Teuchos::RCP<std::vector<double>> global_incr_eps_mfscross;
    global_incr_eps_mfscross = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    Teuchos::RCP<std::vector<double>> global_incr_eps_mfsrey;
    global_incr_eps_mfsrey = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation by viscous forces

    Teuchos::RCP<std::vector<double>> global_incr_eps_visc;
    global_incr_eps_visc = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of dissipation/production by convection

    Teuchos::RCP<std::vector<double>> global_incr_eps_conv;
    global_incr_eps_conv = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    // (in plane) averaged values of subgrid stresses resulting from supg and cross term

    Teuchos::RCP<std::vector<double>> global_incrcrossstress;
    global_incrcrossstress = Teuchos::RCP(new std::vector<double>(stresssize, 0.0));

    // (in plane) averaged values of subgrid stresses resulting from reynolds stresses
    Teuchos::RCP<std::vector<double>> global_incrreystress;
    global_incrreystress = Teuchos::RCP(new std::vector<double>(stresssize, 0.0));


    //--------------------------------------------------
    // global sums

    // compute global sum, volume
    discret_->get_comm().SumAll(local_vol->data(), global_vol->data(), presize);

    // compute global sum, element sizes
    discret_->get_comm().SumAll(local_incrhk->data(), global_incrhk->data(), presize);

    // compute global sum, element sizes in viscous regime, Bazilevs parameter
    discret_->get_comm().SumAll(local_incrhbazilevs->data(), global_incrhbazilevs->data(), presize);

    // compute global sum, element sizes
    discret_->get_comm().SumAll(local_incrstrle->data(), global_incrstrle->data(), presize);

    // compute global sum, gradient based element sizes
    discret_->get_comm().SumAll(local_incrgradle->data(), global_incrgradle->data(), presize);

    // compute global sums, stabilisation parameters
    discret_->get_comm().SumAll(local_incrtauM->data(), global_incrtauM->data(), presize);
    discret_->get_comm().SumAll(local_incrtauC->data(), global_incrtauC->data(), presize);

    // compute global sum, mk
    discret_->get_comm().SumAll(local_incrmk->data(), global_incrmk->data(), presize);

    // compute global sums, momentum equation residuals
    discret_->get_comm().SumAll(local_incrres->data(), global_incrres->data(), velsize);
    discret_->get_comm().SumAll(local_incrres_sq->data(), global_incrres_sq->data(), velsize);
    discret_->get_comm().SumAll(
        local_incrtauinvsvel->data(), global_incrtauinvsvel->data(), velsize);
    discret_->get_comm().SumAll(local_incrabsres->data(), global_incrabsres->data(), presize);

    discret_->get_comm().SumAll(local_incrsvelaf->data(), global_incrsvelaf->data(), velsize);
    discret_->get_comm().SumAll(local_incrsvelaf_sq->data(), global_incrsvelaf_sq->data(), velsize);
    discret_->get_comm().SumAll(local_incrabssvelaf->data(), global_incrabssvelaf->data(), presize);

    // compute global sums, incompressibility residuals
    discret_->get_comm().SumAll(local_incrresC->data(), global_incrresC->data(), presize);
    discret_->get_comm().SumAll(local_incrresC_sq->data(), global_incrresC_sq->data(), presize);

    discret_->get_comm().SumAll(local_incrspressnp->data(), global_incrspressnp->data(), presize);
    discret_->get_comm().SumAll(
        local_incrspressnp_sq->data(), global_incrspressnp_sq->data(), presize);

    // compute global sums, disspiation rates

    discret_->get_comm().SumAll(local_incr_eps_pspg->data(), global_incr_eps_pspg->data(), presize);
    discret_->get_comm().SumAll(local_incr_eps_supg->data(), global_incr_eps_supg->data(), presize);
    discret_->get_comm().SumAll(
        local_incr_eps_cross->data(), global_incr_eps_cross->data(), presize);
    discret_->get_comm().SumAll(local_incr_eps_rey->data(), global_incr_eps_rey->data(), presize);
    discret_->get_comm().SumAll(
        local_incr_eps_graddiv->data(), global_incr_eps_graddiv->data(), presize);
    discret_->get_comm().SumAll(
        local_incr_eps_eddyvisc->data(), global_incr_eps_eddyvisc->data(), presize);
    discret_->get_comm().SumAll(local_incr_eps_visc->data(), global_incr_eps_visc->data(), presize);
    discret_->get_comm().SumAll(local_incr_eps_conv->data(), global_incr_eps_conv->data(), presize);
    discret_->get_comm().SumAll(local_incr_eps_avm3->data(), global_incr_eps_avm3->data(), presize);
    discret_->get_comm().SumAll(local_incr_eps_mfs->data(), global_incr_eps_mfs->data(), presize);
    discret_->get_comm().SumAll(
        local_incr_eps_mfscross->data(), global_incr_eps_mfscross->data(), presize);
    discret_->get_comm().SumAll(
        local_incr_eps_mfsrey->data(), global_incr_eps_mfsrey->data(), presize);

    // compute global sums, subgrid stresses
    discret_->get_comm().SumAll(
        local_incrcrossstress->data(), global_incrcrossstress->data(), stresssize);
    discret_->get_comm().SumAll(
        local_incrreystress->data(), global_incrreystress->data(), stresssize);


    for (int rr = 0; rr < velsize; ++rr)
    {
      (*sumres_)[rr] += (*global_incrres)[rr];
      (*sumres_sq_)[rr] += (*global_incrres_sq)[rr];
      (*sumsvelaf_)[rr] += (*global_incrsvelaf)[rr];
      (*sumsvelaf_sq_)[rr] += (*global_incrsvelaf_sq)[rr];

      (*sumtauinvsvel_)[rr] += (*global_incrtauinvsvel)[rr];
    }
    for (int rr = 0; rr < presize; ++rr)
    {
      (*sumabsres_)[rr] += (*global_incrabsres)[rr];
      (*sumabssvelaf_)[rr] += (*global_incrabssvelaf)[rr];

      (*sumhk_)[rr] += (*global_incrhk)[rr];
      (*sumhbazilevs_)[rr] += (*global_incrhbazilevs)[rr];
      (*sumstrle_)[rr] += (*global_incrstrle)[rr];
      (*sumgradle_)[rr] += (*global_incrgradle)[rr];

      (*sumtau_m_)[rr] += (*global_incrtauM)[rr];
      (*sumtau_c_)[rr] += (*global_incrtauC)[rr];

      (*summk_)[rr] += (*global_incrmk)[rr];

      (*sumres_c_)[rr] += (*global_incrresC)[rr];
      (*sumres_c_sq_)[rr] += (*global_incrresC_sq)[rr];
      (*sumspressnp_)[rr] += (*global_incrspressnp)[rr];
      (*sumspressnp_sq_)[rr] += (*global_incrspressnp_sq)[rr];

      (*sum_eps_pspg_)[rr] += (*global_incr_eps_pspg)[rr];
      (*sum_eps_supg_)[rr] += (*global_incr_eps_supg)[rr];
      (*sum_eps_cross_)[rr] += (*global_incr_eps_cross)[rr];
      (*sum_eps_rey_)[rr] += (*global_incr_eps_rey)[rr];
      (*sum_eps_graddiv_)[rr] += (*global_incr_eps_graddiv)[rr];
      (*sum_eps_eddyvisc_)[rr] += (*global_incr_eps_eddyvisc)[rr];
      (*sum_eps_visc_)[rr] += (*global_incr_eps_visc)[rr];
      (*sum_eps_conv_)[rr] += (*global_incr_eps_conv)[rr];
      (*sum_eps_avm3_)[rr] += (*global_incr_eps_avm3)[rr];
      (*sum_eps_mfs_)[rr] += (*global_incr_eps_mfs)[rr];
      (*sum_eps_mfscross_)[rr] += (*global_incr_eps_mfscross)[rr];
      (*sum_eps_mfsrey_)[rr] += (*global_incr_eps_mfsrey)[rr];
    }

    for (int rr = 0; rr < stresssize; ++rr)
    {
      (*sum_crossstress_)[rr] += (*global_incrcrossstress)[rr];
      (*sum_reystress_)[rr] += (*global_incrreystress)[rr];
    }

    // reset working arrays
    local_vol = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    local_incrhk = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrhbazilevs = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrstrle = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrgradle = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    local_incrtauC = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrtauM = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    local_incrmk = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    local_incrres = Teuchos::RCP(new std::vector<double>(velsize, 0.0));
    local_incrres_sq = Teuchos::RCP(new std::vector<double>(velsize, 0.0));
    local_incrsvelaf = Teuchos::RCP(new std::vector<double>(velsize, 0.0));
    local_incrsvelaf_sq = Teuchos::RCP(new std::vector<double>(velsize, 0.0));
    local_incrtauinvsvel = Teuchos::RCP(new std::vector<double>(velsize, 0.0));

    local_incrabsres = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrabssvelaf = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    local_incrresC = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrresC_sq = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrspressnp = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incrspressnp_sq = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    local_incr_eps_pspg = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_supg = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_cross = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_rey = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_graddiv = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_eddyvisc = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_visc = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_conv = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_avm3 = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_mfs = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_mfscross = Teuchos::RCP(new std::vector<double>(presize, 0.0));
    local_incr_eps_mfsrey = Teuchos::RCP(new std::vector<double>(presize, 0.0));

    local_incrcrossstress = Teuchos::RCP(new std::vector<double>(stresssize, 0.0));
    local_incrreystress = Teuchos::RCP(new std::vector<double>(stresssize, 0.0));

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrvol", local_vol);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhk", local_incrhk);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrhbazilevs", local_incrhbazilevs);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrstrle", local_incrstrle);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrgradle", local_incrgradle);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauC", local_incrtauC);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauM", local_incrtauM);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrmk", local_incrmk);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres", local_incrres);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrres_sq", local_incrres_sq);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabsres", local_incrabsres);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauinvsvel", local_incrtauinvsvel);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf", local_incrsvelaf);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrsvelaf_sq", local_incrsvelaf_sq);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrabssvelaf", local_incrabssvelaf);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC", local_incrresC);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresC_sq", local_incrresC_sq);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp", local_incrspressnp);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrspressnp_sq", local_incrspressnp_sq);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_pspg", local_incr_eps_pspg);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_supg", local_incr_eps_supg);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_cross", local_incr_eps_cross);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_rey", local_incr_eps_rey);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_graddiv", local_incr_eps_graddiv);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_eddyvisc", local_incr_eps_eddyvisc);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_visc", local_incr_eps_visc);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_conv", local_incr_eps_conv);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_avm3", local_incr_eps_avm3);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfs", local_incr_eps_mfs);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfscross", local_incr_eps_mfscross);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incr_eps_mfsrey", local_incr_eps_mfsrey);

    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrcrossstress", local_incrcrossstress);
    eleparams_.set<Teuchos::RCP<std::vector<double>>>("incrreystress", local_incrreystress);

    if (scatradiscret_ != Teuchos::null)
    {
      // ------------------------------------------------
      // get results from element call via parameter list
      Teuchos::RCP<std::vector<double>> local_scatra_vol =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incrvol");

      Teuchos::RCP<std::vector<double>> local_scatra_incrtauS =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incrtauS");
      Teuchos::RCP<std::vector<double>> local_scatra_incrresS =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incrresS");
      Teuchos::RCP<std::vector<double>> local_scatra_incrresS_sq =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incrresS_sq");

      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_visc =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_visc");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_conv =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_conv");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_eddyvisc =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_eddyvisc");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_avm3 =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_avm3");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_mfs =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_mfs");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_mfscross =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_mfscross");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_mfsrey =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_mfsrey");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_supg =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_supg");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_cross =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_cross");
      Teuchos::RCP<std::vector<double>> local_scatra_incr_eps_rey =
          scatraeleparams_.get<Teuchos::RCP<std::vector<double>>>("incr_scatra_eps_rey");

      int phisize = local_scatra_incrresS->size();

      //--------------------------------------------------
      // vectors to sum over all procs

      // volume of element layers
      Teuchos::RCP<std::vector<double>> global_scatra_vol;
      global_scatra_vol = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      // (in plane) averaged values of tauM/tauC
      Teuchos::RCP<std::vector<double>> global_scatra_incrtauS;
      global_scatra_incrtauS = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      // (in plane) averaged values of resS (^2)
      Teuchos::RCP<std::vector<double>> global_scatra_incrresS;
      global_scatra_incrresS = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      Teuchos::RCP<std::vector<double>> global_scatra_incrresS_sq;
      global_scatra_incrresS_sq = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      // (in plane) averaged values of dissipation by supg term
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_supg;
      global_scatra_incr_eps_supg = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      // (in plane) averaged values of dissipation by cross term
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_cross;
      global_scatra_incr_eps_cross = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      // (in plane) averaged values of dissipation by reynolds term
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_rey;
      global_scatra_incr_eps_rey = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      // (in plane) averaged values of dissipation by eddy viscosity
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_eddyvisc;
      global_scatra_incr_eps_eddyvisc = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_avm3;
      global_scatra_incr_eps_avm3 = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      // (in plane) averaged values of dissipation by mfs model
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_mfs;
      global_scatra_incr_eps_mfs = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_mfscross;
      global_scatra_incr_eps_mfscross = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_mfsrey;
      global_scatra_incr_eps_mfsrey = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      // (in plane) averaged values of dissipation by viscous forces
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_visc;
      global_scatra_incr_eps_visc = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      // (in plane) averaged values of dissipation/production by convection
      Teuchos::RCP<std::vector<double>> global_scatra_incr_eps_conv;
      global_scatra_incr_eps_conv = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      //--------------------------------------------------
      // global sums

      // compute global sum, volume
      discret_->get_comm().SumAll(local_scatra_vol->data(), global_scatra_vol->data(), phisize);

      // compute global sums, stabilisation parameters
      discret_->get_comm().SumAll(
          local_scatra_incrtauS->data(), global_scatra_incrtauS->data(), phisize);

      // compute global sums, incompressibility residuals
      discret_->get_comm().SumAll(
          local_scatra_incrresS->data(), global_scatra_incrresS->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incrresS_sq->data(), global_scatra_incrresS_sq->data(), phisize);

      // compute global sums, disspiation rates

      discret_->get_comm().SumAll(
          local_scatra_incr_eps_supg->data(), global_scatra_incr_eps_supg->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_cross->data(), global_scatra_incr_eps_cross->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_rey->data(), global_scatra_incr_eps_rey->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_eddyvisc->data(), global_scatra_incr_eps_eddyvisc->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_visc->data(), global_scatra_incr_eps_visc->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_conv->data(), global_scatra_incr_eps_conv->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_avm3->data(), global_scatra_incr_eps_avm3->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_mfs->data(), global_scatra_incr_eps_mfs->data(), phisize);
      discret_->get_comm().SumAll(
          local_scatra_incr_eps_mfsrey->data(), global_scatra_incr_eps_mfsrey->data(), phisize);

      for (int rr = 0; rr < presize; ++rr)
      {
        (*sumtau_s_)[rr] += (*global_scatra_incrtauS)[rr];

        (*sumres_s_)[rr] += (*global_scatra_incrresS)[rr];
        (*sumres_s_sq_)[rr] += (*global_scatra_incrresS_sq)[rr];

        (*sum_scatra_eps_supg_)[rr] += (*global_scatra_incr_eps_supg)[rr];
        (*sum_scatra_eps_cross_)[rr] += (*global_scatra_incr_eps_cross)[rr];
        (*sum_scatra_eps_rey_)[rr] += (*global_scatra_incr_eps_rey)[rr];
        (*sum_scatra_eps_eddyvisc_)[rr] += (*global_scatra_incr_eps_eddyvisc)[rr];
        (*sum_scatra_eps_visc_)[rr] += (*global_scatra_incr_eps_visc)[rr];
        (*sum_scatra_eps_conv_)[rr] += (*global_scatra_incr_eps_conv)[rr];
        (*sum_scatra_eps_avm3_)[rr] += (*global_scatra_incr_eps_avm3)[rr];
        (*sum_scatra_eps_mfs_)[rr] += (*global_scatra_incr_eps_mfs)[rr];
        (*sum_scatra_eps_mfscross_)[rr] += (*global_scatra_incr_eps_mfscross)[rr];
        (*sum_scatra_eps_mfsrey_)[rr] += (*global_scatra_incr_eps_mfsrey)[rr];
      }


      // reset working arrays
      local_scatra_vol = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      local_scatra_incrtauS = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incrresS = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incrresS_sq = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      local_scatra_incr_eps_supg = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_cross = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_rey = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_eddyvisc = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_visc = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_conv = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_avm3 = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_mfs = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_mfscross = Teuchos::RCP(new std::vector<double>(phisize, 0.0));
      local_scatra_incr_eps_mfsrey = Teuchos::RCP(new std::vector<double>(phisize, 0.0));

      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>("incrvol", local_scatra_vol);

      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>("incrtauS", local_scatra_incrtauS);

      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>("incrresS", local_scatra_incrresS);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incrresS_sq", local_scatra_incrresS_sq);

      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_supg", local_scatra_incr_eps_supg);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_cross", local_scatra_incr_eps_cross);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_rey", local_scatra_incr_eps_rey);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_eddyvisc", local_scatra_incr_eps_eddyvisc);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_visc", local_scatra_incr_eps_visc);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_conv", local_scatra_incr_eps_conv);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_avm3", local_scatra_incr_eps_avm3);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_mfs", local_scatra_incr_eps_mfs);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_mfscross", local_scatra_incr_eps_mfscross);
      scatraeleparams_.set<Teuchos::RCP<std::vector<double>>>(
          "incr_scatra_eps_mfsrey", local_scatra_incr_eps_mfsrey);
    }  // end if scatra
  }    // end if dissipation

  return;
}  // FLD::TurbulenceStatisticsCha::EvaluateResiduals



/*----------------------------------------------------------------------*

       Compute a time average of the mean values over all steps
          since the last output. Dump the result to file.

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::time_average_means_and_output_of_statistics(const int step)
{
  if (numsamp_ == 0)
  {
    FOUR_C_THROW("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average
  int aux = numele_ * numsamp_;
  if (aux < 1) FOUR_C_THROW("Prevent division by zero.");

  for (unsigned i = 0; i < planecoordinates_->size(); ++i)
  {
    (*sumu_)[i] /= aux;
    (*sumv_)[i] /= aux;
    (*sumw_)[i] /= aux;
    (*sump_)[i] /= aux;

    (*sumuv_)[i] /= aux;
    (*sumuw_)[i] /= aux;
    (*sumvw_)[i] /= aux;

    (*sumsqu_)[i] /= aux;
    (*sumsqv_)[i] /= aux;
    (*sumsqw_)[i] /= aux;
    (*sumsqp_)[i] /= aux;
  }

  Core::FE::Nurbs::NurbsDiscretization* nurbsdis =
      dynamic_cast<Core::FE::Nurbs::NurbsDiscretization*>(&(*discret_));

  if (nurbsdis == nullptr)
  {
    for (unsigned i = 0; i < planecoordinates_->size(); ++i)
    {
      // the pointwise values have already been normalised by
      // "countnodesinplaneonallprocs", so we just divide by
      // the number of time samples
      (*pointsumu_)[i] /= numsamp_;
      (*pointsumv_)[i] /= numsamp_;
      (*pointsumw_)[i] /= numsamp_;
      (*pointsump_)[i] /= numsamp_;

      (*pointsumsqu_)[i] /= numsamp_;
      (*pointsumsqv_)[i] /= numsamp_;
      (*pointsumsqw_)[i] /= numsamp_;
      (*pointsumsqp_)[i] /= numsamp_;
    }
  }

  sumforceu_ /= numsamp_;
  sumforcev_ /= numsamp_;
  sumforcew_ /= numsamp_;

  //----------------------------------------------------------------------
  // evaluate area to calculate u_tau, l_tau (and tau_W)
  double area = 1.0;
  for (int i = 0; i < 3; i++)
  {
    if (i != dim_)
    {
      area *= ((*boundingbox_)(1, i) - (*boundingbox_)(0, i));
    }
  }
  // there are two Dirichlet boundaries
  area *= 2.0;

  //----------------------------------------------------------------------
  // we expect nonzero forces (tractions) only in flow direction

  // ltau is used to compute y+
  double ltau = 0.0;
  if (sumforceu_ > sumforcev_ && sumforceu_ > sumforcew_)
  {
    if (abs(sumforceu_) < 1.0e-12)
    {
      FOUR_C_THROW("zero force during computation of wall shear stress\n");
    }

    ltau = visc_ / std::sqrt(sumforceu_ / dens_ / area);
  }
  else if (sumforcev_ > sumforceu_ && sumforcev_ > sumforcew_)
  {
    ltau = visc_ / std::sqrt(sumforcev_ / dens_ / area);
  }
  else if (sumforcew_ > sumforceu_ && sumforcew_ > sumforcev_)
  {
    ltau = visc_ / std::sqrt(sumforcew_ / dens_ / area);
  }
  else
  {
    FOUR_C_THROW("Cannot determine flow direction by traction (seems to be not unique)");
  }
  if (abs(ltau) < 1.0E-14) FOUR_C_THROW("ltau is zero!");

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->get_comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    if (inflowchannel_)
      s.append(".inflow.flow_statistics");
    else
      s.append(".flow_statistics");

    log = Teuchos::RCP(new std::ofstream(s.c_str(), std::ios::app));
    (*log) << "\n\n\n";
    (*log) << "# Statistics record " << countrecord_;
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";

    (*log) << "# (u_tau)^2 = tau_W/rho : ";
    (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforceu_ / dens_ / area;
    (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforcev_ / dens_ / area;
    (*log) << "   " << std::setw(11) << std::setprecision(4) << sumforcew_ / dens_ / area;
    (*log) << &std::endl;


    (*log) << "#|-------------------";
    (*log) << "----------------------------------------------------------";
    (*log) << "--integration based-------------------------";
    (*log) << "----------------------------------------------------------|";
    (*log) << "-------------------------------------------------point";
    (*log) << "wise---------------------------------------";
    (*log) << "------------|\n";

    (*log) << "#     y            y+";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "        mean u^2      mean v^2      mean w^2";
    (*log) << "      mean u*v      mean u*w      mean v*w      mean p^2";
    (*log) << "       umean         vmean         wmean         pmean";
    (*log) << "        mean u^2      mean v^2      mean w^2";
    (*log) << "       mean p^2 \n";
    (*log) << std::scientific;
    for (unsigned i = 0; i < planecoordinates_->size(); ++i)
    {
      // y and y+
      (*log) << " " << std::setw(11) << std::setprecision(4) << (*planecoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*planecoordinates_)[i] / ltau;

      // integral mean values
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumu_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumv_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumw_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sump_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqu_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqv_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqw_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumuv_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumuw_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumvw_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqp_)[i];

      // pointwise means
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsumu_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsumv_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsumw_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsump_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsumsqu_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsumsqv_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsumsqw_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*pointsumsqp_)[i];
      (*log) << "   \n";
    }
    log->flush();

    // ------------------------------------------------------------------
    // additional output for dynamic Smagorinsky model
    if (smagorinsky_)
    {
      // get the outfile
      Teuchos::RCP<std::ofstream> log_Cs;

      std::string s_smag(statistics_outfilename_);
      s_smag.append(".Cs_statistics");

      log_Cs = Teuchos::RCP(new std::ofstream(s_smag.c_str(), std::ios::app));

      (*log_Cs) << "\n\n\n";
      (*log_Cs) << "# Statistics record " << countrecord_;
      (*log_Cs) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";


      (*log_Cs) << "#     y      ";
      (*log_Cs) << "     Cs      ";
      (*log_Cs) << "   (Cs*hk)^2 ";
      (*log_Cs) << "    visceff  ";
      (*log_Cs) << "    Prt      ";
      (*log_Cs) << "(Cs*hk)^2/Prt";
      (*log_Cs) << "    diffeff  ";
      (*log_Cs) << "     Ci      ";
      (*log_Cs) << "   (Ci*hk)^2 ";
      (*log_Cs) << &std::endl;
      (*log_Cs) << std::scientific;
      for (unsigned rr = 0; rr < sum_cs_->size(); ++rr)
      {
        // we associate the value with the midpoint of the element layer
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        // the five values to be visualized
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumvisceff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumdiffeff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_delta_sq_)[rr]) / (numele_ * numsamp_) << &std::endl;
      }
      log_Cs->flush();
    }  // end smagorinsky_

    // ------------------------------------------------------------------
    // additional output for multifractal subgrid-scale modeling
    if (multifractal_)
    {
      // get the outfile
      Teuchos::RCP<std::ofstream> log_mf;

      std::string s_mf(statistics_outfilename_);
      s_mf.append(".MF_statistics");

      log_mf = Teuchos::RCP(new std::ofstream(s_mf.c_str(), std::ios::app));

      (*log_mf) << "\n\n\n";
      (*log_mf) << "# Statistics record " << countrecord_;
      (*log_mf) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";


      (*log_mf) << "#     y      ";
      (*log_mf) << "  N_stream   ";
      (*log_mf) << "  N_normal   ";
      (*log_mf) << "  N_span     ", (*log_mf) << "  B_stream   ";
      (*log_mf) << "  B_normal   ";
      (*log_mf) << "  B_span     ";
      (*log_mf) << "    Csgs     ";
      (*log_mf) << "    sgvisc   ";
      (*log_mf) << &std::endl;
      (*log_mf) << std::scientific;
      for (unsigned rr = 0; rr < sum_n_stream_->size(); ++rr)
      {
        // we associate the value with the midpoint of the element layer
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        // the three values to be visualised
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sum_n_stream_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sum_n_normal_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sum_n_span_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sum_b_stream_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sum_b_normal_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sum_b_span_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sum_csgs_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << ((*sumsgvisc_)[rr]) / (numele_ * numsamp_) << &std::endl;
      }
      log_mf->flush();
    }  // end multifractal_

    if (subgrid_dissipation_)
    {
      Teuchos::RCP<std::ofstream> log_res;

      // output of residuals and subscale quantities
      std::string s_res(statistics_outfilename_);
      s_res.append(".res_statistics");

      log_res = Teuchos::RCP(new std::ofstream(s_res.c_str(), std::ios::app));

      (*log_res) << "\n\n\n";
      (*log_res) << "# Statistics record " << countrecord_;
      (*log_res) << " ( Steps " << step - numsamp_ + 1 << " -- " << step << " )   ";
      (*log_res) << " (dt " << params_.get<double>("time step size") << ")\n";

      (*log_res) << "#       y    ";
      (*log_res) << "    res_x   ";
      (*log_res) << "      res_y  ";
      (*log_res) << "      res_z  ";
      (*log_res) << "     svel_x  ";
      (*log_res) << "     svel_y  ";
      (*log_res) << "     svel_z  ";

      (*log_res) << "   res_sq_x  ";
      (*log_res) << "   res_sq_y  ";
      (*log_res) << "   res_sq_z  ";
      (*log_res) << "   svel_sq_x ";
      (*log_res) << "   svel_sq_y ";
      (*log_res) << "   svel_sq_z ";

      (*log_res) << " tauinvsvel_x";
      (*log_res) << " tauinvsvel_y";
      (*log_res) << " tauinvsvel_z";

      (*log_res) << "    ||res||  ";
      (*log_res) << "   ||svel||  ";

      (*log_res) << "      resC   ";
      (*log_res) << "    spresnp  ";

      (*log_res) << "    resC_sq  ";
      (*log_res) << "  spresnp_sq ";

      (*log_res) << "    tauM     ";
      (*log_res) << "    tauC     ";

      (*log_res) << "  eps_pspg   ";
      (*log_res) << "  eps_supg   ";
      (*log_res) << "  eps_cross  ";
      (*log_res) << "   eps_rey   ";
      (*log_res) << "  eps_graddiv  ";
      (*log_res) << " eps_eddyvisc";
      (*log_res) << "   eps_visc  ";
      (*log_res) << "   eps_conv  ";
      (*log_res) << "   eps_avm3  ";
      (*log_res) << "   eps_mfs   ";
      (*log_res) << " eps_mfscross";
      (*log_res) << " eps_mfsrey  ";

      (*log_res) << "     hk      ";
      (*log_res) << "   strle     ";
      (*log_res) << "   gradle    ";
      (*log_res) << " h_bazilevs  ";
      (*log_res) << "     Dy      ";
      (*log_res) << " tau_cross_11";
      (*log_res) << " tau_cross_22";
      (*log_res) << " tau_cross_33";
      (*log_res) << " tau_cross_12";
      (*log_res) << " tau_cross_23";
      (*log_res) << " tau_cross_31";
      (*log_res) << " tau_rey_11  ";
      (*log_res) << " tau_rey_22  ";
      (*log_res) << " tau_rey_33  ";
      (*log_res) << " tau_rey_12  ";
      (*log_res) << " tau_rey_23  ";
      (*log_res) << " tau_rey_31  ";
      (*log_res) << " mk          ";
      (*log_res) << "\n";

      (*log_res) << std::scientific;
      for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
      {
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabsres_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabssvelaf_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_sq_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_sq_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_m_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_c_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_pspg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_supg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_cross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_rey_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_graddiv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_eddyvisc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_visc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_conv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_avm3_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfscross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfsrey_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*sumhk_)[rr] / (numele_ * numsamp_)
                   << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumstrle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumgradle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumhbazilevs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*nodeplanes_)[rr + 1] - (*nodeplanes_)[rr] << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*summk_)[rr] / (numele_ * numsamp_)
                   << "  ";


        (*log_res) << &std::endl;
      }
      log_res->flush();
    }  // end subgrid_dissipation_
  }    // end myrank 0


  // log was written, so increase counter for records
  countrecord_++;

  return;

}  // TurbulenceStatisticsCha::time_average_means_and_output_of_statistics


/*----------------------------------------------------------------------*

      Compute a time average of the mean values over all steps
       of the sampling period so far. Dump the result to file.

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::dump_statistics(const int step)
{
  if (numsamp_ == 0)
  {
    FOUR_C_THROW("No samples to do time average");
  }

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average
  int aux = numele_ * numsamp_;

  //----------------------------------------------------------------------
  // evaluate area to calculate u_tau, l_tau (and tau_W)
  double area = 1.0;
  for (int i = 0; i < 3; i++)
  {
    if (i != dim_)
    {
      area *= ((*boundingbox_)(1, i) - (*boundingbox_)(0, i));
    }
  }
  // there are two Dirichlet boundaries
  area *= 2;

  //----------------------------------------------------------------------
  // we expect nonzero forces (tractions) only in flow direction

  // ltau is used to compute y+
  double ltau = 0;
  if (sumforceu_ > sumforcev_ && sumforceu_ > sumforcew_)
  {
    ltau = visc_ / std::sqrt(sumforceu_ / dens_ / (area * numsamp_));
  }
  else if (sumforcev_ > sumforceu_ && sumforcev_ > sumforcew_)
  {
    ltau = visc_ / std::sqrt(sumforcev_ / dens_ / (area * numsamp_));
  }
  else if (sumforcew_ > sumforceu_ && sumforcew_ > sumforcev_)
  {
    ltau = visc_ / std::sqrt(sumforcew_ / dens_ / (area * numsamp_));
  }
  else
  {
    FOUR_C_THROW("Cannot determine flow direction by traction (seems to be not unique)");
  }

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->get_comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    if (inflowchannel_)
      s.append(".inflow.flow_statistics");
    else
      s.append(".flow_statistics");

    log = Teuchos::RCP(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent incompressible channel flow (first- and second-order "
              "moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";

    (*log) << "# (u_tau)^2 = tau_W/rho : ";
    (*log) << "   " << std::setw(11) << std::setprecision(4)
           << sumforceu_ / (area * numsamp_) / dens_;
    (*log) << "   " << std::setw(11) << std::setprecision(4)
           << sumforcev_ / (area * numsamp_) / dens_;
    (*log) << "   " << std::setw(11) << std::setprecision(4)
           << sumforcew_ / (area * numsamp_) / dens_;
    (*log) << &std::endl;

    (*log) << "#     y            y+";
    (*log) << "           umean         vmean         wmean         pmean";
    (*log) << "        mean u^2      mean v^2      mean w^2     mean p^2";
    (*log) << "      mean u*v      mean u*w      mean v*w\n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < planecoordinates_->size(); ++i)
    {
      (*log) << " " << std::setw(11) << std::setprecision(4) << (*planecoordinates_)[i];
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*planecoordinates_)[i] / ltau;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumu_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumv_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumw_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sump_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqu_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqv_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqw_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumsqp_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumuv_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumuw_)[i] / aux;
      (*log) << "   " << std::setw(11) << std::setprecision(4) << (*sumvw_)[i] / aux;
      (*log) << "\n";
    }
    log->flush();
  }

  if (discret_->get_comm().MyPID() == 0)
  {
    // ------------------------------------------------------------------
    // additional output for dynamic Smagorinsky model
    if (smagorinsky_)
    {
      // get the outfile
      Teuchos::RCP<std::ofstream> log_Cs;

      std::string s_smag(statistics_outfilename_);
      s_smag.append(".Cs_statistics");

      log_Cs = Teuchos::RCP(new std::ofstream(s_smag.c_str(), std::ios::out));
      (*log_Cs)
          << "# Statistics for turbulent incompressible channel flow (Smagorinsky constant)\n\n";

      (*log_Cs) << "\n\n\n";
      (*log_Cs) << "# Statistics record " << countrecord_;
      (*log_Cs) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";


      (*log_Cs) << "#     y      ";
      (*log_Cs) << "     Cs      ";
      (*log_Cs) << "   (Cs*hk)^2 ";
      (*log_Cs) << "    visceff  ";
      (*log_Cs) << "    Prt      ";
      (*log_Cs) << "(Cs*hk)^2/Prt";
      (*log_Cs) << "    diffeff  ";
      (*log_Cs) << "     Ci      ";
      (*log_Cs) << "   (Ci*hk)^2 ";
      (*log_Cs) << &std::endl;
      (*log_Cs) << std::scientific;
      for (unsigned rr = 0; rr < sum_cs_->size(); ++rr)
      {
        // we associate the value with the midpoint of the element layer
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        // the five values to be visualized
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumvisceff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumdiffeff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_delta_sq_)[rr]) / (numele_ * numsamp_) << &std::endl;
      }
      log_Cs->flush();
    }  // end smagorinsky_

    if (subgrid_dissipation_)
    {
      Teuchos::RCP<std::ofstream> log_res;

      // output of residuals and subscale quantities
      std::string s_res(statistics_outfilename_);
      s_res.append(".res_statistics");

      log_res = Teuchos::RCP(new std::ofstream(s_res.c_str(), std::ios::out));

      (*log_res) << "\n\n\n";
      (*log_res) << "# Statistics record " << countrecord_;
      (*log_res) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")   ";
      (*log_res) << " (dt " << params_.get<double>("time step size") << ")\n";

      (*log_res) << "#       y    ";
      (*log_res) << "    res_x   ";
      (*log_res) << "      res_y  ";
      (*log_res) << "      res_z  ";
      (*log_res) << "     svel_x  ";
      (*log_res) << "     svel_y  ";
      (*log_res) << "     svel_z  ";

      (*log_res) << "   res_sq_x  ";
      (*log_res) << "   res_sq_y  ";
      (*log_res) << "   res_sq_z  ";
      (*log_res) << "   svel_sq_x ";
      (*log_res) << "   svel_sq_y ";
      (*log_res) << "   svel_sq_z ";

      (*log_res) << " tauinvsvel_x";
      (*log_res) << " tauinvsvel_y";
      (*log_res) << " tauinvsvel_z";

      (*log_res) << "    ||res||  ";
      (*log_res) << "   ||svel||  ";

      (*log_res) << "      resC   ";
      (*log_res) << "    spresnp  ";

      (*log_res) << "    resC_sq  ";
      (*log_res) << "  spresnp_sq ";

      (*log_res) << "    tauM     ";
      (*log_res) << "    tauC     ";

      (*log_res) << "  eps_pspg   ";
      (*log_res) << "  eps_supg   ";
      (*log_res) << "  eps_cross  ";
      (*log_res) << "   eps_rey   ";
      (*log_res) << "  eps_graddiv  ";
      (*log_res) << " eps_eddyvisc";
      (*log_res) << "   eps_visc  ";
      (*log_res) << "   eps_conv  ";
      (*log_res) << "   eps_avm3  ";
      (*log_res) << "   eps_mfs   ";
      (*log_res) << " eps_mfscross";
      (*log_res) << " eps_mfsrey  ";

      (*log_res) << "     hk      ";
      (*log_res) << "   strle     ";
      (*log_res) << "   gradle    ";
      (*log_res) << " h_bazilevs  ";
      (*log_res) << "     Dy      ";
      (*log_res) << " tau_cross_11";
      (*log_res) << " tau_cross_22";
      (*log_res) << " tau_cross_33";
      (*log_res) << " tau_cross_12";
      (*log_res) << " tau_cross_23";
      (*log_res) << " tau_cross_31";
      (*log_res) << " tau_rey_11  ";
      (*log_res) << " tau_rey_22  ";
      (*log_res) << " tau_rey_33  ";
      (*log_res) << " tau_rey_12  ";
      (*log_res) << " tau_rey_23  ";
      (*log_res) << " tau_rey_31  ";
      (*log_res) << " mk          ";
      (*log_res) << "\n";

      (*log_res) << std::scientific;
      for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
      {
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabsres_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabssvelaf_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_sq_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_sq_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_m_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_c_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_pspg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_supg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_cross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_rey_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_graddiv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_eddyvisc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_visc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_conv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_avm3_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfscross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfsrey_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*sumhk_)[rr] / (numele_ * numsamp_)
                   << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumstrle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumgradle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumhbazilevs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*nodeplanes_)[rr + 1] - (*nodeplanes_)[rr] << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*summk_)[rr] / (numele_ * numsamp_)
                   << "  ";

        (*log_res) << &std::endl;
      }
      log_res->flush();
    }  // end subgrid_dissipation_
  }    // end myrank 0

  return;

}  // TurbulenceStatisticsCha::DumpStatistics


/*----------------------------------------------------------------------*

     Compute a time average of the mean values for low-Mach-number
          flow over all steps of the sampling period so far.
                      Dump the result to file.

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::dump_loma_statistics(const int step)
{
  if (numsamp_ == 0) FOUR_C_THROW("No samples to do time average");

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average
  int aux = numele_ * numsamp_;

  //----------------------------------------------------------------------
  // evaluate area of bottom and top wall to calculate u_tau, l_tau (and tau_W)
  // and assume that area of bottom and top wall are equal
  double area = 1.0;
  for (int i = 0; i < 3; i++)
  {
    if (i != dim_) area *= ((*boundingbox_)(1, i) - (*boundingbox_)(0, i));
  }
  const double areanumsamp = area * numsamp_;

  //----------------------------------------------------------------------
  // we expect nonzero forces (tractions) only in flow direction

  // rho_w and tau_w at bottom and top wall
  const double rhowb = (*sumrho_)[0] / aux;
  const double rhowt = (*sumrho_)[planecoordinates_->size() - 1] / aux;
  double tauwb = 0.0;
  double tauwt = 0.0;
  if (sumforcebu_ > sumforcebv_ && sumforcebu_ > sumforcebw_)
  {
    tauwb = sumforcebu_ / areanumsamp;
    tauwt = sumforcetu_ / areanumsamp;
  }
  else if (sumforcebv_ > sumforcebu_ && sumforcebv_ > sumforcebw_)
  {
    tauwb = sumforcebv_ / areanumsamp;
    tauwt = sumforcetv_ / areanumsamp;
  }
  else if (sumforcebw_ > sumforcebu_ && sumforcebw_ > sumforcebv_)
  {
    tauwb = sumforcebw_ / areanumsamp;
    tauwt = sumforcetw_ / areanumsamp;
  }
  else
    FOUR_C_THROW("Cannot determine flow direction by traction (appears not unique)");

  // heat flux at the wall is trueresidual of energy equation
  // multiplied by the specific heat
  const double qwb = sumqwb_ * shc_ / areanumsamp;
  const double qwt = sumqwt_ * shc_ / areanumsamp;

  // u_tau and l_tau at bottom and top wall as well as mean values
  const double utaub = std::sqrt(tauwb / rhowb);
  const double utaut = std::sqrt(tauwt / rhowt);
  double Ttaub = 0.0;
  if (rhowb * utaub < -2e-9 or rhowb * utaub > 2e-9) Ttaub = qwb / (rhowb * shc_ * utaub);
  double Ttaut = 0.0;
  if (rhowt * utaut < -2e-9 or rhowt * utaut > 2e-9) Ttaut = qwt / (rhowt * shc_ * utaut);

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->get_comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    if (inflowchannel_)
      s.append(".inflow.loma_statistics");
    else
      s.append(".loma_statistics");

    log = Teuchos::RCP(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent variable-density channel flow at low Mach number (first- "
              "and second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";

    (*log) << "# bottom wall: tauwb, rhowb, u_taub, qwb, Ttaub : ";
    (*log) << "   " << std::setw(17) << std::setprecision(10) << tauwb;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << rhowb;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << utaub;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << qwb;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << Ttaub;
    (*log) << &std::endl;

    (*log) << "# top wall:    tauwt, rhowt, u_taut, qwt, Ttaut : ";
    (*log) << "   " << std::setw(17) << std::setprecision(10) << tauwt;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << rhowt;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << utaut;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << qwt;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << Ttaut;
    (*log) << &std::endl;

    (*log) << "#        y";
    (*log) << "                  umean               vmean               wmean               pmean "
              "            rhomean               Tmean             mommean           rhouTmean";
    (*log) << "              mean u^2            mean v^2            mean w^2            mean p^2  "
              "        mean rho^2            mean T^2";
    (*log) << "            mean u*v            mean u*w            mean v*w            mean u*T    "
              "        mean v*T            mean w*T\n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < planecoordinates_->size(); ++i)
    {
      (*log) << " " << std::setw(17) << std::setprecision(10) << (*planecoordinates_)[i];
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumu_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumv_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sump_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumrho_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sum_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumrhou_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumrhou_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqu_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqv_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqp_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqrho_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsq_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumuv_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumuw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumvw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumu_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumv_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumw_t_)[i] / aux;
      (*log) << "\n";
    }
    log->flush();

    if (subgrid_dissipation_)
    {
      Teuchos::RCP<std::ofstream> log_res;

      // output of residuals and subscale quantities
      std::string s_res(statistics_outfilename_);
      s_res.append(".res_statistics");

      log_res = Teuchos::RCP(new std::ofstream(s_res.c_str(), std::ios::out));

      (*log_res) << "# Statistics for turbulent incompressible channel flow (residuals and "
                    "subscale quantities)\n";
      (*log_res) << "# All values are first averaged over the integration points in an element \n";
      (*log_res)
          << "# and after that averaged over a whole element layer in the homogeneous plane\n\n";

      (*log_res) << "\n\n\n";
      (*log_res) << "# Statistics record " << countrecord_;
      (*log_res) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")   ";
      (*log_res) << " (dt " << params_.get<double>("time step size") << ")\n";

      (*log_res) << "#       y    ";
      (*log_res) << "    res_x   ";
      (*log_res) << "      res_y  ";
      (*log_res) << "      res_z  ";
      (*log_res) << "     svel_x  ";
      (*log_res) << "     svel_y  ";
      (*log_res) << "     svel_z  ";

      (*log_res) << "   res_sq_x  ";
      (*log_res) << "   res_sq_y  ";
      (*log_res) << "   res_sq_z  ";
      (*log_res) << "   svel_sq_x ";
      (*log_res) << "   svel_sq_y ";
      (*log_res) << "   svel_sq_z ";

      (*log_res) << " tauinvsvel_x";
      (*log_res) << " tauinvsvel_y";
      (*log_res) << " tauinvsvel_z";

      (*log_res) << "    ||res||  ";
      (*log_res) << "   ||svel||  ";

      (*log_res) << "      resC   ";
      (*log_res) << "    spresnp  ";

      (*log_res) << "    resC_sq  ";
      (*log_res) << "  spresnp_sq ";

      (*log_res) << "    tauM     ";
      (*log_res) << "    tauC     ";

      (*log_res) << "  eps_pspg   ";
      (*log_res) << "  eps_supg   ";
      (*log_res) << "  eps_cross  ";
      (*log_res) << "   eps_rey   ";
      (*log_res) << "  eps_graddiv  ";
      (*log_res) << " eps_eddyvisc";
      (*log_res) << "   eps_visc  ";
      (*log_res) << "   eps_conv  ";
      (*log_res) << "   eps_avm3  ";
      (*log_res) << "   eps_mfs   ";
      (*log_res) << " eps_mfscross";
      (*log_res) << " eps_mfsrey  ";

      (*log_res) << "     hk      ";
      (*log_res) << "   strle     ";
      (*log_res) << "   gradle    ";
      (*log_res) << " h_bazilevs  ";
      (*log_res) << "     Dy      ";
      (*log_res) << " tau_cross_11";
      (*log_res) << " tau_cross_22";
      (*log_res) << " tau_cross_33";
      (*log_res) << " tau_cross_12";
      (*log_res) << " tau_cross_23";
      (*log_res) << " tau_cross_31";
      (*log_res) << " tau_rey_11  ";
      (*log_res) << " tau_rey_22  ";
      (*log_res) << " tau_rey_33  ";
      (*log_res) << " tau_rey_12  ";
      (*log_res) << " tau_rey_23  ";
      (*log_res) << " tau_rey_31  ";
      (*log_res) << " mk          ";
      (*log_res) << "\n";

      (*log_res) << std::scientific;
      for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
      {
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabsres_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabssvelaf_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_sq_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_sq_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_m_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_c_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_pspg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_supg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_cross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_rey_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_graddiv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_eddyvisc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_visc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_conv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_avm3_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfscross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfsrey_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*sumhk_)[rr] / (numele_ * numsamp_)
                   << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumstrle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumgradle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumhbazilevs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*nodeplanes_)[rr + 1] - (*nodeplanes_)[rr] << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*summk_)[rr] / (numele_ * numsamp_)
                   << "  ";

        (*log_res) << &std::endl;
      }
      log_res->flush();


      Teuchos::RCP<std::ofstream> log_res_scatra;

      // output of residuals and subscale quantities
      std::string s_res_scatra(statistics_outfilename_);
      s_res_scatra.append(".res_scatra_statistics");

      log_res_scatra = Teuchos::RCP(new std::ofstream(s_res_scatra.c_str(), std::ios::out));

      (*log_res_scatra) << "# Statistics for turbulent incompressible channel flow with scalar "
                           "transport (residuals and subscale quantities)\n";
      (*log_res_scatra)
          << "# All values are first averaged over the integration points in an element \n";
      (*log_res_scatra)
          << "# and after that averaged over a whole element layer in the homogeneous plane\n\n";
      (*log_res_scatra)
          << "#                           THIS IS THE SCATRA FILE                          \n\n";

      (*log_res_scatra) << "\n\n\n";
      (*log_res_scatra) << "# Statistics record " << countrecord_;
      (*log_res_scatra) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")   ";
      (*log_res_scatra) << " (dt " << params_.get<double>("time step size") << ")\n";

      (*log_res_scatra) << "#       y    ";

      (*log_res_scatra) << "      resS   ";
      (*log_res_scatra) << "    resS_sq  ";
      (*log_res_scatra) << "    tauS     ";

      (*log_res_scatra) << "  eps_supg   ";
      (*log_res_scatra) << "  eps_cross  ";
      (*log_res_scatra) << "   eps_rey   ";
      (*log_res_scatra) << " eps_eddyvisc";
      (*log_res_scatra) << "   eps_visc  ";
      (*log_res_scatra) << "   eps_conv  ";
      (*log_res_scatra) << "   eps_avm3  ";
      (*log_res_scatra) << "   eps_mfs   ";
      (*log_res_scatra) << " eps_mfscross";
      (*log_res_scatra) << " eps_mfsrey  ";

      (*log_res_scatra) << "\n";

      (*log_res_scatra) << std::scientific;
      for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
      {
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sumres_s_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sumres_s_sq_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sumtau_s_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_supg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_cross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_rey_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_eddyvisc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_visc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_conv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_avm3_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_mfs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_mfscross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_mfsrey_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res_scatra) << &std::endl;
      }
      log_res_scatra->flush();


    }  // end subgrid_dissipation_

    // ------------------------------------------------------------------
    // additional output for dynamic Smagorinsky model
    if (smagorinsky_)
    {
      // get the outfile
      Teuchos::RCP<std::ofstream> log_Cs;

      std::string s_smag(statistics_outfilename_);
      s_smag.append(".Cs_statistics");

      log_Cs = Teuchos::RCP(new std::ofstream(s_smag.c_str(), std::ios::out));
      (*log_Cs)
          << "# Statistics for turbulent incompressible channel flow (Smagorinsky constant)\n\n";

      (*log_Cs) << "\n\n\n";
      (*log_Cs) << "# Statistics record " << countrecord_;
      (*log_Cs) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";


      (*log_Cs) << "#     y      ";
      (*log_Cs) << "     Cs      ";
      (*log_Cs) << "   (Cs*hk)^2 ";
      (*log_Cs) << "    visceff  ";
      (*log_Cs) << "    Prt      ";
      (*log_Cs) << "(Cs*hk)^2/Prt";
      (*log_Cs) << "    diffeff  ";
      (*log_Cs) << "     Ci      ";
      (*log_Cs) << "   (Ci*hk)^2 ";
      (*log_Cs) << &std::endl;
      (*log_Cs) << std::scientific;
      for (unsigned rr = 0; rr < sum_cs_->size(); ++rr)
      {
        // we associate the value with the midpoint of the element layer
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        // the five values to be visualized
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumvisceff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumdiffeff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_delta_sq_)[rr]) / (numele_ * numsamp_) << &std::endl;
      }
      log_Cs->flush();
    }  // end smagorinsky_
  }

  return;

}  // TurbulenceStatisticsCha::DumpLomaStatistics


/*----------------------------------------------------------------------*

     Compute a time average of the mean values for turbulent passive
     scalar transport over all steps of the sampling period so far.
                      Dump the result to file.

 -----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::dump_scatra_statistics(const int step)
{
  if (numsamp_ == 0) FOUR_C_THROW("No samples to do time average");

  //----------------------------------------------------------------------
  // the sums are divided by the number of samples to get the time average
  int aux = numele_ * numsamp_;

  //----------------------------------------------------------------------
  // evaluate area of bottom and top wall to calculate u_tau, l_tau (and tau_W)
  // and assume that area of bottom and top wall are equal
  double area = 1.0;
  for (int i = 0; i < 3; i++)
  {
    if (i != dim_) area *= ((*boundingbox_)(1, i) - (*boundingbox_)(0, i));
  }
  const double areanumsamp = area * numsamp_;

  //----------------------------------------------------------------------
  // we expect nonzero forces (tractions) only in flow direction

  // tau_w at bottom and top wall
  double tauwb = 0.0;
  double tauwt = 0.0;
  if (sumforcebu_ > sumforcebv_ && sumforcebu_ > sumforcebw_)
  {
    tauwb = sumforcebu_ / areanumsamp;
    tauwt = sumforcetu_ / areanumsamp;
  }
  else if (sumforcebv_ > sumforcebu_ && sumforcebv_ > sumforcebw_)
  {
    tauwb = sumforcebv_ / areanumsamp;
    tauwt = sumforcetv_ / areanumsamp;
  }
  else if (sumforcebw_ > sumforcebu_ && sumforcebw_ > sumforcebv_)
  {
    tauwb = sumforcebw_ / areanumsamp;
    tauwt = sumforcetw_ / areanumsamp;
  }
  else
    FOUR_C_THROW("Cannot determine flow direction by traction (appears not unique)");

  // flux at the wall is trueresidual of conv-diff equation
  const double qwb = sumqwb_ / areanumsamp;
  const double qwt = sumqwt_ / areanumsamp;

  // u_tau and l_tau at bottom and top wall as well as mean values
  const double utaub = std::sqrt(tauwb / dens_);
  const double utaut = std::sqrt(tauwt / dens_);
  double Ttaub = 0.0;
  if (utaub < -2e-9 or utaub > 2e-9) Ttaub = qwb / utaub;
  double Ttaut = 0.0;
  if (utaut < -2e-9 or utaut > 2e-9) Ttaut = qwt / utaut;

  //----------------------------------------------------------------------
  // output to log-file
  Teuchos::RCP<std::ofstream> log;
  if (discret_->get_comm().MyPID() == 0)
  {
    std::string s(statistics_outfilename_);
    if (inflowchannel_)
      s.append(".inflow.flow_statistics");
    else
      s.append(".flow_statistics");

    log = Teuchos::RCP(new std::ofstream(s.c_str(), std::ios::out));
    (*log) << "# Statistics for turbulent passiv scalar transport in channel (first- and "
              "second-order moments)";
    (*log) << "\n\n\n";
    (*log) << "# Statistics record ";
    (*log) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";

    (*log) << "# bottom wall: tauwb, u_taub, qwb, Ttaub : ";
    (*log) << "   " << std::setw(17) << std::setprecision(10) << tauwb;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << utaub;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << qwb;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << Ttaub;
    (*log) << &std::endl;

    (*log) << "# top wall:    tauwt, u_taut, qwt, Ttaut : ";
    (*log) << "   " << std::setw(17) << std::setprecision(10) << tauwt;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << utaut;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << qwt;
    (*log) << "   " << std::setw(17) << std::setprecision(10) << Ttaut;
    (*log) << &std::endl;

    (*log) << "#        y";
    (*log) << "                  umean               vmean               wmean               pmean "
              "              Tmean";
    (*log) << "              mean u^2            mean v^2            mean w^2            mean p^2  "
              "          mean T^2";
    (*log) << "            mean u*v            mean u*w            mean v*w            mean u*T    "
              "        mean v*T            mean w*T\n";

    (*log) << std::scientific;
    for (unsigned i = 0; i < planecoordinates_->size(); ++i)
    {
      (*log) << " " << std::setw(17) << std::setprecision(10) << (*planecoordinates_)[i];
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumu_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumv_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sump_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sum_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqu_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqv_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsqp_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumsq_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumuv_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumuw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumvw_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumu_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumv_t_)[i] / aux;
      (*log) << "   " << std::setw(17) << std::setprecision(10) << (*sumw_t_)[i] / aux;
      (*log) << "\n";
    }
    log->flush();

    if (subgrid_dissipation_)
    {
      Teuchos::RCP<std::ofstream> log_res;

      // output of residuals and subscale quantities
      std::string s_res(statistics_outfilename_);
      s_res.append(".res_statistics");

      log_res = Teuchos::RCP(new std::ofstream(s_res.c_str(), std::ios::out));

      (*log_res) << "# Statistics for turbulent incompressible channel flow (residuals and "
                    "subscale quantities)\n";
      (*log_res) << "# All values are first averaged over the integration points in an element \n";
      (*log_res)
          << "# and after that averaged over a whole element layer in the homogeneous plane\n\n";

      (*log_res) << "\n\n\n";
      (*log_res) << "# Statistics record " << countrecord_;
      (*log_res) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")   ";
      (*log_res) << " (dt " << params_.get<double>("time step size") << ")\n";

      (*log_res) << "#       y    ";
      (*log_res) << "    res_x   ";
      (*log_res) << "      res_y  ";
      (*log_res) << "      res_z  ";
      (*log_res) << "     svel_x  ";
      (*log_res) << "     svel_y  ";
      (*log_res) << "     svel_z  ";

      (*log_res) << "   res_sq_x  ";
      (*log_res) << "   res_sq_y  ";
      (*log_res) << "   res_sq_z  ";
      (*log_res) << "   svel_sq_x ";
      (*log_res) << "   svel_sq_y ";
      (*log_res) << "   svel_sq_z ";

      (*log_res) << " tauinvsvel_x";
      (*log_res) << " tauinvsvel_y";
      (*log_res) << " tauinvsvel_z";

      (*log_res) << "    ||res||  ";
      (*log_res) << "   ||svel||  ";

      (*log_res) << "      resC   ";
      (*log_res) << "    spresnp  ";

      (*log_res) << "    resC_sq  ";
      (*log_res) << "  spresnp_sq ";

      (*log_res) << "    tauM     ";
      (*log_res) << "    tauC     ";

      (*log_res) << "  eps_pspg   ";
      (*log_res) << "  eps_supg   ";
      (*log_res) << "  eps_cross  ";
      (*log_res) << "   eps_rey   ";
      (*log_res) << "  eps_graddiv  ";
      (*log_res) << " eps_eddyvisc";
      (*log_res) << "   eps_visc  ";
      (*log_res) << "   eps_conv  ";
      (*log_res) << "   eps_avm3  ";
      (*log_res) << "   eps_mfs   ";
      (*log_res) << " eps_mfscross";
      (*log_res) << " eps_mfsrey  ";

      (*log_res) << "     hk      ";
      (*log_res) << "   strle     ";
      (*log_res) << "   gradle    ";
      (*log_res) << " h_bazilevs  ";
      (*log_res) << "     Dy      ";
      (*log_res) << " tau_cross_11";
      (*log_res) << " tau_cross_22";
      (*log_res) << " tau_cross_33";
      (*log_res) << " tau_cross_12";
      (*log_res) << " tau_cross_23";
      (*log_res) << " tau_cross_31";
      (*log_res) << " tau_rey_11  ";
      (*log_res) << " tau_rey_22  ";
      (*log_res) << " tau_rey_33  ";
      (*log_res) << " tau_rey_12  ";
      (*log_res) << " tau_rey_23  ";
      (*log_res) << " tau_rey_31  ";
      (*log_res) << " mk          ";
      (*log_res) << "\n";

      (*log_res) << std::scientific;
      for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
      {
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumsvelaf_sq_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtauinvsvel_)[3 * rr + 2] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabsres_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumabssvelaf_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumres_c_sq_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumspressnp_sq_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_m_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumtau_c_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_pspg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_supg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_cross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_rey_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_graddiv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_eddyvisc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_visc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_conv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_avm3_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfscross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_eps_mfsrey_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*sumhk_)[rr] / (numele_ * numsamp_)
                   << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumstrle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumgradle_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sumhbazilevs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*nodeplanes_)[rr + 1] - (*nodeplanes_)[rr] << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_crossstress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 1] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 2] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 3] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 4] / (numele_ * numsamp_) << "  ";
        (*log_res) << std::setw(11) << std::setprecision(4)
                   << (*sum_reystress_)[6 * rr + 5] / (numele_ * numsamp_) << "  ";

        (*log_res) << std::setw(11) << std::setprecision(4) << (*summk_)[rr] / (numele_ * numsamp_)
                   << "  ";

        (*log_res) << &std::endl;
      }
      log_res->flush();


      Teuchos::RCP<std::ofstream> log_res_scatra;

      // output of residuals and subscale quantities
      std::string s_res_scatra(statistics_outfilename_);
      s_res_scatra.append(".res_scatra_statistics");

      log_res_scatra = Teuchos::RCP(new std::ofstream(s_res_scatra.c_str(), std::ios::out));

      (*log_res_scatra) << "# Statistics for turbulent incompressible channel flow with scalar "
                           "transport (residuals and subscale quantities)\n";
      (*log_res_scatra)
          << "# All values are first averaged over the integration points in an element \n";
      (*log_res_scatra)
          << "# and after that averaged over a whole element layer in the homogeneous plane\n\n";
      (*log_res_scatra)
          << "#                           THIS IS THE SCATRA FILE                          \n\n";

      (*log_res_scatra) << "\n\n\n";
      (*log_res_scatra) << "# Statistics record " << countrecord_;
      (*log_res_scatra) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")   ";
      (*log_res_scatra) << " (dt " << params_.get<double>("time step size") << ")\n";

      (*log_res_scatra) << "#       y    ";

      (*log_res_scatra) << "      resS   ";
      (*log_res_scatra) << "    resS_sq  ";
      (*log_res_scatra) << "    tauS     ";

      (*log_res_scatra) << "  eps_supg   ";
      (*log_res_scatra) << "  eps_cross  ";
      (*log_res_scatra) << "   eps_rey   ";
      (*log_res_scatra) << " eps_eddyvisc";
      (*log_res_scatra) << "   eps_visc  ";
      (*log_res_scatra) << "   eps_conv  ";
      (*log_res_scatra) << "   eps_avm3  ";
      (*log_res_scatra) << "   eps_mfs   ";
      (*log_res_scatra) << " eps_mfscross";
      (*log_res_scatra) << " eps_mfsrey  ";

      (*log_res_scatra) << "\n";

      (*log_res_scatra) << std::scientific;
      for (unsigned rr = 0; rr < nodeplanes_->size() - 1; ++rr)
      {
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sumres_s_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sumres_s_sq_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sumtau_s_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_supg_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_cross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_rey_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_eddyvisc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_visc_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_conv_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_avm3_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_mfs_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_mfscross_)[rr] / (numele_ * numsamp_) << "  ";
        (*log_res_scatra) << std::setw(11) << std::setprecision(4)
                          << (*sum_scatra_eps_mfsrey_)[rr] / (numele_ * numsamp_) << "  ";

        (*log_res_scatra) << &std::endl;
      }
      log_res_scatra->flush();


    }  // end subgrid_dissipation_

    // ------------------------------------------------------------------
    // additional output for multifractal subgrid-scale modeling
    if ((not inflowchannel_) and multifractal_)
    {
      // get the outfile
      Teuchos::RCP<std::ofstream> log_mf;

      std::string s_mf(statistics_outfilename_);
      s_mf.append(".MF_statistics");

      log_mf = Teuchos::RCP(new std::ofstream(s_mf.c_str(), std::ios::out));

      (*log_mf) << "# Statistics for turbulent passiv scalar transport in channel (multifractal "
                   "subgrid-scales parameters)";
      (*log_mf) << "\n\n\n";
      (*log_mf) << "# Statistics record ";
      (*log_mf) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";


      (*log_mf) << "#     y      ";
      (*log_mf) << "  N_stream   ";
      (*log_mf) << "  N_normal   ";
      (*log_mf) << "  N_span     ", (*log_mf) << "  B_stream   ";
      (*log_mf) << "  B_normal   ";
      (*log_mf) << "  B_span     ";
      (*log_mf) << "    Csgs     ";
      (*log_mf) << "    Nphi     ";
      (*log_mf) << "    Dphi     ";
      (*log_mf) << "  Csgs_phi   ";
      (*log_mf) << "    sgvisc   ";
      (*log_mf) << &std::endl;
      (*log_mf) << std::scientific;
      for (unsigned rr = 0; rr < sum_n_stream_->size(); ++rr)
      {
        // we associate the value with the midpoint of the element layer
        (*log_mf) << std::setw(11) << std::setprecision(4)
                  << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        // the values to be visualised
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_n_stream_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_n_normal_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_n_span_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_b_stream_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_b_normal_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_b_span_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_csgs_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_nphi_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_dphi_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sum_csgs_phi_)[rr]) / aux << "  ";
        (*log_mf) << std::setw(11) << std::setprecision(4) << ((*sumsgvisc_)[rr]) / aux
                  << &std::endl;
      }
      log_mf->flush();
    }  // end multifractal_

    // ------------------------------------------------------------------
    // additional output for dynamic Smagorinsky model
    if (smagorinsky_)
    {
      // get the outfile
      Teuchos::RCP<std::ofstream> log_Cs;

      std::string s_smag(statistics_outfilename_);
      s_smag.append(".Cs_statistics");

      log_Cs = Teuchos::RCP(new std::ofstream(s_smag.c_str(), std::ios::out));
      (*log_Cs)
          << "# Statistics for turbulent incompressible channel flow (Smagorinsky constant)\n\n";

      (*log_Cs) << "\n\n\n";
      (*log_Cs) << "# Statistics record " << countrecord_;
      (*log_Cs) << " (Steps " << step - numsamp_ + 1 << "--" << step << ")\n";


      (*log_Cs) << "#     y      ";
      (*log_Cs) << "     Cs      ";
      (*log_Cs) << "   (Cs*hk)^2 ";
      (*log_Cs) << "    visceff  ";
      (*log_Cs) << "    Prt      ";
      (*log_Cs) << "(Cs*hk)^2/Prt";
      (*log_Cs) << "    diffeff  ";
      (*log_Cs) << "     Ci      ";
      (*log_Cs) << "   (Ci*hk)^2 ";
      (*log_Cs) << &std::endl;
      (*log_Cs) << std::scientific;
      for (unsigned rr = 0; rr < sum_cs_->size(); ++rr)
      {
        // we associate the value with the midpoint of the element layer
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << 0.5 * ((*nodeplanes_)[rr + 1] + (*nodeplanes_)[rr]) << "  ";

        // the five values to be visualized
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumvisceff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_cs_delta_sq_prt_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sumdiffeff_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_)[rr]) / (numele_ * numsamp_) << "  ";
        (*log_Cs) << std::setw(11) << std::setprecision(4)
                  << ((*sum_ci_delta_sq_)[rr]) / (numele_ * numsamp_) << &std::endl;
      }
      log_Cs->flush();
    }  // end smagorinsky_
  }

  return;

}  // TurbulenceStatisticsCha::dump_scatra_statistics


/*----------------------------------------------------------------------*

                  Reset sums and number of samples to 0

  ----------------------------------------------------------------------*/
void FLD::TurbulenceStatisticsCha::clear_statistics()
{
  // reset the number of samples
  numsamp_ = 0;

  // reset forces (mean values and values at bottom and top wall)
  sumforceu_ = 0.0;
  sumforcev_ = 0.0;
  sumforcew_ = 0.0;
  sumforcebu_ = 0.0;
  sumforcebv_ = 0.0;
  sumforcebw_ = 0.0;
  sumforcetu_ = 0.0;
  sumforcetv_ = 0.0;
  sumforcetw_ = 0.0;

  // reset integral and pointwise averages
  for (unsigned i = 0; i < planecoordinates_->size(); ++i)
  {
    (*sumu_)[i] = 0.0;
    (*sumv_)[i] = 0.0;
    (*sumw_)[i] = 0.0;
    (*sump_)[i] = 0.0;
    (*sumrho_)[i] = 0.0;
    (*sum_t_)[i] = 0.0;

    (*sumsqu_)[i] = 0.0;
    (*sumsqv_)[i] = 0.0;
    (*sumsqw_)[i] = 0.0;
    (*sumsqp_)[i] = 0.0;
    (*sumsqrho_)[i] = 0.0;
    (*sumsq_t_)[i] = 0.0;

    (*sumuv_)[i] = 0.0;
    (*sumuw_)[i] = 0.0;
    (*sumvw_)[i] = 0.0;
    (*sumu_t_)[i] = 0.0;
    (*sumv_t_)[i] = 0.0;
    (*sumw_t_)[i] = 0.0;

    (*pointsumu_)[i] = 0.0;
    (*pointsumv_)[i] = 0.0;
    (*pointsumw_)[i] = 0.0;
    (*pointsump_)[i] = 0.0;

    (*pointsumsqu_)[i] = 0.0;
    (*pointsumsqv_)[i] = 0.0;
    (*pointsumsqw_)[i] = 0.0;
    (*pointsumsqp_)[i] = 0.0;
  }

  meanvelnp_->PutScalar(0.0);
  if (physicaltype_ == Inpar::FLUID::loma) meanscanp_->PutScalar(0.0);

  // reset sampling for dynamic Smagorinsky model
  if (smagorinsky_)
  {
    for (unsigned rr = 0; rr < sum_cs_->size(); ++rr)
    {
      // reset value
      (*sum_cs_)[rr] = 0.0;
      (*sum_cs_delta_sq_)[rr] = 0.0;
      (*sumvisceff_)[rr] = 0.0;
      (*sum_prt_)[rr] = 0.0;
      (*sum_cs_delta_sq_prt_)[rr] = 0.0;
      (*sumdiffeff_)[rr] = 0.0;
      (*sum_ci_)[rr] = 0.0;
      (*sum_ci_delta_sq_)[rr] = 0.0;
    }
  }  // end smagorinsky_

  // reset sampling for multifractal subgrid scales
  if (multifractal_)
  {
    for (unsigned rr = 0; rr < sum_n_stream_->size(); ++rr)
    {
      // reset value
      (*sum_n_stream_)[rr] = 0.0;
      (*sum_n_normal_)[rr] = 0.0;
      (*sum_n_span_)[rr] = 0.0;
      (*sum_b_stream_)[rr] = 0.0;
      (*sum_b_normal_)[rr] = 0.0;
      (*sum_b_span_)[rr] = 0.0;
      (*sum_csgs_)[rr] = 0.0;
      (*sumsgvisc_)[rr] = 0.0;
      (*sum_nphi_)[rr] = 0.0;
      (*sum_dphi_)[rr] = 0.0;
      (*sum_csgs_phi_)[rr] = 0.0;
    }
  }  // end multifractal_

  // reset residuals and subscale averages
  if (subgrid_dissipation_)
  {
    for (unsigned rr = 0; rr < sumres_->size() / 3; ++rr)
    {
      (*sumres_)[3 * rr] = 0.0;
      (*sumres_)[3 * rr + 1] = 0.0;
      (*sumres_)[3 * rr + 2] = 0.0;

      (*sumsvelaf_)[3 * rr] = 0.0;
      (*sumsvelaf_)[3 * rr + 1] = 0.0;
      (*sumsvelaf_)[3 * rr + 2] = 0.0;

      (*sumres_sq_)[3 * rr] = 0.0;
      (*sumres_sq_)[3 * rr + 1] = 0.0;
      (*sumres_sq_)[3 * rr + 2] = 0.0;

      (*sumsvelaf_sq_)[3 * rr] = 0.0;
      (*sumsvelaf_sq_)[3 * rr + 1] = 0.0;
      (*sumsvelaf_sq_)[3 * rr + 2] = 0.0;

      (*sumtauinvsvel_)[3 * rr] = 0.0;
      (*sumtauinvsvel_)[3 * rr + 1] = 0.0;
      (*sumtauinvsvel_)[3 * rr + 2] = 0.0;

      for (int mm = 0; mm < 6; ++mm)
      {
        (*sum_crossstress_)[6 * rr + mm] = 0.0;
        (*sum_reystress_)[6 * rr + mm] = 0.0;
      }
    }
    for (unsigned rr = 0; rr < sumres_c_->size(); ++rr)
    {
      (*sumabsres_)[rr] = 0.0;
      (*sumabssvelaf_)[rr] = 0.0;

      (*sumhk_)[rr] = 0.0;
      (*sumhbazilevs_)[rr] = 0.0;
      (*sumstrle_)[rr] = 0.0;
      (*sumgradle_)[rr] = 0.0;

      (*sumtau_m_)[rr] = 0.0;
      (*sumtau_c_)[rr] = 0.0;

      (*summk_)[rr] = 0.0;

      (*sum_eps_pspg_)[rr] = 0.0;
      (*sum_eps_supg_)[rr] = 0.0;
      (*sum_eps_cross_)[rr] = 0.0;
      (*sum_eps_rey_)[rr] = 0.0;
      (*sum_eps_graddiv_)[rr] = 0.0;
      (*sum_eps_eddyvisc_)[rr] = 0.0;
      (*sum_eps_visc_)[rr] = 0.0;
      (*sum_eps_conv_)[rr] = 0.0;
      (*sum_eps_mfs_)[rr] = 0.0;
      (*sum_eps_mfscross_)[rr] = 0.0;
      (*sum_eps_mfsrey_)[rr] = 0.0;
      (*sum_eps_avm3_)[rr] = 0.0;

      (*sumres_c_)[rr] = 0.0;
      (*sumspressnp_)[rr] = 0.0;

      (*sumres_c_sq_)[rr] = 0.0;
      (*sumspressnp_sq_)[rr] = 0.0;
    }
    for (unsigned rr = 0; rr < sumres_c_->size(); ++rr)
    {
      (*sumtau_s_)[rr] = 0.0;

      (*sum_scatra_eps_supg_)[rr] = 0.0;
      (*sum_scatra_eps_cross_)[rr] = 0.0;
      (*sum_scatra_eps_rey_)[rr] = 0.0;
      (*sum_scatra_eps_eddyvisc_)[rr] = 0.0;
      (*sum_scatra_eps_visc_)[rr] = 0.0;
      (*sum_scatra_eps_conv_)[rr] = 0.0;
      (*sum_scatra_eps_mfs_)[rr] = 0.0;
      (*sum_scatra_eps_mfscross_)[rr] = 0.0;
      (*sum_scatra_eps_mfsrey_)[rr] = 0.0;
      (*sum_scatra_eps_avm3_)[rr] = 0.0;

      (*sumres_s_)[rr] = 0.0;
      (*sumres_s_sq_)[rr] = 0.0;
    }
  }  // end subgrid_dissipation_

  return;
}  // TurbulenceStatisticsCha::ClearStatistics


void FLD::TurbulenceStatisticsCha::store_scatra_discret_and_params(
    Teuchos::RCP<Core::FE::Discretization> scatradis,
    Teuchos::RCP<Teuchos::ParameterList> scatraparams,
    Teuchos::RCP<Teuchos::ParameterList> scatraextraparams,
    Teuchos::RCP<Teuchos::ParameterList> scatratimeparams)
{
  scatradiscret_ = scatradis;
  scatraparams_ = scatraparams;
  scatraextraparams_ = scatraextraparams;
  scatratimeparams_ = scatratimeparams;

  if (discret_->get_comm().MyPID() == 0)
  {
    std::cout << "Additional information:" << std::endl;
    std::cout << "-> added ScaTra discretization to channel-flow-statistics manager\n" << std::endl;
  }

  if (physicaltype_ == Inpar::FLUID::incompressible)  // not required for loma
  {
    // get diffusivity from material definition --- for computation
    // of additional mfs-statistics
    int id = Global::Problem::instance()->materials()->first_id_by_type(Core::Materials::m_scatra);
    if (id == -1)
      FOUR_C_THROW("Could not find scatra material");
    else
    {
      const Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance()->materials()->parameter_by_id(id);
      const Mat::PAR::ScatraMat* actmat = static_cast<const Mat::PAR::ScatraMat*>(mat);

      double diffus =
          Mat::PAR::ScatraMat(*actmat).get_parameter(actmat->diff, -1);  // actmat->diffusivity_;
      // calculate Schmidt number
      // visc is the kinematic viscosity here
      scnum_ = visc_ / diffus;
      if (dens_ != 1.0) FOUR_C_THROW("Kinematic quantities assumed!");
    }
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
