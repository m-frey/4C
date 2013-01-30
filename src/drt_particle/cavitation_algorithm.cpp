/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.cpp

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "particle_timint_centrdiff.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_calc.H"
#include "../drt_mat/newtonianfluid.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_cut/cut_position.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
CAVITATION::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : PARTICLE::Algorithm(comm,params),
  fluiddis_(Teuchos::null)
{

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the cavitation algorithm                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // transfer particles into their correct bins
    TransferParticles();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished

}


/*----------------------------------------------------------------------*
 | setup of the system                                      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupSystem()
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Init()
{
  particledis_ = DRT::Problem::Instance()->GetDis("particle");

  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(),false));
  fluid_ = fluid->FluidFieldrcp();

  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  CreateBins();

  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcs();

  if(binrowmap->NumGlobalElements() > fluiddis_->NumGlobalElements() / 4.0)
    IO::cout << "\n\n\n WARNING: Reduction of number of bins recommended!! Increase cutoff radius. \n\n\n" << IO::endl;

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>,PARTICLE::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBins(homelessparticles);

  // ghost bins and particles according to the bins
  SetupGhosting(binrowmap);

  // some output
  IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);
  DRT::UTILS::PrintParallelDistribution(*fluiddis_);

  // access structural dynamic params list which will be possibly modified while creating the time integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // create time integrator based on structural time integration
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> particles =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(), const_cast<Teuchos::ParameterList&>(sdyn), particledis_));
  particles_ = particles->StructureFieldrcp();

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrepareTimeStep()
{
  PARTICLE::Algorithm::PrepareTimeStep();
  fluid_->PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Integrate()
{
  CalculateAndApplyFluidForcesToParticles();
  PARTICLE::Algorithm::Integrate();
//  fluid_->NonlinearSolve(); // after discussion with Ursula
  fluid_->MultiCorrector();
  return;
}


/*----------------------------------------------------------------------*
 | calculate fluid forces on particle and apply it         ghamm 01/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateAndApplyFluidForcesToParticles()
{
  fluiddis_->ClearState();
  particledis_->ClearState();

  //TODO 1: Test whether these states are updated enough
  // at the beginning of the time step: veln := velnp(previous step)
  fluiddis_->SetState("veln",fluid_->Veln());
  fluiddis_->SetState("velnm",fluid_->Velnm());

  particledis_->SetState("bubblepos", particles_->ExtractDispn());
  particledis_->SetState("bubblevel", particles_->ExtractVeln());
  particledis_->SetState("bubbleacc", particles_->ExtractAccn());
  // miraculous transformation from row to col layout ...
  Teuchos::RCP<const Epetra_Vector> bubblepos = particledis_->GetState("bubblepos");
  Teuchos::RCP<const Epetra_Vector> bubblevel = particledis_->GetState("bubblevel");
  Teuchos::RCP<const Epetra_Vector> bubbleacc = particledis_->GetState("bubbleacc");

  // bubble radius of layout nodal col map --> SetState not possible
  Teuchos::RCP<Epetra_Vector> bubbleradius = LINALG::CreateVector(*particledis_->NodeColMap(),false);
  LINALG::Export(*Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp(),*bubbleradius);

  // vectors to be filled with forces
  Teuchos::RCP<Epetra_Vector> bubbleforces = LINALG::CreateVector(*particledis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> fluidforces = LINALG::CreateVector(*fluiddis_->DofRowMap(),true);

  if(fluiddis_->NumMyColElements() <= 0)
    dserror("there is no fluid element to ask for material parameters");
  Teuchos::RCP<MAT::NewtonianFluid> actmat = rcp_dynamic_cast<MAT::NewtonianFluid>(fluiddis_->lColElement(0)->Material());
  if(actmat == Teuchos::null)
    dserror("type cast of fluid material failed");

  double fluiddensity = actmat->Density();
  double fluiddynviscosity = actmat->Viscosity();
  double particledensity = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ParticleDensity();

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // all particles (incl ghost) are evaluated
  const int numcolbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numcolbin; ++ibin)
  {
    DRT::Element* actbin = particledis_->lColElement(ibin);
    DRT::Node** particlesinbin = actbin->Nodes();

    // calculate forces for all particles in this bin
    for(int iparticle=0; iparticle<actbin->NumNode();iparticle++)
    {
      const DRT::Node* currparticle = particlesinbin[iparticle];
      // fill particle position
      LINALG::Matrix<3,1> particleposition;
      std::vector<int> lm_b = particledis_->Dof(currparticle);
      double posx = bubblepos->Map().LID(lm_b[0]);
      for (int dim=0; dim<3; dim++)
      {
        particleposition(dim) = (*bubblepos)[posx+dim];
      }


      // 1st step: element coordinates of particle position in fluid element

      // variables to store information about element in which the particle is located
      DRT::Element* targetfluidele = NULL;
      LINALG::Matrix<3,1> elecoord(true);

      // find out in which fluid element the current particle is located
      std::set<int> fluideleIdsinbin = extendedfluidghosting_[actbin->Id()];
      std::set<int>::const_iterator eleiter;
      for(eleiter=fluideleIdsinbin.begin();eleiter!=fluideleIdsinbin.end();++eleiter)
      {
        DRT::Element* fluidele = fluiddis_->gElement(*eleiter);
        DRT::Node** fluidnodes = fluidele->Nodes();
        const int numnode = fluidele->NumNode();

        // fill currentpositions of fluid element
        std::map<int,LINALG::Matrix<3,1> > currentfluidpositions;
        for (int j=0; j < numnode; j++)
        {
          const DRT::Node* node = fluidnodes[j];
          LINALG::Matrix<3,1> currpos;
          for (int dim=0; dim<3; dim++)
          {
            currpos(dim,0) = node->X()[dim];
          }
          currentfluidpositions[node->Id()] = currpos;
        }
        const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(fluidele, currentfluidpositions));

        //get coordinates of the particle position in parameter space of the element
        bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, elecoord);

        if(insideele == true)
        {
          targetfluidele = fluidele;
          // leave loop over all fluid eles
          break;
        }
      }


      // 2nd step: forces on this bubble are calculated
      if(targetfluidele == NULL)
      {
        cout << "INFO: currparticle with Id: " << currparticle->Id() << " and position: " << particleposition(0) << " "
            << particleposition(1) << " " << particleposition(2) << " " << " does not have an underlying fluid element -> no forces calculated" << endl;
        // do not assemble forces for this bubble and continue with next bubble
        continue;
      }

      // get element location vector and ownerships
      std::vector<int> lm_f;
      std::vector<int> lmowner_f;
      std::vector<int> lmstride;
      targetfluidele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

      // Reshape element matrices and vectors and initialize to zero
      elevector1.Size(lm_f.size());
      elevector2.Size(lm_f.size());

      // set action in order to calculate the velocity and material derivative of the velocity
      Teuchos::ParameterList params;
      params.set<int>("action",FLD::calc_mat_derivative_u);
      params.set<double>("timestep",Dt());
      double tmpelecoords[3];
      for (int dim=0; dim<3; dim++)
      {
        tmpelecoords[dim] = elecoord(dim);
      }
      params.set<double*>("elecoords",tmpelecoords);

      // call the element specific evaluate method (elevec1 = fluid vel u; elevec2 = mat deriv of fluid vel, elevec3 = dummy)
      targetfluidele->Evaluate(params,*fluiddis_,lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

//      cout << "fluid vel at bubble position: " << elevector1(0) << "  " << elevector1(1) << "  " << elevector1(2) << "  " << endl;
//      cout << "mat deriv of fluid at bubble position: " << elevector2(0) << "  " << elevector2(1) << "  " << elevector2(2) << "  " << endl;

      // get bubble velocity and acceleration
      std::vector<double> v_bub(lm_b.size());
      DRT::UTILS::ExtractMyValues(*bubblevel,v_bub,lm_b);
      std::vector<double> acc_bub(lm_b.size());
      DRT::UTILS::ExtractMyValues(*bubbleacc,acc_bub,lm_b);

      // get bubble radius
      double r_bub = (*bubbleradius)[ particledis_->NodeColMap()->LID(currparticle->Id()) ];

      // bubble Reynolds number
      LINALG::Matrix<3,1> v_rel(false);
      for (int dim=0; dim<3; dim++)
      {
        v_rel(dim) = elevector1[dim] - v_bub[dim];
      }
      double v_relabs = v_rel.Norm2();
      double Re_b = 2.0 * r_bub * v_relabs * fluiddensity / fluiddynviscosity;

//      cout << "v_rel: " << v_rel(0) << "  " << v_rel(1) << "  " << v_rel(2) << "  " << endl;
//      cout << "v_relabs: " << v_relabs << endl;
//      cout << "bubble Reynolds number: " << Re_b << endl;

      // variable to sum forces for the current bubble under observation
      Epetra_SerialDenseVector forcecurrbubble(3);
      /*------------------------------------------------------------------*/
      //// 2.1) drag force = 0.5 * c_d * rho_l * Pi * r_b^2 * |u-v| * (u-v) or
      //// Stokes law for very small Re: drag force = 6.0 * Pi * mu_l * r_b * (u-v)
      double coeff1 = 0.0;
      if(Re_b < 0.1)
      {
        coeff1 = 6.0 * M_PI * fluiddynviscosity * r_bub;
      }
      else
      {
        double c_d = 0.0;
        if(Re_b < 1000.0)
          c_d = 24.0 * (1.0 + 0.15 * pow(Re_b,0.687)) / Re_b;
        else
          c_d = 0.44;

        coeff1 = 0.5 * c_d * fluiddensity * M_PI * r_bub * r_bub * v_relabs;
      }

      LINALG::Matrix<3,1> dragforce(true);
      dragforce.Update(coeff1, v_rel);
      //assemble
      for(int dim=0; dim<3; dim++)
        forcecurrbubble[dim] += dragforce(dim);
      /*------------------------------------------------------------------*/

//      /*------------------------------------------------------------------*/
//      //// 2.2) virtual/added mass = c_VM * rho_l * volume_b * ( Du/Dt - Dv/Dt )
//      // OUTDATED: see 2.4)
//      double c_VM = 0.5;
//      double vol_bub = 4.0 / 3.0 * M_PI * r_bub * r_bub* r_bub;
//      double coeff2 = c_VM * fluiddensity * vol_bub;
//      // material derivative of fluid velocity
//      LINALG::Matrix<3,1> Du_Dt(false);
//      for (int dim=0; dim<3; dim++)
//      {
//        Du_Dt(dim) = elevector2[dim];
//      }
//      cout << "Du_Dt: " << Du_Dt(0) << "  " << Du_Dt(1) << "  " << Du_Dt(2) << "  " << endl;
//      // material derivative of bubble velocity
//      LINALG::Matrix<3,1> Dv_Dt(false);
//      for (int dim=0; dim<3; dim++)
//      {
//        Dv_Dt(dim) = acc_bub[dim];
//      }
//      LINALG::Matrix<3,1> addedmassforce(true);
//        addedmassforce.Update(coeff2, Du_Dt, -coeff2, Dv_Dt);
//      //assemble
//      for(int dim=0; dim<3; dim++)
//        forcecurrbubble[dim] += addedmassforce(dim);
//      /*------------------------------------------------------------------*/

      /*------------------------------------------------------------------*/
      //// 2.3) buoyancy and gravity forces = volume_b * ( rho_bub - rho_l ) * g
      double vol_bub = 4.0 / 3.0 * M_PI * r_bub * r_bub* r_bub;
      double coeff3 = vol_bub * ( particledensity - fluiddensity);
      LINALG::Matrix<3,1> gravityforce(true);
      gravityforce.Update(coeff3, gravity_acc_);
      //assemble
      for(int dim=0; dim<3; dim++)
        forcecurrbubble[dim] += gravityforce(dim);
      /*------------------------------------------------------------------*/

      /*------------------------------------------------------------------*/
      //// 2.4) implicit treatment of bubble acceleration in added mass, other forces explicit
      //// final force = \frac{ sum all forces (2.1, ..) + c_VM * rho_l * volume_b * Du/Dt }{ 1 + c_VM * rho_l / rho_b }
      double c_VM = 0.5;
      double coeff4 = c_VM * fluiddensity * vol_bub;
      double coeff5 = 1.0 + c_VM * fluiddensity / particledensity;
      // material derivative of fluid velocity
      LINALG::Matrix<3,1> Du_Dt(false);
      LINALG::Matrix<3,1> sumcurrforces(false);
      for (int dim=0; dim<3; dim++)
      {
        Du_Dt(dim) = elevector2[dim];
        sumcurrforces(dim) = forcecurrbubble[dim];
      }

      LINALG::Matrix<3,1> finalforce(true);
      finalforce.Update(1.0/coeff5, sumcurrforces, coeff4/coeff5, Du_Dt);
      //assemble
      for(int dim=0; dim<3; dim++)
        forcecurrbubble[dim] = finalforce(dim);
      /*------------------------------------------------------------------*/


      // 3rd step: assemble bubble forces

      // assemble of bubble forces can be done without further need of LINALG::Export (ghost bubbles also evaluated)
      std::vector<int> lmowner_b(lm_b.size(), myrank_);
      LINALG::Assemble(*bubbleforces,forcecurrbubble,lm_b,lmowner_b);

      // assemble of fluid forces can only be done on row elements (all bubbles in this element must be considered: ghost near boundary)
      const int numnode = targetfluidele->NumNode();
      Epetra_SerialDenseVector funct(numnode);
      // get shape functions of the element; evaluated at the bubble position --> distribution
      DRT::UTILS::shape_function_3D(funct,elecoord(0,0),elecoord(1,0),elecoord(2,0),targetfluidele->Shape());
      // prepare assembly for fluid forces (pressure degrees do not have to be filled); actio = reactio --> minus sign
      int dim = 3;
      Epetra_SerialDenseVector val(numnode*(dim+1));
      for(int iter=0; iter<numnode; iter++)
      {
        for(int k=0; k<dim; k++)
        {
          val[iter*(dim+1) + k] = -1.0 * funct[iter] * forcecurrbubble(k);
        }
      }
      // do assembly of bubble forces on fluid
      LINALG::Assemble(*fluidforces,val,lm_f,lmowner_f);
    } // end iparticle

  } // end ibin


  // 4th step: apply forces to bubbles and fluid field
  particledis_->SetState("bubbleforces", bubbleforces);

  return;
}
/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Update()
{
  PARTICLE::Algorithm::Update();
  fluid_->Update();
  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ReadRestart(int restart)
{
  PARTICLE::Algorithm::ReadRestart(restart);
  fluid_->ReadRestart(restart);
  return;
}


/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CreateBins()
{
  // if not yet specified, get XAABB_ from underlying discretization
  if( XAABB_(2,1) > 0.9e12  and  XAABB_(2,1) < 1.1e12 )
  {
    IO::cout << "XAABB is computed based on the underlying fluid discretization" << IO::endl;
    XAABB_ = GEO::getXAABBofDis(*fluiddis_);
    // local bounding box
    double locmin[3] = {XAABB_(0,0), XAABB_(1,0), XAABB_(2,0)};
    double locmax[3] = {XAABB_(0,1), XAABB_(1,1), XAABB_(2,1)};
    // global bounding box
    double globmin[3];
    double globmax[3];
    // do the necessary communication
    Comm().MinAll(&locmin[0], &globmin[0], 3);
    Comm().MaxAll(&locmax[0], &globmax[0], 3);

    for(int dim=0; dim<3; dim++)
    {
      XAABB_(dim,0) = globmin[dim];
      XAABB_(dim,1) = globmax[dim];
    }
  }

  // divide global bounding box into bins
  for (int dim = 0; dim < 3; dim++)
  {
    // std::floor leads to bins that are at least of size cutoff_radius
    bin_per_dir_[dim] = (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_);
    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
  }

  IO::cout << "Global bounding box size: " << XAABB_;
  IO::cout << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << IO::endl;

  return;
}


/*----------------------------------------------------------------------*
| bins are distributed to the processors                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CAVITATION::Algorithm::DistributeBinsToProcs()
{
  std::vector<int> rowbins;

  // NOTE: This part of the setup can be the bottleneck because vectors of all bins
  // are needed on each proc (memory issue!!); std::map could perhaps help when gathering
  // num fluid nodes in each bin, then block wise communication after copying data to vector

  int numbins = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];

  // loop over all fluid nodes and determine their number in each bin (init with 0)
  std::vector<int> mynumnodes_per_bin(numbins,0);
  for(int inode=0; inode<fluiddis_->NumMyRowNodes(); inode++)
  {
    int ijk[3] = {0,0,0};
    const double* currpos = fluiddis_->lRowNode(inode)->X();
    for(int dim=0; dim < 3; dim++)
    {
      ijk[dim] = (int)((currpos[dim]-XAABB_(dim,0)) / bin_size_[dim]);
    }

    int binId = ConvertijkToGid(&ijk[0]);
    ++mynumnodes_per_bin[binId];
  }

  // find maximum number of nodes in each bin over all procs (init with -1)
  std::vector<int> maxnumnodes_per_bin(numbins,-1);
  fluiddis_->Comm().MaxAll(&mynumnodes_per_bin[0], &maxnumnodes_per_bin[0], numbins);

  // it is possible that several procs have the same number of nodes in a bin
  // only proc which has maximum number of nodes in a bin writes its rank
  std::vector<int> myrank_per_bin(numbins,-1);
  for(int i=0; i<numbins; i++)
  {
    if(mynumnodes_per_bin[i] == maxnumnodes_per_bin[i])
      myrank_per_bin[i] = myrank_;
  }

  mynumnodes_per_bin.clear();
  maxnumnodes_per_bin.clear();

  // find maximum myrank for each bin over all procs (init with -1)
  std::vector<int> maxmyrank_per_bin(numbins,-1);
  fluiddis_->Comm().MaxAll(&myrank_per_bin[0], &maxmyrank_per_bin[0], numbins);

  // distribute bins to proc with highest rank
  for(int gid=0; gid<numbins; gid++)
  {
    if(myrank_ == maxmyrank_per_bin[gid])
    {
      Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEBIN","dummy", gid, myrank_);
      particledis_->AddElement(bin);
      rowbins.push_back(gid);
    }
  }

  myrank_per_bin.clear();
  maxmyrank_per_bin.clear();

  // return binrowmap
  return Teuchos::rcp(new Epetra_Map(-1,(int)rowbins.size(),&rowbins[0],0,Comm()));
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap)
{
  // 1st and 2nd step
  PARTICLE::Algorithm::SetupGhosting(binrowmap);

  // 3st step: extend ghosting of underlying fluid discretization according to bin distribution
  {
    // gather all fluid roweles /*and their nodes*/ in each bin
    std::map<int, std::set<int> > fluideles;
    /*std::map<int, std::set<int> > fluidnodes;*/
    for (int lid=0;lid<fluiddis_->NumMyRowElements();++lid)
    {
      DRT::Element* fluidele = fluiddis_->lRowElement(lid);
      DRT::Node** nodes = fluidele->Nodes();

      /*const int numnode = fluidele->NumNode();*/
      for(int inode=0; inode<fluidele->NumNode(); inode++)
      {
        DRT::Node* currnode = nodes[inode];
        const double* pos = currnode->X();
        int ijk[3] = {-1,-1,-1};
        for(int dim=0; dim < 3; dim++)
        {
          ijk[dim] = (int)(((pos)[dim]-XAABB_(dim,0)) / bin_size_[dim]);
        }
        int gidofbin = ConvertijkToGid(&ijk[0]);
        /*
        const int* nodeids = fluidele->NodeIds();
        // insert node ids into set
        fluidnodes[gidofbin].insert(nodeids,nodeids+numnode);
        */
        // insert ele id into set; one element can be part of several bins
        fluideles[gidofbin].insert(fluidele->Id());
      }
    }

    // do communication to gather all elements for extended ghosting
    const int numproc = fluiddis_->Comm().NumProc();

    for (int iproc = 0; iproc < numproc; ++iproc)
    {
      // first: proc i tells all procs how many col bins it has
      int numbin = bincolmap_->NumMyElements();
      fluiddis_->Comm().Broadcast(&numbin, 1, iproc);
      // second: proc i tells all procs which col bins it has
      std::vector<int> binid(numbin,0);
      if(iproc == myrank_)
      {
        int* bincolmap = bincolmap_->MyGlobalElements();
        for (int i=0; i<numbin; ++i)
          binid[i] = bincolmap[i];
      }
      fluiddis_->Comm().Broadcast(&binid[0], numbin, iproc);

      // loop over all own bins and find requested ones
      std::map<int, std::set<int> > sdata;
      std::map<int, std::set<int> > rdata;

      for(int i=0; i<numbin; i++)
      {
        sdata[binid[i]].insert(fluideles[binid[i]].begin(),fluideles[binid[i]].end());
      }

      LINALG::Gather<int>(sdata, rdata, 1, &iproc, fluiddis_->Comm());

      // proc i has to store the received data
      if(iproc == myrank_)
      {
        extendedfluidghosting_ = rdata;
      }
    }

    //reduce map of sets to one set and copy to a vector to create fluidcolmap
    std::set<int> redufluideleset;
    std::map<int, std::set<int> >::iterator iter;
    for(iter=extendedfluidghosting_.begin(); iter!= extendedfluidghosting_.end(); ++iter)
    {
      redufluideleset.insert(iter->second.begin(),iter->second.end());
    }
    std::vector<int> fluidcolgids(redufluideleset.begin(),redufluideleset.end());
    Teuchos::RCP<Epetra_Map> fluidcolmap = Teuchos::rcp(new Epetra_Map(-1,(int)fluidcolgids.size(),&fluidcolgids[0],0,Comm()));

    // create ghosting for fluid eles (each knowing its node ids)
    fluiddis_->ExportColumnElements(*fluidcolmap);

    // create a set of node IDs for each proc (row + ghost)
    std::set<int> nodes;
    for (int lid=0;lid<fluidcolmap->NumMyElements();++lid)
    {
      DRT::Element* ele = fluiddis_->gElement(fluidcolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0;inode<ele->NumNode();inode++)
        nodes.insert(nodeids[inode]);
    }

    // copy nodegids to a vector and create nodecolmap
    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,Comm()));

    // create ghosting for nodes
    fluiddis_->ExportColumnNodes(*nodecolmap);

    // do a final fillcomplete to build connectivity
    fluiddis_->FillComplete(true,true,true);

  }

#ifdef DEBUG
  // check whether each particle has an underlying fluid element
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for(int i=0; i<fluiddis_->NumMyColNodes(); i++)
  {
    DRT::Node* node = fluiddis_->lColNode(i);
    LINALG::Matrix<3,1> currpos;

    for (int a=0; a<3; a++)
    {
      currpos(a) = node->X()[a];
    }
    currentpositions.insert(std::pair<int,LINALG::Matrix<3,1> >(node->Id(),currpos));
  }
  // start loop over all particles
  for(int k=0; k<particledis_->NumMyColNodes(); k++)
  {
    DRT::Node* particle = particledis_->lColNode(k);
    const double* pos = particle->X();
    LINALG::Matrix<3,1> projpoint;
    for(int dim=0; dim<3; dim++)
      projpoint(dim) = pos[dim];
    bool foundele = false;
    for(int i=0; i<fluiddis_->NumMyColElements(); i++)
    {
      DRT::Element* fluidele = fluiddis_->lColElement(i);

      LINALG::Matrix<3,1> elecoord(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(Teuchos::rcp(fluidele,false), currentpositions));

      //get coordinates of the particle position in parameter space of the element
      foundele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, projpoint, elecoord);

      // THIS IS JUST TO CHECK WHETHER GEO::currentToVolumeElementCoordinates delivers correct results
      const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      LINALG::Matrix<3,numnode> xyze_linalg;
      for(int dim=0;dim<3;dim++)
        for(size_t n=0;n<numnode;n++)
          xyze_linalg(dim,n) = xyze(dim,n);
      GEO::CUT::Position<DRT::Element::hex8> pos(xyze_linalg, projpoint);
      bool withinele = pos.ComputeTol(GEO::TOL7);
      LINALG::Matrix<3,1> elecoordCut = pos.LocalCoordinates();

      if(abs(elecoordCut(0) - elecoord(0)) > GEO::TOL12 or abs(elecoordCut(1) - elecoord(1)) > GEO::TOL12
          or abs(elecoordCut(2) - elecoord(2)) > GEO::TOL12 or withinele != foundele)
        dserror("GEO::currentToVolumeElementCoordinates delivers different results compared to GEO::CUT::Position");
      // END: THIS IS JUST TO CHECK WHETHER GEO::currentToVolumeElementCoordinates delivers correct results

      if(foundele == true)
        break;
    }
    if(foundele == false)
      dserror("particle (Id:%d) was found which does not have fluid support", particle->Id());
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fluid_->CreateFieldTest());
  PARTICLE::Algorithm::TestResults(comm);
  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Output()
{
  fluid_->Output();
  PARTICLE::Algorithm::Output();
  return;
}

