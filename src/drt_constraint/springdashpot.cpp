/*!----------------------------------------------------------------------
\file patspec.cpp

\brief Methods for spring and dashpot constraints / boundary conditions:

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/


#include "springdashpot.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <iostream>
#include "../drt_io/io_pstream.H" // has to go before io.H
#include "../drt_io/io.H"
#include <Epetra_Time.h>



/*----------------------------------------------------------------------*
 |                                                             mhv 01/14|
 *----------------------------------------------------------------------*/
void SPRINGDASHPOT::SpringDashpot(Teuchos::RCP<DRT::Discretization> dis)
{

  //------------test discretization for presence of spring dashpot tissue condition
  std::vector<DRT::Condition*> springdashpotcond;

  // this is a copy of the discretization passed in to hand it to the element level.
  DRT::Discretization& actdis = *(dis);
  dis->GetCondition("SpringDashpot",springdashpotcond);
  if ((int)springdashpotcond.size())
  {
    //params->set("havespringdashpot", true); // so the condition is evaluated in STR::ApplyForceStiffSpringDashpot
    if (!dis->Comm().MyPID()) IO::cout << "Computing area for spring dashpot condition...\n";

    // loop over all spring dashpot conditions
    for (int cond=0; cond<(int)springdashpotcond.size(); ++cond)
    {

      const std::string* dir = springdashpotcond[cond]->Get<std::string>("direction");

      // a vector for all row nodes to hold element area contributions
      Epetra_Vector nodalarea(*dis->NodeRowMap(),true);
      // a vector for all row dofs to hold normals interpolated to the nodes
      Epetra_Vector refnodalnormals(*dis->DofRowMap(),true);

      std::map<int,Teuchos::RCP<DRT::Element> >& geom = springdashpotcond[cond]->Geometry();
      std::map<int,Teuchos::RCP<DRT::Element> >::iterator ele;
      for (ele=geom.begin(); ele != geom.end(); ++ele)
      {
        DRT::Element* element = ele->second.get();

        Teuchos::ParameterList eparams;
        Teuchos::ParameterList eparams2;
        eparams.set("action","calc_struct_area");
        eparams.set("area",0.0);
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        element->LocationVector(actdis,lm,lmowner,lmstride);
        const int eledim = (int)lm.size();
        Epetra_SerialDenseMatrix dummat(0,0);
        Epetra_SerialDenseVector dumvec(0);
        Epetra_SerialDenseVector elevector;
        elevector.Size(eledim);
        element->Evaluate(eparams,actdis,lm,dummat,dummat,dumvec,dumvec,dumvec);

        // when refsurfnormal direction is chosen, call the respective element evaluation routine and assemble
        if (*dir == "refsurfnormal")
        {
          eparams2.set("action","calc_ref_nodal_normals");
          element->Evaluate(eparams2,actdis,lm,dummat,dummat,elevector,dumvec,dumvec);
          LINALG::Assemble(refnodalnormals,elevector,lm,lmowner);
        }

        DRT::Element::DiscretizationType shape = element->Shape();

        double a = eparams.get("area",-1.0);

        // loop over all nodes of the element that share the area
        // do only contribute to my own row nodes
        double apernode = 0.;
        for (int i=0; i<element->NumNode(); ++i)
        {
          /* here we have to take care to assemble the right stiffness to the nodes!!! (mhv 05/2014):
          we do some sort of "manual" gauss integration here since we have to pay attention to assemble
          the correct stiffness in case of quadratic surface elements*/

          switch(shape)
          {
          case DRT::Element::tri3:
            apernode = a / element->NumNode();
          break;
          case DRT::Element::tri6:
          {
            //integration of shape functions over parameter element surface
            double int_N_cornernode = 0.;
            double int_N_edgemidnode = 1./6.;

            int numcornernode = 3;
            int numedgemidnode = 3;

            double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;

            //corner nodes
            if (i==0) apernode = int_N_cornernode * a / weight;
            if (i==1) apernode = int_N_cornernode * a / weight;
            if (i==2) apernode = int_N_cornernode * a / weight;
            //edge mid nodes
            if (i==3) apernode = int_N_edgemidnode * a / weight;
            if (i==4) apernode = int_N_edgemidnode * a / weight;
            if (i==5) apernode = int_N_edgemidnode * a / weight;
          }
          break;
          case DRT::Element::quad4:
            apernode = a / element->NumNode();
          break;
          case DRT::Element::quad8:
          {
            //integration of shape functions over parameter element surface
            double int_N_cornernode = -1./3.;
            double int_N_edgemidnode = 4./3.;

            int numcornernode = 4;
            int numedgemidnode = 4;

            double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode;

            //corner nodes
            if (i==0) apernode = int_N_cornernode * a / weight;
            if (i==1) apernode = int_N_cornernode * a / weight;
            if (i==2) apernode = int_N_cornernode * a / weight;
            if (i==3) apernode = int_N_cornernode * a / weight;
            //edge mid nodes
            if (i==4) apernode = int_N_edgemidnode * a / weight;
            if (i==5) apernode = int_N_edgemidnode * a / weight;
            if (i==6) apernode = int_N_edgemidnode * a / weight;
            if (i==7) apernode = int_N_edgemidnode * a / weight;
          }
          break;
          case DRT::Element::quad9:
          {
            //integration of shape functions over parameter element surface
            double int_N_cornernode = 1./9.;
            double int_N_edgemidnode = 4./9.;
            double int_N_centermidnode = 16./9.;

            int numcornernode = 4;
            int numedgemidnode = 4;
            int numcentermidnode = 1;

            double weight = numcornernode * int_N_cornernode + numedgemidnode * int_N_edgemidnode + numcentermidnode * int_N_centermidnode;

            //corner nodes
            if (i==0) apernode = int_N_cornernode * a / weight;
            if (i==1) apernode = int_N_cornernode * a / weight;
            if (i==2) apernode = int_N_cornernode * a / weight;
            if (i==3) apernode = int_N_cornernode * a / weight;
            //edge mid nodes
            if (i==4) apernode = int_N_edgemidnode * a / weight;
            if (i==5) apernode = int_N_edgemidnode * a / weight;
            if (i==6) apernode = int_N_edgemidnode * a / weight;
            if (i==7) apernode = int_N_edgemidnode * a / weight;
            //center mid node
            if (i==8) apernode = int_N_centermidnode * a / weight;
          }
          break;
          case DRT::Element::nurbs9:
            dserror("Not yet implemented for Nurbs! To do: Apply the correct weighting of the area per node!");
          break;
          default:
              dserror("shape type unknown!\n");
          }

          int gid = element->Nodes()[i]->Id();
          if (!dis->NodeRowMap()->MyGID(gid)) continue;
          nodalarea[dis->NodeRowMap()->LID(gid)] += apernode;
        }
      } // for (ele=geom.begin(); ele != geom.end(); ++ele)

      // now we have the area per node, put it in a vector that is equal to the nodes vector
      // consider only my row nodes
      const std::vector<int>* nodes = springdashpotcond[cond]->Nodes();

      std::vector<double> apern(nodes->size(),0.0);
      for (int i=0; i<(int)nodes->size(); ++i)
      {
        int gid = (*nodes)[i];
        if (!nodalarea.Map().MyGID(gid)) continue;
        apern[i] = nodalarea[nodalarea.Map().LID(gid)];
      }
      // set vector to the condition
      (*springdashpotcond[cond]).Add("areapernode",apern);


      if (*dir == "refsurfnormal")
      {
        std::vector<double> refndnorms(3 * nodes->size(),0.0);

        const std::vector<int>& nds = *nodes;
        for (int j=0; j<(int)nds.size(); ++j)
        {
          if (dis->NodeRowMap()->MyGID(nds[j]))
          {
            int gid = nds[j];

            DRT::Node* node = dis->gNode(gid);
            if (!node) dserror("Cannot find global node %d",gid);

            int numdof = dis->NumDof(node);
            std::vector<int> dofs = dis->Dof(node);

            assert (numdof==3);

            for (int k=0; k<numdof; ++k)
            {
              //if (!refnodalnormals.Map().MyGID(dofs[k])) continue;
              refndnorms[numdof*j+k] = refnodalnormals[refnodalnormals.Map().LID(dofs[k])];
            }
          }
        }
        // set vector to the condition
        (*springdashpotcond[cond]).Add("refnodalnormals",refndnorms);
      }

    } // for (int cond=0; cond<(int)springdashpotcond.size(); ++cond)

  }// if ((int)springdashpotcond.size())
  //-------------------------------------------------------------------------


  return;
}


/*----------------------------------------------------------------------*
 |                                                             mhv 01/14|
 *----------------------------------------------------------------------*/
void SPRINGDASHPOT::EvaluateSpringDashpot(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<LINALG::SparseOperator> stiff,
    Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> velo,
    Teuchos::ParameterList parlist)
{

  if (disp==Teuchos::null) dserror("Cannot find displacement state in discretization");

  double gamma = parlist.get("scale_gamma",0.0);
  double beta = parlist.get("scale_beta",1.0);
  double ts_size = parlist.get("time_step_size",1.0);

  std::vector<DRT::Condition*> springdashpotcond;
  discret->GetCondition("SpringDashpot",springdashpotcond);


  for (int i=0; i<(int)springdashpotcond.size(); ++i)
  {
    const std::vector<int>* nodes = springdashpotcond[i]->Nodes();
    double springstiff_tens = springdashpotcond[i]->GetDouble("SPRING_STIFF_TENS");
    double springstiff_comp = springdashpotcond[i]->GetDouble("SPRING_STIFF_COMP");
    double springoffset = springdashpotcond[i]->GetDouble("SPRING_OFFSET");
    double dashpotvisc = springdashpotcond[i]->GetDouble("DASHPOT_VISCOSITY");
    const std::string* dir = springdashpotcond[i]->Get<std::string>("direction");

    const std::vector<double>* areapernode = springdashpotcond[i]->Get< std::vector<double> >("areapernode");
    const std::vector<double>* refnodalnormals = springdashpotcond[i]->Get< std::vector<double> >("refnodalnormals");


    const std::vector<int>& nds = *nodes;
    for (int j=0; j<(int)nds.size(); ++j)
    {

      if (discret->NodeRowMap()->MyGID(nds[j]))
      {
        int gid = nds[j];
        double nodalarea = (*areapernode)[j];

        DRT::Node* node = discret->gNode(gid);

        if (!node) dserror("Cannot find global node %d",gid);


        int numdof = discret->NumDof(node);
        std::vector<int> dofs = discret->Dof(node);

        assert (numdof==3);

        // displacement vector of condition nodes
        std::vector<double> u(numdof);
        for (int k=0; k<numdof; ++k)
        {
          u[k] = (*disp)[disp->Map().LID(dofs[k])];
        }

        // velocity vector of condition nodes
        std::vector<double> v(numdof);
        for (int k=0; k<numdof; ++k)
        {
          v[k] = (*velo)[velo->Map().LID(dofs[k])];
        }

        // assemble into residual and stiffness matrix for case that spring / dashpot acts in every surface dof direction
        if (*dir == "all")
        {
          for (int k=0; k<numdof; ++k)
          {
            double val = nodalarea*(springstiff_tens*(u[k]-springoffset) + dashpotvisc*v[k]);
            double dval = nodalarea*(springstiff_tens + dashpotvisc*gamma/(beta*ts_size));

            if (springstiff_tens != springstiff_comp)
            {
              dserror("SPRING_STIFF_TENS != SPRING_STIFF_COMP: Different spring moduli for tension and compression not supported "
                  "when specifying 'all' as DIRECTION (no ref surface normal information is calculated for that case)! "
                  "Only possible for DIRECTION 'refsurfnormal'.");
            }

            int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
            if (err) dserror("SumIntoGlobalValues failed!");
            stiff->Assemble(dval,dofs[k],dofs[k]);
          }
        }

        // assemble into residual and stiffness matrix for case that spring / dashpot acts in reference surface normal direction
        else if (*dir == "refsurfnormal")
        {
          // extract averaged nodal ref normal and compute its absolute value
          std::vector<double> unitrefnormal(numdof);
          double temp_ref = 0.;
          double proj = 0.;
          for (int k=0; k<numdof; ++k)
          {
            unitrefnormal[k] = (*refnodalnormals)[numdof*j+k];
            temp_ref += unitrefnormal[k]*unitrefnormal[k];
          }

          double unitrefnormalabsval = sqrt(temp_ref);

          for(int k=0; k<numdof; ++k)
          {
            unitrefnormal[k] /= unitrefnormalabsval;
          }

          // projection of displacement vector to refsurfnormal (u \cdot N)
          for(int k=0; k<numdof; ++k)
          {
            proj += u[k]*unitrefnormal[k];
          }

          //dyadic product of ref normal with itself (N \otimes N)
          Epetra_SerialDenseMatrix N_x_N(numdof,numdof);
          for (int l=0; l<numdof; ++l)
            for (int m=0; m<numdof; ++m)
              N_x_N(l,m)=unitrefnormal[l]*unitrefnormal[m];

          for (int k=0; k<numdof; ++k)
          {
            for (int m=0; m<numdof; ++m)
            {
              double val = 0.;
              double dval = 0.;

              // projection of displacement to refsurfnormal negative: we have a tensile spring
              if (proj < 0)
              {
                val = nodalarea*N_x_N(k,m)*(springstiff_tens*(u[m]-springoffset) + dashpotvisc*v[m]);
                dval = nodalarea*(springstiff_tens + dashpotvisc*gamma/(beta*ts_size))*N_x_N(k,m);
              }

              // projection of displacement to refsurfnormal positive: we have a compressive spring
              if (proj >= 0)
              {
                val = nodalarea*N_x_N(k,m)*(springstiff_comp*(u[m]-springoffset) + dashpotvisc*v[m]);
                dval = nodalarea*(springstiff_comp + dashpotvisc*gamma/(beta*ts_size))*N_x_N(k,m);
              }

              int err = fint->SumIntoGlobalValues(1,&val,&dofs[k]);
              if (err) dserror("SumIntoGlobalValues failed!");
              stiff->Assemble(dval,dofs[k],dofs[m]);
            }
          }
        }
        else dserror("Invalid direction option! Choose DIRECTION all or DIRECTION refsurfnormal!");

      } //node owned by processor

    } //loop over nodes

  } //loop over conditions

  return;
}



