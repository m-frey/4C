/*!------------------------------------------------------------------------------------------------*
\file xfluid_timeinnt_std_SemiLagrange.cpp

\brief provides the SemiLagrangean class

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241

</pre>
 *------------------------------------------------------------------------------------------------*/

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_integrationcell.H"
#include "../drt_cut/cut_volumecell.H"

#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"

#include "dofkey.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_xfem/xfem_fluiddofset.H"
#include "../drt_xfem/xfem_fluidwizard.H"

#include "xfluid_timeInt_std_SemiLagrange.H"



//#define DEBUG_SEMILAGRANGE

/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrange Back-Tracking algorithm constructor                                 schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_SemiLagrange::XFLUID_SemiLagrange(
    XFEM::XFLUID_TIMEINT_BASE&                                     timeInt,          /// time integration base class object
    const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,  /// reconstruction map for nodes and its dofsets
    INPAR::XFEM::XFluidTimeInt&                                    timeIntType,      /// type of time integration
    const RCP<Epetra_Vector>                                       veln,             /// velocity at time t^n in col map
    const double&                                                  dt,               /// time step size
    const double&                                                  theta,            /// OST theta
    bool                                                           initialize        /// is initialization?
) :
XFLUID_STD(
    timeInt,
    reconstr_method,
    timeIntType,
    veln,
    dt,
    initialize),
theta_default_(theta),
relTolIncr_(1.0e-10),
relTolRes_(1.0e-10)
{
  return;
} // end constructor




/*------------------------------------------------------------------------------------------------*
 * Semi-Lagrangean Back-Tracking main algorithm                                      schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::compute(
    std::vector<RCP<Epetra_Vector> >& newRowVectorsn
)
{
  const int nsd = 3; // 3 dimensions for a 3d fluid element
  handleVectors(newRowVectorsn);

  // REMARK: in case of a new FGI iteration we have values! at new position
  newIteration_prepare(newVectors_);

  switch (FGIType_)
  {
  case FRS1FGI1_:
  {
    IO::cout << "\nXFLUID_SemiLagrange::compute: case FRS1FGI1_" << IO::endl;
    resetState(TimeIntData::basicStd_,TimeIntData::currSL_);
    break;
  }
  case FRSNot1_:
  {
    IO::cout << "\nXFLUID_SemiLagrange::compute: case FRSNot1_" << IO::endl;
    resetState(TimeIntData::doneStd_,TimeIntData::currSL_);
    break;
  }
  case FRS1FGINot1_:
  {
    IO::cout << "\nXFLUID_SemiLagrange::compute: case FRS1FGINot1_" << IO::endl;
    reinitializeData();
    resetState(TimeIntData::basicStd_,TimeIntData::currSL_);
    resetState(TimeIntData::doneStd_,TimeIntData::currSL_);
    break;
  }
  default: dserror("not implemented"); break;
  } // end switch


#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n----------------------------------------------------------------------------------------- " ;
  IO::cout << "\nReconstruct data with SEMILAGRANGEAN algorithm for " << timeIntData_->size() << " dofsets " ;
  IO::cout << "\n----------------------------------------------------------------------------------------- " << IO::endl;
#endif

  /*----------------------------------------------------*
   * first part: get the correct origin for the node    *
   * in a lagrangian point of view using a newton loop  *
   *----------------------------------------------------*/
  int counter = 0; // loop counter to avoid infinite loops

  // loop over nodes which still don't have and may get a good startvalue
  while (true)
  {
    counter += 1;

    // counter limit because maximal max_iter Newton iterations with maximal
    // numproc processor changes per iteration (avoids infinite loop)
    if (!globalNewtonFinished(counter))
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n==============================================";
      IO::cout << "\n CONTINUE GLOBAL NEWTON (" << counter << ") on proc " << myrank_;
      IO::cout << "\n==============================================" << IO::endl;
#endif

      // loop over all nodes (their std-dofsets) that have been chosen for SEMI-Lagrangean reconstruction
      for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
          data!=timeIntData_->end(); data++)
      {

#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t * STD-SL algorithm for node " << data->node_.Id();
#endif

        //------------------------------------
        // find element the initial startpoint lies in, if not done yet
        //------------------------------------
        if(data->initial_eid_ == -1) // initial point not found yet
        {
          bool initial_elefound = false;              // true if the element for a point was found on the processor
          DRT::Element* initial_ele = NULL;           // pointer to the element where start point lies in
          LINALG::Matrix<nsd,1> initial_xi(true);     // local transformed coordinates of x w.r.t found ele

          // set the element pointer where the initial point lies in!
          elementSearch(initial_ele,data->initialpoint_,initial_xi,initial_elefound);

          if(!initial_elefound)
          {
            if (data->searchedProcs_ < numproc_)
            {
              // set state to nextSL to proceed with these data on the next proc
              data->state_ = TimeIntData::nextSL_;
              data->searchedProcs_ += 1; // increment counter that the element the point lies in has not been found on this processor
              data->initial_eid_ = -1;
            }
            else // all procs searched -> initial point not in domain
            {
              data->state_ = TimeIntData::failedSL_;
              dserror("<<< WARNING! Initial point for node %d for finding the Lagrangean origin not in domain! >>>", data->node_.Id());
            }
          }
          else
          {
            data->initial_eid_ = initial_ele->Id();
            data->startpoint_ = data->initialpoint_; // start with the initial point as startpoint approximation
            data->initial_ele_owner_ = myrank_;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t\t -> Initial point found in element " << data->initial_eid_;
#endif

          }
        }

#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t\t -> start with start point approximation: " << data->startpoint_;
#endif

        if (data->state_ == TimeIntData::currSL_) // do not proceed when nextSL is set for current data
        {
          //------------------------------------
          // find element the initial startpoint lies in
          // if found, then get the element information
          //------------------------------------

          // Initialization
          bool elefound = false;              // true if the element for a point was found on the processor
          DRT::Element* ele = NULL;           // pointer to the element where start point lies in
          LINALG::Matrix<nsd,1> xi(true);     // local transformed coordinates of x w.r.t found ele
          LINALG::Matrix<nsd,1> vel(true);    // velocity of the start point approximation

          // search for an element where the current startpoint lies in
          elementSearch(ele,data->startpoint_,xi,elefound);

          //------------------------------------
          // if element is found on this proc, the newton iteration to find a better startpoint can start
          //------------------------------------
          if(elefound)
          {
#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t\t\t ... start point approximation found in element: " << ele->Id();
#endif
            //----------------------------------------------
            DRT::Element * initial_ele = discret_->gElement(data->initial_eid_);

            if(initial_ele == NULL)
            {
              dserror("initial element %d not available on proc %d! -> This issue can be solved, see code!", initial_ele->Id(), myrank_);
              //TODO: modify the ChangedSide check for intersections with all sides in the boundary-dis
              // this is not so efficient but should not be called not so often
              // the check itself does not need information about the background elements
            }

            data->changedside_ = ChangedSide(ele, data->startpoint_,  false,
                                             initial_ele, data->initialpoint_, false);

            //----------------------------------------------
            // get dofset w.r.t to old interface position
            bool step_np = false;
            std::vector<int> nds_curr;
            getNodalDofSet(ele, data->startpoint_,nds_curr, data->last_valid_vc_, step_np);


            //=========================================================
            // how to continue if no side changing comparison possible...
            //=========================================================
            if(data->changedside_) // how to continue in newton loop, or stop the newton loop
            {
              if (!continueForChangingSide(&*data, ele, nds_curr)) continue; // continue with next TimintData
            }
            else // point did not change the side
            {
              // this set is a valid fluid set on the right side of the interface
              data->last_valid_nds_ = nds_curr;
              data->last_valid_ele_ = ele->Id();
              data->nds_ = nds_curr;
            }

            //-------------------------------------------------
            // Newton loop just for sensible points
            //-------------------------------------------------

            //----------------------------------------------
            // compute the velocity at startpoint
            LINALG::Matrix<nsd,nsd> vel_deriv_tmp(true); // dummy matrix for velocity derivatives
            getGPValues(ele, xi, data->nds_, step_np, vel, vel_deriv_tmp, false);

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t\t\t ... computed velocity at start point approximation: " << vel;
#endif

            //----------------------------------------------------------------------------------------
            // call the Newton loop to get the right Lagrangean origin
            // REMARK: if newton loop is converged then return the element,
            //         local coordinates and velocity at lagrangean origin
            NewtonLoop(ele,&*data,xi,vel,elefound);
            //----------------------------------------------------------------------------------------


            // if iteration can go on (that is when startpoint is on
            // correct interface side and iter < max_iter)
            if ((data->counter_<newton_max_iter_) and (data->state_==TimeIntData::currSL_))
            {
              // if element is not found in a newton step, look at another processor and so add
              // all according data to the vectors which will be sent to the next processor
              if (!elefound)
              {
                // here, a point has not been found on this proc for a second time
                data->searchedProcs_ = 2;
                data->state_ = TimeIntData::nextSL_;
              }
              else
              {
                //----------------------------------------------------------------------------------------
                // newton iteration converged to a good startpoint and so the data can be used to go on
                if(data->accepted_)
                {
                  callBackTracking(ele,&*data,xi,"standard");
                }
                else // converged Lagrangean origin, however not accepted
                {
                  //----------------------------------------------------------------------
                  // a Lagrangean origin has been found but it does not lie in the fluid
                  //----------------------------------------------------------------------

                  LINALG::Matrix<nsd,1> proj_x(true);

                  FindNearestSurfPoint(data->startpoint_, proj_x, data->last_valid_vc_, "idispn");

                  proj_x = data->initialpoint_;

                  elementSearch(ele,proj_x,xi,elefound);

                  if( elefound )
                  {
                    // check if the projected point still remains in the same element
                    if(ele->Id() == data->last_valid_ele_ )
                    {
                      // TODO: check if the point lies on Boundary of the last valid vc

                      //set the new startpoint
                      data->startpoint_ = proj_x;

                      callBackTracking(ele,&*data,xi,"standard");
                    }
                    else
                    {
                      dserror("projection of startpoint lies in another element compared to the point to be projected");
                      data->state_ = TimeIntData::failedSL_;
                    }
                  }
                  else dserror("element where the projection point lies in not available on this proc");

                }
                //----------------------------------------------------------------------------------------
              }
            } // end if
            // maximum number of iterations reached or converged origin on wrong interface side
            else if (data->counter_ == newton_max_iter_ or (!data->accepted_))
            {
              // do not use the lagrangian origin since this case is strange and potential dangerous
              data->state_ = TimeIntData::failedSL_;

#ifdef DEBUG_SEMILAGRANGE
              IO::cout << " <<< WARNING: newton iteration to find start value did not converge! >>>" << IO::endl;
#endif
            } // not converged in max_iter
          } // if(elefound)
          // if element is not found, look at another processor and so add all
          // according data to the vectors which will be sent to the next processor
          else // (!elefound)
          {
            if (data->searchedProcs_ < numproc_)
            {
              data->state_ = TimeIntData::nextSL_;
              data->searchedProcs_ += 1;
            }
            else // all procs searched -> point not in domain
            {
              data->state_ = TimeIntData::failedSL_;
              IO::cout << " <<< WARNING! Lagrangian start point not in domain! >>>" << IO::endl;
            }
          } // end if elefound

        }
      } // end loop over all nodes stored in TiminitData marked for Semilagrangean algorithm
    } // !globalNewtonFinished(counter)
    else
    {
      // reset the state to failed
      resetState(TimeIntData::currSL_,TimeIntData::failedSL_);
    }

    //===================================================================
    //                     PARALLEL COMMUNICATION

    // export nodes and according data for which the startpoint isn't still found (next_ vector) to next proc
    bool procDone = globalNewtonFinished();

#ifdef DEBUG_SEMILAGRANGE
    if(procDone)
    {
      IO::cout << "\n==============================================";
      IO::cout << "\n FINISHED GLOBAL NEWTON on proc " << myrank_;
      IO::cout << "\n==============================================" << IO::endl;
    }
#endif

    exportIterData(procDone);

    // convergencecheck: procfinished == 1 just if all procs have finished
    if (procDone)
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n-------------------------------------------------";
      IO::cout << "\n\t\t\t !!!!!!!!!! procDone!!!!!!!!";
      IO::cout << "\n-------------------------------------------------" << IO::endl;
#endif
      break;
    }

    //===================================================================

  } // end while loop over searched nodes


  /*-----------------------------------------------------------------------------*
   * second part: get sensible startvalues for nodes where the algorithm failed, *
   * using another algorithm, and combine the "Done" and the "Failed" - vectors  *
   *-----------------------------------------------------------------------------*/
  if (FGIType_==FRSNot1_) // failed nodes stay equal after the first computation
    clearState(TimeIntData::failedSL_);
  else
  {

    //===================================================================
    //                     PARALLEL COMMUNICATION
    exportAlternativAlgoData(); // export data of failed nodes
    //===================================================================

    getDataForNotConvergedNodes(); // compute final data for failed nodes
  }



  /*-----------------------------------------------------------*
   * third part: set the computed values into the state vector *
   *-----------------------------------------------------------*/
  // send the computed startvalues for every node which needs
  // new start data to the processor where the node is
  exportFinalData();

  // now every proc has the whole data for the nodes and so the data can be set to the right place now
  setFinalData();

#ifdef DEBUG
  if (counter > 8*numproc_) // too much loops shouldnt be if all this works
    cout << "WARNING: semiLagrangeExtrapolation seems to run an infinite loop!" << endl;
#endif
} // end semiLagrangeExtrapolation



/*------------------------------------------------------------------------------------------------*
 * Main Newton loop of the Semi-Lagrangian Back-Tracking algorithm                   schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::NewtonLoop(
    DRT::Element*&          ele,              /// pointer to element
    TimeIntData*            data,             /// current data
    LINALG::Matrix<3,1>&    xi,               /// local coordinates of point
    LINALG::Matrix<3,1>&    vel,              /// velocity at current point
    bool&                   elefound          /// is element found ?
)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t -> XFLUID_SemiLagrange::NewtonLoop" << IO::endl;
#endif

  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization
  LINALG::Matrix<nsd,1> residuum(true);             // residuum of the newton iteration
  LINALG::Matrix<nsd,1> incr(true);                 // increment of the newton system

  // coordinates of endpoint of Lagrangian characteristics
  LINALG::Matrix<nsd,1> origNodeCoords(true);
  for (int i=0;i<nsd;i++)
    origNodeCoords(i) = data->node_.X()[i];

  if(data->node_.Id() == 1656) cout << "node coords: " << origNodeCoords << endl;

  //-------------------------------------------------------
  // initialize residual (Theta = 0 at predictor step)
  residuum.Clear();

  // data->vel_ = vel^(n+1) for FGI>1, vel = vel^n
  residuum.Update((1.0-Theta(data)),vel,Theta(data),data->vel_);   // dt*v(data->startpoint_)
  residuum.Update(1.0,data->startpoint_,-1.0,origNodeCoords,dt_);  // R = data->startpoint_ - data->node_ + dt*v(data->startpoint_)

  //==========================================================
  // (re-)start the Newton-loop on this processor
  //==========================================================
  while(data->counter_ < newton_max_iter_)  // newton loop
  {

#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t\t NewtonLoop("<< data->counter_ <<"): residuum " << residuum << IO::endl;
#endif

    data->counter_ += 1;

    //-------------------------------------
    // compute a new Newton iteration
    //-------------------------------------
    switch (ele->Shape())
    {
    case DRT::Element::hex8:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
      NewtonIter<numnode,DRT::Element::hex8>(ele,data,xi,residuum,incr,elefound);
    }
    break;
    case DRT::Element::hex20:
    {
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
      NewtonIter<numnode,DRT::Element::hex20>(ele,data,xi,residuum,incr,elefound);
    }
    break;
    default:
      dserror("element type not yet implemented in time integration"); break;
    }; // end switch element type

    //=========================================================
    // continue on this proc if the new startpoint approximation is also on this processor
    //=========================================================
    if (elefound) // element of data->startpoint_ at this processor
    {
      IO::cout << "\n\t\t\t\t ... elefound " << ele->Id();

      DRT::Element* initial_ele = discret_->gElement(data->initial_eid_);

      if(initial_ele == NULL) dserror("element where initial point lies in not available on proc %d, no ChangedSide comparison possible", myrank_);

      data->changedside_ = ChangedSide(ele, data->startpoint_,false, initial_ele, data->initialpoint_, false);


      bool step_np = false; // new timestep or old timestep
      std::vector<int> nds_curr;
      getNodalDofSet(ele, data->startpoint_, nds_curr, data->last_valid_vc_, step_np);

      //=========================================================
      // how to continue if no side changing comparison possible...
      //=========================================================
      if(data->changedside_) // how to continue in newton loop, or stop the newton loop
      {
        if (!continueForChangingSide(data, ele, nds_curr)) break;
      }
      else // point did not change the sie
      {
        // this set is a valid fluid set on the right side of the interface
        data->last_valid_nds_ = nds_curr;
        data->last_valid_ele_ = ele->Id();
        data->nds_ = nds_curr;
      }

      //=========================================================
      // continue if the semi-lagrangean origin approximation did not change the side w.r.t the initial start point
      //=========================================================

      //-------------------------------------------------------
      // compute the velocity at startpoint
      LINALG::Matrix<nsd,nsd> vel_deriv_tmp(true); // dummy matrix
      getGPValues(ele, xi, data->nds_, step_np, vel, vel_deriv_tmp, false);

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t\t\t ... computed velocity at start point approximation: " << vel;
#endif

      //-------------------------------------------------------
      //reset residual
      residuum.Clear();
      residuum.Update((1.0-Theta(data)),vel,Theta(data),data->vel_);   // dt*v(data->startpoint_)
      residuum.Update(1.0,data->startpoint_,-1.0,origNodeCoords,dt_);  // R = data->startpoint_ - data->movNode_ + dt*v(data->startpoint_)

      //-------------------------------------------------------
      // convergence criterion
      if (data->startpoint_.Norm2()>1e-3)
      {
        if (incr.Norm2()/data->startpoint_.Norm2() < relTolIncr_ && residuum.Norm2()/data->startpoint_.Norm2() < relTolRes_)
        {

          if(data->changedside_ == false)
          {
            data->accepted_ = true;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!"      ;
            IO::cout << "\n\t  LAGRANGEAN ORIGIN ACCEPTED"    ;
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }
          else
          {
            data->accepted_ = false;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!"     ;
            IO::cout << "\n\t  LAGRANGEAN ORIGIN NOT (!!!) ACCEPTED";
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }

          break;
        }
      }
      else
      {
        if (incr.Norm2() < relTolIncr_ && residuum.Norm2() < relTolRes_)
        {

          if(data->changedside_ == false)
          {
            data->accepted_ = true;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!"     ;
            IO::cout << "\n\t  LAGRANGEAN ORIGIN ACCEPTED"   ;
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }
          else
          {
            data->accepted_ = false;

#ifdef DEBUG_SEMILAGRANGE
            IO::cout << "\n\t*******************************";
            IO::cout << "\n\t    NewtonLoop: converged!"     ;
            IO::cout << "\n\t  LAGRANGEAN ORIGIN NOT (!!!) ACCEPTED";
            IO::cout << "\n\t*******************************" << IO::endl;
#endif
          }

          break;
        }
      }

    } // end if elefound is true
    //=========================================================
    // stop Newton loop on this proc since the new startpoint approximation is not on this processor
    //=========================================================
    else // element of data->startpoint_ not at this processor
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\t <<< !!! element not found on this proc -> stop Newton loop on this proc !!! >>>" << IO::endl;
#endif
      break; // stop newton loop on this proc
    }
  } // end while Newton loop



#ifdef DEBUG_SEMILAGRANGE
  // did newton iteration converge?
  if(data->counter_ == newton_max_iter_)
  {
    IO::cout << "\t <<< WARNING: newton iteration for finding start value not converged for point !!! >>>" << IO::endl;
  }
#endif


} // end function NewtonLoop



/*------------------------------------------------------------------------------------------------*
 * One Newton iteration of the Semi-Lagrangian Back-Tracking algorithm               schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::NewtonIter(
    DRT::Element*&          ele,              /// pointer to element to be updated
    TimeIntData*            data,             /// current data to be updated
    LINALG::Matrix<3,1>&    xi,               /// local coordinates w.r.t ele to be updated
    LINALG::Matrix<3,1>&    residuum,         /// residual for semilagrangean backtracking
    LINALG::Matrix<3,1>&    incr,             /// computed increment for lagrangean origin to be updated
    bool&                   elefound          /// element found ?
)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t\t\t ... new iteration";
#endif

  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization
  LINALG::Matrix<nsd,1>   vel_dummy(true); // dummy matrix for the velocity
  LINALG::Matrix<nsd,nsd> vel_deriv(true); // matrix for the velocity derivatives
  LINALG::Matrix<nsd,nsd> sysmat(true);    // matrix for the newton system

  int step_np = false;

  // compute the velocity derivatives at startpoint
  getGPValues(ele, xi, data->nds_, step_np, vel_dummy, vel_deriv, true);

  // build sysmat
  // JAC = I + dt(1-theta)*velDerivXY
  sysmat.Update((1.0-Theta(data))*dt_,vel_deriv); // dt*(1-theta)dN/dx

  if(data->node_.Id() == 1656) cout << "velderiv " << vel_deriv << endl;

  for (int i=0;i<nsd;i++)
    sysmat(i,i) += 1.0; // I + dt*velDerivXY

  if(data->node_.Id() == 1656) cout << "sysmat " << sysmat << endl;

  // invers system Matrix built
  sysmat.Invert();

  if(data->node_.Id() == 1656) cout << "sysmat_invert " << sysmat << endl;
  if(data->node_.Id() == 1656) cout << "residuum " << residuum << endl;


  //solve Newton iteration
  incr.Clear();
  incr.Multiply(-1.0,sysmat,residuum); // incr = -Systemmatrix^-1 * residuum

  if(data->node_.Id() == 1656)  cout << "incr " << incr << endl;

  // update iteration
  for (int i=0;i<nsd;i++)
    data->startpoint_(i) += incr(i);

#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t\t\t ... new approximate startvalue is " << data->startpoint_(0) << " "
                                                              << data->startpoint_(1) << " "
                                                              << data->startpoint_(2) << IO::endl;
#endif

  // find the element the new approximation lies in
  elementSearch(ele, data->startpoint_, xi,elefound);


  return;
} // end function NewtonIter


/*------------------------------------------------------------------------------------------------*
 * check if newton iteration searching for the Lagrangian origin has finished        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_SemiLagrange::globalNewtonFinished( int counter ) const
{
  if (counter == newton_max_iter_*numproc_)
    return true; // maximal number of iterations reached
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if ((data->state_==TimeIntData::currSL_) or
        (data->state_==TimeIntData::nextSL_))
    {
      return false; // one node requires more data
    }
  }
  return true; // if no more node requires data, we are done
}


/*------------------------------------------------------------------------------------------------*
 * Decide how or if to continue when the startpoint approximation changed the side  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_SemiLagrange::continueForChangingSide(
    TimeIntData*            data,             ///< current data to be updated
    DRT::Element*           ele,              ///< pointer to element the current point lies in
    std::vector<int>&       nds_curr          ///< nds-vector of current volumecell the current startpoint approximation lies in
)
{
#if(1)
  //--------------------------------------------------------------------------------------
  //ALTERNATIVE: CONTINUE NEWTON-ALGO when startvalue changed side during newton
  // maybe the newton turns back to the right interface side


  if(nds_curr == data->last_valid_nds_ and ele->Id() == data->last_valid_ele_)
  {
    // the new newton step is within the same element and has the same nds-vector(same cell set)
    // but changed the side, then we are at the tip of a thin structure -> failed
#ifdef DEBUG_SEMILAGRANGE
    IO::cout << "\n -----------------------------------------------------------------------------------------------";
    IO::cout << "\n <<< Startpoint approximation moved within one fld-vc, but the trace intersects the side >>>";
    IO::cout << "\n                          CHANGED SIDE ";
    IO::cout << "\n Newton stopped! We are at the tip of a thin structure! -> leave newton loop >>>";
    IO::cout << "\n -----------------------------------------------------------------------------------------------" << IO::endl;
#endif
    data->state_ = TimeIntData::failedSL_;

    return false;
  }
  else if(nds_curr != data->last_valid_nds_ and ele->Id() == data->last_valid_ele_)
  {
    // the new newton step is within the same element but has a different nds-vector
    // we are within the structure or changed the side completely
    // it can be that the newton iterations goes back on the right side
    // -> continue within this element using the last valid nds-vector

    nds_curr = data->last_valid_nds_;

    return true;
  }
  else if( ele->Id() != data->last_valid_ele_)
  {
    // within the newton the element and the side have changed
#ifdef DEBUG_SEMILAGRANGE
    IO::cout << " <<< Newton for lagrangian origin can not be continued, iteration changed the side and the element! -> leave newton loop >>>" << IO::endl;
#endif
    data->state_ = TimeIntData::failedSL_;

    return false;
  }
  else dserror("case not possible");


#else
  //--------------------------------------------------------------------------------------
  //ALTERNATIVE: STOP NEWTON-ALGO when startvalue changed side during newton
#ifdef DEBUG_SEMILAGRANGE
  cout << "!!! CHANGED SIDE within Newton loop !!!! -> leave newton loop" << endl;
#endif
  data->state_ = TimeIntData::failedSL_;
  return false; // leave newton loop if point is on wrong domain side
  //--------------------------------------------------------------------------------------
#endif

  return false;
}


/*------------------------------------------------------------------------------------------------*
 * Computing final data where Semi-Lagrangian approach failed                        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::getDataForNotConvergedNodes()
{

  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // remark: all data have to be sent to the processor where
  //         the startpoint lies before calling this function
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::failedSL_)
    {
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "WARNING: failedSL -> alternative algo!" << IO::endl;
      IO::cout << "node " << data->node_.Id() << IO::endl;
      IO::cout << "use initial point: " << data->initialpoint_ << IO::endl;
#endif

      // Initialization
      DRT::Element* ele = NULL;          // pointer to the element where pseudo-Lagrangian origin lies in
      LINALG::Matrix<nsd,1> xi(true);    // local coordinates of pseudo-Lagrangian origin
      LINALG::Matrix<nsd,1> vel(true);   // velocity at pseudo-Lagrangian origin
      bool elefound = false;             // true if an element for a point was found on the processor

      // search for an element where the current startpoint lies in
      elementSearch(ele,data->initialpoint_,xi,elefound);

      // if found, give out all data at the startpoint
      if(elefound)
      {
        bool step_np = false; // data w.r.t old interface position
        getNodalDofSet(ele, data->initialpoint_,data->nds_, data->last_valid_vc_, step_np);
      }
      else // possibly slave node looked for element of master node or vice versa
      {
        dserror("element not found");
      }

      //-------------------------------------------------------------------
      // call the back Tracking computation based on the initial point
      // which is a rough approximation of the lagrangian origin
      callBackTracking(ele,&*data,xi,static_cast<const char*>("failing"));
      //-------------------------------------------------------------------

    } // if(failedSL_)
  } // end loop over nodes

  return;
} // end getDataForNotConvergedNodes


/*------------------------------------------------------------------------------------------------*
 * rewrite data for new computation                                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::newIteration_prepare(
    std::vector<RCP<Epetra_Vector> > newRowVectors
)
{
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    data->searchedProcs_ = 1;
    data->counter_ = 0;
    data->velValues_.clear();
    data->presValues_.clear();
  }

  //TODO: this has to be still done
  newIteration_nodalData(newRowVectors); // data at t^n+1 not used in predictor
  newRowVectors.clear(); // no more needed
}



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::newIteration_nodalData(
    std::vector<RCP<Epetra_Vector> > newRowVectors
)
{
//  dserror("adapt implementation of this function");

  IO::cout << "newIteration_nodalData not implemented yet"  << IO::endl;

//  const int nsd = 3;
//
//  // data about column vectors required
//  const Epetra_Map& newdofcolmap = *discret_->DofColMap();
//  std::map<XFEM::DofKey,XFEM::DofGID> newNodalDofColDistrib;
//  newdofman_->fillNodalDofColDistributionMap(newNodalDofColDistrib);
//
//  std::vector<RCP<Epetra_Vector> > newColVectors;
//
//  for (size_t index=0;index<newRowVectors.size();index++)
//  {
//    RCP<Epetra_Vector> tmpColVector = Teuchos::rcp(new Epetra_Vector(newdofcolmap,true));
//    newColVectors.push_back(tmpColVector);
//    LINALG::Export(*newRowVectors[index],*newColVectors[index]);
//  }
//
//  // computed data
//  std::vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1(static_cast<int>(oldVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
//  std::vector<LINALG::Matrix<1,nsd> > presnpDeriv1(static_cast<int>(oldVectors_.size()),LINALG::Matrix<1,nsd>(true));
//
//  std::vector<LINALG::Matrix<nsd,nsd> > velnpDeriv1Tmp(static_cast<int>(oldVectors_.size()),LINALG::Matrix<nsd,nsd>(true));
//  std::vector<LINALG::Matrix<1,nsd> > presnpDeriv1Tmp(static_cast<int>(oldVectors_.size()),LINALG::Matrix<1,nsd>(true));
//
//  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
//      data!=timeIntData_->end(); data++)
//  {
//    DRT::Node& node = data->node_;
//
//    std::vector<const DRT::Element*> eles;
//    addPBCelements(&node,eles);
//    const int numeles=eles.size();
//
//    for (size_t i=0;i<newColVectors.size();i++)
//    {
//      velnpDeriv1[i].Clear();
//      presnpDeriv1[i].Clear();
//    }
//
//    for (int iele=0;iele<numeles;iele++)
//    {
//      const DRT::Element* currele = eles[iele]; // current element
//
//      switch (currele->Shape())
//      {
//      case DRT::Element::hex8:
//      {
//        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
//        computeNodalGradient<numnode,DRT::Element::hex8>(newColVectors,newdofcolmap,newNodalDofColDistrib,currele,&node,velnpDeriv1Tmp,presnpDeriv1Tmp);
//      }
//      break;
//      case DRT::Element::hex20:
//      {
//        const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
//        computeNodalGradient<numnode,DRT::Element::hex20>(newColVectors,newdofcolmap,newNodalDofColDistrib,currele,&node,velnpDeriv1Tmp,presnpDeriv1Tmp);
//      }
//      break;
//      default:
//        dserror("xfem assembly type not yet implemented in time integration");
//      };
//
//      for (size_t i=0;i<newColVectors.size();i++)
//      {
//        velnpDeriv1[i]+=velnpDeriv1Tmp[i];
//        presnpDeriv1[i]+=presnpDeriv1Tmp[i];
//      }
//    } // end loop over elements around node
//
//    // set transport velocity at this node
//    const int gid = node.Id();
//    const std::set<XFEM::FieldEnr>& fieldenrset(newdofman_->getNodeDofSet(gid));
//    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldenrset.begin();
//        fieldenr != fieldenrset.end();++fieldenr)
//    {
//      const DofKey newdofkey(gid, *fieldenr);
//
//      if (fieldenr->getEnrichment().Type() == XFEM::Enrichment::typeStandard)
//      {
//        if (fieldenr->getField() == XFEM::PHYSICS::Velx)
//        {
//          const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
//          data->vel_(0) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
//        }
//        else if (fieldenr->getField() == XFEM::PHYSICS::Vely)
//        {
//          const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
//          data->vel_(1) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
//        }
//        else if (fieldenr->getField() == XFEM::PHYSICS::Velz)
//        {
//          const int newdofpos = newNodalDofColDistrib.find(newdofkey)->second;
//          data->vel_(2) = (*newColVectors[0])[newdofcolmap.LID(newdofpos)];
//        }
//      }
//    } // end loop over fieldenr
//
//    for (size_t i=0;i<newColVectors.size();i++)
//    {
//      velnpDeriv1[i].Scale(1.0/numeles);
//      presnpDeriv1[i].Scale(1.0/numeles);
//    }
//
//    data->velDeriv_ = velnpDeriv1;
//    data->presDeriv_ = presnpDeriv1;
//    //    cout << "after setting transportvel is " << data->vel_ << ", velderiv is " << velnpDeriv1[0]
//    //         << " and presderiv is " << presnpDeriv1[0] << endl;
//  }
}



/*------------------------------------------------------------------------------------------------*
 * reinitialize data for new computation                                         winklmaier 11/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::reinitializeData()
{
  dserror("adapt implementation of this function");
  dserror("adapt, how to get nds_np?");

//  int nds_np = -1;
//
//  cout << "in SemiLagrange::reinitializeData" << endl;
//  const int nsd = 3; // dimension
//  LINALG::Matrix<nsd,1> dummyStartpoint; // dummy startpoint for comparison
//  for (int i=0;i<nsd;i++) dummyStartpoint(i) = 777.777;
//
//  // fill curr_ structure with the data for the nodes which changed interface side
//  for (int lnodeid=0; lnodeid<discret_->NumMyColNodes(); lnodeid++)  // loop over processor nodes
//  {
//    DRT::Node* currnode = discret_->lColNode(lnodeid);
//    // node on current processor which changed interface side
//    if ((currnode->Owner() == myrank_) &&
//        (interfaceSideCompare((*phinp_)[lnodeid],(*phinpi_)[lnodeid])==false))
//    {
//      if (interfaceSideCompare((*phinp_)[lnodeid],(*phin_)[lnodeid]) == false) // real new side
//        timeIntData_->push_back(TimeIntData(
//            *currnode,
//            nds_np,
//            LINALG::Matrix<nsd,1>(true),
//            std::vector<LINALG::Matrix<nsd,nsd> >(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true)),
//            std::vector<LINALG::Matrix<1,nsd> >(oldVectors_.size(),LINALG::Matrix<1,nsd>(true)),
//            dummyStartpoint,
////            (*phinp_)[lnodeid],
//            1,
//            0,
//            std::vector<int>(1,-1),
//            std::vector<int>(1,-1),
//            INFINITY,
//            TimeIntData::predictor_));
//      else // other side than last FSI, but same side as old solution at last time step
//      {
//        for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
//            data!=timeIntData_->end(); data++)
//        {
//          const int nodeid = currnode->Id();
//
//          // 1) delete data
//          if (data->node_.Id()==nodeid)
//            timeIntData_->erase(data);
//
//          // 2) reset value of old solution
//          // get nodal velocities and pressures with help of the field set of node
//          const std::set<XFEM::FieldEnr>& fieldEnrSet(newdofman_->getNodeDofSet(nodeid));
//          for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fieldEnrSet.begin();
//              fieldenr != fieldEnrSet.end();++fieldenr)
//          {
//            const DofKey dofkey(nodeid, *fieldenr);
//            const int newdofpos = newNodalDofRowDistrib_.find(dofkey)->second;
//            const int olddofpos = oldNodalDofColDistrib_.find(dofkey)->second;
//            switch (fieldenr->getEnrichment().Type())
//            {
//            case XFEM::Enrichment::typeJump :
//            case XFEM::Enrichment::typeKink : break; // just standard dofs
//            case XFEM::Enrichment::typeStandard :
//            case XFEM::Enrichment::typeVoid :
//            {
//              for (size_t index=0;index<newVectors_.size();index++) // reset standard dofs due to old solution
//                (*newVectors_[index])[newdofrowmap_.LID(newdofpos)] =
//                    (*oldVectors_[index])[olddofcolmap_.LID(olddofpos)];
//              break;
//            }
//            case XFEM::Enrichment::typeUndefined : break;
//            default :
//            {
//              cout << fieldenr->getEnrichment().enrTypeToString(fieldenr->getEnrichment().Type()) << endl;
//              dserror("unknown enrichment type");
//              break;
//            }
//            } // end switch enrichment
//          } // end loop over fieldenr
//        } // end loop over nodes
//      }
//    }
//  } // end loop over processor nodes
//
//  startpoints();
//
//  // test loop if all initial startpoints have been computed
//  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
//      data!=timeIntData_->end(); data++)
//  {
//    if (data->startpoint_==dummyStartpoint)
//      dserror("WARNING! No enriched node on one interface side found!\nThis "
//          "indicates that the whole area is at one side of the interface!");
//  } // end loop over nodes
} // end function reinitializeData


/*------------------------------------------------------------------------------------------------*
 * call back-tracking of data at final Lagrangian origin of a point                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::callBackTracking(
    DRT::Element*&        ele,                 /// pointer to element
    TimeIntData*          data,                /// data
    LINALG::Matrix<3,1>&  xi,                  /// local coordinates
    const char*           backTrackingType     /// type of backTracking
)
{

  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    backTracking<numnode,DRT::Element::hex8>(ele,data,xi,backTrackingType);
  }
  break;
  case DRT::Element::hex20:
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex20>::numNodePerElement;
    backTracking<numnode,DRT::Element::hex20>(ele,data,xi,backTrackingType);
  }
  break;
  default:
    dserror("xfem assembly type not yet implemented in time integration"); break;
  };
} // end backTracking


/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final Lagrangian origin of a point                       schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::backTracking(
    DRT::Element*&        fittingele,          /// pointer to element
    TimeIntData*          data,                /// data
    LINALG::Matrix<3,1>&  xi,                  /// local coordinates
    const char*           backTrackingType     /// type of backTracking
)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n==============================================";
  IO::cout << "\n BACK-TRACKING on proc " << myrank_;
  IO::cout << "\n==============================================" << IO::endl;
#endif


  const int nsd = 3; // dimension

  if ((strcmp(backTrackingType,static_cast<const char*>("standard"))!=0) and
      (strcmp(backTrackingType,static_cast<const char*>("failing"))!=0))
    dserror("backTrackingType not implemented");

#ifdef DEBUG_SEMILAGRANGE
  if (strcmp(backTrackingType,static_cast<const char*>("standard"))==0)
  {
        cout << "\n--------------------------------------------------\n"
             << "\nnode: " << data->node_
             << "\ncomputed LAGRANGEAN ORIGIN  (startpoint) "
             << data->startpoint_
             << "with xi-coord. " << xi
             << "in element " << *fittingele
             << "\n--------------------------------------------------"<< endl;
  }
  if(strcmp(backTrackingType,static_cast<const char*>("failing"))==0)
  {
        cout << "\n--------------------------------------------------\n"
             << "\nnode: " << data->node_
             << "\nused <<<PSEUDO>>> LAGRANGEAN ORIGIN (initialpoint) "
             << data->initialpoint_
             << "with xi-coord. " << xi
             << "in element " << *fittingele
             << "\n--------------------------------------------------"<< endl;
  }
#endif

  //---------------------------------------------------------------------------------
  // Initialization

  LINALG::Matrix<3,1> lagrangeanOrigin(true);    // the applied Lagrangean origin (the real computed or an approximated)

  if(strcmp(backTrackingType,static_cast<const char*>("standard"))==0)
    lagrangeanOrigin = data->startpoint_; // use the computed start point approximation or the projected start point
  else if(strcmp(backTrackingType,static_cast<const char*>("failing"))==0)
    lagrangeanOrigin = data->initialpoint_; // use the initial guess for the Lagrangean origin
  else dserror("backTrackingType not implemented");


  LINALG::Matrix<numnode,1> shapeFcn(true);      // shape function
  LINALG::Matrix<3,numnode> shapeFcnDeriv(true); // shape function derivatives w.r.t xyz
  LINALG::Matrix<nsd,nsd> xji(true);             // invers of jacobian

  double deltaT = 0; // pseudo time-step size, used when the initial point is used instead of the computed lagrangean startpoint

  // data for the final back-tracking
  LINALG::Matrix<nsd,1> vel(true);                                                                // velocity data
  std::vector<LINALG::Matrix<nsd,nsd> > velnDeriv1(oldVectors_.size(),LINALG::Matrix<nsd,nsd>(true));  // first derivation of velocity data

  LINALG::Matrix<1,1> pres(true);                                                                 // pressure data
  std::vector<LINALG::Matrix<1,nsd> > presnDeriv1(oldVectors_.size(),LINALG::Matrix<1,nsd>(true));     // first derivation of pressure data

  std::vector<LINALG::Matrix<nsd,1> > veln(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));            // velocity at t^n
  LINALG::Matrix<nsd,1> transportVeln(true);                                                      // transport velocity at Lagrangian origin (x_Lagr(t^n))


  //---------------------------------------------------------------------------------
  // fill velocity and pressure data at nodes of element ...

  // node velocities of the element nodes for transport velocity
  LINALG::Matrix<nsd,numnode> nodevel(true);
  LINALG::Matrix<numnode,1>   nodepre(true);

  // node velocities of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<nsd,numnode> > nodeveldata(oldVectors_.size(),LINALG::Matrix<nsd,numnode>(true));
  // node pressures of the element nodes for the data that should be changed
  std::vector<LINALG::Matrix<numnode,1> > nodepresdata(oldVectors_.size(),LINALG::Matrix<numnode,1>(true));

  // velocity of the data that shall be changed
  std::vector<LINALG::Matrix<nsd,1> > velValues(oldVectors_.size(),LINALG::Matrix<nsd,1>(true));
  // pressures of the data that shall be changed
  std::vector<double> presValues(oldVectors_.size(),0);


  for (size_t index=0;index<oldVectors_.size();index++)
  {
    nodeveldata[index].Clear();
    nodepresdata[index].Clear();
  }

  DRT::Element* ele = fittingele; // current element

  //---------------------------------------------------------------------------------
  // get shape functions and derivatives at local coordinates

  bool compute_deriv = true;

  evalShapeAndDeriv<numnode,DISTYPE>(
      ele,
      xi,
      xji,
      shapeFcn,
      shapeFcnDeriv,
      compute_deriv
  );

  //-------------------------------------------------------
  // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
  std::vector<int> lm;

  for(int inode=0; inode< numnode; inode++)
  {
    DRT::Node* node = ele->Nodes()[inode];
    std::vector<int> dofs;
    dofset_old_->Dof(*node, data->nds_[inode], dofs );

    int size = dofs.size();

    for (int j=0; j< size; ++j)
    {
      lm.push_back(dofs[j]);
    }
  }

  //-------------------------------------------------------
  // all vectors are based on the same map

  extractNodalValuesFromVector<numnode,DISTYPE>(nodevel,nodepre, veln_,lm);

  for (size_t index=0;index<oldVectors_.size();index++)
    extractNodalValuesFromVector<numnode,DISTYPE>(nodeveldata[index],nodepresdata[index],oldVectors_[index],lm);


  //---------------------------------------------------------------------------------
  // interpolate velocity and pressure values at starting point
  transportVeln.Multiply(nodevel, shapeFcn);

#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\t transportVeln\t" <<  transportVeln << IO::endl;
#endif

  //---------------------------------------------------------------------------------
  // computing pseudo time-step deltaT
  // remark: if x is the Lagrange-origin of node, deltaT = dt with respect to small errors.
  // if it is not, deltaT estimates the time x needs to move to node)
  if (data->type_==TimeIntData::predictor_)
  {
    LINALG::Matrix<nsd,1> diff(data->node_.X());
    diff -=  lagrangeanOrigin; // diff = x_Node - x_Appr

    double numerator   = transportVeln.Dot(diff);          // numerator = v^T*(x_Node-x_Appr)
    double denominator = transportVeln.Dot(transportVeln); // denominator = v^T*v

    if (denominator>1e-15) deltaT = numerator/denominator; // else deltaT = 0 as initialized

#ifdef DEBUG_SEMILAGRANGE
    IO::cout << " \t recomputed modified pseudo time-step size: " << deltaT << IO::endl;
#endif
  }
  else
    deltaT = dt_;


  // interpolate velocity and pressure gradients for all fields at starting point and get final values
  for (size_t index=0;index<oldVectors_.size();index++)
  {
    veln[index].Multiply(nodeveldata[index],shapeFcn);

    velnDeriv1[index].MultiplyNT(1.0,nodeveldata[index],shapeFcnDeriv,1.0);
    presnDeriv1[index].MultiplyTT(1.0,nodepresdata[index],shapeFcnDeriv,1.0);
  } // end loop over vectors to be read from


  for (size_t index=0;index<oldVectors_.size();index++)
  {

    vel.Multiply(1.0-Theta(data),velnDeriv1[index],transportVeln);    // v = (1-theta)*Dv^n/Dx*v^n
    vel.Multiply(Theta(data),data->velDeriv_[index],data->vel_,1.0);  // v = theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n
    vel.Update(1.0,veln[index],deltaT);                               // v = v_n + dt*(theta*Dv^n+1/Dx*v^n+1+(1-theta)*Dv^n/Dx*v^n)
    velValues[index]=vel;

    pres.Multiply(1.0-Theta(data),presnDeriv1[index],transportVeln);   // p = (1-theta)*Dp^n/Dx*v^n
    pres.Multiply(Theta(data),data->presDeriv_[index],data->vel_,1.0); // p = theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n
    pres.MultiplyTN(1.0,nodepresdata[index],shapeFcn,deltaT);          // p = p_n + dt*(theta*Dp^n+1/Dx*v^n+1+(1-theta)*Dp^n/Dx*v^n)
    presValues[index] = pres(0);

#ifdef DEBUG_SEMILAGRANGE
    IO::cout << "\n***********************************************";
    IO::cout << "\n           RECONSTRUCTED VALUES for node " << (data->node_).Id();
    IO::cout << "\nvelocity entry in vector \t" << index << "\t " << vel;
    IO::cout << "pressure entry in vector \t" << index << "\t " << pres(0);
    IO::cout << "\n***********************************************" << IO::endl;
#endif
  } // loop over vectors to be set

  data->velValues_  = velValues;
  data->presValues_ = presValues;
  data->state_      = TimeIntData::doneStd_;


  return;
} // end backTracking


/*------------------------------------------------------------------------------------------------*
 * determine point's dofset in element ele w.r.t old or new interface position       schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::getNodalDofSet(
    DRT::Element*           ele,    /// pointer to element
    LINALG::Matrix<3,1>&    x ,     /// global coordinates of point
    std::vector<int>&       nds,    /// determine the points dofset w.r.t old/new interface position
    GEO::CUT::VolumeCell*&  vc,     /// valid fluid volumecell the point x lies in
    bool                    step_np /// computation w.r.t old or new interface position?
    )
{

  nds.clear();


#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t\t\t ... getNodalDofSet";
#endif


  RCP<XFEM::FluidWizard> wizard = step_np ? wizard_new_ : wizard_old_;

  GEO::CUT::ElementHandle* e = wizard->GetElement(ele);

  bool inside_structure = false;

  if ( e!=NULL ) // element in cut involved
  {
    GEO::CUT::plain_volumecell_set cells;
    e->VolumeCells(cells);

    if(cells.size() == 0) dserror("GEO::CUT::Element %d does not contain any volume cell", ele->Id());

    for(GEO::CUT::plain_volumecell_set::iterator cell_it=cells.begin(); cell_it!=cells.end(); cell_it++)
    {
      GEO::CUT::VolumeCell* cell = *cell_it;
//      if(cell->Contains(x))
      // cell contains the point inside or on one of its boundaries and the cell is an outside (fluid) cell
      if( ((cell->IsThisPointInside(x) == "inside") or (cell->IsThisPointInside(x) == "onBoundary") )
            and cell->Position() == GEO::CUT::Point::outside)
      {
#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t\t\t -> Position of point w.r.t volumecell is " << cell->IsThisPointInside(x)
                 << " \t cell pos = " << cell->Position() << IO::endl;
#endif
        nds = cell->NodalDofSet();

        vc = cell;

        IO::cout << "nds-vector " ;
        for(int i=0; i< (int)nds.size(); i++)
        {
          IO::cout << " " << nds[i];
        }
        IO::cout << "\n";

        return;
      }
      // point lies within the structure or on Boundary
      else if( (cell->IsThisPointInside(x) == "inside" or cell->IsThisPointInside(x) == "onBoundary")
               and (cell->Position() == GEO::CUT::Point::inside))
      {
#ifdef DEBUG_SEMILAGRANGE
        IO::cout << "\n\t\t\t -> Position of point w.r.t volumecell is " << cell->IsThisPointInside(x)
                 << " \t cell pos = " << cell->Position() << IO::endl;
#endif
        // do not return before all the other vcs have been tested, maybe a fluid-cell with onBoundary can be found
        inside_structure = true;
      }
      else
      {
        // the point does not lie inside this vc !
//#ifdef DEBUG_SEMILAGRANGE
//        IO::cout << "\n\t\t\t -> Position of cell " << cell->Position() << " and IsThisPointInside "<< cell->IsThisPointInside(x) << IO::endl;
//#endif
      }
    }

    // return if the structural volume cell is the only one which was found
    if(inside_structure)
    {
      nds.clear();
#ifdef DEBUG_SEMILAGRANGE
      IO::cout << "\n\t\t\t -> Position of point inside structure and not onBoundary of other fluid-vcs -> reset nds to empty vector" << IO::endl;
#endif

      return;
    }

    IO::cout << "error: coordinates of point x " << x << " number of volumecells: " << cells.size() << IO::endl;
    dserror("there is no volume cell in element %d which contains point with coordinates (%f,%f,%f) -> void element???", ele->Id(), x(0), x(1), x(2));

  }
  else // standard element, all its nodes have dofset 0
  {
    int numnode = ele->NumNode();

    for(int inode=0; inode <numnode; inode++)
    {
      nds.push_back(0);
    }
  }

  return;
}


/*------------------------------------------------------------------------------------------------*
 * back-tracking of data at final Lagrangian origin of a point                   winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::extractNodalValuesFromVector(
    LINALG::Matrix<3,numnode>& evel,
    LINALG::Matrix<numnode,1>& epre,
    RCP<Epetra_Vector> vel_vec,
    std::vector<int>& lm
    )
{
  const int nsd = 3;
  const int numdofpernode = nsd +1;

  evel.Clear();
  epre.Clear();

  if(vel_vec == Teuchos::null)
    dserror("vector is null");

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*vel_vec,mymatrix,lm);

  for (int inode=0; inode<numnode; ++inode)  // number of nodes
  {
    for(int idim=0; idim<nsd; ++idim) // number of dimensions
    {
      (evel)(idim,inode) = mymatrix[idim+(inode*numdofpernode)];
    }
    (epre)(inode,0) = mymatrix[nsd+(inode*numdofpernode)];
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 * compute Gradients at side-changing nodes                                      winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
template<const int numnode,DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_SemiLagrange::computeNodalGradient(
    std::vector<RCP<Epetra_Vector> >& newColVectors,
    const Epetra_Map& newdofcolmap,
    std::map<XFEM::DofKey,XFEM::DofGID>& newNodalDofColDistrib,
    const DRT::Element* ele,
    DRT::Node* node,
    std::vector<LINALG::Matrix<3,3> >& velnpDeriv1,
    std::vector<LINALG::Matrix<1,3> >& presnpDeriv1
) const
{ dserror("fix computeNodalGradient");
//  const int nsd = 3;
//
//  for (size_t i=0;i<newColVectors.size();i++)
//  {
//    velnpDeriv1[i].Clear();
//    presnpDeriv1[i].Clear();
//  }
//
//  // shape fcn data
//  LINALG::Matrix<nsd,1> xi(true);
//  LINALG::Matrix<nsd,2*numnode> enrShapeXYPresDeriv1(true);
//
//  // nodal data
//  LINALG::Matrix<nsd,2*numnode> nodevel(true);
//  LINALG::Matrix<1,2*numnode> nodepres(true);
//  LINALG::Matrix<nsd,1> nodecoords(true);
//
//  // dummies for function call
//  LINALG::Matrix<nsd,nsd> jacobiDummy(true);
//  LINALG::Matrix<numnode,1> shapeFcnDummy(true);
//  LINALG::Matrix<2*numnode,1> enrShapeFcnDummy(true);
//
//#ifdef COMBUST_NORMAL_ENRICHMENT
//  LINALG::Matrix<1,numnode> nodevelenr(true);
//  ApproxFuncNormalVector<2,2*numnode> shp;
//#else
//  LINALG::Matrix<nsd,2*numnode> enrShapeXYVelDeriv1(true);
//  LINALG::Matrix<2*numnode,1> enrShapeFcnVel(true);
//#endif
//
//  { // get local coordinates
//    bool elefound = false;
//    LINALG::Matrix<nsd,1> coords(node->X());
//    callXToXiCoords(ele,coords,xi,elefound);
//
//    if (!elefound) // possibly slave node looked for element of master node or vice versa
//    {
//      // get pbcnode
//      bool pbcnodefound = false; // boolean indicating whether this node is a pbc node
//      DRT::Node* pbcnode = NULL;
//      findPBCNode(node,pbcnode,pbcnodefound);
//
//      // get local coordinates
//      LINALG::Matrix<nsd,1> pbccoords(pbcnode->X());
//      callXToXiCoords(ele,pbccoords,xi,elefound);
//
//      if (!elefound) // now something is really wrong...
//        dserror("element of a row node not on same processor as node?! BUG!");
//    }
//  }
//
//#ifdef COMBUST_NORMAL_ENRICHMENT
//  pointdataXFEMNormal<numnode,DISTYPE>(
//      ele,
//#ifdef COLLAPSE_FLAME_NORMAL
//      coords,
//#endif
//      xi,
//      jacobiDummy,
//      shapeFcnDummy,
//      enrShapeFcnDummy,
//      enrShapeXYPresDeriv1,
//      shp,
//      true);
//#else
//  pointdataXFEM<numnode,DISTYPE>(
//      (DRT::Element*)ele,
//      xi,
//      jacobiDummy,
//      shapeFcnDummy,
//      enrShapeFcnVel,
//      enrShapeFcnDummy,
//      enrShapeXYVelDeriv1,
//      enrShapeXYPresDeriv1,
//      true);
//#endif
//
//  //cout << "shapefcnvel is " << enrShapeFcnVel << ", velderiv is " << enrShapeXYVelDeriv1 << " and presderiv is " << enrShapeXYPresDeriv1 << endl;
//  for (size_t i=0;i<newColVectors.size();i++)
//  {
//#ifdef COMBUST_NORMAL_ENRICHMENT
//    elementsNodalData<numnode>(
//        ele,
//        newColVectors[i],
//        dofman_,
//        newdofcolmap,
//        newNodalDofColDistrib,
//        nodevel,
//        nodevelenr,
//        nodepres);
//
//    const int* nodeids = currele->NodeIds();
//    size_t velncounter = 0;
//
//    // vderxy = enr_derxy(j,k)*evelnp(i,k);
//    for (int inode = 0; inode < numnode; ++inode)
//    {
//      // standard shape functions are identical for all vector components
//      // e.g. shp.velx.dx.s == shp.vely.dx.s == shp.velz.dx.s
//      velnpDeriv1[i](0,0) += nodevel(0,inode)*shp.velx.dx.s(inode);
//      velnpDeriv1[i](0,1) += nodevel(0,inode)*shp.velx.dy.s(inode);
//      velnpDeriv1[i](0,2) += nodevel(0,inode)*shp.velx.dz.s(inode);
//
//      velnpDeriv1[i](1,0) += nodevel(1,inode)*shp.vely.dx.s(inode);
//      velnpDeriv1[i](1,1) += nodevel(1,inode)*shp.vely.dy.s(inode);
//      velnpDeriv1[i](1,2) += nodevel(1,inode)*shp.vely.dz.s(inode);
//
//      velnpDeriv1[i](2,0) += nodevel(2,inode)*shp.velz.dx.s(inode);
//      velnpDeriv1[i](2,1) += nodevel(2,inode)*shp.velz.dy.s(inode);
//      velnpDeriv1[i](2,2) += nodevel(2,inode)*shp.velz.dz.s(inode);
//
//      const int gid = nodeids[inode];
//      const std::set<XFEM::FieldEnr>& enrfieldset = olddofman_->getNodeDofSet(gid);
//
//      for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
//          enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
//      {
//        if (enrfield->getField() == XFEM::PHYSICS::Veln)
//        {
//          velnpDeriv1[i](0,0) += nodevelenr(0,velncounter)*shp.velx.dx.n(velncounter);
//          velnpDeriv1[i](0,1) += nodevelenr(0,velncounter)*shp.velx.dy.n(velncounter);
//          velnpDeriv1[i](0,2) += nodevelenr(0,velncounter)*shp.velx.dz.n(velncounter);
//
//          velnpDeriv1[i](1,0) += nodevelenr(0,velncounter)*shp.vely.dx.n(velncounter);
//          velnpDeriv1[i](1,1) += nodevelenr(0,velncounter)*shp.vely.dy.n(velncounter);
//          velnpDeriv1[i](1,2) += nodevelenr(0,velncounter)*shp.vely.dz.n(velncounter);
//
//          velnpDeriv1[i](2,0) += nodevelenr(0,velncounter)*shp.velz.dx.n(velncounter);
//          velnpDeriv1[i](2,1) += nodevelenr(0,velncounter)*shp.velz.dy.n(velncounter);
//          velnpDeriv1[i](2,2) += nodevelenr(0,velncounter)*shp.velz.dz.n(velncounter);
//
//          velncounter += 1;
//        }
//      }
//    }
//
//#else
//    elementsNodalData<numnode>(
//        (DRT::Element*)ele,
//        newColVectors[i],
//        newdofman_,
//        newdofcolmap,
//        newNodalDofColDistrib,
//        nodevel,
//        nodepres);
//
//    velnpDeriv1[i].MultiplyNT(1.0,nodevel,enrShapeXYVelDeriv1,1.0);
//    presnpDeriv1[i].MultiplyNT(1.0,nodepres,enrShapeXYPresDeriv1,1.0);
//#endif
//  }
} // end function compute nodal gradient



/*------------------------------------------------------------------------------------------------*
 * get the time integration factor theta fitting to the computation type             schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
double XFEM::XFLUID_SemiLagrange::Theta(TimeIntData* data) const
{
  double theta = -1.0;

  switch (data->type_)
  {
  case TimeIntData::predictor_: theta = 0.0;            break;
  case TimeIntData::standard_ : theta = theta_default_; break;
  default: dserror("type not implemented"); break;
  }

  if (theta < 0.0) dserror("something wrong");

  return theta;
} // end function theta



/*------------------------------------------------------------------------------------------------*
 * export alternative algo data to neighbour proc                                    schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::exportAlternativAlgoData()
{

  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  std::vector<std::vector<TimeIntData> > dataVec(numproc_);

  // fill vectors with the data
  for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
      data!=timeIntData_->end(); data++)
  {
    if (data->state_==TimeIntData::failedSL_)
    {
      dataVec[data->initial_ele_owner_].push_back(*data);
    }
  }

  clearState(TimeIntData::failedSL_);
  timeIntData_->insert(timeIntData_->end(),
      dataVec[myrank_].begin(),
      dataVec[myrank_].end());

  dataVec[myrank_].clear(); // clear the set data from the vector

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest higher neighbour...)
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    DRT::PackBuffer dataSend; // vector including all data that has to be send to dest proc

    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;

    // pack data to be sent
    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      if (data->state_==TimeIntData::failedSL_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data=dataVec[dest].begin();
        data!=dataVec[dest].end(); data++)
    {
      if (data->state_==TimeIntData::failedSL_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    // clear the no more needed data
    dataVec[dest].clear();

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0,(double*)coords,0);
      int nds_np;
      LINALG::Matrix<nsd,1> vel;
      std::vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      std::vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> initialpoint;
      int initial_eid;
      int initial_ele_owner;
      int newtype;

      unpackNode(posinData,dataRecv,node);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,nds_np);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,initialpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,initial_eid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,initial_ele_owner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);

      timeIntData_->push_back(TimeIntData(
          node,
          nds_np,
          vel,
          velDeriv,
          presDeriv,
          initialpoint,
          initial_eid,
          initial_ele_owner,
          (TimeIntData::type)newtype)); // startOwner is current proc
    } // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end loop over processors
} // end exportAlternativAlgoData



/*------------------------------------------------------------------------------------------------*
 * export data while Newton loop to neighbour proc                                   schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_SemiLagrange::exportIterData(
    bool& procDone
)
{
#ifdef DEBUG_SEMILAGRANGE
  IO::cout << "\n\t=============================";
  IO::cout << "\n\t  export Iteration Data  ";
  IO::cout << "\n\t=============================" << IO::endl;
#endif

  const int nsd = 3; // 3 dimensions for a 3d fluid element

  // Initialization
  int dest = myrank_+1; // destination proc (the "next" one)
  if(myrank_ == (numproc_-1))
    dest = 0;

  int source = myrank_-1; // source proc (the "last" one)
  if(myrank_ == 0)
    source = numproc_-1;


  /*-------------------------------------------*
   * first part: send procfinished in order to *
   * check whether all procs have finished     *
   *-------------------------------------------*/
  for (int iproc=0;iproc<numproc_-1;iproc++)
  {
    DRT::PackBuffer dataSend;

    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));
    dataSend.StartPacking();
    DRT::ParObject::AddtoPack(dataSend,static_cast<int>(procDone));

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
    size_t posinData = 0;
    int allProcsDone;

    //unpack received data
    DRT::ParObject::ExtractfromPack(posinData,dataRecv,allProcsDone);

    if (allProcsDone==0)
      procDone = 0;

    // processors wait for each other
    discret_->Comm().Barrier();
  }

  /*--------------------------------------*
   * second part: if not all procs have   *
   * finished send data to neighbour proc *
   *--------------------------------------*/
  if (!procDone)
  {
    DRT::PackBuffer dataSend;

    // fill vectors with the data
    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::nextSL_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend,data->counter_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data=timeIntData_->begin();
        data!=timeIntData_->end(); data++)
    {
      if (data->state_==TimeIntData::nextSL_)
      {
        packNode(dataSend,data->node_);
        DRT::ParObject::AddtoPack(dataSend,data->nds_np_);
        DRT::ParObject::AddtoPack(dataSend,data->vel_);
        DRT::ParObject::AddtoPack(dataSend,data->velDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->presDeriv_);
        DRT::ParObject::AddtoPack(dataSend,data->initialpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_eid_);
        DRT::ParObject::AddtoPack(dataSend,data->initial_ele_owner_);
        DRT::ParObject::AddtoPack(dataSend,data->startpoint_);
        DRT::ParObject::AddtoPack(dataSend,data->searchedProcs_);
        DRT::ParObject::AddtoPack(dataSend,data->counter_);
        DRT::ParObject::AddtoPack(dataSend,(int)data->type_);
      }
    }

    clearState(TimeIntData::nextSL_);

    std::vector<char> dataRecv;
    sendData(dataSend,dest,source,dataRecv);

    // pointer to current position of group of cells in global string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      double coords[nsd] = {0.0};
      DRT::Node node(0,(double*)coords,0);
      int nds_np;
      LINALG::Matrix<nsd,1> vel;
      std::vector<LINALG::Matrix<nsd,nsd> > velDeriv;
      std::vector<LINALG::Matrix<1,nsd> > presDeriv;
      LINALG::Matrix<nsd,1> initialpoint;
      int initial_eid;
      int initial_ele_owner;
      LINALG::Matrix<nsd,1> startpoint;
      int searchedProcs;
      int iter;
      int newtype;

      unpackNode(posinData,dataRecv,node);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,nds_np);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,vel);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,velDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,presDeriv);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,initialpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,initial_eid);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,initial_ele_owner);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,startpoint);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,searchedProcs);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,iter);
      DRT::ParObject::ExtractfromPack(posinData,dataRecv,newtype);

      timeIntData_->push_back(TimeIntData(
          node,
          nds_np,
          vel,
          velDeriv,
          presDeriv,
          initialpoint,
          initial_eid,
          initial_ele_owner,
          startpoint,
          searchedProcs,
          iter,
          (TimeIntData::type)newtype));
    } // end loop over number of points to get

    // processors wait for each other
    discret_->Comm().Barrier();
  } // end if procfinished == false
} // end exportIterData


