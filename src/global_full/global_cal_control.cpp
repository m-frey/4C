/*!----------------------------------------------------------------------
\file global_cal_control.cpp
\brief routine to control execution phase

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_structure/stru_dyn_nln_drt.H"
#include "../drt_fluid/fluid_dyn_nln_drt.H"
#include "../drt_scatra/scatra_dyn.H"
#include "../drt_ale/ale_dyn.H"
#include "../drt_fsi/fsi_dyn.H"
#include "../drt_fs3i/fs3i_dyn.H"
#include "../drt_loma/loma_dyn.H"
#include "../drt_elch/elch_dyn.H"
#include "../drt_combust/combust_dyn.H"
#include "../drt_opti/topopt_dyn.H"
#include "../drt_thermo/thr_dyn.H"
#include "../drt_tsi/tsi_dyn.H"
#include "../drt_art_net/art_net_dyn_drt.H"
#include "../drt_red_airways/red_airways_dyn_drt.H"
#include "../drt_stru_ale/stru_ale_dyn.H"
#include "../drt_poroelast/poro_dyn.H"
#include "../drt_ssi/ssi_dyn.H"
#include "../drt_particle/particle_dyn.H"
#include "../drt_stru_multi/microstatic_npsupport.H"


/*----------------------------------------------------------------------*
 |  routine to control execution phase                   m.gee 6/01     |
 *----------------------------------------------------------------------*/
void ntacal()
{
  int restart = DRT::Problem::Instance()->Restart();

  // choose the entry-routine depending on the problem type
  switch (DRT::Problem::Instance()->ProblemType())
  {
    case prb_structure:
      caldyn_drt();
      break;
    case prb_fluid:
      dyn_fluid_drt(restart);
      break;
    case prb_scatra:
      scatra_dyn(restart);
      break;
    case prb_fluid_xfem2:
      fluid_xfem2_drt();
      break;
    case prb_fluid_fluid_ale:
      fluid_fluid_ale_drt();
      break;
    case prb_fluid_fluid_fsi:
      fluid_fluid_fsi_drt();
    break;
    case prb_fluid_fluid:
      fluid_fluid_drt(restart);
      break;
    case prb_fluid_ale:
      fluid_ale_drt();
      break;
    case prb_freesurf:
      fluid_freesurf_drt();
      break;

    case prb_fsi:
    case prb_fsi_lung:
      fsi_ale_drt();
      break;
    case prb_fsi_xfem:
      xfsi_drt();
      break;

    case prb_gas_fsi:
    case prb_biofilm_fsi:
    case prb_thermo_fsi:
    case prb_tfsi_aero:
      fs3i_dyn();
      break;

    case prb_ale:
      dyn_ale_drt();

    case prb_thermo:
      thr_dyn_drt();
      break;

    case prb_tsi:
      tsi_dyn_drt();
      break;

    case prb_loma:
      loma_dyn(restart);
      break;

    case prb_elch:
      elch_dyn(restart);
      break;

    case prb_combust:
      combust_dyn();
      break;

    case prb_fluid_topopt:
      fluid_topopt_dyn();
      break;

    case prb_art_net:
      dyn_art_net_drt();
      break;

    case prb_red_airways:
      dyn_red_airways_drt();
      break;

    case prb_struct_ale:
      stru_ale_dyn_drt(restart);
      break;

    case prb_poroelast:
      poroelast_drt();
      break;
    case prb_poroscatra:
      poro_scatra_drt();
      break;
    case prb_ssi:
      ssi_drt();
      break;
    case prb_redairways_tissue:
      redairway_tissue_dyn();
      break;
    case prb_particle:
      particle_drt();
      break;

    case prb_np_support:
      STRUMULTI::np_support_drt();
      break;

    default:
      dserror("solution of unknown problemtyp %d requested", DRT::Problem::Instance()->ProblemType());
      break;
  }

}
