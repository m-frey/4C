/*!----------------------------------------------------------------------
\file so_hex20_multiscale.cpp
\brief

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#include "so_hex20.H"
#include "../drt_lib/drt_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/micromaterial.H"

using namespace std; // cout etc.


/*----------------------------------------------------------------------*
 |  homogenize material density (public)                                |
 *----------------------------------------------------------------------*/
// this routine is intended to determine a homogenized material
// density for multi-scale analyses by averaging over the initial volume

void DRT::ELEMENTS::So_hex20::soh20_homog(ParameterList&  params)
{
  double homogdens = 0.;
  const static std::vector<double> weights = soh20_weights();
  const double density = Material()->Density();

  for (int gp=0; gp<NUMGPT_SOH20; ++gp)
  {
    homogdens += detJ_[gp] * weights[gp] * density;
  }

  double homogdensity = params.get<double>("homogdens", 0.0);
  params.set("homogdens", homogdensity+homogdens);

  return;
}


/*----------------------------------------------------------------------*
 |  Read restart on the microscale                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_read_restart_multi()
{
  const int ele_ID = Id();
  RefCountPtr<MAT::Material> mat = Material();

  for (int gp=0; gp<NUMGPT_SOH20; ++gp)
  {

    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

    micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_readrestart");
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Initialization in case of inverse analyses                  lw 08/11|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_multi_invana_init()
{
  const int ele_ID = Id();
  RefCountPtr<MAT::Material> mat = Material();

  for (int gp=0; gp<NUMGPT_SOH20; ++gp)
  {
    MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

    micro->Evaluate(NULL, NULL, NULL, NULL, gp, ele_ID, 0., 0., "multi_invana_init");
  }
  return;
}

#endif
#endif
