/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch_NP.cpp

\brief evaluation of scatra elements for elch

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_elch_NP.H"
#include "scatra_ele_parameter_elch.H"


/*----------------------------------------------------------------------*
 * Add dummy mass matrix to sysmat
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::PrepMatAndRhsInitialTimeDerivative(
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra
)
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  /*----------------------------------------------------------------------*/
  // element integration loop
  /*----------------------------------------------------------------------*/
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // loop starts at k=numscal_ !!
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = fac*my::funct_(vi); // no density required here
      const int fvi = vi*my::numdofpernode_+my::numscal_;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+my::numscal_;

        elemat1_epetra(fvi,fui) += v*my::funct_(ui);
      }
    }
  }

  // set zero for the rhs of the potential
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+my::numscal_;

    elevec1_epetra[fvi] = 0.0; // zero out!
  }

  return;
}

/*----------------------------------------------------------------------*
 * Get Conductivity
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::GetConductivity(
  const enum INPAR::ELCH::EquPot    equpot,
  double&                           sigma_all,
  Epetra_SerialDenseVector&         sigma
)
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerElch> dme = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerElch>(my::diffmanager_);

  // calculate conductivity of electrolyte solution
  const double frt = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->FRT();
  const double factor = frt*INPAR::ELCH::faraday_const; // = F^2/RT

  // get concentration of transported scalar k at integration point
  std::vector<double> conint(my::numscal_);
  for (int k = 0;k<my::numscal_;++k)
    conint[k] = my::funct_.Dot(my::ephinp_[k]);

  // Dilute solution theory:
  // Conductivity is computed by
  // sigma = F^2/RT*Sum(z_k^2 D_k c_k)
  for(int k=0; k < my::numscal_; k++)
  {
    double sigma_k = factor*dme->GetValence(k)*dme->GetIsotropicDiff(k)*dme->GetValence(k)*conint[k];
    sigma[k] += sigma_k; // insert value for this ionic species
    sigma_all += sigma_k;

    // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
    if(equpot==INPAR::ELCH::equpot_enc_pde_elim)
    {
      sigma_all += factor*dme->GetIsotropicDiff(my::numscal_)*dme->GetValence(my::numscal_)*dme->GetValence(k)*(-conint[k]);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 * Calculate Mat and Rhs for electric potential field        ehrl 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CalcMatAndRhsElectricPotentialField(
  Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
  const enum INPAR::ELCH::EquPot    equpot,
  Epetra_SerialDenseMatrix&         emat,
  Epetra_SerialDenseVector&         erhs,
  const double                      fac,
  Teuchos::RCP<ScaTraEleDiffManagerElch>& dme
)
{
  // calculate conductivity of electrolyte solution
  const double frt = dynamic_cast<DRT::ELEMENTS::ScaTraEleParameterElch*>(my::scatrapara_)->FRT();

  double sigmaint(0.0);
  for (int k=0; k<my::numscal_; ++k)
  {
    double sigma_k = frt*dme->GetValence(k)*dme->GetIsotropicDiff(k)*dme->GetValence(k)*vm->ConInt(k);
    sigmaint += sigma_k;

    // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
    if(equpot==INPAR::ELCH::equpot_enc_pde_elim)
      sigmaint += frt*dme->GetValence(k)*dme->GetIsotropicDiff(my::numscal_)*dme->GetValence(my::numscal_)*(-vm->ConInt(k));

    // diffusive terms on rhs
    const double vrhs = fac*dme->GetIsotropicDiff(k)*dme->GetValence(k);
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+my::numscal_;
      double laplawf(0.0);
      my::GetLaplacianWeakFormRHS(laplawf,vm->GradPhi(k),vi);
      erhs[fvi] -= vrhs*laplawf;
      // effect of eliminated species c_m has to be added (c_m = - 1/z_m \sum_{k=1}^{m-1} z_k c_k)
      if(equpot==INPAR::ELCH::equpot_enc_pde_elim)
        erhs[fvi] -= -fac*dme->GetValence(k)*dme->GetIsotropicDiff(my::numscal_)*laplawf;
    }

    // provide something for conc. dofs: a standard mass matrix
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int    fvi = vi*my::numdofpernode_+k;
      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;
        emat(fvi,fui) += fac*my::funct_(vi)*my::funct_(ui);
      }
    }
  } // for k

  // ----------------------------------------matrix entries
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int    fvi = vi*my::numdofpernode_+my::numscal_;
    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+my::numscal_;
      double laplawf(0.0);
      my::GetLaplacianWeakForm(laplawf,ui,vi);
      emat(fvi,fui) += fac*sigmaint*laplawf;
    }

    double laplawf(0.0);
    my::GetLaplacianWeakFormRHS(laplawf,vm->GradPot(),vi);
    erhs[fvi] -= fac*sigmaint*laplawf;
  }

  return;
}

/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     gjb 06/08|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchNP<distype>::CalculateFlux(
    LINALG::Matrix<my::nsd_,1>&     q,      //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const int                       k,          //!< index of current scalar
    const double                    fac,         //!< integration factor
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,  //!< variable manager
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme                                 //!< diffusion manager
  )
{
  // dynamic cast to elch-specific diffusion manager
  Teuchos::RCP<ScaTraEleInternalVariableManagerElchNP <my::nsd_,my::nen_> > vmnp
    = Teuchos::rcp_dynamic_cast<ScaTraEleInternalVariableManagerElchNP <my::nsd_,my::nen_> >(vm);

  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // add different flux contributions as specified by user input
  switch (fluxtype)
  {
  case INPAR::SCATRA::flux_total_domain:
    // convective flux contribution
    q.Update(vmnp->ConInt(k),vmnp->ConVelInt());

    // no break statement here!
  case INPAR::SCATRA::flux_diffusive_domain:
    // diffusive flux contribution
    q.Update(-(dme->GetIsotropicDiff(k)),vmnp->GradPhi(k),1.0);

    q.Update(-myelch::elchpara_->FRT()*dme->GetIsotropicDiff(k)*dme->GetValence(k)*vmnp->ConInt(k),vmnp->GradPot(),1.0);

    break;
  default:
    dserror("received illegal flag inside flux evaluation for whole domain"); break;
  };

  return;
} // ScaTraCalc::CalculateFlux


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchNP<DRT::Element::nurbs27>;
