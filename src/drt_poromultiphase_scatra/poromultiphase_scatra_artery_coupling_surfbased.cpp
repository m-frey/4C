/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for surface-based (non-conforming) coupling between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_artery_coupling_surfbased.H"
#include "drt_globalproblem.H"
#include "linalg_utils_densematrix_communication.H"

#include "poromultiphase_scatra_artery_coupling_pair.H"

#include <Epetra_MultiVector.h>

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::PoroMultiPhaseScaTraArtCouplSurfBased(
    Teuchos::RCP<DRT::Discretization> arterydis, Teuchos::RCP<DRT::Discretization> contdis,
    const Teuchos::ParameterList& couplingparams, const std::string& condname,
    const std::string& artcoupleddofname, const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplNonConforming(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname)
{
  // user info
  if (myrank_ == 0)
  {
    std::cout << "<                                                  >" << std::endl;
    PrintOutCouplingMethod();
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "\n";
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::PreEvaluateCouplingPairs()
{
  const int numpatch_axi = DRT::Problem::Instance()
                               ->PoroFluidMultiPhaseDynamicParams()
                               .sublist("ARTERY COUPLING")
                               .get<int>("NUMPATCH_AXI");
  const int numpatch_rad = DRT::Problem::Instance()
                               ->PoroFluidMultiPhaseDynamicParams()
                               .sublist("ARTERY COUPLING")
                               .get<int>("NUMPATCH_RAD");
  const int numartele = arterydis_->NumGlobalElements();
  const int numgp_per_artele = numpatch_axi * numpatch_rad * 25;
  const int numgp_desired = numgp_per_artele * numartele;

  // this vector keeps track of evaluation of GPs
  Teuchos::RCP<Epetra_MultiVector> gp_vector =
      Teuchos::rcp(new Epetra_MultiVector(*arterydis_->ElementColMap(), numgp_per_artele));

  // pre-evaluate
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++) coupl_elepairs_[i]->PreEvaluate(gp_vector);

  // delete the inactive pairs
  coupl_elepairs_.erase(std::remove_if(coupl_elepairs_.begin(), coupl_elepairs_.end(), IsNotActive),
      coupl_elepairs_.end());

  // the following takes care of a very special case, namely, if a GP on the lateral surface lies
  // exactly in between two or more 3D elements owned by different procs.
  // In that case, the GP is duplicated across all owning procs.
  // We detect such cases and communicate them to all procs. below by adapting the
  // "gp_vector", to contain the multiplicity of the respective GP. Later the GP weight inside the
  // coupling pair is scaled by the inverse of the multiplicity

  int duplicates = 0;
  if (Comm().NumProc() > 1)
  {
    int mygpvec[numgp_per_artele];
    int sumgpvec[numgp_per_artele];
    std::fill(mygpvec, mygpvec + numgp_per_artele, 0);
    std::fill(sumgpvec, sumgpvec + numgp_per_artele, 0);
    // loop over all GIDs
    for (int gid = gp_vector->Map().MinAllGID(); gid <= gp_vector->Map().MaxAllGID(); gid++)
    {
      // reset
      std::fill(sumgpvec, sumgpvec + numgp_per_artele, 0);

      const int mylid = gp_vector->Map().LID(gid);
      // if not owned or ghosted fill with zeros
      if (mylid < 0) std::fill(mygpvec, mygpvec + numgp_per_artele, 0);
      // else get the GP vector
      else
        for (int igp = 0; igp < numgp_per_artele; igp++)
          mygpvec[igp] = static_cast<int>(((*gp_vector)[igp])[mylid]);

      // communicate to all via summation
      Comm().SumAll(mygpvec, sumgpvec, numgp_per_artele);

      // this is ok for now, either the GID does not exist or the entire element protrudes.
      // Inform user and continue
      if (*std::max_element(sumgpvec, sumgpvec + numgp_per_artele) < 1)
      {
        std::cout << "WARNING! No GP of element  " << gid + 1 << " could be projected!"
                  << std::endl;
        continue;
      }

      // if one entry is equal to zero, this GP could not be projected
      if (*std::min_element(sumgpvec, sumgpvec + numgp_per_artele) < 1)
        dserror("It seems as if one GP could not be projected");

      // find number of duplicates
      int sum = 0;
      sum = std::accumulate(sumgpvec, sumgpvec + numgp_per_artele, sum);
      duplicates += sum - numgp_per_artele;

      // if owned or ghosted by this proc. and if duplicates have been detected, replace entry in
      // gp_vector
      if (mylid >= 0 && sum > numgp_per_artele)
      {
        for (int igp = 0; igp < numgp_per_artele; igp++)
        {
          int err = gp_vector->ReplaceMyValue(mylid, igp, static_cast<double>(sumgpvec[igp]));
          if (err != 0) dserror("ReplaceMyValue failed with error code %d!", err);
        }
      }
    }
  }

  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
    coupl_elepairs_[i]->DeleteUnnecessaryGPs(gp_vector);

  int total_num_gp = 0;
  int numgp = 0;

  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    // segment ID not needed in this case, just set to zero
    coupl_elepairs_[i]->SetSegmentID(0);
    numgp = numgp + coupl_elepairs_[i]->NumGP();
  }
  // safety check
  Comm().SumAll(&numgp, &total_num_gp, 1);
  if (numgp_desired != total_num_gp - duplicates)
    dserror("It seems as if some GPs could not be projected");

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);
  if (contdis_->Name() == "porofluid" && myrank_ == 0)
    std::cout << "Only " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs are active" << std::endl;

  // print out summary of pairs
  if (contdis_->Name() == "porofluid" &&
      (DRT::INPUT::IntegralValue<int>(couplingparams_, "PRINT_OUT_SUMMARY_PAIRS")))
  {
    if (myrank_ == 0)
      std::cout << "In total " << numgp_desired << " GPs (" << numgp_per_artele
                << " per artery element) required for lateral surface coupling" << std::endl;
    std::cout << "Proc. " << myrank_ << " evaluates " << numgp << " GPs "
              << "(" << (double)(numgp) / (double)(total_num_gp)*100.0 << "% of all GPs)"
              << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::IsNotActive(
    const Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
        coupling_pair)
{
  return not coupling_pair->IsActive();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::Setup()
{
  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Setup();

  // error-checks
  if (has_varying_diam_) dserror("Varying diameter not yet possible for surface-based coupling");
  if (!evaluate_in_ref_config_)
    dserror("Evaluation in current configuration not yet possible for surface-based coupling");

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::Evaluate(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  if (!issetup_) dserror("Setup() has not been called");

  if (!porofluidmanagersset_)
  {
    // pre-evaluate the pairs --> has to be done here since radius inside the material is required
    PreEvaluateCouplingPairs();
  }

  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Evaluate(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::SetupSystem(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art,
    Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_art)
{
  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetupSystem(sysmat, rhs,
      sysmat_cont, sysmat_art, rhs_cont, rhs_art, dbcmap_cont, dbcmap_art->CondMap(),
      dbcmap_art->CondMap());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::ApplyMeshMovement()
{
  if (!evaluate_in_ref_config_)
    dserror("Evaluation in current configuration not possible for surface-based coupling");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::BloodVesselVolumeFraction()
{
  dserror("Output of vessel volume fraction not possible for surface-based coupling");

  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplSurfBased::PrintOutCouplingMethod() const
{
  std::cout << "<   surface-based formulation                      >" << std::endl;
  PoroMultiPhaseScaTraArtCouplNonConforming::PrintOutCouplingMethod();
}
