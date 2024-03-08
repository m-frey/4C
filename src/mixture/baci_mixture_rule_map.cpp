/*----------------------------------------------------------------------*/
/*! \file

\brief Mixture rule for homogenized constrained mixtures with mass fractions defined as discrete
values per element

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_rule_map.hpp"

#include "baci_global_data.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_mat_par_material.hpp"
#include "baci_mixture_constituent.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_ConfigDefs.h>
#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>


BACI_NAMESPACE_OPEN

namespace
{
  std::vector<double> GetValidateMassFractions(
      const std::unordered_map<int, std::vector<double>>& mass_fractions_map, const int ele_id_key,
      const std::size_t num_constituents)
  {
    auto it = mass_fractions_map.find(ele_id_key);
    if (it == mass_fractions_map.end())
    {
      dserror("Element id %d not found in the mass fraction map supplied by csv file.", ele_id_key);
    }

    if (it->second.size() != num_constituents)
    {
      dserror(
          "Number of mass fractions for element id %d does not match the number of constituents "
          "%d.",
          ele_id_key, num_constituents);
    }
    const std::vector<double> massfracs = it->second;

    // check, whether the mass frac sums up to 1
    double sum = 0.0;
    for (double massfrac : massfracs) sum += massfrac;
    if (std::abs(1.0 - sum) > 1e-8)
      dserror(
          "Mass fractions for element id %d don't sum up to 1, which is unphysical.", ele_id_key);

    return massfracs;
  }
}  // namespace

MIXTURE::PAR::MapMixtureRule::MapMixtureRule(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureRule(matdata),
      initial_reference_density_(*matdata->Get<double>("DENS")),
      num_constituents_(*matdata->Get<int>("NUMCONST")),
      mass_fractions_map_(
          *matdata->Get<std::unordered_map<int, std::vector<double>>>("MASSFRACMAPFILE")){};

std::unique_ptr<MIXTURE::MixtureRule> MIXTURE::PAR::MapMixtureRule::CreateRule()
{
  return std::make_unique<MIXTURE::MapMixtureRule>(this);
}

MIXTURE::MapMixtureRule::MapMixtureRule(MIXTURE::PAR::MapMixtureRule* params)
    : MixtureRule(params), params_(params)
{
}

void MIXTURE::MapMixtureRule::Setup(Teuchos::ParameterList& params, const int eleGID)
{
  MixtureRule::Setup(params, eleGID);
}

void MIXTURE::MapMixtureRule::UnpackMixtureRule(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureRule::UnpackMixtureRule(position, data);
}

void MIXTURE::MapMixtureRule::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, const int gp,
    const int eleGID)
{
  // define temporary matrices
  static CORE::LINALG::Matrix<6, 1> cstress;
  static CORE::LINALG::Matrix<6, 6> ccmat;

  // evaluate the mass fractions at the given element id (one based entires in the csv file)
  auto massfracs =
      GetValidateMassFractions(params_->mass_fractions_map_, eleGID + 1, Constituents().size());

  // Iterate over all constituents and add all stress/cmat contributions
  for (std::size_t i = 0; i < Constituents().size(); ++i)
  {
    MixtureConstituent& constituent = *Constituents()[i];
    cstress.Clear();
    ccmat.Clear();
    constituent.Evaluate(F, E_strain, params, cstress, ccmat, gp, eleGID);

    // add stress contribution to global stress
    double constituent_density = params_->initial_reference_density_ * massfracs[i];
    S_stress.Update(constituent_density, cstress, 1.0);
    cmat.Update(constituent_density, ccmat, 1.0);
  }
}


BACI_NAMESPACE_CLOSE