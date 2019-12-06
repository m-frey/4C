/*----------------------------------------------------------------------*/
/*! \file

\brief This holds the implementation of the non-abstract methods of the Mixture constituents
 interface

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "mixture_constituent.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "mixture_constituent_elasthyper.H"

// Constructor of the mixture constituent parameters
MIXTURE::PAR::MixtureConstituent::MixtureConstituent(
    const Teuchos::RCP<MAT::PAR::Material>& matdata, const double ref_mass_fraction)
    : MAT::PAR::Parameter(matdata), ref_mass_fraction_(ref_mass_fraction)
{
}

// Create instance of the constituent from the parameters
Teuchos::RCP<MAT::Material> MIXTURE::PAR::MixtureConstituent::CreateMaterial()
{
  dserror("Cannot create mixture constituent from this method. Use CreateConstituent() instead.");
  return Teuchos::null;
}

// Create the parameters of the constituents from the material number and the reference mass
// fraction
MIXTURE::PAR::MixtureConstituent* MIXTURE::PAR::MixtureConstituent::Factory(
    int matnum, const double ref_mass_fraction)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
  {
    dserror("List of materials cannot be accessed in the global problem instance.");
  }

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
  {
    dserror("List of materials in the global problem instance is empty.");
  }

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  switch (curmat->Type())
  {
    case INPAR::MAT::mix_elasthyper:
    {
      if (curmat->Parameter() == nullptr)
      {
        curmat->SetParameter(
            new MIXTURE::PAR::MixtureConstituent_ElastHyper(curmat, ref_mass_fraction));
      }
      auto* params =
          dynamic_cast<MIXTURE::PAR::MixtureConstituent_ElastHyper*>(curmat->Parameter());
      return params;
    }
    default:
      break;
  }
  dserror("The referenced material with id %d is not registered as a Mixture Constituent!", matnum);
  return nullptr;
}

// Empty constructor
MIXTURE::MixtureConstituent::MixtureConstituent()
    : initialReferenceDensity_(0.0), numgp_(0), has_read_element_(false), is_setup_(false)
{
}

//! Init is called once at the beginning to setup the number of GPs and the Parameter List
void MIXTURE::MixtureConstituent::ReadElement(const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Init must only be called once
  if (has_read_element_) dserror("ReadElement() is called multiple times. Just once allowed.");
  has_read_element_ = true;
  numgp_ = numgp;
}

// Setup of the mixture constituents and all its subparts
void MIXTURE::MixtureConstituent::Setup(Teuchos::ParameterList& params)
{
  // Setup must be called after Init()
  if (!has_read_element_) dserror("ReadElement() must be called before Setup()");

  // Setup must only be called once
  if (is_setup_) dserror("Setup() is called multiple times. Just once allowed.");
  is_setup_ = true;
}

// Pack everything for distribution to other processors
void MIXTURE::MixtureConstituent::PackConstituent(DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data, numgp_);
  DRT::ParObject::AddtoPack(data, static_cast<const int>(has_read_element_));
  DRT::ParObject::AddtoPack(data, static_cast<const int>(is_setup_));
}

// Unpack base constituent data, need to be called by every derived class
void MIXTURE::MixtureConstituent::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // make sure we have a pristine material
  has_read_element_ = false;
  numgp_ = 0;
  is_setup_ = false;

  DRT::ParObject::ExtractfromPack(position, data, numgp_);

  has_read_element_ = (bool)DRT::ParObject::ExtractInt(position, data);
  is_setup_ = (bool)DRT::ParObject::ExtractInt(position, data);
}