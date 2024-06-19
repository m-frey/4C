/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material for SPH fluid

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "4C_mat_particle_sph_fluid.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHFluidType Mat::ParticleMaterialSPHFluidType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      ParticleMaterialBase(matdata),
      ParticleMaterialThermo(matdata),
      refDensFac_(matdata.parameters.get<double>("REFDENSFAC")),
      exponent_(matdata.parameters.get<double>("EXPONENT")),
      backgroundPressure_(matdata.parameters.get<double>("BACKGROUNDPRESSURE")),
      bulkModulus_(matdata.parameters.get<double>("BULK_MODULUS")),
      dynamicViscosity_(matdata.parameters.get<double>("DYNAMIC_VISCOSITY")),
      bulkViscosity_(matdata.parameters.get<double>("BULK_VISCOSITY")),
      artificialViscosity_(matdata.parameters.get<double>("ARTIFICIAL_VISCOSITY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::ParticleMaterialSPHFluid::create_material()
{
  return Teuchos::rcp(new Mat::ParticleMaterialSPHFluid(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::ParticleMaterialSPHFluidType::Create(
    const std::vector<char>& data)
{
  Mat::ParticleMaterialSPHFluid* particlematsph = new Mat::ParticleMaterialSPHFluid();
  particlematsph->unpack(data);
  return particlematsph;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid(Mat::PAR::ParticleMaterialSPHFluid* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialSPHFluid::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*---------------------------------------------------------------------------*
 | unpack                                                     sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialSPHFluid::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      // note: dynamic_cast needed due diamond inheritance structure
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<Mat::PAR::ParticleMaterialSPHFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
