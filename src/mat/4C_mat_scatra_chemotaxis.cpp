/*----------------------------------------------------------------------*/
/*! \file
 \brief
This file contains the base material for chemotactic scalars.

\level 3

*----------------------------------------------------------------------*/

#include "4C_mat_scatra_chemotaxis.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ScatraChemotaxisMat::ScatraChemotaxisMat(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      numscal_(matdata.parameters.get<int>("NUMSCAL")),
      pair_(matdata.parameters.get<std::vector<int>>("PAIR")),
      chemocoeff_(matdata.parameters.get<double>("CHEMOCOEFF"))
{
  // Some checks for more safety
  if (numscal_ != (int)pair_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d", numscal_,
        pair_.size());

  // is there exactly one '1' (i.e. attractant) and at least one '-1' (i.e. chemotractant)?
  int numpos = 0;
  int numneg = 0;
  for (int i = 0; i < numscal_; i++)
  {
    if (pair_.at(i) > 1e-10)
      numpos++;
    else if (pair_.at(i) < -1e-10)
      numneg++;
  }
  if (numpos != 1 or numneg != 1)
    FOUR_C_THROW(
        "Each PAIR vector must contain exactly one '-1' (i.e. chemotractant) and exactly one '1' "
        "(i.e. attractant)!");

  return;
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraChemotaxisMat::create_material()
{
  return Teuchos::rcp(new Mat::ScatraChemotaxisMat(this));
}


Mat::ScatraChemotaxisMatType Mat::ScatraChemotaxisMatType::instance_;


Core::Communication::ParObject* Mat::ScatraChemotaxisMatType::Create(const std::vector<char>& data)
{
  Mat::ScatraChemotaxisMat* scatra_chemotaxis_mat = new Mat::ScatraChemotaxisMat();
  scatra_chemotaxis_mat->unpack(data);
  return scatra_chemotaxis_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraChemotaxisMat::ScatraChemotaxisMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraChemotaxisMat::ScatraChemotaxisMat(Mat::PAR::ScatraChemotaxisMat* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraChemotaxisMat::pack(Core::Communication::PackBuffer& data) const
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraChemotaxisMat::unpack(const std::vector<char>& data)
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
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ScatraChemotaxisMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
