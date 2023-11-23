/*----------------------------------------------------------------------------*/
/*! \file

\brief One beam contact segment living on an element pair

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "baci_beamcontact_beam3contactvariables.H"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  constructor (public)                                     meier 01/14|
 *----------------------------------------------------------------------*/
template <const int numnodes, const int numnodalvalues>
CONTACT::Beam3contactvariables<numnodes, numnodalvalues>::Beam3contactvariables(
    std::pair<TYPE, TYPE>& closestpoint, std::pair<int, int>& segids, std::pair<int, int>& intids,
    const double& pp, TYPE jacobi)
    : closestpoint_(closestpoint),
      segids_(segids),
      intids_(intids),
      jacobi_(jacobi),
      gap_(0.0),
      normal_(CORE::LINALG::Matrix<3, 1, TYPE>(true)),
      pp_(pp),
      ppfac_(0.0),
      dppfac_(0.0),
      fp_(0.0),
      dfp_(0.0),
      energy_(0.0),
      integratedenergy_(0.0),
      angle_(0.0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  end: constructor
 *----------------------------------------------------------------------*/

// Possible template cases: this is necessary for the compiler
template class CONTACT::Beam3contactvariables<2, 1>;
template class CONTACT::Beam3contactvariables<3, 1>;
template class CONTACT::Beam3contactvariables<4, 1>;
template class CONTACT::Beam3contactvariables<5, 1>;
template class CONTACT::Beam3contactvariables<2, 2>;
