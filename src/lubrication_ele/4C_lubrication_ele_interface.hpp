/*--------------------------------------------------------------------------*/
/*! \file

\brief Interface of Lubrication element evaluation routine

\level 3

*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_LUBRICATION_ELE_INTERFACE_HPP
#define FOUR_C_LUBRICATION_ELE_INTERFACE_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_ParameterList.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::Elements
{
  class Element;
}

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    /// Interface base class for LubricationEleCalc
    /*!
      This class exists to provide a common interface for all template
      versions of LubricationEleCalc.
     */
    class LubricationEleInterface
    {
     public:
      /**
       * Virtual destructor.
       */
      virtual ~LubricationEleInterface() = default;

      /// Default constructor.
      LubricationEleInterface() = default;

      /// Setup element evaluation
      virtual int SetupCalc(CORE::Elements::Element* ele, DRT::Discretization& discretization) = 0;

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function; it
        must be defined in LubricationEleCalc.
        The Evaluate() method is meant only for the assembling of the
        linearized matrix and the right hand side
       */
      virtual int Evaluate(CORE::Elements::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;

      virtual int EvaluateEHLMon(CORE::Elements::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;

      /*!
        This class does not provide a definition for this function; it
        must be defined in LubricationEleCalc.
        The EvaluateService() method is meant for everything not related to
        the assembling of the linearized matrix and the right hand side
      */
      virtual int EvaluateService(CORE::Elements::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;
    };
  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
