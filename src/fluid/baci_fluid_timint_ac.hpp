/*-----------------------------------------------------------*/
/*! \file

\brief Fluid time integrator for FS3I-AC problems


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_AC_HPP
#define BACI_FLUID_TIMINT_AC_HPP


#include "baci_config.hpp"

#include "baci_fluid_implicit_integration.hpp"

BACI_NAMESPACE_OPEN

namespace FLD
{
  class TimIntAC : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntAC(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /// read restart from step
    void ReadRestart(int step) override;

    /// write output
    void Output() override;

   protected:
   private:
  };  // class TimIntAC

}  // namespace FLD


BACI_NAMESPACE_CLOSE

#endif  // FLUID_TIMINT_AC_H
