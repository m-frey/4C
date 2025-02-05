// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FLUID_ELE_CALC_IMMERSED_HPP
#define FOUR_C_FLUID_ELE_CALC_IMMERSED_HPP

#include "4C_config.hpp"

#include "4C_fluid_ele_calc.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  namespace Elements
  {
    class FluidImmersedBase;

    template <Core::FE::CellType distype>
    class FluidEleCalcImmersed : public FluidEleCalc<distype>
    {
      typedef Discret::Elements::FluidEleCalc<distype> my;
      using my::nen_;
      using my::nsd_;

     protected:
      /// private Constructor since we are a Singleton.
      FluidEleCalcImmersed();

     public:
      /// Singleton access method
      static FluidEleCalcImmersed<distype>* instance(
          Core::Utils::SingletonAction action = Core::Utils::SingletonAction::create);


     protected:
      /*!
        Evaluate

        \param eid              (i) element id
        \param discretization   (i) fluid discretization the element belongs to
        \param lm               (i) location matrix of element
        \param params           (i) element parameter list
        \param mat              (i) material
        \param elemat1_epetra   (o) element matrix to calculate
        \param elemat2_epetra   (o) element matrix to calculate
        \param elevec1_epetra   (o) element vector to calculate
        \param elevec2_epetra   (o) element vector to calculate
        \param elevec3_epetra   (o) element vector to calculate
        \param offdiag          (i) flag indicating whether diagonal or off diagonal blocks are to
        be calculated

       */
      int evaluate(Discret::Elements::Fluid* ele, Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          std::shared_ptr<Core::Mat::Material>& mat,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra, bool offdiag = false) override;

      //! compute residual of momentum equation and subgrid-scale velocity
      void compute_subgrid_scale_velocity(
          const Core::LinAlg::Matrix<nsd_, nen_>& eaccam,  ///< acceleration at time n+alpha_M
          double& fac1,                                    ///< factor for old s.-s. velocities
          double& fac2,                                    ///< factor for old s.-s. accelerations
          double& fac3,     ///< factor for residual in current s.-s. velocities
          double& facMtau,  ///< facMtau = modified tau_M (see code)
          int iquad,        ///< integration point
          double* saccn,    ///< s.-s. acceleration at time n+alpha_a / n
          double* sveln,    ///< s.-s. velocity at time n+alpha_a / n
          double* svelnp    ///< s.-s. velocity at time n+alpha_f / n+1
          ) override;

      //! Provide linearization of Garlerkin momentum residual with respect to the velocities
      void lin_gal_mom_res_u(Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
                                 lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          const double& timefacfac             ///< = timefac x fac
          ) override;

      //! Compute element matrix and rhs entries: inertia, convective andyn
      //! reactive terms of the Galerkin part
      void inertia_convection_reaction_gal_part(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                                    estif_u,  ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,         ///< rhs forces velocity
          Core::LinAlg::Matrix<nsd_ * nsd_, nen_>&
              lin_resM_Du,  ///< linearisation of the Garlerkin momentum residual
          Core::LinAlg::Matrix<nsd_, 1>&
              resM_Du,          ///< linearisation of the Garlerkin momentum residual
          const double& rhsfac  ///< right-hand-side factor
          ) override;

      //! Compute element matrix entries: continuity terms of the Garlerkin part and rhs
      void continuity_gal_part(
          Core::LinAlg::Matrix<nen_, nen_ * nsd_>& estif_q_u,  ///< block (weighting function q x u)
          Core::LinAlg::Matrix<nen_, 1>& preforce,             ///< rhs forces pressure
          const double& timefacfac,                            ///< = timefac x fac
          const double& timefacfacpre,
          const double& rhsfac  ///< right-hand-side factor
          ) override;

      //! Compute element matrix entries: conservative formulation
      void conservative_formulation(Core::LinAlg::Matrix<nen_ * nsd_, nen_ * nsd_>&
                                        estif_u,       ///< block (weighting function v x u)
          Core::LinAlg::Matrix<nsd_, nen_>& velforce,  ///< rhs forces velocity
          const double& timefacfac,                    ///< = timefac x fac
          const double& rhsfac                         ///< right-hand-side factor
          ) override;

      // current element
      Discret::Elements::FluidImmersedBase* immersedele_;
      // number of current gp
      int gp_iquad_;

    };  // class FluidEleCalcImmersed
  }     // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
