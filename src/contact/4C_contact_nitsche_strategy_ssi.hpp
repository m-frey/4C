/*----------------------------------------------------------------------------*/
/*! \file
\brief Nitsche ssi contact solving strategy

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_NITSCHE_STRATEGY_SSI_HPP
#define FOUR_C_CONTACT_NITSCHE_STRATEGY_SSI_HPP

#include "4C_config.hpp"

#include "4C_contact_nitsche_strategy.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
   * @brief Contact solving strategy with Nitsche's method.
   *
   * This is a specialization of the abstract contact algorithm as defined in AbstractStrategy.
   * For a more general documentation of the involved functions refer to CONTACT::AbstractStrategy.
   */
  class NitscheStrategySsi : public NitscheStrategy
  {
   public:
    //! Standard constructor
    NitscheStrategySsi(const Epetra_Map* dofRowMap, const Epetra_Map* nodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        const Teuchos::RCP<Epetra_Comm>& comm, double alphaf, int maxdof)
        : NitscheStrategy(
              dofRowMap, nodeRowMap, params, std::move(interface), dim, comm, alphaf, maxdof)
    {
    }

    //! Shared data constructor
    NitscheStrategySsi(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
        const Epetra_Map* dofRowMap, const Epetra_Map* nodeRowMap,
        const Teuchos::ParameterList& params,
        std::vector<Teuchos::RCP<CONTACT::Interface>> interface, int dim,
        const Teuchos::RCP<const Epetra_Comm>& comm, double alphaf, int maxdof)
        : NitscheStrategy(data_ptr, dofRowMap, nodeRowMap, params, std::move(interface), dim, comm,
              alphaf, maxdof)
    {
    }

    void apply_force_stiff_cmt(Teuchos::RCP<Core::LinAlg::Vector> dis,
        Teuchos::RCP<Core::LinAlg::SparseOperator>& kt, Teuchos::RCP<Core::LinAlg::Vector>& f,
        const int step, const int iter, bool predictor) override
    {
      if (kt != Teuchos::null && f != Teuchos::null)
      {
        NitscheStrategy::apply_force_stiff_cmt(dis, kt, f, step, iter, predictor);
      }
      curr_state_eval_ = true;
    }

    /*!
     * @note currently there is no support for friction; furthermore we can not call
     * 'update_trace_ineq_etimates' at the moment since at this stage we do not have the gp
     * concentrations -> so far no adaptive nitsche penalty possible!
     */
    void evaluate_reference_state() override;

    /*!
     * @brief initialize TraceHE parameter of the mortar elements
     *
     * @note: this approximation is done here since update_trace_ineq_etimates() can not be called
     * at the moment
     */
    void init_trace_he();

    void integrate(const CONTACT::ParamsInterface& cparams) override;

    void set_state(
        const enum Mortar::StateType& statename, const Core::LinAlg::Vector& vec) override;

    void set_parent_state(const enum Mortar::StateType& statename, const Core::LinAlg::Vector& vec,
        const Core::FE::Discretization& dis) override;

    Teuchos::RCP<const Core::LinAlg::Vector> get_rhs_block_ptr(
        const enum CONTACT::VecBlockType& bp) const override;

    /*!
     * @brief get the pointer to the matrix block
     *
     * @param[in] bt   block type of requested matrix block
     * @return pointer to matrix block
     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bp) const;

    //! don't want = operator
    NitscheStrategySsi operator=(const NitscheStrategySsi& old) = delete;
    //! don't want copy constructor
    NitscheStrategySsi(const NitscheStrategySsi& old) = delete;

   protected:
    Teuchos::RCP<Epetra_FEVector> setup_rhs_block_vec(
        const enum CONTACT::VecBlockType& bt) const override;

    Teuchos::RCP<Core::LinAlg::SparseMatrix> setup_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bt) override;

    void complete_matrix_block_ptr(
        const enum CONTACT::MatBlockType& bt, Teuchos::RCP<Core::LinAlg::SparseMatrix> kc) override;

    //! current scalar state vector
    Teuchos::RCP<Core::LinAlg::Vector> curr_state_scalar_;
    //! ScaTra residual
    Teuchos::RCP<Epetra_FEVector> fs_;
    //! linearization of ScaTra residual w.r.t ScaTra dofs
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kss_;
    //! linearization of ScaTra residual w.r.t displacement dofs
    Teuchos::RCP<Core::LinAlg::SparseMatrix> ksd_;
    //! linearization of displacement residual w.r.t ScaTra dofs
    Teuchos::RCP<Core::LinAlg::SparseMatrix> kds_;
  };
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
