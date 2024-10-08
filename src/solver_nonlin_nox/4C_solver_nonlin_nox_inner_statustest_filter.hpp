/*----------------------------------------------------------------------------*/
/*! \file

\brief Inner status test class for constraint problems. Filter
       techniques are based on ideas from multi-objective optimization:

       - Control of the two distinct goals of minimization of the objective
         function and satisfaction of the constraints.

       - Unlike merit functions, filter methods keep these two goals separate



\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_INNER_STATUSTEST_FILTER_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_INNER_STATUSTEST_FILTER_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensevector.hpp"
#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"
#include "4C_solver_nonlin_nox_inner_statustest_generic.hpp"
#include "4C_solver_nonlin_nox_meritfunction_infeasibility.hpp"

#include <Teuchos_RCP.hpp>

#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    class Group;
    namespace LineSearch
    {
      class Generic;
    }  // namespace LineSearch
    namespace Inner
    {
      namespace StatusTest
      {
        /** \brief binary class to compare two RCPs */
        template <typename T>
        struct RcpComp
        {
          inline bool operator()(const Teuchos::RCP<T>& i, const Teuchos::RCP<T>& j) const
          {
            return i.get() < j.get();
          }
        };

        /// status types of the filter acceptability test
        enum class FilterStatusType : char
        {
          passed_point_by_point,  ///< passed the filter acceptance test
          rejected,               ///< rejected by the filter
          unevaluated             ///< unevaluated
        };

        /*--------------------------------------------------------------------*/
        /** \brief container for all filter input variables
         *
         *  \author hiermeier */
        struct FilterParams
        {
          /** \brief Armijo status test.
           *
           *  Necessary for the L-type/f-type check or if no active constraints
           *  are left, thus the problem degenerates to an unconstrained problem. */
          Teuchos::RCP<Generic> armijo_ = Teuchos::null;

          /// set of the used infeasibility measures
          std::vector<Teuchos::RCP<::NOX::MeritFunction::Generic>> infeasibility_vec_;

          /// weight for the objective function value
          double weight_objective_func_ = 0.0;

          /// weight for the infeasibility function values
          double weight_infeasibility_func_ = 0.0;

          /// exponent needed for the L-/f-type switching condition
          double sf_ = 0.0;

          /// exponent needed for the L-/f-type switching condition
          double st_ = 0.0;

          /** This factor is used to address possible errors in the minimal step
           *  length estimation. These estimates are based on models and might be
           *  too large and, thus, this parameter can be used as an additional
           *  reduction factor. */
          double gamma_alpha_ = 0.0;

          /** \brief True if Second Order Correction shall be used
           *
           *  It is highly recommended to active it. */
          bool use_soc_ = false;

          /// second order correction strategy
          NOX::Nln::CorrectionType soc_type_ = NOX::Nln::CorrectionType::vague;

          /** If a good iterate is in N consecutive Newton iterates blocked by
           *  the filter set, while it would have been accepted by the remaining
           *  tests and, furthermore, it shows a sufficient reduction wrt the
           *  constraint violation, then the filter is reinitialized. The number
           *  N is specified by this parameter. */
          unsigned consecutive_blocking_iterates_ = 0;

          /** If during one Newton iterate a good iterate is in N consecutive line
           *  search steps blocked by the filter set, while it would have been
           *  accepted by the remaining tests and, furthermore, it shows a sufficient
           *  reduction wrt the constraint violation, then the filter is reinitialized.
           *  The number N is specified by this parameter. */
          unsigned consecutive_blocking_ls_steps_ = 0;

          /// after each reinitialization the max theta value is reduced by this factor
          double max_theta_blocking_red_ = 0.0;

          /// initial scaling for the max theta values
          double init_max_theta_blocking_scaling_ = 0.0;
        };

        /*--------------------------------------------------------------------*/
        /** \brief Filter inner status test
         *
         *  \author hiermeier */
        class Filter : public Generic
        {
          class Point;
          typedef std::pair<Teuchos::RCP<Point>, Teuchos::RCP<Point>> plain_point_pair;
          typedef std::pair<Teuchos::RCP<const Point>, Teuchos::RCP<const Point>>
              plain_const_point_pair;
          typedef std::vector<Teuchos::RCP<Point>> plain_point_set;
          typedef std::vector<Teuchos::RCP<::NOX::MeritFunction::Generic>> plain_merit_func_set;

         public:
          /// \brief constructor
          /** \param armijo                    (in) : internal armijo check if the f-type switching
           *                                          condition is fulfilled
           *  \param infeasibility_vec         (in) : vector of all defined infeasibility merit
           * functions \param weight_objective_func     (in) : weight for the objective function
           * values in the filter set \param weight_infeasibility_func (in) : weight for the
           * infeasibility function values in the filter set \param sf                        (in) :
           * exponent of the objective contributions in the f-type switching condition \param st
           * (in) : exponent of the infeasibility contributions in the f-type switching condition
           *  \param gamma_alpha               (in) : safety factor for the minimal step length
           * estimates ( must be < 1.0 ) \param utils                     (in) : in/output stream
           * manager object
           *
           *  \author hiermeier \date 04/17 */
          Filter(const FilterParams& fparams, const ::NOX::Utils& utils);

          /** \brief %Test the inner stopping criterion
           *
           *  The test can (and should, if possible) be skipped if
           *  checkType is NOX::StatusType::None.  If the test is skipped, then
           *  the status should be set to ::NOX::StatusTest::Unevaluated. */
          StatusType check_status(const Interface::Required& interface,
              const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
              ::NOX::StatusTest::CheckType checkType) override;

          //! Return the result of the most recent inner checkStatus call
          StatusType get_status() const override;

          //! Output formatted description of inner stopping test to output stream.
          std::ostream& print(std::ostream& stream, int indent = 0) const override;

         private:
          /// initialize the (filter) points
          void init_points(const Interface::Required& interface,
              const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp);

          /// evaluate and set new trial point
          void set_trial_point(
              const ::NOX::MeritFunction::Generic& merit_func, const ::NOX::Abstract::Group& grp);

          /// reset the internal state at the beginning of a new Newton iteration
          void reset();

          enum NOX::Nln::Inner::StatusTest::FilterStatusType acceptability_check(
              const Point& trial_fp);

          /** \brief Perform a pre-selection to avoid the unnecessary point-by-point
           *  comparison in obvious cases
           *
           * \author hiermeier */
          unsigned prefiltering(const Point& trial_fp);

          /** identify all points in the current filter set which are not
           * dominated by the new trial point */
          void identify_non_dominated_filter_points(
              const Point& trial_fp, const unsigned non_dominated_index);

          /// set-up all model terms
          void setup_model_terms(const ::NOX::Abstract::Vector& dir,
              const ::NOX::Abstract::Group& grp, const Interface::Required& interface);

          /// execute the sufficient reduction check
          enum NOX::Nln::Inner::StatusTest::StatusType sufficient_reduction_check(
              const Point& trial_fp) const;

          /// Is the step still larger than the minimal step length estimate?
          NOX::Nln::Inner::StatusTest::StatusType is_admissible_step(
              const ::NOX::Solver::Generic& solver, const double& step) const;

          /// access the active set status
          enum ::NOX::StatusTest::StatusType get_active_set_status(
              const ::NOX::Solver::Generic& solver) const;

          /// get the specified constraint tolerance
          double get_constraint_tolerance(const ::NOX::Solver::Generic& solver) const;

          /** \brief Augment the current filter
           *
           *  The filter is only augmented if its a non-f-type step and the
           *  filter criterion in respect to the previous filter point is fulfilled.
           *
           *  \author hiermeier \date 04/17 */
          void augment_filter();

          /** \brief Compute the minimal step length estimates based on the different models
           *
           *  This method initiates the calculation of the step length estimates based
           *  on the objective function model, the infeasibility measures and the
           *  f-type switching condition. In the end a final minimal step length estimate
           *  is set.
           *
           *  \author hiermeier \date 04/17 */
          void compute_minimal_step_length_estimates();

          /// Compute the minimal step length estimate based on the objective function model
          double minimal_step_length_estimate_of_obj_func_filter_check() const;

          /** \brief Check the F-Type switching condition for the given step-length
           *
           *  \param step (in): current step-length
           *
           *  \return TRUE if the condition is fulfilled, otherwise false.
           *
           *  \author hiermeier \date 04/17 */
          bool check_f_type_switching_condition(const double step) const;

          /** \brief Evaluate the F-Type switching condition
           *
           *  \param step (in) : current step-length
           *  \param d    (in) : internal scaling factor (must be larger than zero)
           *
           *  \author hiermeier \date 04/17 */
          double compute_f_type_switching_condition(const double step, const double d) const;

          /** \brief Compute the minimal step length estimate based on the f-type switching
           *         condition
           *
           *  This routine uses a local Newton scheme for the calculation of the estimate.
           *
           *  \author hiermeier \date 04/17 */
          double minimal_step_length_estimate_of_f_type_condition() const;

          /// evaluate the objective model based on the given step length
          double get_obj_model(const double step) const;

          /** \brief Translate the filter acceptability status
           *
           *  The following completes the sentence "The filter ...".
           *
           *  \author hiermeier \date 08/17  */
          inline std::string filter_status_to_string(enum FilterStatusType filter_status) const
          {
            switch (filter_status)
            {
              case FilterStatusType::passed_point_by_point:
                return "accepted the trial filter point via point by point comparison";
              case FilterStatusType::rejected:
                return "rejected the trial filter point";
              case FilterStatusType::unevaluated:
                return "acceptability test is unevaluated";
              default:
                return "acceptability test has an undefined status";
            }
          }

          /// executed in the end of the check status test
          StatusType post_check_status(const NOX::Nln::LineSearch::Generic& linesearch,
              const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
              ::NOX::StatusTest::CheckType checkType);

          /// actual test
          void execute_check_status(const NOX::Nln::LineSearch::Generic& linesearch,
              const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
              ::NOX::StatusTest::CheckType checkType);

          /// recover from a back up state if the SOC fails
          void recover_from_backup(::NOX::Abstract::Group& grp) const;

          /// throw an error if all strategies fail and the step is too short
          void throw_if_step_too_short(const NOX::Nln::LineSearch::Generic& linesearch,
              const ::NOX::Solver::Generic& solver) const;

         protected:
          //! Status of the inner filter status test
          NOX::Nln::Inner::StatusTest::StatusType status_;

         private:
          /*------------------------------------------------------------------------*/
          /// nested backup state class
          class BackupState
          {
           public:
            /** \brief Create a backup of the lastly accepted state at the very
             *  beginning of each new Newton step
             *
             *  \param(in) grp: group containing the state which shall be considered
             *                  for the back-up.
             *  \param(in) dir: current search direction
             *
             *  \author hiermeier \date 12/17 */
            void create(const ::NOX::Abstract::Group& grp, const ::NOX::Abstract::Vector& dir);

            /// recover from the back-up
            void recover(::NOX::Abstract::Group& grp) const;

           private:
            /** \brief check the recovered state
             *
             *  If the L2-norm of the recovered rhs differs more than machine
             *  precision from the backup state rhs, an error will be thrown. */
            void check_recovered_state(const ::NOX::Abstract::Vector& f) const;

           private:
            Teuchos::RCP<::NOX::Epetra::Vector> xvector_ = Teuchos::null;
            double normf_ = 0.0;
          };

          /// \brief Second Order Correction base class
          /** The base class is a empty dummy class which is going to be built, if
           *  no SOC steps shall be considered.
           *
           *  \author hiermeier \date 12/17 */
          class SOCBase
          {
           public:
            /// create the second order correction object
            static Teuchos::RCP<SOCBase> create(
                Filter& filter, const bool use_soc, const CorrectionType user_type)
            {
              if (use_soc)
                return Teuchos::RCP<SOCBase>(new SecondOrderCorrection(filter, user_type));
              else
                return Teuchos::RCP<SOCBase>(new SOCBase(filter, NOX::Nln::CorrectionType::vague));
            }

            /// base class constructor
            SOCBase(Filter& filter, NOX::Nln::CorrectionType user_type)
                : filter_(filter), user_type_(user_type){/* empty */};

            /// delete default constructor
            SOCBase() = delete;

            /// use default destructor
            virtual ~SOCBase() = default;

            /// The base class does not perform a SOC step
            virtual StatusType execute(const NOX::Nln::LineSearch::Generic& linesearch,
                const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
                ::NOX::StatusTest::CheckType checkType)
            {
              return filter_.get_status();
            };

           protected:
            /// call back to the filter
            Filter& filter_;

            /// user defined SOC tye
            const NOX::Nln::CorrectionType user_type_;
          };

          /// \brief Concrete implementation of a Second Order Correction class
          class SecondOrderCorrection : public SOCBase
          {
           public:
            SecondOrderCorrection(Filter& filter, NOX::Nln::CorrectionType user_type)
                : SOCBase(filter, user_type){/* empty */};

            SecondOrderCorrection() = delete;

            /// compute the SOC step
            StatusType execute(const NOX::Nln::LineSearch::Generic& linesearch,
                const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
                ::NOX::StatusTest::CheckType checkType) override;

           private:
            /// compute the SOC system
            void compute_system(NOX::Nln::Group& grp, const ::NOX::Solver::Generic& solver) const;

            /// solve the SOC system
            void solve(const NOX::Nln::LineSearch::Generic& linesearch,
                const ::NOX::Solver::Generic& solver, ::NOX::Abstract::Group& grp) const;

            /// postprocess the SOC step
            void postprocess(const NOX::Nln::LineSearch::Generic& linesearch,
                const ::NOX::Solver::Generic& solver, ::NOX::Abstract::Group& grp,
                ::NOX::StatusTest::CheckType checkType);

            /// Which SOC system shall be used?
            CorrectionType which_type(const ::NOX::Solver::Generic& solver) const;

            /// Use an automatic type choice (recommended).
            CorrectionType automatic_type_choice(const ::NOX::Solver::Generic& solver) const;

            /// print infos about the SOC step
            void print(std::ostream& os) const;

           private:
            /// currently chosen SOC type
            CorrectionType curr_type_ = CorrectionType::vague;

            /// execution time for the SOC step
            double time_exe_ = 0.0;

            /// needed time for the recovery
            double time_recover_ = 0.0;

            /// Are we currently in a SOC step? Avoids recursive calls.
            bool issoc_ = false;

            /// What is the result of the SOC step?
            StatusType soc_status_ = status_unevaluated;
          };

          /*------------------------------------------------------------------------*/
          /** \brief Helps to detect a blocking filter set due to old historic information
           *
           * \author hiermeier */
          class Blocking
          {
           public:
            /// simple create method
            static Teuchos::RCP<Blocking> create(Filter& filter, const FilterParams& fparams)
            {
              return Teuchos::RCP(new Blocking(filter, fparams));
            }

            /// constructor
            Blocking(Filter& filter, const FilterParams& fparams)
                : consecutive_iter_(fparams.consecutive_blocking_iterates_),
                  consecutive_ls_steps_(fparams.consecutive_blocking_ls_steps_),
                  max_theta_red_(fparams.max_theta_blocking_red_),
                  init_max_theta_scaling_(fparams.init_max_theta_blocking_scaling_),
                  filter_(filter){};

            /* \brief check if the current filter set blocks good iterates
             *
             * If the filter rejected the trial point but the inner test would accept it
             * we have an indicator for a blocking filter set. This can happen due to old
             * historic information which is not reliable for the current neighborhood. */
            void check(const NOX::Nln::LineSearch::Generic& linesearch,
                const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
                const Point& rejected_fp);

            /** add the current iterate to the blocking set if the sufficient
             *  reduction criterion with respect to the user-specified max theta
             *  value is fulfilled. */
            void add_filter_iterate(const int newton_iter, const Point& rejected_fp);

            /// initialize a reinitialization of the filter
            void reinitialize_filter();

            /** print some information about the blocking scenario and a possible
             *  reinitialization */
            void print_info(std::ostream& os) const;

            // vector containing consecutive blocking filter iteration numbers
            std::vector<std::pair<unsigned, unsigned>> filter_iterates_;

            /// number of consecutive allowed blocking Newton iterates
            const unsigned consecutive_iter_;

            /// number of consecutive allowed blocking line search steps
            const unsigned consecutive_ls_steps_;

            /** max theta reduction value applied in case of an reinitialization
             *  and used for the additional sufficient reducition check */
            const double max_theta_red_;

            /// initial max theta scaling
            const double init_max_theta_scaling_;

           private:
            /// filter call-back
            Filter& filter_;
          };

          /*------------------------------------------------------------------------*/
          /// nested structure representing a set of infeasibility measures
          struct Infeasibility
          {
            Infeasibility(const plain_merit_func_set& infeasibility_vec)
                : vector_(infeasibility_vec), number_(infeasibility_vec.size())
            { /* empty */
            }

            /// evaluate the function values of the infeasibility merit functions
            void computef(double* theta_values, const ::NOX::Abstract::Group& grp) const;

            /// find the maximal infeasibility measure in a set of theta values
            unsigned find_max_theta_id(double* theta_values) const;

            /// compute the slope of all infeasibility merit functions
            void compute_slope(const ::NOX::Abstract::Vector& dir,
                const ::NOX::Abstract::Group& grp, double* theta_slope_values) const;

            /// compute mixed 2-nd order terms of all infeasibility merit functions
            void compute_mixed2nd_order_terms(const ::NOX::Abstract::Vector& dir,
                const ::NOX::Abstract::Group& grp, double* theta_mixed_values) const;

            /** \brief compute the over-all minimal step length estimate based on
             *  all infeasibility merit functions */
            double minimal_step_length_estimate(
                const double* accepted_theta, const double* theta_slope) const;

            /// set of infeasibility merit functions
            const plain_merit_func_set vector_;

            /// total number of all infeasibility measures
            const unsigned number_;

          };  // struct Infeasibility

          /*------------------------------------------------------------------------*/
          /// nested class representing the filter point
          class Point
          {
           public:
            /// create a new point ( NOT a filter point)
            static Teuchos::RCP<Point> create(const ::NOX::MeritFunction::Generic& merit_func,
                const Infeasibility& infeasibility_func, const ::NOX::Abstract::Group& grp);

            /// create a new filter point from an existing point
            static Teuchos::RCP<Point> make_filter_point(const Point& p, const bool do_scaling);

            /// (re)set all global point member variables
            static void reset_static_members(const unsigned num_obj_coords,
                const unsigned num_theta_coords, const double weight_objective_func,
                const double weight_infeasibility_func, const double init_max_theta_scale);

            /** (re)set global point member variables at the beginning of each new
             *  Newton iteration */
            static void reset_static_members();

            /// (re)set the global margin safety factors
            static void set_margin_safety_factors();

            static void reinit_filter(plain_point_set& filter,
                const Infeasibility& infeasibility_func, const double& downscale_fac);

           private:
            static void clear_filter_point_register();

            /** \brief Add a new filter point to the register
             *
             *  Before we add a new filter point all unused filter points in the
             *  register are removed. */
            static void add_filter_point_to_register(const Teuchos::RCP<Point>& fp_ptr);

            /** scale the coordinate with the given %id as soon as the corresponding
             * scaling changes its state from invalid to valid */
            static void scale_coordinate_of_all_registered_filter_points(const int id);

            static void set_initial_scaled_max_theta_value(const int id, const double& val);

            static void scale_max_theta_values(const double& fac);

            /// constructor
            Point()
                : is_filter_point_(false),
                  is_feasible_(false),
                  norm_(-1.0),
                  max_theta_id_(0),
                  coords_(num_coords_, true),
                  margin_(num_coords_, true)
            { /*nothing to do here */
            }

            /// copy constructor
            Point(const Point& point)
                : is_filter_point_(point.is_filter_point_),
                  is_feasible_(point.is_feasible_),
                  norm_(point.norm_),
                  max_theta_id_(point.max_theta_id_),
                  coords_(point.coords_),
                  margin_(point.margin_)
            { /* nothing to do here */
            }

           public:
            /// calculate and set the point norm value
            void set_norm();

            /// scale point coordinates
            void scale();

            /// set margin values for each filter coordinate
            void set_margin();

            // return true if the current point is feasible with respect to one
            // infeasibility measure
            bool is_feasible(const double tol) const;

            //
            bool is_sufficiently_reduced_compared_to_max_theta(const double& red_fac) const;

            /** \brief access the point coordinate %index
             *
             *  \param index (in) : id of the point coordinate entry */
            inline double& operator()(unsigned index) { return coords_(index); }

            /** \brief access the point coordinate %index (read-only)
             *
             *  \param index (in) : id of the point coordinate entry */
            inline const double& operator()(unsigned index) const { return coords_(index); }

            /// access the data pointer of the point coordinates
            inline double* data() { return coords_.values(); }

            /// access the data pointer of the point coordinates (read-only)
            inline const double* data() const { return coords_.values(); }

            /// return the maximal infeasibility measure of this point
            inline double max_theta() const { return coords_(num_obj_coords_ + max_theta_id_); }

            /// return the scaling factor of the maximal infeasibility measure of this point
            inline double scale_of_max_theta() const
            {
              return scale_(num_obj_coords_ + max_theta_id_);
            }

            /// print the current point
            std::ostream& print(
                std::ostream& stream, int par_indent_length, const ::NOX::Utils* u) const;
            inline std::ostream& print(
                std::ostream& stream, int par_indent_length, const ::NOX::Utils& u) const
            {
              return print(stream, par_indent_length, &u);
            };
            inline std::ostream& print(std::ostream& stream, const ::NOX::Utils& u) const
            {
              return print(stream, 0, u);
            };
            std::ostream& print(std::ostream& stream) const { return print(stream, 0, nullptr); };

            /// the point is a filter point
            bool is_filter_point_;

            /// the point is feasible
            mutable bool is_feasible_;

            /// norm of the filter point
            double norm_;

            /// id of the maximal infeasibility measure coordinate of this point
            int max_theta_id_;

            /// filter point coordinates
            Core::LinAlg::SerialDenseVector coords_;

            /// margin of each filter point coordinate
            Core::LinAlg::SerialDenseVector margin_;

            /// global number of coordinates per filter point
            static unsigned num_coords_;

            /// global number of objective coordinates per filter point
            static unsigned num_obj_coords_;

            /// global scaling factor of the objective function margin
            static double gamma_obj_;

            /// global scaling factor of the infeasibility function margin
            static double gamma_theta_;

            /// global max theta scale (defined via the blocking class)
            static double global_init_max_theta_scale_;

            static std::vector<bool> isvalid_scaling_;

            /// global scaling of each coordinate
            static Core::LinAlg::SerialDenseVector scale_;

            /// global weights for the filter point scaling
            static Core::LinAlg::SerialDenseVector weights_;

            /// global maximal infeasibility values
            static Core::LinAlg::SerialDenseVector global_scaled_max_thetas_;

           private:
            /// this vector contains all registered filter points
            static std::set<Teuchos::RCP<Point>, RcpComp<Point>> filter_point_register_;

          };  // struct Point

          Infeasibility theta_;

          /// pair of the current trial point (second) and the previous accepted point (first)
          plain_const_point_pair curr_points_;

          /// pair of the current trial filter point (second) and the previous accepted point
          /// (first)
          plain_point_pair curr_fpoints_;

          /// ordered set of filter points
          plain_point_set filter_;

          /** set of non-dominated filter points, these points won't be removed during
           *  the filter augmentation */
          plain_point_set non_dominated_filter_points_;

          /* backup state object. E.g. useful for recovery of the last accepted step
           * if the second order correction step fails to achieve a better solution. */
          BackupState backup_;

          Teuchos::RCP<SOCBase> soc_;

          /// blocking object. Helps to detect a blocking filter due to historic information
          Teuchos::RCP<Blocking> blocking_ptr_;
          Blocking& blocking_;

          /// safety factor for the minimal step length check
          const double gamma_alpha_;

          /// minimal step length estimate derived from the objective function model
          double amin_obj_;

          /// minimal step length estimate derived from the infeasibility function models
          double amin_theta_;

          /// minimal step length estimate derived from the f-type switching condition
          double amin_ftype_;

          /// over-all minimal step length estimate
          double amin_;

          /// exponent of the objective merit function contributions in the f-type condition
          const double sf_;

          /// exponent of the theta/infeasibility merit function contributions in the f-type
          /// condition
          const double st_;

          /// linear model terms / slopes of the objective and infeasibility merit-functions
          Core::LinAlg::SerialDenseVector model_lin_terms_;

          /// mixed 2-nd order terms of the objective and infeasibility merit-functions
          Core::LinAlg::SerialDenseVector model_mixed_terms_;

          /// armijo inner status test object
          Teuchos::RCP<Generic> armijo_test_;

          /** \brief Does the current step fulfill the f-type switching condition?
           *
           *  If this variable is TRUE, the inner armijo test will be checked. */
          bool is_ftype_step_;

          // theta_min value used to skip the f-type condition in pre-asymptotic phase
          double theta_min_ftype_;

          /// internal status of the filter method
          enum FilterStatusType filter_status_;

          /// nox output management object
          const ::NOX::Utils& utils_;

          /// output precision of the scientific numbers in the print methods
          static constexpr int OUTPUT_PRECISION = 15;
        };  // class Filter
      }     // namespace StatusTest
    }       // namespace Inner
  }         // namespace Nln
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
