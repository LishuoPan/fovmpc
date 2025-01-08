//
// Created by lishuo on 8/31/24.
//

#include <cbf/optimization/CBFQPOperations.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    CBFQPOperations<T, DIM>::CBFQPOperations(std::shared_ptr<FovCBF> cbf) : cbf_(cbf) {}

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::CostAddition
    CBFQPOperations<T, DIM>::desiredControlCost(const VectorDIM &desired_u) {
        Matrix quadratic_term(DIM, DIM);
        quadratic_term.setZero();
        Vector linear_term(DIM);
        linear_term.setZero();

        // set the quadratic term
        quadratic_term.setIdentity();
        // set the linear term
        linear_term = -2.0 * desired_u;
        // set the constant term
        T constant_term = desired_u.transpose() * desired_u;

        return CostAddition(quadratic_term, linear_term, constant_term);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::safetyConstraint(const Vector &state,
                                              const Vector &target_state) {
        Vector coefficients = -1.0 * cbf_->getSafetyConstraints(state, target_state);
        T bound = cbf_->getSafetyBound(state, target_state);
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::leftBorderConstraint(const Vector &state,
                                                  const Vector &target_state) {
        Vector coefficients = -1.0 * cbf_->getLBConstraints(state, target_state);
        T bound = cbf_->getLBBound(state, target_state);
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rightBorderConstraint(const Vector &state,
                                                   const Vector &target_state) {
        Vector coefficients = -1.0 * cbf_->getRBConstraints(state, target_state);
        T bound = cbf_->getRBBound(state, target_state);
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rangeConstraint(const Vector &state,
                                             const Vector &target_state) {
        Vector coefficients = -1.0 * cbf_->getRangeConstraints(state, target_state);
        T bound = cbf_->getRangeBound(state, target_state);
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::leftBorderConstraintWithSlackVar(const Vector &state,
                                                  const Vector &target_state,
                                                  const T &slack) {
        Vector coefficients = -1.0 * cbf_->getLBConstraints(state, target_state);
        T bound = cbf_->getLBBound(state, target_state);
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound+slack);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rightBorderConstraintWithSlackVar(const Vector &state,
                                                   const Vector &target_state,
                                                   const T &slack) {
        Vector coefficients = -1.0 * cbf_->getRBConstraints(state, target_state);
        T bound = cbf_->getRBBound(state, target_state);
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound+slack);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::LinearConstraint
    CBFQPOperations<T, DIM>::rangeConstraintWithSlackVar(const Vector &state,
                                             const Vector &target_state,
                                             const T &slack) {
        Vector coefficients = -1.0 * cbf_->getRangeConstraints(state, target_state);
        T bound = cbf_->getRangeBound(state, target_state);
        return LinearConstraint(coefficients, std::numeric_limits<T>::lowest(), bound);
    }

    template <typename T, unsigned int DIM>
    typename CBFQPOperations<T, DIM>::DecisionVariableBounds
    CBFQPOperations<T, DIM>::controlBoundConstraint(const VectorDIM &u_min,
                                                    const VectorDIM &u_max) {
        VectorDIM lower_bounds;
        VectorDIM upper_bounds;
        for (size_t d = 0; d < DIM; ++d) {
            lower_bounds(d) = u_min(d);
            upper_bounds(d) = u_max(d);
        }
        return DecisionVariableBounds(lower_bounds, upper_bounds);
    }

    template class CBFQPOperations<double, 3U>;

} // cbf