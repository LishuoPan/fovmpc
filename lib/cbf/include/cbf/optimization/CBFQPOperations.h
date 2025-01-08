//
// Created by lishuo on 8/31/24.
//

#ifndef CBF_OPTIMIZATION_CBFQPOPERATIONS_H
#define CBF_OPTIMIZATION_CBFQPOPERATIONS_H

#include <qpcpp/QPOperations.h>
#include <cbf/detail/cbf.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    class CBFQPOperations {
    public:
        using CostAddition = typename qpcpp::QPOperations<T>::CostAddition;
        using LinearConstraint = typename qpcpp::QPOperations<T>::LinearConstraint;
        using DecisionVariableBounds = typename qpcpp::QPOperations<T>::DecisionVariableBounds;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using Vector = math::Vector<T>;
        using Matrix = math::Matrix<T>;

        CBFQPOperations(std::shared_ptr<FovCBF> cbf);
        size_t numDecisionVariables() const {return DIM;}

        CostAddition desiredControlCost(const VectorDIM& desired_u);
        LinearConstraint safetyConstraint(const Vector& state, const Vector& target_state);
        LinearConstraint leftBorderConstraint(const Vector& state, const Vector& target_state);
        LinearConstraint rightBorderConstraint(const Vector& state, const Vector& target_state);
        LinearConstraint rangeConstraint(const Vector& state, const Vector& target_state);
        LinearConstraint leftBorderConstraintWithSlackVar(const Vector& state, const Vector& target_state, const T &slack);
        LinearConstraint rightBorderConstraintWithSlackVar(const Vector& state, const Vector& target_state, const T &slack);
        LinearConstraint rangeConstraintWithSlackVar(const Vector& state, const Vector& target_state, const T &slack);
        DecisionVariableBounds controlBoundConstraint(const VectorDIM& u_min, const VectorDIM& u_max);

    private:
        std::shared_ptr<FovCBF> cbf_;
    };

} // cbf

#endif //CBF_OPTIMIZATION_CBFQPOPERATIONS_H
