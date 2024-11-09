//
// Created by lishuo on 2/13/24.
//

#ifndef SPLINES_BEZIERQPOPERATIONS_H
#define SPLINES_BEZIERQPOPERATIONS_H

#include <math/Types.h>
#include <math/Helpers.h>
#include <splines/optimization/SingleParameterCurveQPOperations.h>

namespace splines {
    template <typename T, unsigned int DIM>
    class BezierQPOperations : public SingleParameterCurveQPOperations<T, DIM> {
    public:
        using Base = SingleParameterCurveQPOperations<T, DIM>;
        using Hyperplane = typename Base::Hyperplane;
        using AlignedBox = typename Base::AlignedBox;
        using Vector = math::Vector<T>;
        using Matrix = math::Matrix<T>;
        using Row = math::Row<T>;
        using VectorDIM = math::VectorDIM<T, DIM>;
        using CostAddition = typename Base::CostAddition;
        using LinearConstraint = typename Base::LinearConstraint;
        using DecisionVariableBounds = typename Base::DecisionVariableBounds;
        using SingleParameterCurve = splines::SingleParameterCurve<T, DIM>;

        struct Params {
            size_t num_control_points_;
            T max_parameter_;
        };

        BezierQPOperations(Params &p);

        virtual std::size_t numDecisionVariables() const override;
        virtual void set_max_parameter(T new_max_parameter) override;
        T max_parameter() const override;

        /// recover curve representation from solution
        virtual std::unique_ptr<SingleParameterCurve> generateCurveFromSolution(const Vector& decision_variables) const override;
        virtual Row evalBasisRow(unsigned int dimension, T parameter, uint64_t derivative_degree) const override;

        /// objective functions
        virtual CostAddition integratedSquaredDerivativeCost(uint64_t derivative_degree, T lambda) const override;

        virtual CostAddition evalCost(T parameter, uint64_t derivative_degree,
                                      const VectorDIM& target, T lambda) const override;

        /// constraints
        virtual std::vector<LinearConstraint> evalConstraint(
                T parameter, uint64_t derivative_degree,
                const VectorDIM& target) const override;

        virtual std::vector<LinearConstraint> evalBound(
                T parameter, uint64_t derivative_degree,
                const VectorDIM& LB,
                const VectorDIM& UB) const override;

        virtual DecisionVariableBounds boundingBoxConstraint(
                const AlignedBox& bounding_box) const override;

        virtual DecisionVariableBounds positionConstraint(const VectorDIM& position) const override;

        /**
         * @brief Constaint the given derivative of the curve to be inside the given
         * bounding box at all points
         *
         * @details Any derivative of a Bezier curve is also a Bezier curve, control
         * points of which are constant multiplies of the control points of the
         * original curve. This function returns constrains for the control points
         * of *this curve so that the control points of the derivative curve is
         * contained within the bounding box.
         *
         * @param bounding_box Bounding box that the given derivative of the curve
         * must be contained in
         * @param derivative_degree Derivative degree of the curve
         * @return std::vector<LinearConstraint> Set of constraints that enforces
         * the given derivative of the curve to be inside the given bounding box.
         */
        virtual std::vector<LinearConstraint>
        boundingBoxConstraintAll(const AlignedBox& bounding_box,
                                 uint64_t derivative_degree) const override;

        // epsilon is used to make the <= to < for the hyperplane constraints
        virtual std::vector<LinearConstraint> hyperplaneConstraintAll(
                const Hyperplane& hyperplane, T epsilon) const override;

        virtual std::vector<LinearConstraint>
        hyperplaneConstraintAt(T parameter, const Hyperplane& hyperplane, T epsilon) const override;

    private:
        std::size_t num_control_points_;
        T max_parameter_;

    };

} // splines

#endif //SPLINES_BEZIERQPOPERATIONS_H
