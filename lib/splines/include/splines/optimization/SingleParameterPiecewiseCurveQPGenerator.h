//
// Created by lishuo on 2/13/24.
//

#ifndef SPLINES_SINGLEPARAMETERPIECEWISECURVEQPGENERATOR_H
#define SPLINES_SINGLEPARAMETERPIECEWISECURVEQPGENERATOR_H

#include <qpcpp/Problem.h>
#include <math/Helpers.h>
#include <splines/curves/SingleParameterPiecewiseCurve.h>
#include <splines/optimization/SingleParameterCurveQPOperations.h>

namespace splines {
    template <typename T, unsigned int DIM>
    class SingleParameterPiecewiseCurveQPGenerator {
    public:
        // types for QPOperations
        using SingleParameterCurveQPOperations = splines::SingleParameterCurveQPOperations<T, DIM>;
        using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<T, DIM>;
        using CostAddition = typename SingleParameterCurveQPOperations::CostAddition;
        using LinearConstraint = typename SingleParameterCurveQPOperations::LinearConstraint;
        using DecisionVariableBounds =
                typename SingleParameterCurveQPOperations::DecisionVariableBounds;
        // types for math
        using AlignedBox = typename SingleParameterCurveQPOperations::AlignedBox;
        using Hyperplane = typename SingleParameterCurveQPOperations::Hyperplane;
        using VectorDIM = typename SingleParameterCurveQPOperations::VectorDIM;

        // adds operations for the next piece
        void addPiece(std::unique_ptr<SingleParameterCurveQPOperations>&& piece_operations_ptr);

        // number of pieces that this generated contains
        std::size_t numPieces() const;

        // returns the maximum parameter of the piecewise curve. return status is
        // not ok if there is no pieces.
        T max_parameter() const;

        // return a reference to the problem instance generated by this generator
        qpcpp::Problem<T>& problem();

        // generates the piecewise curve from the solution to the problem(). it is
        // assumed that variables of the problem() instance is filled by a solver
        // before a call to this function
        SingleParameterPiecewiseCurve generateCurveFromSolution() const;

        // add
        // lambda*\int_{0}^{max_parameter}||df^{derivative_degree}(t)/dt^{derivative_degree}||^2dt
        // cost to the generated problem
        void addIntegratedSquaredDerivativeCost(uint64_t derivative_degree,
                                                T lambda);

        // add lambda*||df^{derivative_degree}/dt^{deriative_degree}(parameter) -
        // target||^2 cost to the generated qp
        void addEvalCost(T parameter, uint64_t derivative_degree,
                         const VectorDIM& target, T lambda);

        // add df^{derivative_degree}/dt^{derivative_degree}(paramter) = target
        // constraint to the generated qp
        void addEvalConstraint(T parameter, uint64_t derivative_degree,
                               const VectorDIM& target);

        // adds a continuity constraint between pieces with indices piece_idx and
        // piece_idx + 1 on the derivative_degree^th derivative. return status is
        // not ok if piece_idx or piece_idx+1 is out of range
        void addContinuityConstraint(std::size_t piece_idx,
                                     uint64_t derivative_degree);

//        // constraints the given piece to be in the negative side of the given
//        // hyperplane
//        void addHyperplaneConstraintForPiece(std::size_t piece_idx,
//                                             const Hyperplane& hyperplane);

        // constraints the curve to be inside the given bounding box
        void addBoundingBoxConstraint(const AlignedBox& bounding_box);

        /**
         * @brief Adds a constraint that the derivative of the curve is contained in
         * the given box at all points.
         *
         * @param aligned_box box that the derivative is constrained to be inside
         * @param derivative_degree degree of the derivative
         * @return absl::Status status of the operation
         */
        void addBoundingBoxConstraintAll(const AlignedBox& aligned_box,
                                         uint64_t derivative_degree);

        // constraints the given piece to be on the given position
        void addPositionConstraintForPiece(std::size_t piece_idx,
                                           const VectorDIM& position);

        // constraints the given piece to be in the negative side of the given
        // hyperplane
        void addHyperplaneConstraintForPiece(std::size_t piece_idx,
                                             const Hyperplane& hyperplane,
                                             T epsilon = 1e-8);

        // adds given decision_variable_bounds for a piece to the generated qp. if
        // there are already decision variable bounds for a piece, minimum of the
        // upper bounds and the maximum of the lower bounds between old bounds and
        // the new bounds are taken
        void addDecisionVariableBoundsForPiece(std::size_t piece_idx,
                                               const DecisionVariableBounds& decision_variable_bounds);

    private:
        // Container class for piece index and parameter pair. return type of
        // getPieceIndexAndParameter
        class PieceIndexAndParameter {
        public:
            PieceIndexAndParameter(std::size_t piece_idx, T parameter);
            std::size_t piece_idx() const;
            T parameter() const;

        private:
            std::size_t piece_idx_;
            T parameter_;
        };

        // adds given cost_addition for a piece to the generated qp. return
        // status is not OK if internal calls to problem_ instance are faulty,
        // which should not logically happen
        void addCostAdditionForPiece(std::size_t piece_idx,
                                     const CostAddition& cost_addition);

        // adds given linear_constraint for a piece to the generated qp. return
        // status is not OK if internal calls to problem_ instance are faulty, which
        // should not logically happen
        void addLinearConstraintForPiece(std::size_t piece_idx,
                                         const LinearConstraint& linear_constraint);

        // returns the piece index and parameter within the piece with the piece
        // index that the `parameter` argument corresponds to.
        //
        // e.g. if there are 4 pieces with max parameters [1.0,2.0,3.0,4.0] in the
        // piecewise curve, getPieceIndexAndParameter(1.2) returns
        // PieceIndexAndParameter(1, 0.2), because parameter 1.2 corresponds to
        // piece with index 1 and parameter 0.2. getPieceIndexAndParameter(1.0)
        // returns PieceIndexAndParameter(1, 0.0) and getPieceIndexAndParameter(0.0)
        // returns PieceIndexAndParameter(0, 0.0).
        PieceIndexAndParameter getPieceIndexAndParameter(T parameter) const;

        // contains pointers to operations for pieces
        std::vector<std::unique_ptr<SingleParameterCurveQPOperations>> piece_operations_ptrs_;

        // contains cumulative max parameters of pieces.
        // cumulative_max_parameters_[i] gives the total max parameters of pieces
        // with indices [0, ..., i]
        std::vector<T> cumulative_max_parameters_;

        // quadratic program generated by this generator
        qpcpp::Problem<T> problem_;

        // [piece-idx][decision-variable-idx] => Variable<T>*
        std::vector<std::vector<qpcpp::Variable<T>*>> variable_map_;

    };

} // splines

#endif //SPLINES_SINGLEPARAMETERPIECEWISECURVEQPGENERATOR_H
