//
// Created by lishuo on 9/1/24.
//

#ifndef CBF_CBFCONTROL_H
#define CBF_CBFCONTROL_H

#include <cbf/detail/cbf.h>
#include <cbf/optimization/CBFQPOperations.h>
#include <cbf/optimization/CBFQPGenerator.h>
#include <qpcpp/solvers/CPLEX.h>

namespace cbf {
    template <typename T, unsigned int DIM>
    class CBFControl {
    public:
        using CBFQPOperations = cbf::CBFQPOperations<T, DIM>;
        using CBFQPGenerator = cbf::CBFQPGenerator<T, DIM>;
        using Problem = qpcpp::Problem<T>;
        using CPLEXSolver = qpcpp::CPLEXSolver<T>;
        using SolveStatus = qpcpp::SolveStatus;
        using Vector = math::Vector<T>;
        using VectorDIM = math::VectorDIM<T, DIM>;

        CBFControl(std::shared_ptr<FovCBF> cbf);
        ~CBFControl()=default;
        bool optimize(VectorDIM &cbf_u, const VectorDIM &desired_u,
                      const Vector &state, const Vector &target_state,
                      const VectorDIM& u_min, const VectorDIM& u_max);
        bool optimizeWithSlackVariables(VectorDIM &cbf_u, const VectorDIM &desired_u,
                      const Vector &state, const std::vector<VectorDIM> &other_robots_states, const std::vector<T> &slacks,
                      const VectorDIM& u_min, const VectorDIM& u_max);

    private:
        CBFQPGenerator qp_generator_;

    };

} // cbf

#endif //CBF_CBFCONTROL_H
