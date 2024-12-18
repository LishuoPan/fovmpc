//
// Created by lishuo on 9/22/24.
//

#include <mpc_cbf/controller/BezierIMPCCBF.h>

namespace mpc_cbf {
    template <typename T, unsigned int DIM>
    BezierIMPCCBF<T, DIM>::BezierIMPCCBF(Params &p, std::shared_ptr<DoubleIntegrator> model_ptr,
                                         std::shared_ptr<FovCBF> fov_cbf_ptr, uint64_t bezier_continuity_upto_degree,
                                         std::shared_ptr<const CollisionShape> collision_shape_ptr,
                                         int impc_iter)
    : bezier_continuity_upto_degree_(bezier_continuity_upto_degree),
    collision_shape_ptr_(collision_shape_ptr),
    impc_iter_(impc_iter) {
        // Initiate the qp_generator
        std::unique_ptr<PiecewiseBezierMPCCBFQPOperations> piecewise_mpc_cbf_operations_ptr =
                std::make_unique<PiecewiseBezierMPCCBFQPOperations>(p, model_ptr, fov_cbf_ptr);
        qp_generator_.addPiecewise(std::move(piecewise_mpc_cbf_operations_ptr));

        // load mpc tuning params
        mpc_tuning_ = p.mpc_params.tuning_;

        T h = p.mpc_params.h_;
        T Ts = p.mpc_params.Ts_;
        assert(Ts <= h);
        assert((h - (int)(h/Ts) * Ts) == 0);
        // control sequence U control point coefficient
        int u_interp = int(h / Ts);
        ts_samples_ = Vector::LinSpaced(u_interp, 0, h - Ts); // [I,1]

        // for position pred
        h_ = p.mpc_params.h_;
        k_hor_ = p.mpc_params.k_hor_;
        h_samples_ = Vector::LinSpaced(k_hor_, 0, (k_hor_-1)*h_);
    }

    template <typename T, unsigned int DIM>
    bool BezierIMPCCBF<T, DIM>::optimize(std::vector<SingleParameterPiecewiseCurve> &result_curves,
                                         const State &current_state,
                                         const std::vector<VectorDIM> &other_robot_positions,
                                         const std::vector<Matrix> &other_robot_covs,
                                         const Vector &ref_positions) {
        bool success = true;
        assert(impc_iter_ >= 1);

        VectorDIM current_pos = current_state.pos_;
        VectorDIM current_vel = current_state.vel_;

        // sort neighbors based on distance
        std::vector<std::pair<VectorDIM, Matrix>> other_pos_covs(other_robot_positions.size());                 // each elem is a tuple (pos, cov)
        for (int i = 0; i < other_pos_covs.size(); ++i)
        {
            other_pos_covs[i] = std::make_pair(other_robot_positions.at(i), other_robot_covs.at(i));
        }
        std::vector<std::pair<VectorDIM, Matrix>> sorted_pos_covs = other_pos_covs;
        std::sort(sorted_pos_covs.begin(), sorted_pos_covs.end(), [this, current_pos](std::pair<VectorDIM, Matrix> a, std::pair<VectorDIM, Matrix> b){
            return compareDist(current_pos, a, b);
            });

        // Define new sorted vectors of positions and covariances
        std::vector<VectorDIM> sorted_robots(other_robot_positions.size());
        std::vector<Matrix> sorted_covs(other_robot_positions.size());
        std::cout << "Sorted robots: \n"; 
        for (int i = 0; i < sorted_robots.size(); ++i)
        {
            sorted_robots[i] = sorted_pos_covs.at(i).first;
            sorted_covs[i] = sorted_pos_covs.at(i).second;
            std::cout << sorted_robots.at(i).transpose() << std::endl;
        }

        // define slack variables weights (exp decay)
        T w_init = 100.0;
        T decay_factor = 0.25;
        Vector weights(sorted_robots.size());
        for (int i = 0; i < weights.size(); ++i)
        {
            weights(i) = w_init * pow(decay_factor, i);
        }

        for (size_t iter = 0; iter < impc_iter_; ++iter) {
            // reset the problem
            qp_generator_.problem().resetProblem();

            size_t num_pieces = qp_generator_.piecewise_mpc_qp_generator_ptr()->numPieces();
            // add the position error cost
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addPositionErrorPenaltyCost(current_state, ref_positions);
//            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalPositionErrorPenaltyCost(ref_positions);
            // minimize the control effort for the curve up to specified degree
            for (size_t d = 1; d <= bezier_continuity_upto_degree_; ++d) {
                qp_generator_.piecewise_mpc_qp_generator_ptr()->addIntegratedSquaredDerivativeCost(d,
                                                                                                   mpc_tuning_.w_u_eff_);
            }

            // add the current state constraints
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 0, current_pos);
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalConstraint(0, 1, current_vel);

            // add the continuity between pieces
            for (size_t piece_index = 0; piece_index < num_pieces - 1; ++piece_index) {
                for (size_t d = 0; d < bezier_continuity_upto_degree_; ++d) {
                    qp_generator_.piecewise_mpc_qp_generator_ptr()->addContinuityConstraint(piece_index, d);
                }
            }

            // add the collision avoidance constraints
            AlignedBox robot_bbox_at_zero = collision_shape_ptr_->boundingBox(VectorDIM::Zero());
            for (size_t i = 0; i < other_robot_positions.size(); ++i) {
                const VectorDIM& other_robot_pos = other_robot_positions.at(i);

                // TODO this is model specific, here the last dimension is orientation data
                VectorDIM current_xy = current_pos; current_xy(DIM-1) = 0;
                VectorDIM other_robot_xy = other_robot_pos; other_robot_xy(DIM-1) = 0;

                Hyperplane hyperplane = separating_hyperplanes::voronoi<T, DIM>(current_xy, other_robot_xy);
                // shifted hyperplane
                Hyperplane shifted_hyperplane = math::shiftHyperplane<T, DIM>(hyperplane, robot_bbox_at_zero);
                qp_generator_.piecewise_mpc_qp_generator_ptr()->addHyperplaneConstraintForPiece(0, shifted_hyperplane); // add the collision avoidance constraint on the first segment
            }

            // cbf constraints
            if (iter == 0) {
                for (size_t i = 0; i < other_robot_positions.size(); ++i) {
                    Vector other_xy(2);
                    other_xy << other_robot_positions.at(i)(0), other_robot_positions.at(i)(1);

                    Matrix other_xy_cov(2, 2);
                    other_xy_cov = other_robot_covs.at(i).block(0,0,2,2);
                    // compute the distance to target ellipse
                    T distance_to_ellipse = distanceToEllipse(current_pos, other_xy, other_xy_cov);
                    T slack_value = sigmoid(distance_to_ellipse - 3*0.3); // TODO pass in param

                    qp_generator_.addSafetyCBFConstraint(current_state, other_xy, slack_value);
                    qp_generator_.addFovLBConstraint(current_state, other_xy, slack_value);
                    qp_generator_.addFovRBConstraint(current_state, other_xy, slack_value);
                }
            } else if (iter > 0 && success) {
                // pred the robot's position in the horizon, use for the CBF constraints
                std::vector<State> pred_states;
                std::vector<T> slack_values;
                int cbf_horizon = 2; // TODO pass in this argument
                for (size_t k = 0; k < cbf_horizon; ++k) {
                    State pred_state;
                    pred_state.pos_ = result_curves.back().eval(h_samples_(k), 0);
                    pred_state.vel_ = result_curves.back().eval(h_samples_(k), 1);
                    pred_states.push_back(pred_state);
                }
                for (size_t i = 0; i < other_robot_positions.size(); ++i) {
                    Vector other_xy(2);
                    other_xy << other_robot_positions.at(i)(0), other_robot_positions.at(i)(1);

                    Matrix other_xy_cov(2, 2);
                    other_xy_cov = other_robot_covs.at(i).block(0,0,2,2);
                    // compute the slack value
                    for (size_t k = 0; k < cbf_horizon; ++k) {
                        T distance_to_ellipse = distanceToEllipse(pred_states.at(k).pos_, other_xy, other_xy_cov);
                        slack_values.push_back(sigmoid(distance_to_ellipse - 3*0.3)); // TODO pass in param
                    }
//                    qp_generator_.addPredSafetyCBFConstraints(pred_states, other_xy);
                    qp_generator_.addSafetyCBFConstraint(current_state, other_xy, slack_values.at(0));
                    qp_generator_.addPredFovLBConstraints(pred_states, other_xy, slack_values);
                    qp_generator_.addPredFovRBConstraints(pred_states, other_xy, slack_values);
                }
            }

            // dynamics constraints
            T a_max = 10;
            T v_max = 1;
            VectorDIM a_min_vec = {-a_max, -a_max, -a_max};
            VectorDIM a_max_vec = {a_max, a_max, a_max};
            VectorDIM v_min_vec = {-v_max, -v_max, -v_max};
            VectorDIM v_max_vec = {v_max, v_max, v_max};
            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalBoundConstraints(2, a_min_vec, a_max_vec);
//            qp_generator_.piecewise_mpc_qp_generator_ptr()->addEvalBoundConstraints(1, v_min_vec, v_max_vec);

//            AlignedBox acc_derivative_bbox(a_min_vec, a_max_vec);
//            AlignedBox vel_derivative_bbox(v_min_vec, v_max_vec);
//            qp_generator_.piecewise_mpc_qp_generator_ptr()->addBoundingBoxConstraintAll(acc_derivative_bbox, 2);
//            qp_generator_.piecewise_mpc_qp_generator_ptr()->addBoundingBoxConstraintAll(vel_derivative_bbox, 1);

            // solve QP
            Problem &problem = qp_generator_.problem();
            CPLEXSolver cplex_solver;
            SolveStatus solve_status = cplex_solver.solve(problem);
            if (solve_status == SolveStatus::OPTIMAL) {
                success = true;
            } else {
                success = false;
            }

            // generate bezier
            result_curves.push_back(qp_generator_.piecewise_mpc_qp_generator_ptr()->generateCurveFromSolution());
        }

        return success;
    }

    template <typename T, unsigned int DIM>
    bool BezierIMPCCBF<T, DIM>::compareDist(const VectorDIM& p_current,
                                         const std::pair<VectorDIM, Matrix>& a,
                                         const std::pair<VectorDIM, Matrix>& b) {
        return distanceToEllipse(p_current, a.first, a.second) < distanceToEllipse(p_current, b.first, b.second);
    }

    template <typename T, unsigned int DIM>
    T BezierIMPCCBF<T, DIM>::sigmoid(T x) {
        return 1 / (1 + exp(-2*x));
    }

    template <typename T, unsigned int DIM>
    T BezierIMPCCBF<T, DIM>::distanceToEllipse(const VectorDIM &robot_pos,
                                               const Vector &target_mean,
                                               const Matrix &target_cov) {
        assert(DIM == 3); // DIM other than 3 is not implemented yet.

        if (!isinf(target_cov(0, 0))) {
            Eigen::EigenSolver<Matrix> es(target_cov.block(0, 0, DIM-1, DIM-1));
            Vector eigenvalues = es.eigenvalues().real();
            Matrix eigenvectors = es.eigenvectors().real();

            // s = 4.605 for 90% confidence interval
            // s = 5.991 for 95% confidence interval
            // s = 9.210 for 99% confidence interval
            T s = 4.605;
            T a = sqrt(s * eigenvalues(0)); // major axis
            T b = sqrt(s * eigenvalues(1)); // minor axis

            // a could be smaller than b, so swap them
            if (a < b)
            {
                T temp = a;
                a = b;
                b = temp;
            }

            int m = 0; // higher eigenvalue index
            int l = 1; // lower eigenvalue index
            if (eigenvalues(1) > eigenvalues(0))
            {
                m = 1;
                l = 0;
            }

            T theta = atan2(eigenvectors(1, m), eigenvectors(0, m)); // angle of the major axis wrt positive x-asis (ccw rotation)
            if (theta < 0.0) {
                theta += M_PI;
            } // angle in [0, 2pi]

            T slope = atan2(-target_mean(1) + robot_pos(1), -target_mean(0) + robot_pos(0));
            T x_n = target_mean(0) + a * cos(slope - theta) * cos(theta) - b * sin(slope - theta) * sin(theta);
            T y_n = target_mean(1) + a * cos(slope - theta) * sin(theta) + b * sin(slope - theta) * cos(theta);

            Vector p_near(2);
            p_near << x_n, y_n;

            T dist = sqrt(pow(p_near(0) - robot_pos(0), 2) + pow(p_near(1) - robot_pos(1), 2));

            if (isnan(dist)) {
                dist = 5.0;
                return dist;
            }

            // Check if robot is inside ellipse
            T d = sqrt(pow(target_mean(0) - robot_pos(0), 2) + pow(target_mean(1) - robot_pos(1), 2));
            T range = sqrt(pow(target_mean(0) - p_near(0), 2) + pow(target_mean(1) - p_near(1), 2));

            if (d < range) {
                return -dist;
            } else {
                return dist;
            }
        }

        return -5.0;
    }

    template <typename T, unsigned int DIM>
    void BezierIMPCCBF<T, DIM>::resetProblem() {
        qp_generator_.problem().resetProblem();
    }

    template <typename T, unsigned int DIM>
    typename BezierIMPCCBF<T, DIM>::Vector
    BezierIMPCCBF<T, DIM>::generatorDerivativeControlInputs(uint64_t derivative_degree) {
        Matrix U_basis = qp_generator_.piecewise_mpc_qp_generator_ptr()->piecewise_operations_ptr()->evalSamplingBasisMatrix(ts_samples_, derivative_degree); // [3I, num_piece*dim*num_control_pts]
        Vector control_pts_variables_value = qp_generator_.piecewise_mpc_qp_generator_ptr()->getVariablesValue(); // [num_piece*dim*num_control_pts, 1]
        return U_basis * control_pts_variables_value;
    }

    template class BezierIMPCCBF<double, 3U>;
} // mpc_cbf