//
// Created by lishuo on 9/21/24.
//

#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <mpc_cbf/controller/BezierIMPCCBF.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <particle_filter/detail/particle_filter.h>
#include <nlohmann/json.hpp>
#include <cxxopts.hpp>
#include <fstream>
#include <random>

math::Vector<double> closestPointOnEllipse(const math::VectorDIM<double, 3U> &robot_pos,
                         const math::Vector<double> &target_mean,
                         const math::Matrix<double> &target_cov) {
    if (!isinf(target_cov(0, 0))) {
        Eigen::EigenSolver<math::Matrix<double>> es(target_cov.block(0, 0, 3U-1, 3U-1));
        math::Vector<double> eigenvalues = es.eigenvalues().real();
        math::Matrix<double> eigenvectors = es.eigenvectors().real();

        // s = 4.605 for 90% confidence interval
        // s = 5.991 for 95% confidence interval
        // s = 9.210 for 99% confidence interval
        double s = 4.605;
        double a = sqrt(s * eigenvalues(0)); // major axis
        double b = sqrt(s * eigenvalues(1)); // minor axis

        // a could be smaller than b, so swap them
        if (a < b)
        {
            double temp = a;
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

        double theta = atan2(eigenvectors(1, m), eigenvectors(0, m)); // angle of the major axis wrt positive x-asis (ccw rotation)
        if (theta < 0.0) {
            theta += M_PI;
        } // angle in [0, 2pi]

        double slope = atan2(-target_mean(1) + robot_pos(1), -target_mean(0) + robot_pos(0));
        double x_n = target_mean(0) + a * cos(slope - theta) * cos(theta) - b * sin(slope - theta) * sin(theta);
        double y_n = target_mean(1) + a * cos(slope - theta) * sin(theta) + b * sin(slope - theta) * cos(theta);

        math::Vector<double> p_near(2);
        p_near << x_n, y_n;
        return p_near;
//        double dist = sqrt(pow(p_near(0) - robot_pos(0), 2) + pow(p_near(1) - robot_pos(1), 2));
//
//        if (isnan(dist)) {
//            dist = 5.0;
//            return dist;
//        }
//
//        // Check if robot is inside ellipse
//        double d = sqrt(pow(target_mean(0) - robot_pos(0), 2) + pow(target_mean(1) - robot_pos(1), 2));
//        double range = sqrt(pow(target_mean(0) - p_near(0), 2) + pow(target_mean(1) - p_near(1), 2));
//
//        if (d < range) {
//            return -dist;
//        } else {
//            return dist;
//        }
    } else {
        math::Vector<double> p_near(2);
        p_near << 0, 0;
        return p_near;
    }

//    return -5.0;
}


double convertYawInRange(double yaw) {
    assert(yaw < 2 * M_PI && yaw > -2 * M_PI);
    if (yaw > M_PI) {
        return yaw - 2 * M_PI;
    } else if (yaw < -M_PI) {
        return yaw + 2 * M_PI;
    } else {
        return yaw;
    }
}

bool insideFOV(Eigen::VectorXd robot, Eigen::VectorXd target, double fov, double range)
{
    double yaw = robot(2);

    Eigen::Matrix3d R;
    R << cos(yaw), sin(yaw), 0.0,
            -sin(yaw), cos(yaw), 0.0,
            0.0, 0.0, 1.0;
    Eigen::VectorXd t_local = R.block<2,2>(0,0) * (target.head(2) - robot.head(2));
    double dist = t_local.norm();
    double angle = abs(atan2(t_local(1), t_local(0)));
    if (angle <= 0.5*fov && dist <= range)
    {
        return true;
    } else
    {
        return false;
    }
}

math::VectorDIM<double, 3U> convertToClosestYaw(model::State<double, 3U>& state, const math::VectorDIM<double, 3U>& goal) {
    using VectorDIM = math::VectorDIM<double, 3U>;
    using Vector = math::Vector<double>;

//    std::cout << "start get the current_yaw\n";
    double current_yaw = state.pos_(2);
    // generate the candidate desire yaws
    Vector candidate_yaws(3);

//    std::cout << "start build candidate yaw\n";
//    std::cout << goal_(2) << "\n";
//    std::cout << M_PI << "\n";
    candidate_yaws << goal(2), goal(2) + 2 * M_PI, goal(2) - 2 * M_PI;

//    std::cout << "compute the offset\n";
    Vector candidate_yaws_offset(3);
    candidate_yaws_offset << std::abs(candidate_yaws(0) - current_yaw), std::abs(candidate_yaws(1) - current_yaw), std::abs(candidate_yaws(2) - current_yaw);


//    std::cout << "start to find the argmin\n";
    // Find the index of the minimum element
    int argmin_index = 0;
    double min_value = candidate_yaws_offset(0);

    for (int i = 1; i < candidate_yaws_offset.size(); ++i) {
        if (candidate_yaws_offset(i) < min_value) {
            min_value = candidate_yaws_offset(i);
            argmin_index = i;
        }
    }

    VectorDIM converted_goal;
    converted_goal << goal(0), goal(1), candidate_yaws(argmin_index);
//    std::cout << "converted_goal: " << converted_goal.transpose() << "\n";
    return converted_goal;

}

math::Vector<double> addRandomNoise(const math::Vector<double>& xt, double pos_std, double vel_std, unsigned int dim=3U) {
    std::random_device rd;
    std::mt19937 gen(rd());
    math::Vector<double> result_xt = xt;
    double sample = 0.0;
    std::normal_distribution<double> distribution_position(0.0, pos_std);
    std::normal_distribution<double> distribution_velocity(0.0, vel_std);
    for (int i = 0; i < 2*dim; ++i) {
        if (i < dim) {
            sample = distribution_position(gen);
        } else {
            sample = distribution_velocity(gen);
        }
        result_xt(i) += sample;
    }
    return result_xt;
}

int main(int argc, char* argv[]) {
    constexpr unsigned int DIM = 3U;
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using BezierIMPCCBF = mpc_cbf::BezierIMPCCBF<double, DIM>;
    using State = model::State<double, DIM>;
    using StatePropagator = model::StatePropagator<double>;
    using VectorDIM = math::VectorDIM<double, DIM>;
    using Vector = math::Vector<double>;
    using Matrix = math::Matrix<double>;
    using AlignedBox = math::AlignedBox<double, DIM>;
    using AlignedBoxCollisionShape = math::AlignedBoxCollisionShape<double, DIM>;
    using ParticleFilter = pf::ParticleFilter;

    using PiecewiseBezierParams = mpc::PiecewiseBezierParams<double, DIM>;
    using MPCParams = mpc::MPCParams<double>;
    using FoVCBFParams = cbf::FoVCBFParams<double>;
    using BezierMPCCBFParams = mpc_cbf::PiecewiseBezierMPCCBFQPOperations<double, DIM>::Params;
    using IMPCParams = mpc_cbf::BezierIMPCCBF<double, DIM>::IMPCParams;
    using IMPCCBFParams = mpc_cbf::BezierIMPCCBF<double, DIM>::Params;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;

    using json = nlohmann::json;

    // cxxopt
    cxxopts::Options options(
            "simulation",
            "fovmpc simulation");
    options.add_options()
            ("instance_type", "instance type for simulations",
             cxxopts::value<std::string>()->default_value("circle"))
            ("num_robots", "number of robots in the simulation",
             cxxopts::value<int>()->default_value(std::to_string(2)))
            ("fov", "fov degree",
             cxxopts::value<int>()->default_value(std::to_string(120)))
            ("slack_mode", "bool to turn on and off the slack mode",
             cxxopts::value<int>()->default_value(std::to_string(1)))
            ("slack_decay", "slack variable cost decay rate",
             cxxopts::value<double>()->default_value(std::to_string(0.1)))
            ("write_filename", "write to json filename",
             cxxopts::value<std::string>()->default_value("../../../experiments/instances/results/circle2States.json"));
    auto option_parse = options.parse(argc, argv);

    // load experiment config
    std::cout << "loading experiment settings...\n";
//    std::string experiment_config_filename = "../../../config/config.json";
    std::string instance_type = option_parse["instance_type"].as<std::string>();
    std::string instance_path = instance_type+"_instances/";
    const int num_robots_parse = option_parse["num_robots"].as<int>();
    const int fov_beta_parse = option_parse["fov"].as<int>();
    std::string experiment_config_filename = "../../../experiments/instances/"+instance_path+instance_type+std::to_string(num_robots_parse)+"_config.json";
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);
    // piecewise bezier params
    size_t num_pieces = experiment_config_json["bezier_params"]["num_pieces"];
    size_t num_control_points = experiment_config_json["bezier_params"]["num_control_points"];
    double piece_max_parameter = experiment_config_json["bezier_params"]["piece_max_parameter"];
    // mpc params
    double h = experiment_config_json["mpc_params"]["h"];
    double Ts = experiment_config_json["mpc_params"]["Ts"];
    int k_hor = experiment_config_json["mpc_params"]["k_hor"];
    double w_pos_err = experiment_config_json["mpc_params"]["mpc_tuning"]["w_pos_err"];
    double w_u_eff = experiment_config_json["mpc_params"]["mpc_tuning"]["w_u_eff"];
    int spd_f = experiment_config_json["mpc_params"]["mpc_tuning"]["spd_f"];
    Vector p_min = Vector::Zero(2);
    p_min << experiment_config_json["mpc_params"]["physical_limits"]["p_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["p_min"][1];
    Vector p_max = Vector::Zero(2);
    p_max << experiment_config_json["mpc_params"]["physical_limits"]["p_max"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["p_max"][1];
    // fov cbf params
    double fov_beta = double(fov_beta_parse) * M_PI / 180.0;
    std::cout << "fov_beta: " << double(fov_beta_parse) << "\n";
    double fov_Ds = experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0];
    double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];

    // robot physical params
    VectorDIM v_min;
    v_min << experiment_config_json["mpc_params"]["physical_limits"]["v_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["v_min"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["v_min"][2];
    VectorDIM v_max;
    v_max << experiment_config_json["mpc_params"]["physical_limits"]["v_max"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["v_max"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["v_max"][2];
    VectorDIM a_min;
    a_min << experiment_config_json["mpc_params"]["physical_limits"]["a_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][2];
    VectorDIM a_max;
    a_max << experiment_config_json["mpc_params"]["physical_limits"]["a_max"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][2];

    VectorDIM aligned_box_collision_vec;
    aligned_box_collision_vec <<
                              experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][0],
            experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][1],
            experiment_config_json["robot_params"]["collision_shape"]["aligned_box"][2];
    AlignedBox robot_bbox_at_zero = {-aligned_box_collision_vec, aligned_box_collision_vec};
    std::shared_ptr<const AlignedBoxCollisionShape> aligned_box_collision_shape_ptr =
            std::make_shared<const AlignedBoxCollisionShape>(robot_bbox_at_zero);

    double pos_std = experiment_config_json["mpc_params"]["physical_limits"]["pos_std"];
    double vel_std = experiment_config_json["mpc_params"]["physical_limits"]["vel_std"];

    // create params
    PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points, piece_max_parameter};
    MPCParams mpc_params = {h, Ts, k_hor, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, v_min, v_max, a_min, a_max}};
    FoVCBFParams fov_cbf_params = {fov_beta, fov_Ds, fov_Rs};

    // filter params
    int num_particles = 100;
    Matrix initCov = 1.0*Eigen::MatrixXd::Identity(DIM-1, DIM-1);
    Matrix processCov = 0.25*Eigen::MatrixXd::Identity(DIM-1, DIM-1);
    Matrix measCov = 0.05*Eigen::MatrixXd::Identity(DIM-1, DIM-1);

    // json for record
    std::string JSON_FILENAME = option_parse["write_filename"].as<std::string>();
    json states;
    states["dt"] = h;
    states["Ts"] = Ts;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h);
    std::shared_ptr<DoubleIntegratorXYYaw> exe_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts);
    StatePropagator exe_A0 = exe_model_ptr->get_A0(int(h/Ts));
    StatePropagator exe_Lambda = exe_model_ptr->get_lambda(int(h/Ts));
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs, v_min, v_max);
    // init bezier mpc-cbf
    uint64_t bezier_continuity_upto_degree = 4;
    int num_neighbors = experiment_config_json["tasks"]["so"].size() - 1;
    std::cout << "neighbor size: " << num_neighbors << "\n";
    BezierMPCCBFParams bezier_mpc_cbf_params = {piecewise_bezier_params, mpc_params, fov_cbf_params};
    int cbf_horizon = 2;
    int impc_iter = 2;
    double slack_cost = 1000;
    double slack_decay_rate = option_parse["slack_decay"].as<double>();
    std::cout << "slack_decay_rate: " << slack_decay_rate << "\n";
    bool slack_mode = (bool)option_parse["slack_mode"].as<int>();
    std::cout << "slack mode is on: " << slack_mode << "\n";
    IMPCParams impc_params = {cbf_horizon, impc_iter, slack_cost, slack_decay_rate, slack_mode};
    IMPCCBFParams impc_cbf_params = {bezier_mpc_cbf_params, impc_params};
    BezierIMPCCBF bezier_impc_cbf(impc_cbf_params, pred_model_ptr, fov_cbf, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr, num_neighbors);

    // main loop
    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    assert(num_robots_parse==num_robots);
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    for (size_t i = 0; i < num_robots; ++i) {
        // load init states
        State init_state;
        init_state.pos_ << so_json[i][0], so_json[i][1], convertYawInRange(so_json[i][2]);
        init_state.vel_ << VectorDIM::Zero();
        init_states.push_back(init_state);
        Vector current_state(6);
        current_state << init_state.pos_, init_state.vel_;
        current_states.push_back(current_state);
        // load target positions
        VectorDIM target_pos;
        target_pos << sf_json[i][0], sf_json[i][1], sf_json[i][2];
        target_positions.push_back(target_pos);
    }
    // init filters
    std::vector<std::vector<ParticleFilter>> filters;
    std::vector<std::vector<size_t>> neighbor_ids(num_robots, std::vector<size_t>(num_robots-1));
    for (size_t i = 0; i < num_robots; ++i) {
        size_t neighbor_idx = 0;
        for (size_t j = 0; j < num_robots; ++j) {
            if (i == j) {
                continue;
            }
            neighbor_ids[i][neighbor_idx] = j;
            neighbor_idx += 1;
        }
    }

    for (size_t i = 0; i < num_robots; ++i){
        filters.emplace_back(num_robots-1);
        for (size_t j = 0; j < num_robots-1; ++j) {
            size_t neighbor_id = neighbor_ids.at(i).at(j);
            // init particle filter
            Vector neighbor_xy(DIM-1);
            neighbor_xy << init_states.at(neighbor_id).pos_(0), init_states.at(neighbor_id).pos_(1);
            filters[i][j].init(num_particles, neighbor_xy, initCov, processCov, measCov);
        }
        assert(filters[i].size() == num_robots-1);
    }
    assert(filters.size() == num_robots);

    // planning results
    std::vector<bool> init_success(num_robots, false);
    std::vector<std::shared_ptr<SingleParameterPiecewiseCurve>> pred_traj_ptrs(num_robots);
    std::vector<double> traj_eval_ts(num_robots, 0);
    double pred_horizon = num_pieces * piece_max_parameter;

    double sim_runtime = 40;
    double sim_t = 0;
    int loop_idx = 0;
    while (sim_t < sim_runtime) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
//            if (robot_idx == 1) {
//                continue;
//            }
            std::vector<VectorDIM> other_robot_positions;
            std::vector<Matrix> other_robot_covs;
            for (int j = 0; j < num_robots-1; ++j) {
                size_t neighbor_id = neighbor_ids.at(robot_idx).at(j);
                const VectorDIM& ego_pos = init_states.at(robot_idx).pos_;
                const VectorDIM& neighbor_pos = init_states.at(neighbor_id).pos_;
                ParticleFilter &filter = filters.at(robot_idx).at(j);
                // update filter estimate
                filter.predict();
                Vector weights = filter.getWeights();
                Matrix samples = filter.getParticles();
                for (int s = 0; s < num_particles; ++s) {
                    if (insideFOV(ego_pos, samples.col(s), fov_beta, fov_Rs)) {
                        weights[s] /= 10.0;
                    }
                }
                filter.setWeights(weights);

                // emulate vision
                if (insideFOV(ego_pos, neighbor_pos, fov_beta, fov_Rs)) {
                    Vector neighbor_xy(DIM-1);
                    neighbor_xy << neighbor_pos(0), neighbor_pos(1);
                    filter.update(neighbor_xy);
                }

                filter.resample();
                filter.estimateState();

                // update variables
                Vector estimate = filter.getState();
                VectorDIM extended_estimate;
                extended_estimate << estimate(0), estimate(1), 0;
                Matrix cov(DIM-1, DIM-1);
                cov = filter.getDistribution();
                other_robot_positions.push_back(extended_estimate);
                Matrix extended_cov(DIM, DIM);
                extended_cov << cov(0), cov(1), 0,
                                cov(2), cov(3), 0,
                                0,      0,      0;
                other_robot_covs.push_back(extended_cov);
                // log estimate
                states["robots"][std::to_string(robot_idx)]["estimates_mean"][std::to_string(neighbor_id)].push_back({extended_estimate(0), extended_estimate(1), extended_estimate(2)});
                states["robots"][std::to_string(robot_idx)]["estimates_cov"][std::to_string(neighbor_id)].push_back({extended_cov(0), extended_cov(1), extended_cov(2),
                                                                                                                     extended_cov(3), extended_cov(4), extended_cov(5),
                                                                                                                     extended_cov(6), extended_cov(7), extended_cov(8)});

                // for the closest point on ellipse visualization
                Vector p_near = closestPointOnEllipse(ego_pos, estimate, cov);
                states["robots"][std::to_string(robot_idx)]["p_near_ellipse"][std::to_string(neighbor_id)].push_back({p_near(0), p_near(1)});

//                // for debug: fixed estimate
//                other_robot_positions.push_back(init_states.at(neighbor_id).pos_);
//
//                Matrix other_robot_cov(DIM, DIM);
//                other_robot_cov << 0.1, 0, 0,
//                                   0, 0.1, 0,
//                                   0, 0, 0.1;
//                other_robot_covs.push_back(other_robot_cov);
//                states["robots"][std::to_string(robot_idx)]["estimates_mean"][std::to_string(neighbor_id)].push_back({init_states.at(neighbor_id).pos_(0), init_states.at(neighbor_id).pos_(1), init_states.at(neighbor_id).pos_(2)});
//                states["robots"][std::to_string(robot_idx)]["estimates_cov"][std::to_string(neighbor_id)].push_back({other_robot_cov(0), other_robot_cov(1), other_robot_cov(2),
//                                                                                                                     other_robot_cov(3), other_robot_cov(4), other_robot_cov(5),
//                                                                                                                     other_robot_cov(6), other_robot_cov(7), other_robot_cov(8)});
            }
//            BezierIMPCCBF bezier_impc_cbf(bezier_impc_cbf_params, pred_model_ptr, fov_cbf, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr, impc_iter);
            bezier_impc_cbf.resetProblem();
            Vector ref_positions(DIM*k_hor);
            // static target reference
            VectorDIM converted_target_position = convertToClosestYaw(init_states.at(robot_idx), target_positions.at(robot_idx));
            ref_positions = converted_target_position.replicate(k_hor, 1);
//            ref_positions = target_positions.at(robot_idx).replicate(k_hor, 1);

//            std::cout << "ref_positions shape: (" << ref_positions.rows() << ", " << ref_positions.cols() << ")\n";
//            std::cout << "ref_positions: " << ref_positions.transpose() << "\n";
            std::vector<SingleParameterPiecewiseCurve> trajs;
            bool success = bezier_impc_cbf.optimize(trajs, init_states.at(robot_idx),
                                                    other_robot_positions, other_robot_covs,
                                                    ref_positions);
            if (!success) {
                std::cout << "QP failed at ts: " << sim_t;
//                if (!trajs.empty()) {
//                    std::cout << "; But initial planning is successful...";
//                    pred_traj_ptrs.at(robot_idx) = std::make_shared<SingleParameterPiecewiseCurve>(std::move(trajs.back()));
//                    traj_eval_ts.at(robot_idx) = 0;
//                }
                states["robots"][std::to_string(robot_idx)]["QP_success"].push_back(0);
                std::cout << "\n";
            }

//            /* continuous control
            if (success) {
                init_success.at(robot_idx) = true;
//                SingleParameterPiecewiseCurve optim_traj = trajs.back();
                pred_traj_ptrs.at(robot_idx) = std::make_shared<SingleParameterPiecewiseCurve>(std::move(trajs.back()));
                traj_eval_ts.at(robot_idx) = 0;
                states["robots"][std::to_string(robot_idx)]["QP_success"].push_back(1);
            }

//            // log down the prediction
//            double t = 0;
//            while (t <= pred_horizon) {
////                for (size_t traj_index = 0; traj_index < trajs.size(); ++traj_index) {
//                VectorDIM pred_pos;
////                    if (success) {
////                        pred_pos = trajs.at(traj_index).eval(t, 0);
////                    } else {
//                double pred_t = traj_eval_ts.at(robot_idx) + t;
//                if (pred_t > pred_traj_ptrs.at(robot_idx)->max_parameter()) {
//                    pred_t = pred_traj_ptrs.at(robot_idx)->max_parameter();
//                }
//                pred_pos = pred_traj_ptrs.at(robot_idx)->eval(pred_t, 0);
////                    }
//                states["robots"][std::to_string(robot_idx)]["pred_curve"][loop_idx][0].push_back({pred_pos(0), pred_pos(1), pred_pos(2)});
////                }
//                t += 0.05;
//            }

            // log down the trajectory
            double eval_t = 0;
            for (int t_idx = 1; t_idx <= int(h / Ts); ++t_idx) {
                if (init_success.at(robot_idx)) {
                    eval_t = traj_eval_ts.at(robot_idx) + Ts*t_idx;
                    if (eval_t > pred_traj_ptrs.at(robot_idx)->max_parameter()) {
                        eval_t = pred_traj_ptrs.at(robot_idx)->max_parameter();
                    }

                    Vector x_t_pos = pred_traj_ptrs.at(robot_idx)->eval(eval_t, 0);
                    x_t_pos(2) = convertYawInRange(x_t_pos(2));
                    Vector x_t_vel = pred_traj_ptrs.at(robot_idx)->eval(eval_t, 1);
                    Vector x_t(6);
                    x_t << x_t_pos, x_t_vel;
                    x_t = addRandomNoise(x_t, pos_std, vel_std);
                    states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
                    current_states.at(robot_idx) = x_t;
                } else {
                    Vector x_t(6);
                    x_t = current_states.at(robot_idx);
                    x_t = addRandomNoise(x_t, pos_std, vel_std);
                    states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
                    current_states.at(robot_idx) = x_t;
                }
            }
            traj_eval_ts.at(robot_idx) = eval_t;
//             */

            /* discrete control
            Vector U = bezier_impc_cbf.generatorDerivativeControlInputs(2);
//            std::cout << "curve eval: " << traj.eval(1.5, 0) << "\n";

            // log down the optimized curve prediction
            double t = 0;
            while (t <= 1.5) {
                for (size_t traj_index = 0; traj_index < trajs.size(); ++traj_index) {
                    VectorDIM pred_pos = trajs.at(traj_index).eval(t, 0);
                    states["robots"][std::to_string(robot_idx)]["pred_curve"][loop_idx][traj_index].push_back({pred_pos(0), pred_pos(1), pred_pos(2)});
                }
                t += 0.05;
            }
            //

            Matrix x_pos_pred = exe_A0.pos_ * current_states.at(robot_idx) + exe_Lambda.pos_ * U;
            Matrix x_vel_pred = exe_A0.vel_ * current_states.at(robot_idx) + exe_Lambda.vel_ * U;
            assert(int(h/Ts)*DIM == x_pos_pred.rows());
            int ts_idx = 0;
            for (size_t i = 0; i < int(h / Ts); ++i) {
                Vector x_t_pos = x_pos_pred.middleRows(ts_idx, 3);
                x_t_pos(2) = convertYawInRange(x_t_pos(2));
                Vector x_t_vel = x_vel_pred.middleRows(ts_idx, 3);
                Vector x_t(6);
                x_t << x_t_pos, x_t_vel;
                states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
                ts_idx += 3;
                current_states.at(robot_idx) = x_t;
            }
             */

        }
        // update the init_states
        for (size_t robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            init_states.at(robot_idx).pos_(0) = current_states.at(robot_idx)(0);
            init_states.at(robot_idx).pos_(1) = current_states.at(robot_idx)(1);
            init_states.at(robot_idx).pos_(2) = convertYawInRange(current_states.at(robot_idx)(2));
//            init_states.at(robot_idx).pos_(2) = current_states.at(robot_idx)(2);
            init_states.at(robot_idx).vel_(0) = current_states.at(robot_idx)(3);
            init_states.at(robot_idx).vel_(1) = current_states.at(robot_idx)(4);
            init_states.at(robot_idx).vel_(2) = current_states.at(robot_idx)(5);
        }
        sim_t += h;
        loop_idx += 1;
    }

    // write states to json
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;
    return 0;
}
