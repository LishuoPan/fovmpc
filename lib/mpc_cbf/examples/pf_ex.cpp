//
// Created by lishuo on 9/21/24.
//

#include <mpc_cbf/optimization/PiecewiseBezierMPCCBFQPGenerator.h>
#include <model/DoubleIntegratorXYYaw.h>
#include <mpc_cbf/controller/BezierIMPCCBF.h>
#include <math/collision_shapes/AlignedBoxCollisionShape.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <particle_filter/detail/particle_filter.h>
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/utils.h>
#include <tf2_ros/transform_broadcaster.h>

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
    // tf2::Quaternion q(robot(3), robot(4), robot(5), robot(6));
    // tf2::Matrix3x3 m(q);
    // double roll, pitch, yaw;
    // m.getRPY(roll, pitch, yaw);
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

int main() {
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

    using PiecewiseBezierParams = mpc::PiecewiseBezierParams<double, DIM>;
    using MPCParams = mpc::MPCParams<double>;
    using FoVCBFParams = cbf::FoVCBFParams<double>;
    using BezierMPCCBFParams = mpc_cbf::PiecewiseBezierMPCCBFQPOperations<double, DIM>::Params;
    using SingleParameterPiecewiseCurve = splines::SingleParameterPiecewiseCurve<double, DIM>;

    using json = nlohmann::json;
    using ParticleFilter = pf::ParticleFilter;

    // load experiment config
    // std::string experiment_config_filename = "../../../config/config.json";
    std::string experiment_config_filename = "/home/user/catkin_ws/src/fovmpc/config/config.json";
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
    double fov_beta = double(experiment_config_json["fov_cbf_params"]["beta"]) * M_PI / 180.0;
    double fov_Ds = experiment_config_json["fov_cbf_params"]["Ds"];
    double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];

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
    // create params
    PiecewiseBezierParams piecewise_bezier_params = {num_pieces, num_control_points, piece_max_parameter};
    MPCParams mpc_params = {h, Ts, k_hor, {w_pos_err, w_u_eff, spd_f}, {p_min, p_max, a_min, a_max}};
    FoVCBFParams fov_cbf_params = {fov_beta, fov_Ds, fov_Rs};

    // json for record
    std::string JSON_FILENAME = "/home/user/catkin_ws/src/fovmpc/tools/PFCBFXYYawStates.json";
    json states;
    states["dt"] = h;
    states["Ts"] = Ts;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h);
    std::shared_ptr<DoubleIntegratorXYYaw> exe_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(Ts);
    StatePropagator exe_A0 = exe_model_ptr->get_A0(int(h/Ts));
    StatePropagator exe_Lambda = exe_model_ptr->get_lambda(int(h/Ts));
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs);
    // init bezier mpc-cbf
    uint64_t bezier_continuity_upto_degree = 4;
    int impc_iter = 2;
    int num_neighbors = experiment_config_json["tasks"]["so"].size() - 1;
    std::cout << "neighbor size: " << num_neighbors << "\n";
    bool slack_mode = true;
    BezierMPCCBFParams bezier_impc_cbf_params = {piecewise_bezier_params, mpc_params, fov_cbf_params};
    BezierIMPCCBF bezier_impc_cbf(bezier_impc_cbf_params, pred_model_ptr, fov_cbf, bezier_continuity_upto_degree, aligned_box_collision_shape_ptr, impc_iter, num_neighbors, slack_mode);

    

    // main loop
    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];
    
    // filters
    std::vector<ParticleFilter> filters(num_robots);
    int num_particles = 100;
    Matrix initCov = 1.0*Eigen::MatrixXd::Identity(DIM, DIM);
    Matrix processCov = 0.25*Eigen::MatrixXd::Identity(DIM, DIM);
    Matrix measCov = 0.05*Eigen::MatrixXd::Identity(DIM, DIM);
    for (size_t i = 0; i < num_robots; ++i) {
        // load init states
        State init_state;
        init_state.pos_ << so_json[i][0], so_json[i][1], so_json[i][2];
        init_state.vel_ << VectorDIM::Zero();
        init_states.push_back(init_state);
        Vector current_state(6);
        current_state << init_state.pos_, init_state.vel_;
        current_states.push_back(current_state);
        // load target positions
        VectorDIM target_pos;
        target_pos << sf_json[i][0], sf_json[i][1], sf_json[i][2];
        target_positions.push_back(target_pos);
        // init particle filter
        filters[i].init(num_particles, init_state.pos_, initCov, processCov, measCov);
    }

    double sim_runtime = 20;
    double sim_t = 0;
    int loop_idx = 0;
    while (sim_t < sim_runtime) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
//            if (robot_idx == 1) {
//                continue;
//            }
            std::vector<VectorDIM> other_robot_positions;
            std::vector<Matrix> other_robot_covs;
            for (int j = 0; j < num_robots; ++j) {
                if (j==robot_idx) {
                    continue;
                }
                // other_robot_positions.push_back(init_states.at(j).pos_);

                // Matrix other_robot_cov(DIM, DIM);
                // other_robot_cov << 0.1, 0, 0,
                //                    0, 0.1, 0,
                //                    0, 0, 0.1;
                // other_robot_covs.push_back(other_robot_cov);
                filters[j].predict();
                Vector weights = filters[j].getWeights();
                Matrix samples = filters[j].getParticles();
                for (int s = 0; s < num_particles; ++s)
                {
                    if (insideFOV(init_states.at(robot_idx).pos_, samples.col(s), fov_beta, fov_Rs))
                    {
                        weights[s] /= 10.0;
                    }
                }
                filters[j].setWeights(weights);

                // emulate vision
                if (insideFOV(init_states.at(robot_idx).pos_, init_states.at(j).pos_, fov_beta, fov_Rs))
                {
                    filters[j].update(init_states.at(j).pos_);
                }

                filters[j].resample();
                filters[j].estimateState();

                // update variables
                VectorDIM estimate = filters[j].getState();
                Matrix cov(DIM, DIM);
                cov = filters[j].getDistribution();
                other_robot_positions.push_back(estimate);
                other_robot_covs.push_back(cov);
                // log estimate
                states["robots"][std::to_string(robot_idx)]["estimates"][std::to_string(j)].push_back({estimate(0), estimate(1), estimate(2)});
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
                std::cout << "QP failed at ts: " << sim_t << "\n";
            }
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
                Vector x_t_vel = x_vel_pred.middleRows(ts_idx, 3);
                Vector x_t(6);
                x_t << x_t_pos, x_t_vel;
                states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
                ts_idx += 3;
                current_states.at(robot_idx) = x_t;
            }
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
