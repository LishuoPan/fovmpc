//
// Created by lishuo on 9/2/24.
//

#include <model/DoubleIntegratorXYYaw.h>
#include <cbf/detail/cbf.h>
#include <cbf/controller/CBFControl.h>
#include <particle_filter/detail/particle_filter.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <cmath>

constexpr unsigned int DIM = 3U;
using State = model::State<double, DIM>;
using VectorDIM = math::VectorDIM<double, DIM>;
using Vector = math::Vector<double>;
using Matrix = math::Matrix<double>;

math::VectorDIM<double, DIM> proportionalControl(const VectorDIM& current_pos, const VectorDIM& target_pos, const double k_prop) {
    return k_prop * (target_pos - current_pos);
}

math::VectorDIM<double, DIM> criticallyDampedSpringControl(const State& current_state, const VectorDIM& target, const double spring_constant) {
    VectorDIM Fs;
    Fs = spring_constant * target - current_state.pos_;
    VectorDIM Fd;
    Fd = -current_state.vel_ * 2 * sqrt(spring_constant);
    return Fs + Fd;
}

math::VectorDIM<double, DIM> rotateControlInputToBodyFrame(const VectorDIM& control_wf, double radian) {
    math::Matrix<double> R(3,3);
    R << cos(radian), sin(radian), 0,
         -sin(radian), cos(radian), 0,
         0, 0, 1;
    return R * control_wf;
}

math::VectorDIM<double, DIM> rotateControlInputToWorldFrame(const VectorDIM& control_bf, double radian) {
    math::Matrix<double> R(3,3);
    R << cos(radian), sin(radian), 0,
         -sin(radian), cos(radian), 0,
         0, 0, 1;
    return R.transpose() * control_bf;
}

double sigmoid(double x){
    return 1 / (1 + exp(-x));
}

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

double distanceToEllipse(const VectorDIM &robot_pos,
                            const Vector &other_robot_pos,
                            const Matrix &other_robot_cov) {
    if (!isinf(other_robot_cov(0, 0))) {
        Eigen::EigenSolver<Matrix> es(other_robot_cov.block(0, 0, DIM-1, DIM-1));
        Vector eigenvalues = es.eigenvalues().real();
        Matrix eigenvectors = es.eigenvectors().real();

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

        double slope = atan2(-other_robot_pos(1) + robot_pos(1), -other_robot_pos(0) + robot_pos(0));
        double x_n = other_robot_pos(0) + a * cos(slope - theta) * cos(theta) - b * sin(slope - theta) * sin(theta);
        double y_n = other_robot_pos(1) + a * cos(slope - theta) * sin(theta) + b * sin(slope - theta) * cos(theta);

        Vector p_near(2);
        p_near << x_n, y_n;

        double dist = sqrt(pow(p_near(0) - robot_pos(0), 2) + pow(p_near(1) - robot_pos(1), 2));

        if (isnan(dist)) {
            dist = 5.0;
            return dist;
        }

        // Check if robot is inside ellipse
        double d = sqrt(pow(other_robot_pos(0) - robot_pos(0), 2) + pow(other_robot_pos(1) - robot_pos(1), 2));
        double range = sqrt(pow(other_robot_pos(0) - p_near(0), 2) + pow(other_robot_pos(1) - p_near(1), 2));

        if (d < range) {
            return -dist;
        } else {
            return dist;
        }
    }

    return -5.0;
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

int main() {
    using FovCBF = cbf::FovCBF;
    using DoubleIntegratorXYYaw = model::DoubleIntegratorXYYaw<double>;
    using CBFControl = cbf::CBFControl<double, DIM>;
    using json = nlohmann::json;
    using ParticleFilter = pf::ParticleFilter;

    // load experiment config
    // std::string experiment_config_filename = "../../../config/config.json";
    std::string experiment_config_filename = "/home/user/catkin_ws/src/fovmpc/config/config.json";
    std::fstream experiment_config_fc(experiment_config_filename.c_str(), std::ios_base::in);
    json experiment_config_json = json::parse(experiment_config_fc);
    double h = experiment_config_json["mpc_params"]["h"];

    double fov_beta = double(experiment_config_json["fov_cbf_params"]["beta"]) * M_PI / 180.0;
    double fov_Ds = experiment_config_json["fov_cbf_params"]["Ds"];
    double fov_Rs = experiment_config_json["fov_cbf_params"]["Rs"];
    double vmax = 2.0;

    VectorDIM a_min;
    a_min << experiment_config_json["mpc_params"]["physical_limits"]["a_min"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_min"][2];
    VectorDIM a_max;
    a_max << experiment_config_json["mpc_params"]["physical_limits"]["a_max"][0],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][1],
            experiment_config_json["mpc_params"]["physical_limits"]["a_max"][2];
    double VMAX = 1.;
    // json for record
    std::string JSON_FILENAME = "/home/user/catkin_ws/src/fovmpc/tools/CBFXYYawStates.json";
    json states;
    states["dt"] = h;
    // init model
    std::shared_ptr<DoubleIntegratorXYYaw> pred_model_ptr = std::make_shared<DoubleIntegratorXYYaw>(h);
    // init cbf
    std::shared_ptr<FovCBF> fov_cbf = std::make_unique<FovCBF>(fov_beta, fov_Ds, fov_Rs, vmax);

    // load the tasks
    std::vector<State> init_states;
    std::vector<Vector> current_states;
    std::vector<VectorDIM> target_positions;
    size_t num_robots = experiment_config_json["tasks"]["so"].size();
    json so_json = experiment_config_json["tasks"]["so"];
    json sf_json = experiment_config_json["tasks"]["sf"];

    // filters params
    int num_particles = 100;
    Matrix initCov = 1.0*Eigen::MatrixXd::Identity(DIM, DIM);
    Matrix processCov = 0.25*Eigen::MatrixXd::Identity(DIM, DIM);
    Matrix measCov = 0.025*Eigen::MatrixXd::Identity(DIM, DIM);
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
            filters[i][j].init(num_particles, init_states.at(neighbor_id).pos_, initCov, processCov, measCov);
        }
        assert(filters[i].size() == num_robots-1);
    }
    assert(filters.size() == num_robots);

    // control loop
    int loop_idx = 0;
    while (loop_idx < 400) {
        for (int robot_idx = 0; robot_idx < num_robots; ++robot_idx) {
            // estimate other robots positions
            std::vector<VectorDIM> other_robot_positions;
            std::vector<Matrix> other_robot_covs;
            std::vector<double> slack_vars;
            for (int j = 0; j < num_robots-1; ++j)
            {
                size_t neighbor_id = neighbor_ids.at(robot_idx).at(j);
                const VectorDIM& ego_pos = init_states.at(robot_idx).pos_;
                const VectorDIM& neighbor_pos = init_states.at(neighbor_id).pos_;
                ParticleFilter &filter = filters.at(robot_idx).at(j); 
                filter.predict();
                Vector weights = filter.getWeights();
                Matrix samples = filter.getParticles();
                for (int s = 0; s < num_particles; ++s)
                {
                    if (insideFOV(ego_pos, samples.col(s), fov_beta, fov_Rs))
                    {
                        weights[s] /= 10.0;
                    }
                }
                filter.setWeights(weights);

                // emulate vision
                if (insideFOV(ego_pos, neighbor_pos, fov_beta, fov_Rs))
                {
                    filter.update(neighbor_pos);
                }

                filter.resample();
                filter.estimateState();

                // update variables
                VectorDIM estimate = filter.getState();
                Matrix cov(DIM, DIM);
                cov = filter.getDistribution();
                other_robot_positions.push_back(estimate);
                other_robot_covs.push_back(cov);
                // calculate slack var
                double dist_to_ellipse = distanceToEllipse(init_states.at(robot_idx).pos_, estimate, cov);
                // double slack = sigmoid(dist_to_ellipse - 3*fov_Ds);
                double slack = 0.0;
                std::cout << "Neighbor id: " << neighbor_id << " | pos: " << estimate.transpose() << " | real: " << neighbor_pos.transpose() << std::endl;
                slack_vars.push_back(slack);
                // log estimate
                states["robots"][std::to_string(robot_idx)]["estimates_mean"][std::to_string(neighbor_id)].push_back({estimate(0), estimate(1), estimate(2)});
                states["robots"][std::to_string(robot_idx)]["estimates_cov"][std::to_string(neighbor_id)].push_back({cov(0), cov(1), cov(2),
                                                                                                                     cov(3), cov(4), cov(5),
                                                                                                                     cov(6), cov(7), cov(8)});
                // closest point on ellipse
                Vector p_near = closestPointOnEllipse(ego_pos, estimate, cov);
                states["robots"][std::to_string(robot_idx)]["p_near_ellipse"][std::to_string(neighbor_id)].push_back({p_near(0), p_near(1)});

            }
                
            // compute the desired control
            const VectorDIM& target_pos = target_positions.at(robot_idx);
            Vector& current_state = current_states.at(robot_idx);
            // VectorDIM desired_u = proportionalControl(init_states.at(robot_idx).pos_, target_pos, 0.8);
            VectorDIM desired_u = criticallyDampedSpringControl(init_states.at(robot_idx), target_pos, 1.0);
            // cbf control
            CBFControl cbf_control(fov_cbf);
            VectorDIM cbf_u;
            cbf_control.optimizeWithSlackVariables(cbf_u, desired_u, current_state, other_robot_positions, slack_vars, a_min, a_max);
            std::cout << "Desired u: " << desired_u.transpose() << " | optimal u: " << cbf_u.transpose() << std::endl;
            State next_init_state = pred_model_ptr->applyInput(init_states.at(robot_idx), cbf_u);
            init_states.at(robot_idx) = next_init_state;

            current_states.at(robot_idx)(0) = init_states.at(robot_idx).pos_(0);
            current_states.at(robot_idx)(1) = init_states.at(robot_idx).pos_(1);
            current_states.at(robot_idx)(2) = init_states.at(robot_idx).pos_(2);
            current_states.at(robot_idx)(3) = init_states.at(robot_idx).vel_(0);
            current_states.at(robot_idx)(4) = init_states.at(robot_idx).vel_(1);
            current_states.at(robot_idx)(5) = init_states.at(robot_idx).vel_(2);

            // log the state
            Vector x_t(6);
            x_t = current_states.at(robot_idx);
            states["robots"][std::to_string(robot_idx)]["states"].push_back({x_t[0], x_t[1], x_t[2], x_t[3], x_t[4], x_t[5]});
        }
        loop_idx += 1;
    }

    // write states to json
    std::cout << "writing to file " << JSON_FILENAME << ".\n";
    std::ofstream o(JSON_FILENAME, std::ofstream::trunc);
    o << std::setw(4) << states << std::endl;

    return 0;
}