import numpy as np

from Metrics import *
from ComputeCI import *
import colorsys

# from src.fovmpc.experiments.metrics.ComputeCI import percentile_compute_with_inf, percentile_compute

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['DejaVu Serif']
QP_SUCCESS_METRIC = True

def generate_rgb_colors(num_colors):
    output = []
    num_colors += 1 # to avoid the first color
    for index in range(1, num_colors):
        incremented_value = 1.0 * index / num_colors
        output.append(colorsys.hsv_to_rgb(incremented_value, 0.5, 0.7))
    return np.asarray(output)

instance_type = "circle"
# instance_type = "formation"
config_path = "../instances/"+instance_type+"_instances/"
result_path = "/media/lishuo/ssd/RSS2025_results/log"
date = "03222025"
num_exp = 15
min_r = 2
max_r = 8
min_fov=120
max_fov=240
fov_step=120
min_slack_decay=0.2
max_slack_decay=0.3
slack_decay_step=0.1
default_slack_decay=0.2
decay_exp_fovs = [120]
num_robot = range(min_r, max_r+1)
fovs = range(min_fov, max_fov+fov_step, fov_step)
# fovs = []
slack_decays = [1.0]
print("slack decays: ", slack_decays)

experiment_key=[]
experiment_slack_decay_key=[]
experiment_fov_key=[]

for fov in fovs:
    experiment_key.append("SlackOff, "+r"\beta_{H}="+str(fov))
    experiment_slack_decay_key.append(default_slack_decay)
    experiment_fov_key.append(fov)
for fov in fovs:
    experiment_key.append("SlackOn, "+r"\beta_{H}="+str(fov)+r", \gamma_{s}="+str(default_slack_decay))
    experiment_slack_decay_key.append(default_slack_decay)
    experiment_fov_key.append(fov)
# for fov in decay_exp_fovs:
#     for decay in slack_decays:
#         experiment_key.append(r"SlackOn \beta_{H}="+str(fov)+r", \gamma_{s}="+str(decay))
#         experiment_slack_decay_key.append(decay)
#         experiment_fov_key.append(fov)

success_rate_dict = {}
makespan_dict = {}
num_neighbor_in_fov_dict = {}
percent_neighbor_in_fov_dict = {}
QP_success_rate_dict = {}

success_rate_ci_dict={}
success_rate_percentile_dict={}
makespan_ci_dict={}
makespan_percentile_dict={}
num_neighbor_in_fov_ci_dict={}
num_neighbor_in_fov_percentile_dict={}
percent_neighbor_in_fov_ci_dict={}
percent_neighbor_in_fov_percentile_dict={}
QP_success_rate_ci_dict = {}
QP_success_rate_percentile_dict = {}
for key in experiment_key:
    success_rate_dict[key] = [[] for _ in num_robot]  # [entry,], averaged by M, M is the sample
    makespan_dict[key] = [[] for _ in num_robot]  # [entry,], averaged by M, M is the sample
    num_neighbor_in_fov_dict[key] = [[] for _ in num_robot]  # [entry, M]
    percent_neighbor_in_fov_dict[key] = [[] for _ in num_robot]  # [entry, M]
    QP_success_rate_dict[key] = [[] for _ in num_robot]  # [entry, M]

    success_rate_ci_dict[key] = {}
    success_rate_percentile_dict[key]={}
    makespan_ci_dict[key] = {}
    makespan_percentile_dict[key]={}
    num_neighbor_in_fov_ci_dict[key] = {}
    num_neighbor_in_fov_percentile_dict[key]={}
    percent_neighbor_in_fov_ci_dict[key] = {}
    percent_neighbor_in_fov_percentile_dict[key]={}
    QP_success_rate_ci_dict[key] = {}
    QP_success_rate_percentile_dict[key] = {}

# pre-load
for i in range(len(experiment_key)):
    is_slack_off = "SlackOff" in experiment_key[i]
    key = experiment_key[i]
    fov = experiment_fov_key[i]
    slack_decay = experiment_slack_decay_key[i]
    print("key:", key)
    for r_idx, num_r in enumerate(num_robot):
        print("r_idx: ", r_idx, "num_r: ", num_r)
        for exp_idx in range(num_exp):
            if is_slack_off:
                state_json = result_path+date+"/"+instance_type+str(num_r)+"_fov"+str(fov)+"_slack_off"+"_States_"+str(exp_idx)+".json"
            else:
                state_json = result_path+date+"/"+instance_type+str(num_r)+"_fov"+str(fov)+"_decay"+str(slack_decay)+"_States_"+str(exp_idx)+".json"
            config_json = config_path+instance_type+str(num_r)+"_config.json"
            f = open(state_json)
print("pass pre-load, all data exist...")

for i in range(len(experiment_key)):
    is_slack_off = "SlackOff" in experiment_key[i]
    key = experiment_key[i]
    fov = experiment_fov_key[i]
    slack_decay = experiment_slack_decay_key[i]
    for r_idx, num_r in enumerate(num_robot):
        print("r_idx: ", r_idx, "num_r: ", num_r)
        for exp_idx in range(num_exp):
            if is_slack_off:
                state_json = result_path+date+"/"+instance_type+str(num_r)+"_fov"+str(fov)+"_slack_off"+"_States_"+str(exp_idx)+".json"
            else:
                state_json = result_path+date+"/"+instance_type+str(num_r)+"_fov"+str(fov)+"_decay"+str(slack_decay)+"_States_"+str(exp_idx)+".json"
            config_json = config_path+instance_type+str(num_r)+"_config.json"
            states = load_states(state_json)
            goals = np.array(load_states(config_json)["tasks"]["sf"])
            num_robots = len(states["robots"])
            traj = np.array([states["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
            # TODO [optional] downsampling traj
            traj = traj[:,::10,:]
            if QP_SUCCESS_METRIC:
                QP_success = np.array([states["robots"][str(_)]["QP_success"] for _ in range(num_robots)])  # [n_robot, loops]

            collision_shape = np.array(load_states(config_json)["robot_params"]["collision_shape"]["aligned_box"][:2]) - 0.01
            FoV_beta = fov * np.pi/180
            FoV_range = load_states(config_json)["fov_cbf_params"]["Rs"]
            Ts = states["Ts"]
            # TODO [optional] downsampling Ts
            Ts = Ts*10

            # Metric1: success rate
            goal_radius = 1.6
            success, makespan = instance_success(traj, goals, goal_radius, collision_shape)
            print("success: ", success, "makespan: ", makespan*Ts)
            success_rate_dict[key][r_idx].append(float(success))
            makespan_dict[key][r_idx].append(makespan*Ts)

            # Metric2: avg number of neighbor in FoV
            num_neighbors = num_r - 1
            avg_num_neighbor_in_fov = avg_neighbor_in_fov(traj, FoV_beta, makespan)
            mean_avg_num_neighbor_in_fov = np.mean(avg_num_neighbor_in_fov)
            # print("avg_num_neighbor_in_fov: ", avg_num_neighbor_in_fov)
            # print("mean_avg_num_neighbor_in_fov: ", mean_avg_num_neighbor_in_fov)
            num_neighbor_in_fov_dict[key][r_idx].append(mean_avg_num_neighbor_in_fov)
            percent_neighbor_in_fov_dict[key][r_idx].append(mean_avg_num_neighbor_in_fov / num_neighbors)

            # Metric3: QP success rate
            avg_QP_success = avg_QP_success_rate(QP_success)
            mean_avg_QP_success = np.mean(avg_QP_success)
            QP_success_rate_dict[key][r_idx].append(mean_avg_QP_success)


# plots
num_robot_arr = np.array(num_robot)
success_rate_mean_arr_list = []
success_rate_negative_ci_arr_list = []
success_rate_positive_ci_arr_list = []
success_rate_median_arr_list = []
success_rate_q1_arr_list = []
success_rate_q3_arr_list = []
success_rate_label_list = []
success_samples = []

makespan_mean_arr_list = []
makespan_negative_ci_arr_list = []
makespan_positive_ci_arr_list = []
makespan_median_arr_list = []
makespan_q1_arr_list = []
makespan_q3_arr_list = []
makespan_label_list = []
makespan_samples = []

num_neighbor_in_fov_mean_arr_list = []
num_neighbor_in_fov_negative_ci_arr_list = []
num_neighbor_in_fov_positive_ci_arr_list = []
num_neighbor_in_fov_median_arr_list = []
num_neighbor_in_fov_q1_arr_list = []
num_neighbor_in_fov_q3_arr_list = []
num_neighbor_in_fov_label_list = []
num_neighbor_in_fov_samples = []

percent_neighbor_in_fov_mean_arr_list = []
percent_neighbor_in_fov_negative_ci_arr_list = []
percent_neighbor_in_fov_positive_ci_arr_list = []
percent_neighbor_in_fov_median_arr_list = []
percent_neighbor_in_fov_q1_arr_list = []
percent_neighbor_in_fov_q3_arr_list = []
percent_neighbor_in_fov_label_list = []
percent_neighbor_in_fov_samples = []

QP_success_rate_mean_arr_list = []
QP_success_rate_negative_ci_arr_list = []
QP_success_rate_positive_ci_arr_list = []
QP_success_rate_median_arr_list = []
QP_success_rate_q1_arr_list = []
QP_success_rate_q3_arr_list = []
QP_success_rate_label_list = []
QP_success_rate_samples = []

num_colors = 0
for key in experiment_key:
    success_rate_ci_dict[key]["mean"], success_rate_ci_dict[key]["ci"] = CI_compute(np.array(success_rate_dict[key]))
    success_rate_percentile_dict[key]["median"], success_rate_percentile_dict[key]["q1"], success_rate_percentile_dict[key]["q3"] = percentile_compute(np.array(success_rate_dict[key]))
    success_samples.append(np.array(success_rate_dict[key]))
    makespan_ci_dict[key]["mean"], makespan_ci_dict[key]["ci"] = CI_compute_with_inf(np.array(makespan_dict[key]))
    makespan_percentile_dict[key]["median"], makespan_percentile_dict[key]["q1"], makespan_percentile_dict[key]["q3"] = percentile_compute_with_inf(np.array(makespan_dict[key]))
    makespan_samples.append(np.array(makespan_dict[key]))
    print("mean: ", makespan_ci_dict[key]["mean"], "ci: ", makespan_ci_dict[key]["ci"])
    num_neighbor_in_fov_ci_dict[key]["mean"], num_neighbor_in_fov_ci_dict[key]["ci"] = CI_compute_with_inf(np.array(num_neighbor_in_fov_dict[key]))
    num_neighbor_in_fov_percentile_dict[key]["median"], num_neighbor_in_fov_percentile_dict[key]["q1"], num_neighbor_in_fov_percentile_dict[key]["q3"] = percentile_compute_with_inf(np.array(num_neighbor_in_fov_dict[key]))
    num_neighbor_in_fov_samples.append(np.array(num_neighbor_in_fov_dict[key]))
    percent_neighbor_in_fov_ci_dict[key]["mean"], percent_neighbor_in_fov_ci_dict[key]["ci"] = CI_compute(np.array(percent_neighbor_in_fov_dict[key]))
    percent_neighbor_in_fov_percentile_dict[key]["median"], percent_neighbor_in_fov_percentile_dict[key]["q1"], percent_neighbor_in_fov_percentile_dict[key]["q3"] = percentile_compute(np.array(percent_neighbor_in_fov_dict[key]))
    percent_neighbor_in_fov_samples.append(np.array(percent_neighbor_in_fov_dict[key]))
    QP_success_rate_ci_dict[key]["mean"], QP_success_rate_ci_dict[key]["ci"] = CI_compute(np.array(QP_success_rate_dict[key]))
    QP_success_rate_percentile_dict[key]["median"], QP_success_rate_percentile_dict[key]["q1"], QP_success_rate_percentile_dict[key]["q3"] = percentile_compute(np.array(QP_success_rate_dict[key]))
    QP_success_rate_samples.append(np.array(QP_success_rate_dict[key]))

    success_rate_mean_arr_list.append(success_rate_ci_dict[key]["mean"])
    success_rate_negative_ci_arr_list.append(success_rate_ci_dict[key]["ci"])
    success_rate_positive_ci_arr_list.append(success_rate_ci_dict[key]["ci"])
    success_rate_median_arr_list.append(success_rate_percentile_dict[key]["median"])
    success_rate_q1_arr_list.append(success_rate_percentile_dict[key]["q1"])
    success_rate_q3_arr_list.append(success_rate_percentile_dict[key]["q3"])
    success_rate_label_list.append(key)

    makespan_mean_arr_list.append(makespan_ci_dict[key]["mean"])
    makespan_negative_ci_arr_list.append(makespan_ci_dict[key]["ci"])
    makespan_positive_ci_arr_list.append(makespan_ci_dict[key]["ci"])
    makespan_median_arr_list.append(makespan_percentile_dict[key]["median"])
    makespan_q1_arr_list.append(makespan_percentile_dict[key]["q1"])
    makespan_q3_arr_list.append(makespan_percentile_dict[key]["q3"])
    makespan_label_list.append(key)

    num_neighbor_in_fov_mean_arr_list.append(num_neighbor_in_fov_ci_dict[key]["mean"])
    num_neighbor_in_fov_negative_ci_arr_list.append(num_neighbor_in_fov_ci_dict[key]["ci"])
    num_neighbor_in_fov_positive_ci_arr_list.append(num_neighbor_in_fov_ci_dict[key]["ci"])
    num_neighbor_in_fov_median_arr_list.append(num_neighbor_in_fov_percentile_dict[key]["median"])
    num_neighbor_in_fov_q1_arr_list.append(num_neighbor_in_fov_percentile_dict[key]["q1"])
    num_neighbor_in_fov_q3_arr_list.append(num_neighbor_in_fov_percentile_dict[key]["q3"])
    num_neighbor_in_fov_label_list.append(key)

    percent_neighbor_in_fov_mean_arr_list.append(percent_neighbor_in_fov_ci_dict[key]["mean"])
    percent_neighbor_in_fov_negative_ci_arr_list.append(percent_neighbor_in_fov_ci_dict[key]["ci"])
    percent_neighbor_in_fov_positive_ci_arr_list.append(percent_neighbor_in_fov_ci_dict[key]["ci"])
    percent_neighbor_in_fov_median_arr_list.append(percent_neighbor_in_fov_percentile_dict[key]["median"])
    percent_neighbor_in_fov_q1_arr_list.append(percent_neighbor_in_fov_percentile_dict[key]["q1"])
    percent_neighbor_in_fov_q3_arr_list.append(percent_neighbor_in_fov_percentile_dict[key]["q3"])
    percent_neighbor_in_fov_label_list.append(key)

    QP_success_rate_mean_arr_list.append(QP_success_rate_ci_dict[key]["mean"])
    QP_success_rate_negative_ci_arr_list.append(QP_success_rate_ci_dict[key]["ci"])
    QP_success_rate_positive_ci_arr_list.append(QP_success_rate_ci_dict[key]["ci"])
    QP_success_rate_median_arr_list.append(QP_success_rate_percentile_dict[key]["median"])
    QP_success_rate_q1_arr_list.append(QP_success_rate_percentile_dict[key]["q1"])
    QP_success_rate_q3_arr_list.append(QP_success_rate_percentile_dict[key]["q3"])
    QP_success_rate_label_list.append(key)

    num_colors += 1

colors = generate_rgb_colors(num_colors)
assert len(colors) == len(success_rate_mean_arr_list)
list_CI_plot(num_robot_arr, success_rate_mean_arr_list, success_rate_negative_ci_arr_list, success_rate_positive_ci_arr_list, success_rate_label_list, colors, xlabel="Number of Robots", ylabel="Success Rate", save_name="./multi_success_rate_ci_slack_mode.pdf")
list_CI_plot(num_robot_arr, num_neighbor_in_fov_mean_arr_list, num_neighbor_in_fov_negative_ci_arr_list, num_neighbor_in_fov_positive_ci_arr_list, num_neighbor_in_fov_label_list, colors, xlabel="Number of Robots", ylabel="Number of Neighbors in FoV", save_name="./multi_nnif_ci_slack_mode.pdf")
list_CI_plot(num_robot_arr, percent_neighbor_in_fov_mean_arr_list, percent_neighbor_in_fov_negative_ci_arr_list, percent_neighbor_in_fov_positive_ci_arr_list, percent_neighbor_in_fov_label_list, colors, xlabel="Number of Robots", ylabel="Percentage of\n Neighbors in FoV", save_name="./multi_pnif_ci_slack_mode.pdf")
list_CI_plot(num_robot_arr, makespan_mean_arr_list, makespan_negative_ci_arr_list, makespan_positive_ci_arr_list, makespan_label_list, colors, xlabel="Number of Robots", ylabel="Makespan", save_name="./multi_makespan_ci_slack_mode.pdf")
list_CI_plot(num_robot_arr, QP_success_rate_mean_arr_list, QP_success_rate_negative_ci_arr_list, QP_success_rate_positive_ci_arr_list, QP_success_rate_label_list, colors, xlabel="Number of Robots", ylabel="QP Success Rate", save_name="./multi_qp_success_rate_ci_slack_mode.pdf")

histogram_list_CI_plot(num_robot_arr, success_samples, success_rate_mean_arr_list, success_rate_negative_ci_arr_list, success_rate_positive_ci_arr_list, success_rate_label_list, colors, xlabel="Number of Robots", ylabel="Success Rate", legend=True, save_name="./hist_success_rate_ci_slack_mode.pdf")
histogram_list_CI_plot(num_robot_arr, makespan_samples, makespan_mean_arr_list, makespan_negative_ci_arr_list, makespan_positive_ci_arr_list, makespan_label_list, colors, xlabel="Number of Robots", ylabel="Makespan", legend=False, save_name="./hist_makespan_ci_slack_mode.pdf")
histogram_list_CI_plot(num_robot_arr, success_samples, success_rate_median_arr_list, success_rate_q1_arr_list, success_rate_q3_arr_list, success_rate_label_list, colors, xlabel="Number of Robots", ylabel="Success Rate", legend=True, save_name="./hist_success_rate_percentile_slack_mode.pdf")
histogram_list_CI_plot(num_robot_arr, makespan_samples, makespan_median_arr_list, makespan_q1_arr_list, makespan_q3_arr_list, makespan_label_list, colors, xlabel="Number of Robots", ylabel="Makespan", legend=False, save_name="./hist_makespan_percentile_slack_mode.pdf")