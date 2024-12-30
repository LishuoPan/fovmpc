import numpy as np
import json

def load_states(json_filename):
    f = open(json_filename)
    data = json.load(f)
    return data

def rectangles_collide(xA, yA, widthA, heightA, xB, yB, widthB, heightB):
    """
        here xA, yA is the left top corner of recA, widthA, heightA is the width and height of recA
    """
    return (
            xA < xB + widthB and
            xA + widthA > xB and
            yA < yB + heightB and
            yA + heightA > yB
    )

def collision_check(x1, y1, x2, y2, collision_shape):
    """x1, y1 is the center of pos, collision shape is the half of the collision width and height"""
    collision_x = collision_shape[0]
    collision_y = collision_shape[1]
    xA = x1 - collision_x/2
    yA = y1 - collision_y/2
    xB = x2 - collision_x/2
    yB = y2 - collision_y/2
    widthA = 2*collision_x
    heightA = 2*collision_y
    widthB = 2*collision_x
    heightB = 2*collision_y
    return rectangles_collide(xA, yA, widthA, heightA, xB, yB, widthB, heightB)

def instance_success(traj, collision_shape):
    """traj: [n_robot, ts, dim]"""
    n_robot = traj.shape[0]
    ts = traj.shape[1]
    for t in range(ts):
        for i in range(n_robot):
            for j in range(i+1, n_robot):
                pos_1 = traj[i, t, :3]  # [3, ]
                pos_2 = traj[j, t, :3]  # [3, ]
                if collision_check(pos_1[0], pos_1[1], pos_2[0], pos_2[1], collision_shape):
                    return False
    return True

def in_fov(x1, y1, yaw1, x2, y2, fov_beta):
    R = np.array([[np.cos(yaw1), np.sin(yaw1)],
                  [-np.sin(yaw1), np.cos(yaw1)]])
    target_local = R @ np.array([[x2 - x1],
                                 [y2 - y1]])
    angle = np.abs(np.arctan2(target_local[1], target_local[0]))
    return angle < 0.5*fov_beta

def avg_neighbor_in_fov(traj, FoV_beta):
    """traj: [n_robot, ts, dim]"""
    n_robot = traj.shape[0]
    ts = traj.shape[1]
    avg_in_fov = []
    for i in range(n_robot):
        num_neighbor_in_fov = 0
        for j in range(n_robot):
            if i == j:
                continue
            for t in range(ts):
                pos_1 = traj[i, t, :3]  # [3, ]
                pos_2 = traj[j, t, :3]  # [3, ]
                inFoV = in_fov(pos_1[0], pos_1[1], pos_1[2], pos_2[0], pos_2[1], FoV_beta)
                if inFoV:
                    num_neighbor_in_fov += 1
        avg_in_fov.append(num_neighbor_in_fov / ts)
    return avg_in_fov


if __name__ == '__main__':
    states_json = load_states("../instances/results/CBFXYYawStates.json")
    num_robots = len(states_json["robots"])
    traj = np.array([states_json["robots"][str(_)]["states"] for _ in range(num_robots)])  # [n_robot, ts, dim]
    collision_shape = load_states("../../config/config.json")["robot_params"]["collision_shape"]["aligned_box"][:2]
    # Metric1: success rate
    success = instance_success(traj, collision_shape)
    # Metric2: avg number of neighbor in FoV
    FoV_beta = load_states("../../config/config.json")["fov_cbf_params"]["beta"] * np.pi/180
    FoV_range = load_states("../../config/config.json")["fov_cbf_params"]["Rs"]
    avg_num_neighbor_in_fov = avg_neighbor_in_fov(traj, FoV_beta)
    mean_avg_num_neighbor_in_fov = np.mean(avg_num_neighbor_in_fov)
    print("avg_num_neighbor_in_fov: ", avg_num_neighbor_in_fov)

