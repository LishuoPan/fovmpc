cd ~/mavros_ws;
source devel/setup.bash;
roslaunch mavros px4_multi.launch ROBOT_ID:=0 fcu_url:="udp://:14540@127.0.0.1:14557" & sleep 10;
cd ~/PX4-Autopilot;
DONT_RUN=1 make px4_sitl_default gazebo-classic;
source ~/mavros_ws/devel/setup.bash;
source Tools/simulation/gazebo-classic/setup_gazebo.bash $(pwd) $(pwd)/build/px4_sitl_default;
export ROS_PACKAGE_PATH=$ROS_PACKAGE_PATH:$(pwd);
export ROS_PACKAGE_PATH=$ROS_PACKAGE_PATH:$(pwd)/Tools/simulation/gazebo-classic/sitl_gazebo-classic;
roslaunch px4 posix_sitl.launch gui:=false & sleep 10;
cd ~/mavros_ws;
source source devel/setup.bash;
roslaunch fovmpc_controller sim_rviz.launch
roslaunch drone_system fake_target.launch
