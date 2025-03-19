cd ../../lib/mpc_cbf/cmake-build-release
cmake --build ./ --target mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example
date="03202025"
instances="circle"
#instances="formation"
min_r=2
max_r=3
max_exp=2
default_slack_decay=0.2
for (( exp_idx = 0; exp_idx < ${max_exp}; exp_idx++ )); do
  (for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
    for fov in 120 240; do
        echo experiment instance type:${instances} num_robot:${num_r}, exp_idx:${exp_idx}
        ./mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example --instance_type ${instances} --num_robots ${num_r} --fov ${fov} --slack_mode 0 --slack_decay ${default_slack_decay} --write_filename "/media/lishuo/ssd/RSS2025_results/log${date}/${instances}${num_r}_fov${fov}_slack_off_States_${exp_idx}.json"
    done
  done) &

  (for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
      for fov in 120 240; do
          echo experiment instance type:${instances} num_robot:${num_r}, exp_idx:${exp_idx}
          ./mpc_cbf_examples_BezierIMPCCBFPFXYYaw_example --instance_type ${instances} --num_robots ${num_r} --fov ${fov} --slack_decay ${default_slack_decay} --write_filename "/media/lishuo/ssd/RSS2025_results/log${date}/${instances}${num_r}_fov${fov}_decay${default_slack_decay}_States_${exp_idx}.json"
      done
    done) &

done
wait # Wait for all background processes to complete