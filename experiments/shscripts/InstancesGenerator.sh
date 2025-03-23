cd ../instances/
instance_type="circle"
#instance_type="formation"
min_r=2
max_r=20
for (( num_r = ${min_r}; num_r <= ${max_r}; num_r++ )); do
  python3 CircleInstances.py --instance_type ${instance_type} --num_robots ${num_r}
done