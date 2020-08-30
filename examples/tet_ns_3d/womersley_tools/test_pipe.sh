#!/bin/bash

rm -rf SOL_*
rm -rf dot_SOL_*

# number of cpu for analysis
num_cpu=6
# number of cpu for postprocess
num_pos=6

# degree elevation in x y direction
ndxy=0

# degree elevation in z direction
ndz=1

# elements inserted in x, y knots
nexy=15;

# elements inserted in z knot
nez=5;

# preprocess for initialization run
./preprocess3d_init -elemx $nexy -elemy $nexy -elemz $nez \
  -addSDegree $ndxy -addTDegree $ndxy -addUDegree $ndz \
  -val_a 1 -val_b 0 \
  -geo_file ./geometry_3d_cylinder.txt \
  -cpu_size $num_cpu

# initialization run
mpirun -np $num_cpu ./nsinit -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 50

# preprocess for isogeometric analysis
./preprocess3d -elemx $nexy -elemy $nexy -elemz $nez \
  -addSDegree $ndxy -addTDegree $ndxy -addUDegree $ndz \
  -val_a 1 -val_b 0 \
  -geo_file ./geometry_3d_cylinder.txt \
  -cpu_size $num_cpu

# isogeometric analysis
mpirun -np $num_cpu ./ns2nd -nl_rtol 1.0e-12 -nl_atol 1.0e-13 -nl_refreq 3 \
  -init_step 1.00e-1 -fina_time 1.0 -sol_rec_freq 1 -ttan_freq 10000 \
  -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 50

# mesh decomposition for postprocess
./prepostproc -cpu_size $num_pos

# plot results in vtk
mpirun -np $num_pos ./vis_3d_ns -time_start 0 -time_step 1 -time_end 10 \
  -dt 5.0e-2

# calcualte error
mpirun -np $num_pos ./post_compare_manu -sol_time 1.0 \
  -sol_d_name SOL_velo_900000010 \
  -sol_p_name SOL_pres_900000010 \
  -nqpx 5 -nqpy 5 -nqpz 5

# EOF
