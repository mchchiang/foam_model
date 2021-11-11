#!/bin/bash
# model_init.sh

# Directory of this script
sh_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

if (( $# != 4 )); then
    echo "Usage: model_init.sh kappa lambda run run_dir"
    exit 1
fi

kappa=$1     # 
lambda=$2    # 
run=$3       # Trial number
run_dir=$4   # Run directory

pyx=python3 # Choose between python2 and python3

# Required scripts and programs
cell_shape_py="${sh_dir}/cell_shape.py"
init_rand_py="${sh_dir}/init_random.py"

# Format input params
kappa=$($pyx -c "print('{:.3f}'.format($kappa))")
lambda=$($pyx -c "print('{:.3f}'.format($lambda))")

# A function for generating random numbers
function get_rand(){
    # Generate a 4-byte random integer using urandom
    rand=$(od -vAn -N4 -tu4 < /dev/urandom)
    echo $rand
}

# Set the model parameters
ncells=100 # 100
ncell_x=10 # 10
ncell_y=10 # 10
lx=100
ly=100
ideal_radius=8.0 # 8.0
init_radius=5.0 # 4.0
cell_lx=41 # 41
cell_ly=41 # 41
epsilon=0.01 
xi=1.0 # 1.0 - interfacial thickness
max_attempt=1000
rmin=$ideal_radius
rmax=$($pyx -c "print($ideal_radius*2.0)")
period=100000
phase_shift=0

nsteps=500000 # 20000000 # Add another 10^6 steps for equilibration
nequil=10000 # 10000
delta_t=0.5 # 0.5
dump_cm_freq=1000 # 1000
dump_bulk_cm_freq=1000 # 1000
dump_gyr_freq=1000 # 1000
dump_gyr_field_freq=10000 # 1000
dump_vel_freq=1000 # 1000
dump_vel_field_freq=10000 # 1000
dump_deform_freq=1000 # 1000
dump_deform_field_freq=10000 # 1000
dump_field_freq=10000 # 100000
dump_cell_field_freq=10000 # 100000
dump_index_field_freq=10000 # 100000
dump_shape_freq=1000 # 1000
dump_neighbour_freq=1000 # 1000
dump_energy_freq=1000 # 1000
dump_overlap_freq=1000 # 1000
equildump_cm_freq=1000 # 1000
equildump_gyr_freq=10000 # 10000
equildump_field_freq=10000 # 10000
seed=$(get_rand)

tmp_cell_file="cell_$seed.tmp"
tmp_shape_file="shape_$seed.tmp"

# Set a hexagonal lattice
#size=$($pyx triangle.py $ncell_x $ncell_y $confine_radius $tmp_cm_file)
#size=$($pyx square.py $ncell_x $ncell_y $confine_radius $tmp_cm_file)
#size=$($pyx square.2.py 160 138 $ncell_x $ncell_y $confine_radius $tmp_cm_file)
#cutoff=0.2
#seed_2=$(get_rand)
#size=$($pyx triangle_noise.py $ncell_x $ncell_y $confine_radius $cutoff $seed_2 $tmp_cm_file)
#size=$(python triangle_stretched_noise.py $ncell_x $ncell_y $confine_radius $cutoff $seed_2 $tmp_cm_file)

#lx=$(echo $size | awk '{print $3}')
#ly=$(echo $size | awk '{print $6}')

# Generate random positions for cells
seed_2=$(get_rand)
$pyx $init_rand_py $lx $ly $ncells $init_radius $max_attempt $seed_2 $tmp_cell_file
# Sort the cells by position
sort -k3,3g -k2,2g $tmp_cell_file | awk -v r=$ideal_radius -v rmin=$rmin -v rmax=$rmax -v period=$period -v shift=$phase_shift '{
n = NR-1
if (n == 50) {
  print n,$2,$3,"cosine",rmin,rmax,period,shift
} else {
  print n,$2,$3,"const",r
}}' > $tmp_cell_file.tmp
mv $tmp_cell_file.tmp $tmp_cell_file

# Set cell shape
$pyx $cell_shape_py $cell_lx $cell_ly 1.0 $tmp_shape_file circle $init_radius

# Create run directory and set file names
sim_name="cell_N_${ncells}_k_${kappa}_lam_${lambda}_run_${run}"
run_dir="${run_dir}/${sim_name}/"

if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
fi

cell_file="${sim_name}.in"
shape_file="shape_${sim_name}.in"
params_file="params_${sim_name}.txt"
#equildump_cm_file="pos-equil_${sim_name}.dat"
#equildump_gyr_file="gyr-equil_${sim_name}.dat"
#equildump_field_file="field-equil_${sim_name}.dat"
dump_cm_file="pos_${sim_name}.dat"
dump_gyr_file="gyr_${sim_name}.dat"
dump_gyr_field_file="gyr-field_${sim_name}.dat"
dump_vel_file="vel_${sim_name}.dat"
dump_vel_field_file="vel-field_${sim_name}.dat"
dump_deform_file="deform_${sim_name}.dat"
dump_deform_field_file="deform-field_${sim_name}.dat"
dump_field_file="field_${sim_name}.dat"
dump_cell_field_file="cell-field_${sim_name}"
dump_index_field_file="index-field_${sim_name}.dat"
dump_bulk_cm_file="pos-bulk_${sim_name}.dat"
dump_shape_file="shape_${sim_name}.dat"
dump_neighbour_file="neigh_${sim_name}.dat"
dump_energy_file="energy_${sim_name}.dat"
dump_overlap_file="olap_${sim_name}.dat"

# Copy the template file
params_file=${run_dir}/$params_file
mv $tmp_cell_file ${run_dir}/$cell_file
mv $tmp_shape_file ${run_dir}/$shape_file

# Replace macros in template with input values
echo \
"cahnHilliardCoeff = ${kappa}
volumePenaltyCoeff = ${lambda}
repulsionCoeff = ${epsilon}
thickness = ${xi}
cellLx = ${cell_lx}
cellLy = ${cell_ly}
lx = ${lx}
ly = ${ly}
nsteps = ${nsteps}
nequil = ${nequil}
ncells = ${ncells}
dt = ${delta_t}
cell_file = ${cell_file}
shape_file = ${shape_file}
seed = ${seed}
" > $params_file

# Set dumps
function add_dump() {
    params=$1; file=$2;
    if [ "$file" ]; then
	echo "$params $file" >> $params_file 
    fi
}
#add_dump "dump_cm $equildump_cm_freq 0 equil" $equildump_cm_file
#add_dump "dump_gyr $equildump_gyr_freq 0 equil" $equildump_gyr_file
#add_dump "dump_field $equildump_field_freq 0 equil" $equildump_field_file
add_dump "dump_cm $dump_cm_freq 0 main" $dump_cm_file
#add_dump "dump_gyr $dump_gyr_freq 0 main" $dump_gyr_file
#add_dump "dump_gyr_field $dump_gyr_field_freq 0 main" $dump_gyr_field_file
#add_dump "dump_vel $dump_vel_freq 0 main" $dump_vel_file
#add_dump "dump_vel_field $dump_vel_field_freq 0 main" $dump_vel_field_file
#add_dump "dump_deform $dump_deform_freq 0 main" $dump_deform_file
#add_dump "dump_deform_field $dump_deform_field_freq 0 main" $dump_deform_field_file
add_dump "dump_field $dump_field_freq 0 main" $dump_field_file
#add_dump "dump_index_field $dump_index_field_freq 0 main" $dump_index_field_file
add_dump "dump_bulk_cm $dump_bulk_cm_freq main" $dump_bulk_cm_file
add_dump "dump_shape 4 25 4.0 5 31 $dump_shape_freq 0 main" $dump_shape_file
#add_dump "dump_neighbour $dump_neighbour_freq 0 main" $dump_neighbour_file
#add_dump "dump_energy $dump_energy_freq 0 main" $dump_energy_file
#add_dump "dump_overlap $dump_overlap_freq 0 main" $dump_overlap_file

#for (( i=0; $i<$ncells; i++ ))
#do
#    add_dump "dump_cell_field $i $dump_cell_field_freq 0 main" "${dump_cell_field_file}_cell_${i}.dat"
#done
