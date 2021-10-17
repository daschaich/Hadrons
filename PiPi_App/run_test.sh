#!/bin/bash

if [ $# -lt 2 ]; then
  echo "Usage: $0 <config> <time>"
  exit 1
fi

# Input parameters
config=$1
time=$2

# Adjustable parameters
per_node=2
per_task=16
L=24
Nt=64
Lattice="24.24.24.64"
tag=L${L}_l0.035_h0.175

# Strong scaling setup (tasks = tasks_per_node * nodes)
# TODO May be able to make this cleaner
# 1) 32 cores in 1x2=2 tasks with 24x24x24x32 per task
#N=1
#MPIGrid="1.1.1.2"
# 2) 64 cores in 2x2=4 tasks with 24x24x24x16 per task
N=2
MPIGrid="1.1.1.4"
# 3) 128 cores in 4x2=8 tasks with 24x24x12x16 per task
#N=4
#MPIGrid="1.1.2.4"
# 4) 256 cores in 8x2=16 tasks with 24x12x12x16 per task
#N=8
#MPIGrid="1.2.2.4"
# 5) 512 cores in 16x2=32 tasks with 12x12x12x16 per task
#N=16
#MPIGrid="2.2.2.4"
# 6) 1024 cores in 32x2=64 tasks with 12x12x12x8 per task
#N=32
#MPIGrid="2.2.2.8"

tasks=`echo $N | awk -v each="$per_node" '{print($1*each)}'`
echo "$N x $per_node x $per_task"

# Hadrons stuff
#XML=profile_$config.xml
XML=test.xml
bin=test_xml_input

temp=submit$config
echo "#!/bin/bash" > $temp
echo "#SBATCH -A dirac-dp162-CPU" >> $temp
echo "#SBATCH -p skylake" >> $temp
echo "#SBATCH -N $N" >> $temp
echo "#SBATCH --ntasks-per-node=$per_node" >> $temp
echo "#SBATCH --cpus-per-task=$per_task" >> $temp
echo "#SBATCH -t $time" >> $temp
echo "#SBATCH -J xmltest_$tag" >> $temp
echo "#SBATCH -o out.xmltest.%j" >> $temp

# Diagnostic information
echo "echo '=====================JOB DIAGNOTICS========================'" >> $temp
echo "date" >> $temp
echo "echo -n 'This machine is ';hostname" >> $temp
echo "echo -n 'My jobid is '; echo \$SLURM_JOBID" >> $temp
echo "echo 'My path is:' " >> $temp
echo "echo \$PATH" >> $temp
echo "echo 'My job info:'" >> $temp
echo "squeue -j \$SLURM_JOBID" >> $temp


echo "echo '=====================JOB STARTING=========================='" >> $temp
lat=/home/dc-culv1/rds/rds-dirac-dp162/4+6f/24nt64/b4.03_ml0.035_mh0.175/Configs/ckpoint_lat.$config
out=test.$config-${tasks}x${per_task}
if [ ! -e $lat ]; then
  echo "Error: lattice not found"
  rm -f $temp
  exit 1
fi
if [ -e $out ]; then
  echo "Error: output file $out exists"
  rm -f $temp
  exit 1
fi

# Load modules and linker path for hdf5
echo "LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/dc-culv1/lib/hdf5-Sky/install/lib" >> $temp

# Go, with operf and opreport added...
echo "echo \"=== Running MPI application on ${tasks}x${per_task} cpus ===\" >> $out" >> $temp
echo "echo \"srun --mpi=pmi2 -n $tasks $bin $XML --grid $Lattice --mpi $MPIGrid --threads $per_task --comms-concurrent --comms-overlap --shm 2048\" >> $out" >> $temp
echo "srun --mpi=pmi2 -n $tasks $bin $XML --grid $Lattice --mpi $MPIGrid --threads $per_task --comms-concurrent --comms-overlap --shm 2048 >> $out" >> $temp
echo "echo \"=== MPI application finished at \"\`date\`\" ===\" >> $out" >> $temp
echo "echo '========================ALL DONE==========================='" >> $temp

sbatch $temp
rm -f $temp
