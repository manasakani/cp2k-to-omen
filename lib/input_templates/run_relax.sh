#!/bin/bash -l
#
# CP2K on Piz Daint: 16 nodes, 6 MPI task per node, 2 OpenMP threads per task
#
#SBATCH --mail-user=<mkaniselvan@iis.ee.ethz.ch>
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --account=s1119
#SBATCH --job-name=HfO2_1
#SBATCH --time=4:30:00
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=2
#SBATCH --constraint=gpu

#========================================
# load modules and run simulation
module load daint-gpu
module load CP2K
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1
ulimit -s unlimited
srun cp2k.psmp relax.inp > log_relax.out
