#!/bin/bash --login
#SBATCH --account=s1119
#SBATCH --mail-user=<mkaniselvan@iis.ee.ethz.ch>
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --ntasks=320
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=1
#SBATCH --constraint=gpu
#SBATCH --time=01:00:00
#SBATCH --job-name="negf"

export OMP_NUM_THREADS=1
export CRAY_CUDA_MPS=1
module load daint-gpu

srun -n $SLURM_NTASKS OMEN_Daint_gnu64-XC50 omen.cmd > omen_output.out 
