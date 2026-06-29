#!/bin/bash
#SBATCH -J 
#SBATCH -p wzhcnormal
#SBATCH -N 1
#SBATCH -n 64

module purge
source /work/home/acsppswc7s/apprepo/openmx/3.9-hpcx_gcc7.3.1/scripts/env.sh
srun --mpi=pmix_v3 openmx test_sqs-opt.dat
