#!/bin/bash
#SBATCH --job-name=freq-collisions # create a short name for your job
#SBATCH --output=slurm-%A.%a.out   # stdout file
#SBATCH --error=slurm-%A.%a.err    # stderr file
#SBATCH --nodes=1                  # node count
#SBATCH --ntasks=1                 # total number of tasks across all nodes
#SBATCH --cpus-per-task=12         # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G           # memory per cpu-core (4G is default)
#SBATCH --time=23:59:59            # total run time limit (HH:MM:SS)
#SBATCH --array=350,500,50%3       # job array with index values 1705, 170, 20
#SBATCH --mail-type=all            # send email on job start, end and fault
#SBATCH --mail-user=bherber@princeton.edu

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

module purge
module load matlab/R2021b

matlab -nodisplay -nosplash -r "addpath(genpath(\"..\"));test_RandFreqCollisionsSlurm;exit"
