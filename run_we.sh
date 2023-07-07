#!/bin/bash
#SBATCH --job-name=twodimers_westpa
#SBATCH --account=hagan-lab
#SBATCH --partition=hagan-compute,hagan-gpu
##SBATCH --nodelist=compute-1-9,compute-3-16,compute-5-[3,12]
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=25
#SBATCH --output=job_logs/slurm.out
#SBATCH --error=job_logs/slurm.err
#SBATCH --time=20:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=smritipradhan@brandeis.edu
#
# run.sh
#
# Run the weighted ensemble simulation. Make sure you ran init.sh first!
#

#module purge
#module load intel/2017.1.132 amber/16
#module unload python
#export WEST_PYTHON=($which python2.7)
export OPENBLAS_NUM_THREADS=1
source env.sh

rm -f west.log
w_run --work-manager processes "$@" &> west.log
${SLURM_JOB_ID}
