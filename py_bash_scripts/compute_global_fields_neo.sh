#!/bin/bash

# Job name:
#SBATCH --job-name=Ca3Al7v5
#
# Project:
#SBATCH --account=nn8100k
#
# Wall time limit:
#SBATCH --time=00-06:00:00
#
# Other parameters:
# SBATCH --qos=short

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# Error and log files
#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/error_files/Ca3Al7v10.err
#SBATCH --output=/cluster/projects/nn8100k/harish_workspace/log_files/Ca3Al7v10.out

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

# load the Anaconda3
module load Anaconda3/2022.05

# Set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

# Source the conda environment setup
# The variable ${EBROOTANACONDA3} or ${EBROOTMINICONDA3}
# So use one of the following lines
# comes with the module load command
# source ${EBROOTMINICONDA3}/etc/profile.d/conda.sh
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh

# Deactivate any spill-over environment from the login node
conda deactivate &>/dev/null

# Activate the environment by using the full path (not name)
# to the environment. The full path is listed if you do
# conda info --envs at the command prompt.
conda activate /cluster/projects/nn8100k/condastuff/myenv

# syntax for compute stress or force
# srun python compute_stress_or_force.py PATH A Ca B (a or a_rep) a_adh
# A - stress and energy = 1, always choose 1
# B - old or new potential, old = 1, new = 2
# example srun python compute_stress_or_force.py /cluster/projects/nn8100k/harish_workspace/neoinit107/ 1 20e-2 2 1.0 0.5


srun python compute_phi_field.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set5_In3Ca3aa0ar0D0v10Al7Ga3Init0/
srun python neo_positions.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set5_In3Ca3aa0ar0D0v10Al7Ga3Init0/ 100
srun python compute_stress_or_force.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set5_In3Ca3aa0ar0D0v10Al7Ga3Init0/ 1 20.0e-2 2 1.0 1.0
#different version of python needed to run shapeprop
module load ParaView/5.8.0-foss-2020a-Python-3.8.2-mpi 
srun python3.8 compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set5_In3Ca3aa0ar0D0v10Al7Ga3Init0/
