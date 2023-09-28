#!/bin/bash

# Job name:
#SBATCH --job-name=demojob
#
# Project:
#SBATCH --account=nn8100k
#
# Wall time limit:
#SBATCH --time=00-01:59:00
#
# Other parameters:
#SBATCH --qos=short

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5

# Error and log files
#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/error_files/nemdemo31.err
#SBATCH --output=/cluster/projects/nn8100k/harish_workspace/log_files/nemdemo31.out

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load PETSc/3.12.4-intel-2020a-Python-3.8.2
module load ParMETIS/4.0.3-iimpi-2020a
module load CMake/3.15.3-GCCcore-8.3.0

## Do some work:
# srun ./build/divisionless_model init/neoinit1.2d
# srun ./build/divisionless_nem_neoint_model init/neoinit_nematic_2.2d
#srun ./build/divisionless_nem_neoint_model init/neoinit_nematic_2.2d
srun tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230220/nemdemo25/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230220/nemdemo25/data 
srun tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230220/nemdemo26/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230220/nemdemo26/data 
srun tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230220/nemdemo27/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230220/nemdemo27/data 
