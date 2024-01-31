#!/bin/bash

# Job name:
#SBATCH --job-name=phd20230614

# Project:
#SBATCH --account=nn8100k

# Wall time limit:
#SBATCH --time=01-10:00:00

# Other parameters:
#SBATCH --ntasks=100
# SBATCH --mem-per-cpu=1812M
# SBATCH --ntasks-per-node=25
# SBATCH --cpus-per-task=1

# Error and log files
#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/error_files/out_o20230614_set18SHEAR_In3Ca3aa0ar0D0v8Al0Ga3Init1.err
#SBATCH --output=/cluster/projects/nn8100k/harish_workspace/log_files/out_o20230614_set18SHEAR_In3Ca3aa0ar0D0v8Al0Ga3Init1.out

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

module load SuiteSparse/5.10.1-intel-2022a-METIS-5.1.0

## Do some work:
srun /cluster/home/harishpj/softwares/cell_growth_division/build/divisionless_nem_neoint_model /cluster/home/harishpj/softwares/cell_growth_division/init/neoinit_nematic_shear.2d

