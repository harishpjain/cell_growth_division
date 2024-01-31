#!/bin/bash

# Job name:
#SBATCH --job-name=phd20230614
#
# Project:
#SBATCH --account=nn8100k
#
# Wall time limit:
#SBATCH --time=00-71:30:00
#
# Other parameters:

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25

#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/error_files/easy.err
#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/log_files/easy.out
## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default


module load EasyBuild/4.8.2
my_path=/cluster/projects/nn8100k/easybuild
srun eb --prefix=$my_path --robot --parallel=25 PETSc-3.12.4-intel-2019b-Python-3.7.4.eb

