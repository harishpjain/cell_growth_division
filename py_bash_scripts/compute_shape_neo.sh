#!/bin/bash

# Job name:
#SBATCH --job-name=compute_fields
#
# Project:
#SBATCH --account=nn8100k
#
# Wall time limit:
#SBATCH --time=00-03:00:00
#
# Other parameters:

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# Error and log files
#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/error_files/shape5.err
#SBATCH --output=/cluster/projects/nn8100k/harish_workspace/log_files/shape5.out

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

# load the Anaconda3
module load Anaconda3/2022.05
module load ParaView/5.8.0-foss-2020a-Python-3.8.2-mpi 
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

#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_3_a_0_D_0_v_5_Al_0_Ga_1/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_3_a_0_D_0_v_5_Al_0_Ga_2/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_3_a_0_D_0_v_5_Al_0_Ga_3/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_3_a_0_D_0_v_5_Al_0_Ga_4/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_3_a_0_D_0_v_5_Al_0_Ga_5/

#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_1_a_1_D_0_v_5_Al_0_Ga_1/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_1_a_1_D_0_v_5_Al_0_Ga_2/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_1_a_1_D_0_v_5_Al_0_Ga_3/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_1_a_1_D_0_v_5_Al_0_Ga_4/
#srun python compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230220/out20230220_set2_In_3_Ca_1_a_1_D_0_v_5_Al_0_Ga_5/

#srun python3.8 compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set2_In3Ca3aa0ar0D0v16Al0Ga3Init1/
#srun python3.8 compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set2_In3Ca3aa0ar0D0v14Al0Ga3Init1/
srun python3.8 compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set2_In3Ca3aa0ar0D0v19Al0Ga3Init1/
srun python3.8 compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set2_In3Ca3aa0ar0D0v20Al0Ga3Init1/
#srun python3.8 compute_shapeprop.py /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set2_In3Ca3aa0ar0D0v12Al0Ga3Init1/