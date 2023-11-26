#!/bin/bash

# Job name:
#SBATCH --job-name=savejob
#SBATCH --exclusive
#
# Project:
#SBATCH --account=nn8100k
#
# Wall time limit:
#SBATCH --time=00-02:30:00
#
# Other parameters:
# SBATCH --qos=normal

#SBATCH --nodes=4
# SBATCH --ntasks-per-node=25

# Error and log files
#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/error_files/save2.err
#SBATCH --output=/cluster/projects/nn8100k/harish_workspace/log_files/save2.out

## Set up job environment:
#set -o errexit  # Exit the script on any error
#set -o nounset  # Treat any unset variables as an error

module --quiet purge
#module load EasyBuild/4.8.2
#mypath=/cluster/projects/nn8100k/easybuild
#srun -n64 eb Boost.Python-1.71.0-iimpi-2019b.eb --prefix=$mypath --robot --parallel=64 --ignore-locks


#module --quiet purge  # Reset the modules to the system default
#module load PETSc/3.12.4-intel-2020a-Python-3.8.2 #used before upgrade
#module load ParMETIS/4.0.3-iimpi-2020a #used before upgrade
#module load CMake/3.15.3-GCCcore-8.3.0 #used before upgrade
#module load PETSc/3.17.4-foss-2022a 
#module load intel/2022a
#module load ParMETIS/4.0.3-iimpi-2022a
#module load iimpi/2022a
module load SuiteSparse/5.10.1-foss-2021b-METIS-5.1.0 
module load OpenMPI/4.1.1-GCC-11.2.0

#module unload Boost/1.79.0-GCC-11.3.0
#module load Boost.MPI/1.79.0-gompi-2022a
#module load ParMETIS/4.0.3-iimpi-2022b

#module load impi/2021.6.0-intel-compilers-2022.1.0
#module load Valgrind/3.19.0-gompi-2022a 
#srun -n50 -N2 valgrind --leak-check=full --track-origins=yes --log-file=/cluster/projects/nn8100k/harish_workspace/error_files/valgrind_output.%p ./build/divisionless_nem_neoint_model init/neoinit_nematic_5.2d

#mpirun -np 50 valgrind --leak-check=full --track-origins=yes --log-file=/cluster/projects/nn8100k/harish_workspace/error_files/valgrind_output.%p ./build/divisionless_nem_neoint_model init/neoinit_nematic_5.2d

srun -n100 ./build/divisionless_nem_neoint_model init/neoinit_nematic_5.2d
