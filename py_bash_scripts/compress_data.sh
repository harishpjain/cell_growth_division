#!/bin/bash

# Job name:
#SBATCH --job-name=Ca2_data
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
#SBATCH --error=/cluster/projects/nn8100k/harish_workspace/error_files/dataCa2.err
#SBATCH --output=/cluster/projects/nn8100k/harish_workspace/log_files/dataCa2.out

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default

tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v6Al0Ga3Init1/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v6Al0Ga3Init1/data
rm -r /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v6Al0Ga3Init1/data

tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v7Al0Ga3Init1/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v7Al0Ga3Init1/data
rm -r /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v7Al0Ga3Init1/data

tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v8Al0Ga3Init1/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v8Al0Ga3Init1/data
rm -r /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v8Al0Ga3Init1/data

tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v9Al0Ga3Init1/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v9Al0Ga3Init1/data
rm -r /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v9Al0Ga3Init1/data

tar -cvf /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v10Al0Ga3Init1/data.tar /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v10Al0Ga3Init1/data
rm -r /cluster/projects/nn8100k/harish_workspace/phd20230614/o20230614_set3_In3Ca2aa0ar0D0v10Al0Ga3Init1/data
