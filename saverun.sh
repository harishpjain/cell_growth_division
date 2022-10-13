#!/bin/bash
#SBATCH -J savearh
#SBATCH --ntasks=100
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=16:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-neo_workspace3/savearh5.err
#SBATCH --output=/beegfs/ws/1/haja565a-neo_workspace3/savearh5.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load amdis
srun ./build/divisionless_model init/savearh.2d
