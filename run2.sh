#!/bin/bash
#SBATCH -J out3
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=3:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-my_workspace_neo/out3.err
#SBATCH --output=/beegfs/ws/1/haja565a-my_workspace_neo/out3.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load amdis
srun ./build/divisionless_model init/init2.2d
