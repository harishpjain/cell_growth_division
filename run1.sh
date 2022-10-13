#!/bin/bash
#SBATCH -J elonrect
#SBATCH --ntasks=100
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=16:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-my_workspace_neo/elonrect.err
#SBATCH --output=/beegfs/ws/1/haja565a-my_workspace_neo/elonrect.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load amdis
#srun ./build/modified_elongation_model init/supracellsave.2d
srun ./build/divisionless_model init/init1.2d   