#!/bin/bash
#SBATCH -J sups6
#SBATCH --ntasks=16
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=16:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-my_workspace_neo/supracell/sim_6.err
#SBATCH --output=/beegfs/ws/1/haja565a-my_workspace_neo/supracell/sim_6.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load amdis



#srun ./build/modified_elongation_model init/supracellsave.2d

#srun ./build/supracell_model init/supracellread.2d   
srun ./build/divisionless_tensionless_model init/supracellread.2d   