#!/bin/bash
#SBATCH -J correlation
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=4:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-my_workspace/correlation.err
#SBATCH --output=/beegfs/ws/1/haja565a-my_workspace/correlation.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2

srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_1_Ca_4_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_2_Ca_3_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_1_a_2/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_2_a_2/

srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_3_Ca_4_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_2_Ca_4_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_4_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_4_a_2/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_3_Ca_3_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_1_Ca_3_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_3_a_1/
srun python tissue_correlations.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_3_a_2/
