#!/bin/bash
#SBATCH -J savequant
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=3:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-workspacebeegfs/savequant.err
#SBATCH --output=/beegfs/ws/1/haja565a-workspacebeegfs/savequant.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load amdis
module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2

srun python savequant.py /beegfs/ws/1/haja565a-workspacebeegfs/out10b 150