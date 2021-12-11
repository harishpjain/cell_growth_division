#!/bin/bash
#SBATCH -J save2
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/ws/1/haja565a-workspace2/master_thesis/save.out
#SBATCH --error=/scratch/ws/1/haja565a-workspace2/master_thesis/save.err
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500


module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2

srun python saveFieldSingular.py 700h13 150 25.0
srun python saveFieldSingular.py 700a13b 150 25.0
srun python saveFieldSingular.py 700a13b 150 37.5
