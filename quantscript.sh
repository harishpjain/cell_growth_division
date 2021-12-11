#!/bin/bash
#SBATCH -J Quantd
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=8:00:00
#SBATCH --output=/scratch/ws/1/haja565a-workspace2/master_thesis/Quantd.out
#SBATCH --error=/scratch/ws/1/haja565a-workspace2/master_thesis/Quantd.err
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500


srun python savequant.py 700a13b 150
srun python savequant.py 700h13 150
