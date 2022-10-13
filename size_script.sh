#!/bin/bash
#SBATCH -J size
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=24:00:00
#SBATCH --error=/scratch/ws/1/haja565a-workspace2/collision/size.err
#SBATCH --output=/scratch/ws/1/haja565a-workspace2/collision/size.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

srun cd /scratch/ws/1/haja565a-workspace1/collision && du -sh * | sort -hr >> sizes.txt