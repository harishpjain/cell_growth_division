#!/bin/bash
#SBATCH -J out81b
#SBATCH --ntasks=30
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=6:00:00
#SBATCH --error=/scratch/ws/1/haja565a-workspace1/collision/out81c.err
#SBATCH --output=/scratch/ws/1/haja565a-workspace1/collision/out81c.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load amdis
srun ./build/modified_elongation_model init/init4.2d
