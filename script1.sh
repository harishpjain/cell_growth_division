#!/bin/bash
#SBATCH -J out134read
#SBATCH --ntasks=10
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=3:00:00
#SBATCH --error=/scratch/ws/1/haja565a-workspace1/collision/out134read.err
#SBATCH --output=/scratch/ws/1/haja565a-workspace1/collision/out134read.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load amdis
srun ./build/modified_elongation_model init/taurus_model_comparison_random.2d
