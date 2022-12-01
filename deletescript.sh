#!/bin/bash
#SBATCH -J del11
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=24:00:00
#SBATCH --error=beegfs/ws/1/haja565a-workspace/del11.err
#SBATCH --output=beegfs/ws/1/haja565a-workspace/del11.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500


srun rm -rf /beegfs/ws/1/haja565a-workspace/elongated_rectangles/sim11
srun rm -rf /beegfs/ws/1/haja565a-workspace/elongated_rectangles/sim12
srun rm -rf /beegfs/ws/1/haja565a-workspace/elongated_rectangles/sim13
srun rm -rf /beegfs/ws/1/haja565a-workspace/elongated_rectangles/sim14
