#!/bin/bash
#SBATCH -J QT4
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=7:00:00
#SBATCH --output=/scratch/ws/1/haja565a-workspace2/master_thesis/QTens4.out
#SBATCH --error=/scratch/ws/1/haja565a-workspace2/master_thesis/QTens4.err
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2


srun python agefield.py 700a13b 25.0
srun python QTensorFieldDirect.py 700a13b 25.0


srun python agefield.py 700a13b 10.0
srun python agefield.py 700a13b 20.0
srun python agefield.py 700a13b 30.0
srun python agefield.py 700a13b 40.0
srun python agefield.py 700a13b 5.0
srun python agefield.py 700a13b 15.0
srun python agefield.py 700a13b 25.0
srun python agefield.py 700a13b 35.0
srun python agefield.py 700a13b 2.5
srun python agefield.py 700a13b 7.5
srun python agefield.py 700a13b 12.5
srun python agefield.py 700a13b 17.5
srun python agefield.py 700a13b 22.5
srun python agefield.py 700a13b 27.5
srun python agefield.py 700a13b 32.5
srun python agefield.py 700a13b 37.5
