#!/bin/bash
#SBATCH -J defect_global_fields
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=16:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-workspacebeegfs/defectglobalfield.err
#SBATCH --output=/beegfs/ws/1/haja565a-workspacebeegfs/defectglobalfield.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2


#python defect_global_field.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201/out10000_In_1_Ca_1/global_fields/
python strain_fields.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201/out10000_In_1_Ca_1/global_fields/