#!/bin/bash
#SBATCH -J global_fields
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=16:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-my_workspace/globalfield.err
#SBATCH --output=/beegfs/ws/1/haja565a-my_workspace/globalfield.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2

python compute_phi_field.py /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201_vtu_v1_5/out30000_In_3_Ca_4_a_1/
python compute_phi_field.py /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201_vtu_v1_5/out30000_In_4_Ca_4_a_1/
python compute_phi_field.py /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201_vtu_v1_5/out30000_In_2_Ca_4_a_1/
#python compute_global_fields.py /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201_vtu_v1_5/out30000_In_4_Ca_4_a_1/
#python strain_field_interpolated.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201/b_out10000_In_1_Ca_1/

 
#python strain_field_interpolated.py /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203phd24022201_v1_0/out30000_In_1_Ca_1_a_1/ 1000 30000 1 0
#python strain_field_interpolated.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201_v1_0/out30000_In_1_Ca_4_a_1/ 1000 30000 1 0
#python strain_field_interpolated.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201_v1_0/out30000_In_1_Ca_2_a_1/ 1000 30000 1 0
#python strain_field_interpolated.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201_v1_0/out30000_In_1_Ca_3_a_1/ 1000 30000 1 0

###python compute_global_fields.py ${1} ${2} ${3}



