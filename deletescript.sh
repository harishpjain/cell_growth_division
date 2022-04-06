#!/bin/bash
#SBATCH -J del11
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=24:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-my_workspace/del11.err
#SBATCH --output=/beegfs/ws/1/haja565a-my_workspace/del11.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500



srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_3_Ca_1/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_3_Ca_2/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_3_Ca_3/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_3_Ca_4/

srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_4_Ca_1/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_4_Ca_2/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_4_Ca_3/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_4_Ca_4/

srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_6_Ca_1/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_6_Ca_2/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_6_Ca_3/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_6_Ca_4/


srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_7_Ca_1/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_7_Ca_2/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_7_Ca_3/
srun rm -rf /beegfs/ws/1/haja565a-my_workspace/haja565a-workspacebeegfs-1647745203/phd24022201/b_out10000_In_7_Ca_4/