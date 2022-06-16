#!/bin/bash
#SBATCH -J defect_global_fields
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=16:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-my_workspace_neo/defectglobalfield2.err
#SBATCH --output=/beegfs/ws/1/haja565a-my_workspace_neo/defectglobalfield2.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2
#module load SciPy-bundle/2019.10-foss-2019b-Python-2.7.16 
#module load matplotlib/2.1.2-intel-2018a-Python-2.7.14

srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_5_Al_0
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_5_Al_1
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_5_Al_2
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_10_Al_2
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_10_Al_1
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_1_Ca_4_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_2_Ca_3_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_1_a_2/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_2_a_2/

#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_3_Ca_4_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_2_Ca_4_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_4_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_4_a_2/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_3_Ca_3_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_1_Ca_3_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_3_a_1/
#srun python defect_analyser.py /beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/out30000_In_4_Ca_3_a_2/

#python defect_global_field.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201/out10000_In_1_Ca_1/global_fields/
#python strain_fields.py /beegfs/ws/1/haja565a-workspacebeegfs/phd24022201/out10000_In_1_Ca_1/global_fields/


#cp out30000_In_1_Ca_3_a_1/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i1c3a1.npy
#cp out30000_In_3_Ca_3_a_1/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i3c3a1.npy
#cp out30000_In_4_Ca_3_a_1/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i4c3a1.npy
#cp out30000_In_2_Ca_4_a_1/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i2c4a1.npy
#cp out30000_In_3_Ca_4_a_1/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i3c4a1.npy
#cp out30000_In_4_Ca_4_a_1/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i4c4a1.npy
#cp out30000_In_4_Ca_3_a_2/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i4c3a2.npy
#cp out30000_In_4_Ca_3_a_2/global_fields_200/nematics/defect_nature.npy /home/quasar/mydir/defnat_thres7_i4c4a2.npy
