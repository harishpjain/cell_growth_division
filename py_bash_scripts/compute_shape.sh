#!/bin/bash
#SBATCH -J shape
#SBATCH --ntasks=1
#SBATCH --mail-type=end
#SBATCH --mail-user=harish_p.jain@mailbox.tu-dresden.de
#SBATCH --time=16:00:00
#SBATCH --error=/beegfs/ws/1/haja565a-neo_workspace3/shape5.err
#SBATCH --output=/beegfs/ws/1/haja565a-neo_workspace3/shape5.out
#SBATCH -A wir
#SBATCH -p haswell
#SBATCH --mem-per-cpu=2500

module load VTK/8.2.0-intel-2020a-Python-3.8.2
module load matplotlib/3.2.1-intel-2020a-Python-3.8.2
module load ParaView/5.9.0-RC1-egl-mpi-Python-3.8

#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set17_In_3_Ca_3_a_0_D_0_v_5_Al_2/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set17_In_3_Ca_3_a_0_D_0_v_5_Al_7/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set17_In_3_Ca_3_a_0_D_1_v_5_Al_2/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set17_In_3_Ca_3_a_0_D_1_v_5_Al_7/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set18a_In_3_Ca_3_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set18b_In_3_Ca_3_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set18c_In_3_Ca_3_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set18d_In_3_Ca_3_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set19lowgap_In_3_Ca_3_a_1_D_0_v_5_Al_0/

#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set20_In_3_Ca_0_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set20_In_3_Ca_1_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set20_In_3_Ca_2_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set20_In_3_Ca_4_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set20_In_3_Ca_5_a_1_D_0_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set21_In_3_Ca_3_a_1_D_0_v_5_Al_2/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set21_In_3_Ca_3_a_1_D_0_v_5_Al_7/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set21_In_3_Ca_3_a_1_D_0_v_5_Al_14/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set21_In_3_Ca_3_a_1_D_0_v_5_Al_15/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set22_In_3_Ca_3_a_1_D_1_v_5_Al_0/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-neo_workspace3/phd10052022/out_set22_In_3_Ca_3_a_1_D_1_v_5_Al_7/


#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set15_In_3_Ca_3_a_1_D_0_v_5_Al_7/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set15_In_3_Ca_3_a_1_D_0_v_5_Al_8/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set15_In_3_Ca_3_a_1_D_0_v_5_Al_9/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set15_In_3_Ca_3_a_1_D_0_v_5_Al_10/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set15_In_3_Ca_3_a_1_D_0_v_5_Al_11/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set15_In_3_Ca_3_a_1_D_0_v_5_Al_12/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set15_In_3_Ca_3_a_1_D_0_v_5_Al_13/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set16_In_3_Ca_0_a_1_D_0_v_5_Al_2/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set16_In_3_Ca_2_a_1_D_0_v_5_Al_2/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo2/phd10052022/out_set16_In_3_Ca_4_a_1_D_0_v_5_Al_2/

#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set14_In_3_Ca_3_a_1_D_0_v_5_Al_6/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set8_In_3_Ca_3_a_1_D_0_v_1_Al_2b/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_5_Al_0/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set8_In_3_Ca_3_a_1_D_0_v_1_Al_2b/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set8_In_3_Ca_3_a_1_D_0_v_1_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set8_In_3_Ca_3_a_1_D_0_v_3_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set9_In_1_Ca_1_a_1_D_0_v_5_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set9_In_3_Ca_1_a_1_D_0_v_5_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set10_In_1_Ca_3_a_1_D_0_v_5_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set11_In_3_Ca_3_a_1_D_1_v_5_Al_7/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set12_In_3_Ca_3_a_0_D_0_v_5_Al_2/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set12_In_3_Ca_3_a_0_D_0_v_5_Al_6/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set12_In_3_Ca_3_a_0_D_0_v_5_Al_7/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set13_In_3_Ca_3_a_0_D_1_v_5_Al_7/
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set2_In_3_Ca_3_a_1_D_0_v_5_Al_0/  
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_5_Al_1/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set3_In_3_Ca_3_a_1_D_0_v_5_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set2_In_3_Ca_3_a_1_D_0_v_5_Al_3/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set3_In_3_Ca_3_a_1_D_0_v_5_Al_4/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set2_In_3_Ca_3_a_1_D_0_v_5_Al_5/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_In_3_Ca_3_a_1_D_0_v_5_Al_0/  
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set4_In_3_Ca_3_a_0_D_1_v_5_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set4_In_3_Ca_3_a_1_D_1_v_5_Al_2/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set5_In_3_Ca_3_a_1_D_1_v_5_Al_6/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set5_In_3_Ca_3_a_0_D_1_v_5_Al_6/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set6_In_3_Ca_3_a_1_D_0_v_5_Al_7/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set7_In_3_Ca_3_a_1_D_0_v_5_Al_5/ 
#srun python compute_shapeprop.py /beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/out_set7_In_3_Ca_3_a_1_D_0_v_5_Al_2/