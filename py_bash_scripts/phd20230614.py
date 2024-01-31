#from doctest import script_from_examples
import sys
import os
import subprocess
import time
import numpy as np

def generate_file(in_filename, out_filename, params):
    out_ = open(out_filename,'w')
    with open(in_filename,'r') as in_:
        for line in in_:
            for p,v in params.items():
                line = line.replace(p, v)
            out_.write(line)
# {end def}

setup = 'benchmark2'

base_dir = '/cluster/home/harishpj/softwares/cell_growth_division'
exec_dir = base_dir + '/build'
init_dir = base_dir + '/init' 
out_dir = '/cluster/projects/nn8100k/harish_workspace/phd20230614/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

templ_file = init_dir + '/phd20230614_init.templ'
#templ_file = init_dir + '/phd20230614_init_rectangles.templ'
templ_run_file = base_dir + '/run/phd20230614.sh.templ'
#templ_run_file = base_dir + '/run/phd20230614_rectangles.sh.templ'

a = np.array([1.0, 1.5])
aa = np.array([1.0])
ar = np.array([1.0, 0.75, 0.5, 0.25])
In = np.array([2.5e-2, 5.0e-2,  7.5e-2, 10e-2])
Ca = np.array([5.0e-2, 10.0e-2, 15.0e-2, 20.0e-2, 25.0e-2, 30e-2, 35e-2])
D = np.array([0.01, 0.0])
v = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95])
alpha = np.array([0.1, 0.4, 1.0, 0.025, 0.05, 2.0, 5.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 0.01, 0.001])
gamma_active = np.array([-0.7, -0.5, -0.3, 0.0, 0.3, 0.5, 0.7])
initial_folder = np.array(['/cluster/projects/nn8100k/harish_workspace/saved_states/hexatic_random', 
                '/cluster/projects/nn8100k/harish_workspace/saved_states/hexatic_random_1_0'])
#the first folder has packing of 0.99 and second has packing of 1.0
timestep_size = np.array([5e-3, 2e-3, 1e-3, 10e-3, 7.5e-3])

aa_ind = [0]
ar_ind = [0] #[0]
In_ind = [3]
Ca_ind = [3]
D_ind = [0]
v_ind = [5] #[11, 12, 13, 14, 15, 16, 17, 18, 19, 20]#, 4, 3, 2 ,1]
alpha_ind = [7]
gamma_ind = [3]
ini_ind = [1] #[1]
dt_ind = [0]#, 1]#, 2, 3, 4]


multi_exp = ['a', 'b', 'c', 'd', 'e']
multi_indA = 0

print('aa: ', aa[aa_ind])
print('ar: ', ar[ar_ind])
print('In: ', In[In_ind])
print('Ca: ', Ca[Ca_ind])
print('D: ', D[D_ind])
print('v: ', v[v_ind])
print('alpha : ', alpha[alpha_ind])
print('gamma active: ', gamma_active[gamma_ind])
print('initial condition: ', initial_folder[ini_ind])
print('timestep sizes condition: ', timestep_size[dt_ind])


for indI in In_ind:
    for indC in Ca_ind:
        for indAa in aa_ind:
            for indAr in ar_ind:
                for indD in D_ind:
                    for indV in v_ind:
                        for indAl in alpha_ind:
                            for indGa in gamma_ind:
                                for indInit in ini_ind:
                                    for inddt in dt_ind:
                                        ##### change ELONGATION_VEL PARAMETER FROM 3 to 0 to change from shear to no shear 
                                        out_dir_full = out_dir + 'o20230614_set17_In' + str(indI) + 'Ca' + str(indC) + 'aa' + str(indAa) + 'ar' + str(indAr) + 'D' + str(indD) + 'v' + str(indV) + 'Al' + str(indAl) + 'Ga' +  str(indGa) + 'Init' + str(indInit) + 'dt' + str(inddt)
                                        #out_dir_full = out_dir + 'out_set28' + multi_exp[multi_indA] +  '_In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' +  str(indAl)
                                        print(out_dir_full)
                                        name = 'out20230614_set17_In' + str(indI) + 'Ca' + str(indC) + 'aa' + str(indAa) + 'ar' + str(indAr) + 'D' + str(indD) + 'v' + str(indV) + 'Al' +  str(indAl) + 'Ga' +  str(indGa) + 'Init' + str(indInit) + 'dt' + str(inddt)
                                        #name = multi_exp[multi_indA] + 'In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' + str(indAl)
                                        print(name)
                                        print('Ca>', Ca[indC])
                                        print('In>', In[indI])
                                        print('aa>', aa[indAa])
                                        print('ar>', ar[indAr])
                                        print('D>', D[indD])
                                        print('v>', v[indV])
                                        print('Alpha>', alpha[indAl])
                                        print('Gamma_active>', gamma_active[indGa])
                                        print('initial condition>', initial_folder[indInit])
                                        print('Timestep size', timestep_size[inddt])
                                        if not os.path.exists(out_dir_full):
                                            os.makedirs(out_dir_full)
                                        parameters = { \
                                        '#{OUT_DIR}': out_dir_full,
                                        '#{Aa}': str(aa[indAa]),
                                        '#{Ar}': str(ar[indAr]),
                                        '#{IN}': str(In[indI]),
                                        '#{CA}': str(Ca[indC]),
                                        '#{D}': str(D[indD]),
                                        '#{v}': str(v[indV]),
                                        '#{ALPHA}': str(alpha[indAl]),
                                        '#{GAMMA_A}': str(gamma_active[indGa]),
                                        '#{INI_FOLDER}': str(initial_folder[indInit]),
                                        '#{INDEX}': name,
                                        '#{dt}': str(timestep_size[inddt])}

                                        postfix = name
                                        init_file = out_dir + 'i_' + postfix + '.dat.2d'
                                        generate_file(templ_file, init_file, parameters)

                                        parameters['#{INIT_FILE}'] = init_file
                                        run_file = out_dir + 'r_' + postfix + '.sh'
                                        generate_file(templ_run_file, run_file, parameters)
                                        subprocess.call('sbatch ' + run_file, shell=True)
                    
            #multi_indA+=1
