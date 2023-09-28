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
out_dir = '/cluster/projects/nn8100k/harish_workspace/phd20230220/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

templ_file = init_dir + '/phd20230220_init_2.templ'
templ_run_file = base_dir + '/run/phd20230220.sh.templ'

a = np.array([1.0, 1.5])
In = np.array([2.5e-2, 5.0e-2,  7.5e-2, 10e-2])
Ca = np.array([5.0e-2, 10.0e-2, 15.0e-2, 20.0e-2, 25.0e-2, 30e-2, 35e-2])
D = np.array([0.01, 0.0])
v = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
alpha = np.array([0.1, 0.4, 1.0, 0.025, 0.05, 2.0, 5.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 0.01, 0.001])
gamma_active = np.array([-0.7, -0.5, -0.3, 0.0, 0.3, 0.5, 0.7])


a_ind = [1]#[0, 1]
In_ind = [3]
Ca_ind = [1,5]#[1, 3]
D_ind = [0]
v_ind = [5]
alpha_ind = [0]
gamma_ind = [1, 2, 3, 4, 5]

multi_exp = ['a', 'b', 'c', 'd', 'e']
multi_indA = 0

print('a: ', a[a_ind])
print('In: ', In[In_ind])
print('Ca: ', Ca[Ca_ind])
print('D: ', D[D_ind])
print('v: ', v[v_ind])
print('alpha : ', alpha[alpha_ind])
print('gamma active: ', gamma_active[gamma_ind])

for indI in In_ind:
    for indC in Ca_ind:
        for indA in a_ind:
            for indD in D_ind:
                for indV in v_ind:
                    for indAl in alpha_ind:
                        for indGa in gamma_ind:
                            ##### change ELONGATION_VEL PARAMETER FROM 3 to 0 to change from shear to no shear 
                            out_dir_full = out_dir + 'out20230220_set2_In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' +  str(indAl) + '_Ga_' +  str(indGa)
                            #out_dir_full = out_dir + 'out_set28' + multi_exp[multi_indA] +  '_In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' +  str(indAl)
                            print(out_dir_full)
                            name = 'out20230220_set2_' + 'In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' + str(indAl) + '_Ga_' +  str(indGa)
                            #name = multi_exp[multi_indA] + 'In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' + str(indAl)
                            print(name)
                            print('Ca>', Ca[indC])
                            print('In>', In[indI])
                            print('a>', a[indA])
                            print('D>', D[indD])
                            print('v>', v[indV])
                            print('Alpha>', alpha[indAl])
                            print('Gamma_active>', gamma_active[indGa])

                            if not os.path.exists(out_dir_full):
                                os.makedirs(out_dir_full)
                            parameters = { \
                            '#{OUT_DIR}': out_dir_full,
                            '#{A}': str(a[indA]),
                            '#{IN}': str(In[indI]),
                            '#{CA}': str(Ca[indC]),
                            '#{D}': str(D[indD]),
                            '#{v}': str(v[indV]),
                            '#{ALPHA}': str(alpha[indAl]),
                            '#{GAMMA_A}': str(gamma_active[indGa]),
                            '#{INDEX}': name}

                            postfix = name
                            init_file = out_dir + 'i_' + postfix + '.dat.2d'
                            generate_file(templ_file, init_file, parameters)

                            parameters['#{INIT_FILE}'] = init_file
                            run_file = out_dir + 'r_' + postfix + '.sh'
                            generate_file(templ_run_file, run_file, parameters)
                            subprocess.call('sbatch --exclusive ' + run_file, shell=True)
                
            multi_indA+=1
