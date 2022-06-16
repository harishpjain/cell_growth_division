from doctest import script_from_examples
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

base_dir = '/home/haja565a/softwares/cell_growth_division'
exec_dir = base_dir + '/build'
init_dir = base_dir + '/init' 
out_dir = '/beegfs/ws/1/haja565a-my_workspace_neo/phd10052022/' 
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

templ_file = init_dir + '/phd10052022_init.templ'
templ_run_file = base_dir + '/run/phd10052022.sh.templ'

a = np.array([1.0, 1.5])
In = np.array([2.5e-2, 5.0e-2,  7.5e-2, 10e-2])
Ca = np.array([5.0e-2, 10.0e-2, 15.0e-2, 20.0e-2])
D = np.array([0.01, 0.0])
v = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
alpha = np.array([0.1, 0.4, 1.0, 0.025, 0.05, 2.0, 5.0, 0.0])

a_ind = [1]
In_ind = [3]
Ca_ind = [3]
D_ind = [0]
v_ind = [5]
alpha_ind = [2,5]

for indI in In_ind:
    for indC in Ca_ind:
        for indA in a_ind:
            for indD in D_ind:
                for indV in v_ind:
                    for indAl in alpha_ind:
                        out_dir_full = out_dir + 'out_set7_In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' +  str(indAl)
                        print(out_dir_full)
                        name = 'In_' + str(indI) + '_Ca_' + str(indC) + '_a_' + str(indA) + '_D_' + str(indD) + '_v_' + str(indV) + '_Al_' + str(indAl)
                        print(name)
                        print('Ca>', Ca[indC])
                        print('In>', In[indI])
                        print('a>', a[indA])
                        print('D>', D[indD])
                        print('v>', v[indV])
                        print('Alpha>', alpha[indAl])

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
                        '#{INDEX}': name}

                        postfix = name
                        init_file = out_dir + 'i_' + postfix + '.dat.2d'
                        generate_file(templ_file, init_file, parameters)

                        parameters['#{INIT_FILE}'] = init_file
                        run_file = out_dir + 'r_' + postfix + '.sh'
                        generate_file(templ_run_file, run_file, parameters)
                        subprocess.call('sbatch --exclusive ' + run_file, shell=True)
