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
#command = exec_dir + '/collective3'
#out_dir = '/scratch/ws/1/haja565a-workspace1/collision/' + setup
#out_dir = '/scratch/ws/1/haja565a-workspace1/division_less/' 
out_dir = '/beegfs/ws/1/haja565a-my_workspace/phd24022201_vtu_v1_0/' 
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

templ_file = init_dir + '/phd24022201_extensility_init.templ'
templ_run_file = base_dir + '/run/phd24022201.sh.templ'

#Part 1
#a = 1.5# for In 1 to In 7, and a is 1.0 for In 8
#a = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]) #1.5 #
#In = np.array([5.0e-2, 2.5e-2, 1.0e-2, 0.5e-2]) # from In 1 to In 4
#In = np.array([10.0e-2, 15.0e-2, 20.0e-2]) # from In 5 to In 7
#In = np.array([5.0e-2])#, 15.0e-2, 20.0e-2]) # from In 8
#In = np.array([1.0e-2, 0.5e-2]) # fOR In 3 and 4
#Ca = np.array([20.0e-2, 15.0e-2, 10.0e-2, 5.0e-2])

#Part 2
a = np.array([1.5])# for In 1 to In 7, and a is 1.0 for In 8
#a = np.array([1.0])
#a = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]) #1.5 #
#In = np.array([2.5e-2, 5.0e-2,  7.5e-2, 10e-2]) #In 1 to In 4
In = np.array([2.5e-2])
#Ca = np.array([5.0e-2, 10.0e-2, 15.0e-2, 20.0e-2]) #Ca 1 to Ca 4
#Ca = np.array([15.0e-2, 20.0e-2]) #Ca 3 to Ca 4
#Ca = np.array([5.0e-2, 10.0e-2]) #ca 1 to ca 2
Ca = np.array([20.0e-2])
for indI, In_value in enumerate(In):
    for indC, Ca_value in enumerate(Ca):
        for indA, a_value in enumerate(a):
            out_dir_full = out_dir + 'out30000_In_' + str(indI+1) + '_Ca_' +  str(indC+4) + '_a_' + str(indA+1)
            print(out_dir_full)
            name = 'In_' + str(indI+1) + '_Ca_' +  str(indC+4) + '_a_' + str(indA+1)
            print(name)
            print(Ca_value)
            print(In_value)
            print(a_value)

            if not os.path.exists(out_dir_full):
                os.makedirs(out_dir_full)
            parameters = { \
            '#{OUT_DIR}': out_dir_full,
            '#{A}': str(a_value),
            '#{IN}': str(In_value),
            '#{CA}': str(Ca_value),
            '#{INDEX}': name}

            postfix = name
            init_file = out_dir + 'i_' + postfix + '.dat.2d'
            generate_file(templ_file, init_file, parameters)

            parameters['#{INIT_FILE}'] = init_file
            run_file = out_dir + 'r_' + postfix + '.sh'
            generate_file(templ_run_file, run_file, parameters)
            subprocess.call('sbatch --exclusive ' + run_file, shell=True)
