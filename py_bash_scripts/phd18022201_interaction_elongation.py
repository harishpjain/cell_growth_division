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
command = exec_dir + '/collective3'
#out_dir = '/scratch/ws/1/haja565a-workspace1/collision/' + setup
#out_dir = '/scratch/ws/1/haja565a-workspace1/division_less/' 
out_dir = '/beegfs/ws/1/haja565a-workspacebeegfs/phd18022201/' 
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

templ_file = init_dir + '/phd1802202201_interaction_elongation_init.templ'
templ_run_file = base_dir + '/run/phd1802202201.sh.templ'


a = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]) #1.5 #
In = 2.5e-2 #np.array([10e-2, 7.5e-2, 5.0e-2, 2.5e-2, 1.0e-2, 0.75e-2, 0.5e-2]) #
Ca = 10.0e-2#np.array([20.0e-2, 15.0e-2, 10.0e-2, 5.0e-2])

for ind, a_value in enumerate(a):
    out_dir_full = out_dir +  'cout_3000t_Ca_10_In_2_5_a_' + str(ind+1)
    if not os.path.exists(out_dir_full):
        os.makedirs(out_dir_full)
    name = 'a_varyb'
    parameters = { \
    '#{OUT_DIR}': out_dir_full,
    '#{A}': str(a_value),
    '#{IN}': str(In),
    '#{CA}': str(Ca),
    '#{INDEX}': str(ind+1)+name}

    postfix = '_a_' + str(a_value)
    init_file = out_dir + 'i_' + postfix + '.dat.2d'
    generate_file(templ_file, init_file, parameters)

    parameters['#{INIT_FILE}'] = init_file
    run_file = out_dir + 'r_' + postfix + '.sh'
    generate_file(templ_run_file, run_file, parameters)
    subprocess.call('sbatch --exclusive ' + run_file, shell=True)

""" processes = [];
run = 1

N = 24
gamma = 1.0
beta = 0.01
Ca = 0.025
In = 0.1
v0 = 2.0
Pa = 1.0

parameters = { \
  '#{OUTPUT}': out_dir,
  '#{GAMMA}': str(gamma),
  '#{BETA}': str(beta),
  '#{CA}': str(Ca),
  '#{IN}': str(In),
  '#{V0}': str(v0),
  '#{PA}': str(Pa),
  '#{N}': str(N) }

for run in range(5):
  for p in [1,2,4,8,16,24]: #,32,48]:
    postfix = 'N' + str(N) + '_P' + str(p) + '_run' + str(run)
    parameters['#{P}'] = str(p)
    parameters['#{POSTFIX}'] = postfix
        
    init_file = out_dir + '/init_' + postfix + '.dat.2d'
    generate_file(templ_file, init_file, parameters)

    
    parameters['#{COMMAND}'] = command
    parameters['#{INITFILE}'] = init_file
    run_file = out_dir + '/run_' + postfix + '.sh'
    generate_file(templ_run_file, run_file, parameters)

    subprocess.call('sbatch --exclusive ' + run_file, shell=True)
    """
