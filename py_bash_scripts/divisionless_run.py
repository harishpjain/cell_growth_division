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
out_dir = '/beegfs/ws/1/haja565a-workspacebeegfs/division_less/' 
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

templ_file = init_dir + '/divisionless_init.templ'
templ_run_file = base_dir + '/run/divisionless_run.sh.templ'


v0 = np.array([0.1, 0.25, 0.5, 0.75, 1.0])
arh_dir = '/scratch/ws/1/haja565a-workspace1/collision/out152/time_012000'
num_cells = 90 #check how many nonzero volume cells are there

for ind, vel in enumerate(v0):
    out_dir_full = out_dir +  'out_20000t_90c_In_2_5e2_el_' + str(ind+1)
    if not os.path.exists(out_dir_full):
        os.makedirs(out_dir_full)
    name = 'In_2_5e2'
    parameters = { \
    '#{OUT_DIR}': out_dir_full,
    '#{V0}': str(vel),
    '#{ARH_DIR}': str(arh_dir),
    '#{ARH_NUM_CELLS}': str(num_cells),
    '#{TOTAL_CELLS}': str(num_cells),
    '#{INDEX}': str(ind+1)+name}

    postfix = 'demo' + '_v0_' + str(vel)
    init_file = out_dir + '/init/i_' + postfix + '.dat.2d'
    generate_file(templ_file, init_file, parameters)

    parameters['#{INIT_FILE}'] = init_file
    run_file = out_dir + '/run/r_' + postfix + '.sh'
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
