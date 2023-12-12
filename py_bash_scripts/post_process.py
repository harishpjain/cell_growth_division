#from doctest import script_from_examples
import sys
import os
import subprocess
import time
import numpy as np

longname='o20230614_set13_In3Ca3aa0ar0D0v0Al0Ga3Init1dt0'
shortname='o20230614_set13_In3Ca3aa0ar0D0v0Al0Ga3Init1dt0'
Ca='20e-2'
potential= '2'#1 for old(master thesis) and 2 for new
a_rep_or_a = '1.0'#'a_rep' for new potential and 'a' for old potential
a_adh = '1.0'#'a_ash' for new potential and '' for old potential



def generate_file(in_filename, out_filename, params):
    out_ = open(out_filename,'w')
    with open(in_filename,'r') as in_:
        for line in in_:
            for p,v in params.items():
                line = line.replace(p, v)
            out_.write(line)
# {end def}

setup = 'benchmark2'

base_dir = '/cluster/home/harishpj/softwares/cell_growth_division/py_bash_scripts'
out_dir = '/cluster/projects/nn8100k/harish_workspace/postprocess_runfiles/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

templ_file = base_dir + '/post_process.sh.templ'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
parameters = { \
'#{SHORTNAME}': shortname,
'#{LONGNAME}': longname,
'#{Ca}': Ca,
'#{POTENTIAL}': potential,
'#{A_REP_or_A}': a_rep_or_a,
'#{A_ADH}': a_adh
}
run_file = out_dir + 'postproc_' + longname + '.sh'
generate_file(templ_file, run_file, parameters)
subprocess.call('sbatch ' + run_file, shell=True)

