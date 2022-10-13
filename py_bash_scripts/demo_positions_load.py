import matplotlib.pyplot as plt
import sys
import numpy as np
import os
import re
import glob
import pandas as pd
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata
from scipy import ndimage
from scipy.ndimage import gaussian_filter
import load_files as load_files

"""
This code is valid only for square domains. 

Returns:
    _type_: _description_
"""

print(sys.argv)

quantity = 'force'
if int(sys.argv[2]) == 1:
    quantity = 'stress'
    print('checkpoint 1')
if int(sys.argv[2]) == 2:
    print('checkpoint 1')
    quantity = 'force'
    
pbc=True

# --- Parameters
domain_size = np.array([100.0,100.0])
confined = False

if confined:
    mode = 'nearest'
else:
    mode = 'wrap'

interpolation_steps = 500
# --- Parameters

if len(sys.argv) > 1:
    file_pattern = sys.argv[1] + 'neo_positions_p*.csv'
    out_dir = sys.argv[1] + 'stress_fields_500'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
        


print(file_pattern)

if len(sys.argv) > 3:
    endTimestep = int(sys.argv[3])
else:
    endTimestep = np.inf

if len(sys.argv) > 4:
    stride = int(sys.argv[4])
else:
    stride = 100

positions_raw = []
ranks = []
count = 0
for filename in glob.glob(file_pattern):
    if filename.endswith('.csv'):
        print(count)
        print(filename)
        count += 1
        tmp = pd.read_csv(filename)#, skiprows=2)# engine='python')#, engine='python')#-fwf')#, engine='python')
        #print(tmp)
        #print(tmp[1][0])
        positions_raw.append(tmp[['time','rank']])
        # we also need to extract the rank to build the bridge to vtk files
        ranks.append(int(re.findall(r'\d+', filename)[-1]))
ranks = np.array(ranks)
