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

# --- Parameters
domain_size = np.array([100.0,100.0])
confined = False

if confined:
    mode = 'nearest'
else:
    mode = 'wrap'

interpolation_steps = 200
# --- Parameters

if len(sys.argv) > 1:
    file_pattern = sys.argv[1] + 'positions_p*.csv'
    out_dir = sys.argv[1] + 'global_fields_200'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

print(file_pattern)

if len(sys.argv) > 2:
    endTimestep = int(sys.argv[2])
else:
    endTimestep = np.inf

if len(sys.argv) > 3:
    stride = int(sys.argv[3])
else:
    stride = 100

positions_raw = []
ranks = []

for filename in glob.glob(file_pattern):
    tmp = pd.read_csv(filename)
    positions_raw.append(tmp[['time','rank']])
    # we also need to extract the rank to build the bridge to vtk files
    ranks.append(int(re.findall(r'\d+', filename)[-1]))
ranks = np.array(ranks)

sorted_index = np.argsort(ranks)
def sort_ranks(elem):
    return elem.iloc[0]['rank']
positions_raw.sort(key=sort_ranks)
#positions_raw = positions_raw[sorted_index]
ranks = ranks[sorted_index]

# Interpolation
x = np.linspace(0,domain_size[0],interpolation_steps)
y = np.linspace(0,domain_size[1],interpolation_steps)
xx,yy = np.meshgrid(x,y)
xx = np.reshape(xx,(interpolation_steps**2,1))
yy = np.reshape(yy,(interpolation_steps**2,1))

reader = vtk.vtkXMLUnstructuredGridReader()
print(stride)
#row_indices = np.arange(0,positions_raw[0].index.shape[0],stride,dtype=int) uncomment this if 0 problem is fixed
row_indices = np.arange(stride-1,positions_raw[0].index.shape[0],stride,dtype=int)
times = []


count = 0
for ind in row_indices[:-1]:
    # we now have one particular timepoint
    time_1 = positions_raw[0].iloc[ind]['time']
    time_2 = positions_raw[0].iloc[ind+stride]['time']
    print(time_1)
    print(time_2)
    #if not os.path.exists(sys.argv[1] + 'data/phase_p0_' + '{:06.3f}'.format(time) + '.vtu'):
    #    continue
    times.append(time_1)
    phi_all_diff = []
    phi_all_diff_2 = [] 
    
    for rank_ind,rank in enumerate(ranks):
        filename = sys.argv[1] + 'data/phase_p' + str(rank) + '_' + '{:06.3f}'.format(time_1) + '.vtu'
        reader.SetFileName(filename)
        reader.Update()
        data = reader.GetOutput()

        # grid points
        points = data.GetPoints()
        points = vtk_to_numpy(points.GetData())
        # values 
        phi_1 = vtk_to_numpy(data.GetPointData().GetArray(0))


        phi_interp_1 = griddata(points[:,0:2],phi_1,(xx,yy),method='nearest')

        phi_interp_1 = np.reshape(phi_interp_1,(interpolation_steps,interpolation_steps))
        
        filename2 = sys.argv[1] + 'data/phase_p' + str(rank) + '_' + '{:06.3f}'.format(time_2) + '.vtu'
        reader.SetFileName(filename2)
        reader.Update()
        data = reader.GetOutput()

        # grid points
        points = data.GetPoints()
        points = vtk_to_numpy(points.GetData())
        # values 
        phi_2 = vtk_to_numpy(data.GetPointData().GetArray(0))


        phi_interp_2 = griddata(points[:,0:2],phi_2,(xx,yy),method='nearest')

        phi_interp_2 = np.reshape(phi_interp_2,(interpolation_steps,interpolation_steps))
        
        phi_interp_1 = 0.5*phi_interp_1 + 0.5
        phi_interp_2 = 0.5*phi_interp_2 + 0.5
        
        phi_diff = phi_interp_2 - phi_interp_1
        phi_diff_2 = np.fmax(phi_interp_2, 0.0) - np.fmax(phi_interp_1, 0.0)
        phi_all_diff.append(phi_diff)
        phi_all_diff_2.append(phi_diff_2)
    
    # after this we have a list IN THE SAME ORDER AS ranks with all phi
    # now the axis 0 here is the rank axis which we want to remove
    phi_all_diff = np.array(phi_all_diff)
    phi_all_diff_2 = np.array(phi_all_diff_2)

    # global phasefield, given by a lot of 1s and something in between

    phi_glob_diff = np.sum(phi_all_diff,axis=0)
    phi_glob_diff_2 = np.sum(phi_all_diff_2,axis=0)
    # this is more interesting -> this locally gives the rank index the contributed to the max, ie the cell
    #phi_field = rank_max
    #phi_field[phi_field==rank_max] = ranks[rank_max, :]

    np.save(out_dir + '/phi_diff_field' +  '{:06.3f}'.format(time_1),phi_glob_diff)
    np.save(out_dir + '/phi_diff_2_field' +  '{:06.3f}'.format(time_1),phi_glob_diff_2)
