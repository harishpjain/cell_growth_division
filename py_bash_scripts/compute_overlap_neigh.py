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

interpolation_steps = 500
# --- Parameters

if len(sys.argv) > 1:
    file_pattern = sys.argv[1] + 'positions_p*.csv'
    out_dir = sys.argv[1] + 'global_fields_500'
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
ranks = ranks[sorted_index] #now ranks are sorted in ascending order

# Interpolation
x = np.linspace(0,domain_size[0],interpolation_steps)
y = np.linspace(0,domain_size[1],interpolation_steps)
xx,yy = np.meshgrid(x,y)
xx = np.reshape(xx,(interpolation_steps**2,1))
yy = np.reshape(yy,(interpolation_steps**2,1))
np.save(out_dir + '/grid_x',np.reshape(xx,(interpolation_steps,interpolation_steps)))
np.save(out_dir + '/grid_y',np.reshape(yy,(interpolation_steps,interpolation_steps)))
        
        
reader = vtk.vtkXMLUnstructuredGridReader()
print(stride)
print(positions_raw[0].index.shape)
print(positions_raw[0].index)
#row_indices = np.arange(0,positions_raw[0].index.shape[0],stride,dtype=int) uncomment this if 0 problem is fixed
row_indices = np.arange(stride-1,positions_raw[0].index.shape[0],stride,dtype=int)
times = []
phi_gap_fraction = []

count = 0

overlap_neighbours_1 = np.zeros([len(row_indices), len(ranks), len(ranks)], dtype = bool)
overlap_neighbours_2 = np.zeros([len(row_indices), len(ranks), len(ranks)], dtype = bool)
overlap_neighbours_3 = np.zeros([len(row_indices), len(ranks), len(ranks)], dtype = bool)

for indt, ind in enumerate(row_indices):
    # we now have one particular timepoint 
    time = positions_raw[0].iloc[ind]['time']
    print(time)
    #if not os.path.exists(sys.argv[1] + 'data/phase_p0_' + '{:06.3f}'.format(time) + '.vtu'):
    #    continue
    times.append(time)
    phi_all = []
    
    for rank_ind,rank in enumerate(ranks):
        filename = sys.argv[1] + 'data/phase_p' + str(rank) + '_' + '{:06.3f}'.format(time) + '.vtu'
        reader.SetFileName(filename)
        reader.Update()
        data = reader.GetOutput()

        # grid points
        points = data.GetPoints()
        points = vtk_to_numpy(points.GetData())
        # values 
        phi = vtk_to_numpy(data.GetPointData().GetArray(0))


        phi_interp = griddata(points[:,0:2],phi,(xx,yy),method='nearest')

        phi_interp = np.reshape(phi_interp,(interpolation_steps,interpolation_steps))
        phi_all.append(0.5 * phi_interp + 0.5)
    
    # after this we have a list IN THE SAME ORDER AS ranks with all phi
    # now the axis 0 here is the rank axis which we want to remove
    phi_all = np.array(phi_all)

    
    
    #phi_glob = np.max(phi_all,axis=0)
    
    #finding out if cells are neighbours by phasefield overlap
    #NOTE: this method depends on interpolation step size
    overlap_thres1 = 0.01 #for interp steps = 200, domain = 100, use 0.01
    overlap_thres2 = 0.001 #for interp steps = 200, domain = 100, use 0.01
    overlap_thres3 = 0.0001 #for interp steps = 200, domain = 100, use 0.01
    for indr, rank in enumerate(ranks):
        for indr2, rank2 in enumerate(ranks):
            if rank != rank2:
                overlap = np.sum(np.multiply(phi_all[indr], phi_all[indr2]))
                if overlap > overlap_thres1:
                    overlap_neighbours_1[indt][indr][indr2] = 1
                if overlap > overlap_thres2:
                    overlap_neighbours_2[indt][indr][indr2] = 1
                if overlap > overlap_thres3:
                    overlap_neighbours_3[indt][indr][indr2] = 1
       
    
np.save(out_dir + '/timesteps',np.array(times))
np.save(out_dir + '/overlap_neighbours_1',np.array(overlap_neighbours_1))
np.save(out_dir + '/overlap_neighbours_2',np.array(overlap_neighbours_2))
np.save(out_dir + '/overlap_neighbours_3',np.array(overlap_neighbours_3))

