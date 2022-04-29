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

interpolation_steps = 1000
# --- Parameters

if len(sys.argv) > 1:
    file_pattern = sys.argv[1] + 'positions_p*.csv'
    out_dir = sys.argv[1] + 'global_fields'
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
    positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', 'vel_angle', 'nematic_angle', 'total_interaction', 'neighbours', 'confine_interaction', 'growth_rate', '\'S0full\'', '\'S1full\'']])
    # we also need to extract the rank to build the bridge to vtk files
    ranks.append(int(re.findall(r'\d+', filename)[-1]))

# Interpolation
x = np.linspace(0,domain_size[0],interpolation_steps)
y = np.linspace(0,domain_size[1],interpolation_steps)
xx,yy = np.meshgrid(x,y)
xx = np.reshape(xx,(interpolation_steps**2,1))
yy = np.reshape(yy,(interpolation_steps**2,1))

reader = vtk.vtkXMLUnstructuredGridReader()
print(positions_raw)
print(positions_raw[0].index)
print(positions_raw[0].index.shape)
print(stride)
#row_indices = np.arange(0,positions_raw[0].index.shape[0],stride,dtype=int) uncomment this if 0 problem is fixed
row_indices = np.arange(stride-1,positions_raw[0].index.shape[0],stride,dtype=int)
print(row_indices)
times = []

#test_plot
fig,ax = plt.subplots()
ax.plot(positions_raw[0]['v0'], 'r')

"""
rolling average. the velocity and the nematic field saved is the average 
over a range of time about a given time step. Let's say time = t_m, and rolling average
window is 100 indices wide. Then the avg_velocity(t_m) = average of velocities from time 
t_m - 50*delt to t_m + 50*delt. the average is taken over hundred different times.
This rolling average is dependent on cell advection velocity and on timestep size.
The average is taken in a gaussian way

A T1 transition typically lasts for about 100-400 timesteps if delt = 0.005
"""
rolling_average = True
rolling_window  = 100
if(rolling_average):
    print('rolling average performed over a window: ' + str(rolling_window))
    postfix = '_roll_' + str(rolling_window) + '_'
    for rank_ind,rank in enumerate(ranks):
        positions_raw[rank_ind]['v0'] = positions_raw[rank_ind]['v0'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['v1'] = positions_raw[rank_ind]['v1'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S0'] = positions_raw[rank_ind]['S0'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S1'] = positions_raw[rank_ind]['S1'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['\'S0full\''] = positions_raw[rank_ind]['\'S0full\''].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['\'S1full\''] = positions_raw[rank_ind]['\'S1full\''].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
else:
    postfix = '_'

ax.plot(positions_raw[0]['v0'], 'k')
ax.set_ylim(bottom=0, top = 1)
plt.savefig(out_dir + 'test.png')

for ind in row_indices:
    # we now have one particular timepoint
    time = positions_raw[0].iloc[ind]['time']
    print(time)
    #if not os.path.exists(sys.argv[1] + 'data/phase_p0_' + '{:06.3f}'.format(time) + '.vtu'):
    #    continue
    times.append(time)
    phi_all = []
    v0 = []
    v1 = []
    S0 = []
    S1 = []
    S0full = []
    S1full = []
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

        
        half_roll_window = int(rolling_window/2.0)
        v0_tmp = positions_raw[rank_ind].iloc[ind]['v0']
        v1_tmp = positions_raw[rank_ind].iloc[ind]['v1']
        v0.append(v0_tmp)# / np.sqrt(v0_tmp**2 + v1_tmp**2))
        v1.append(v1_tmp)# / np.sqrt(v0_tmp**2 + v1_tmp**2))
        S0.append(positions_raw[rank_ind].iloc[ind]['S0'])
        S1.append(positions_raw[rank_ind].iloc[ind]['S1'])
        S0full.append(positions_raw[rank_ind].iloc[ind]['\'S0full\''])
        S1full.append(positions_raw[rank_ind].iloc[ind]['\'S1full\''])
    
    # after this we have a list IN THE SAME ORDER AS ranks with all phi
    # now the axis 0 here is the rank axis which we want to remove
    phi_all = np.array(phi_all)
    v0 = np.array(v0)
    v1 = np.array(v1)
    S0 = np.array(S0)
    S1 = np.array(S1)
    S0full = np.array(S0full)
    S1full = np.array(S1full)

    # global phasefield, given by a lot of 1s and something in between

    phi_glob = np.max(phi_all,axis=0)
    # this is more interesting -> this locally gives the rank the contributed to the max, ie the cell
    rank_max = np.argmax(phi_all,axis=0)

    v0_glob = v0[rank_max] * phi_glob
    v1_glob = v1[rank_max] * phi_glob
    norm_v = np.sqrt(v0_glob**2 + v1_glob**2)
    v0_glob = v0_glob# / norm_v
    v1_glob = v1_glob# / norm_v
    S0_glob = S0[rank_max] * phi_glob
    S1_glob = S1[rank_max] * phi_glob
    S0full_glob = S0full[rank_max] * phi_glob
    S1full_glob = S1full[rank_max] * phi_glob
    norm_S = np.sqrt(S0_glob**2 + S1_glob**2)
    S0_glob = S0_glob / norm_S
    S1_glob = S1_glob / norm_S
    np.save(out_dir + '/v0_glob' + postfix + '{:06.3f}'.format(time),v0_glob)
    np.save(out_dir + '/v1_glob' + postfix + '{:06.3f}'.format(time),v1_glob)
    np.save(out_dir + '/S0_glob' + postfix + '{:06.3f}'.format(time),S0_glob)
    np.save(out_dir + '/S1_glob' + postfix + '{:06.3f}'.format(time),S1_glob)
    np.save(out_dir + '/S0full_glob' + postfix + '{:06.3f}'.format(time),S0full_glob)
    np.save(out_dir + '/S1full_glob' + postfix + '{:06.3f}'.format(time),S1full_glob)
    np.save(out_dir + '/phi_glob' + postfix + '{:06.3f}'.format(time),phi_glob)

    np.save(out_dir + '/grid_x_' + '{:06.3f}'.format(time),np.reshape(xx,(interpolation_steps,interpolation_steps)))
    np.save(out_dir + '/grid_y_' + '{:06.3f}'.format(time),np.reshape(yy,(interpolation_steps,interpolation_steps)))
    print('Finished time ' + str(time))
np.save(out_dir + '/timesteps',np.array(times))
