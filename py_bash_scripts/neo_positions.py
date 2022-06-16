"""
This script recalculates the positions and velocities of each cell. This is because while keeping the
center_of_mass parameter in the init file turned off, the center of mass is calculated correctly for 
cells traversing the periodic domain. However, it does so with an accuracy related to the FEM grid 
refinement. So, the update of positions occur only if the cell's COM moves from being close to one grid
point to another. But, the velocities are incorrectly calculated because of the sudden jump in the 
COM when traversing the periodic boundary. In this code, the shape of the cell which is saved in the
VTU files is used to recalculate the COM to arbitrary accuracy without being limited to a grid point
and the velocity is calculated correctly while considering the sudden jumps in position values that 
occur at the boundaries. The limitation of this method is that the position and velocity is saved only
when the VTU data is available.

This code requires: phi_field.npy and timesteps.npy which is generated
                    from the python script "compute_phi_field.py"
"""

import os
import math
import numpy as np
import sys
import glob
import csv
import re
import pandas as pd

pos_input_dir = sys.argv[1]
phi_input_dir = pos_input_dir + '/global_fields_200' #directory of phi_field and timesteps.npy 

result_dir = pos_input_dir

domain_dimension = np.array([100, 100])
delta_t = 0.005 #timestep size in simulations
stride = 100 #number of timesteps between saved vtu files

del_t = delta_t*stride #accounting for both timestep size and the stride

num_cells = int(sys.argv[2])
cells = np.arange(num_cells)

times = np.load(phi_input_dir + '/timesteps.npy')

pos_x = np.zeros([num_cells, len(times)])
pos_y = np.zeros([num_cells, len(times)])
vel_x = np.zeros([num_cells, len(times)])
vel_y = np.zeros([num_cells, len(times)])

for ind_t, time in enumerate(times):
    print(time)
   
    phi_rankwise = np.load(phi_input_dir + '/phi_field' + '{:06.3f}'.format(time) +'.npy')
    try:
        grid_X = np.load(phi_input_dir + '/grid_x_' + '{:06.3f}'.format(time) +'.npy')
        grid_Y = np.load(phi_input_dir + '/grid_y_' + '{:06.3f}'.format(time) +'.npy')
    except:
        grid_X = np.load(phi_input_dir + '/grid_x.npy')
        grid_Y = np.load(phi_input_dir + '/grid_y.npy')
    
    for ind_r, rank in enumerate(cells):
        #creating a phase field 
        phasefield = np.zeros(phi_rankwise.shape)
        phasefield[phi_rankwise==(rank)] = 1
        
        #fixing for periodic BC
        maxX = np.max(grid_X[phasefield>0], initial=-1)
        minX = np.min(grid_X[phasefield>0], initial=101)
        maxY = np.max(grid_Y[phasefield>0], initial=-1)
        minY = np.min(grid_Y[phasefield>0], initial=101)
        #maxX = np.max(grid_X[phasefield>0], initial=0)
        #minX = np.min(grid_X[phasefield>0], initial=0)
        #maxY = np.max(grid_Y[phasefield>0], initial=0)
        #minY = np.min(grid_Y[phasefield>0], initial=0)
        
        grid_Xc = grid_X.copy()
        grid_Yc = grid_Y.copy()
        if((maxX - minX) > 0.5*domain_dimension[0]):
           grid_Xc[grid_Xc<0.5*domain_dimension[0]] += domain_dimension[0]
        if((maxY - minY) > 0.5*domain_dimension[1]):
           grid_Yc[grid_Yc<0.5*domain_dimension[1]] += domain_dimension[1]
        
        #note: below values need to be can go beyond the periodic domain.
        pos_x[ind_r, ind_t] = np.sum((grid_Xc.flatten()*phasefield.flatten()))/np.sum(phasefield)
        pos_y[ind_r, ind_t] = np.sum((grid_Yc.flatten()*phasefield.flatten()))/np.sum(phasefield)

#fixing values such that the positions are within the domain
pos_x[pos_x > domain_dimension[0]] -= domain_dimension[0]
pos_y[pos_y > domain_dimension[0]] -= domain_dimension[0]

vel_x[:, 1:] = np.diff(pos_x, axis = 1)
vel_y[:, 1:] = np.diff(pos_y, axis = 1)

#accounting for jump at boundary for positive values
vel_x[vel_x > 0.25*domain_dimension[0]] -= domain_dimension[0]
vel_y[vel_y > 0.25*domain_dimension[1]] -= domain_dimension[1]

#accounting for jump at boundary for negative values
vel_x[vel_x < -0.25*domain_dimension[0]] += domain_dimension[0]
vel_y[vel_y < -0.25*domain_dimension[1]] += domain_dimension[1]

#fixing for timestep size
vel_x /= del_t
vel_y /= del_t

# Reading positions raw to get other variables
ranks = []
positions_raw = []
file_pattern = pos_input_dir + '/positions_p*.csv'
for filename in glob.glob(file_pattern):
    tmp = pd.read_csv(filename)
    positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', 'vel_angle', 'nematic_angle', 'total_interaction', 'neighbours', 'confine_interaction', 'growth_rate', 'S0full', 'S1full']])
    # we also need to extract the rank to build the bridge to vtk files
    ranks.append(int(re.findall(r'\d+', filename)[-1]))
stride = 100 #time gap between vtu data
ranks = np.array(ranks)
sorted_index = np.argsort(ranks)
def sort_ranks(elem):
    return elem.iloc[0]['rank']
positions_raw.sort(key=sort_ranks)
#positions_raw = positions_raw[sorted_index]
ranks = ranks[sorted_index]

#rolling_average
rolling_average = True
rolling_window  = 100
if(rolling_average):
    print('rolling average performed over a window: ' + str(rolling_window))
    postfix = '_roll_' + str(rolling_window) + '_'
    for rank_ind,rank in enumerate(ranks):
        positions_raw[rank_ind]['S0'] = positions_raw[rank_ind]['S0'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S1'] = positions_raw[rank_ind]['S1'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['growth_rate'] = positions_raw[rank_ind]['growth_rate'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['total_interaction'] = positions_raw[rank_ind]['total_interaction'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['confine_interaction'] = positions_raw[rank_ind]['confine_interaction'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S0full'] = positions_raw[rank_ind]['S0full'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S1full'] = positions_raw[rank_ind]['S1full'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
else:
    postfix = '_'

#creating new database
for ind_c, cell in enumerate(cells):
    time_ = np.array(positions_raw[cell]['time'][stride-1::stride])
    x0_ = pos_x[ind_c]
    x1_ = pos_y[ind_c]
    v0_ = vel_x[ind_c]
    v1_ = vel_y[ind_c]
    S0_ = np.array(positions_raw[cell]['S0'][stride-1::stride])
    S1_ = np.array(positions_raw[cell]['S1'][stride-1::stride])
    S0full_ = np.array(positions_raw[cell]['S0full'][stride-1::stride])
    S1full_ = np.array(positions_raw[cell]['S1full'][stride-1::stride])
    total_interaction_ = np.array(positions_raw[cell]['total_interaction'][stride-1::stride])
    confine_interaction_ = np.array(positions_raw[cell]['confine_interaction'][stride-1::stride])
    growth_rate_ = np.array(positions_raw[cell]['growth_rate'][stride-1::stride])
    r_ = np.array(positions_raw[cell]['r'][stride-1::stride])
    neighbours_ = np.array(positions_raw[cell]['neighbours'][stride-1::stride])
    rank_ = np.ones(len(x0_), dtype = 'int')*cell
    
    data = pd.DataFrame({ "time": time_,
                        "x0": x0_,
                        "x1": x1_,
                        "v0": v0_,
                        "v1": v1_,
                        "S0": S0_,
                        "S1": S1_,
                        "S0full": S0full_,
                        "S1full": S1full_,
                        "total_interaction": total_interaction_,
                        "confine_interaction": confine_interaction_,
                        "growth_rate": growth_rate_,
                        "r": r_,
                        "neighbours": neighbours_,
                        "rank": rank_},
                       index = time_)
    data.to_csv(pos_input_dir + '/neo_positions_p' + str(cell) + '.csv')