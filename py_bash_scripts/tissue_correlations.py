import os
import math
import numpy as np
import pandas as pd
import sys
import glob
import csv
import re
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors 
from matplotlib.collections import LineCollection
plt.rcParams['figure.figsize'] = (12,12)

input_dir_pos = sys.argv[1]#'i3c4a1_v1_0_neo' 
input_dir = input_dir_pos + '/global_fields_200'
result_dir = input_dir + '/correlations'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
else: 
    print('nothing')

times = np.load(input_dir + '/timesteps.npy')
#times = np.array([1.5, 2.5])#, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])#
jump = 1
times = times[2::jump]
stride = 100*jump #time gap between vtu data
postfix = '_stride_' + str(stride)
print(times)

ranks = []
positions_raw = []
file_pattern = input_dir_pos + '/positions_p*.csv'
for filename in glob.glob(file_pattern):
    tmp = pd.read_csv(filename)
    positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', 'vel_angle', 'nematic_angle', 'total_interaction', 'neighbours', 'confine_interaction', 'growth_rate', 'S0full', 'S1full']])
    # we also need to extract the rank to build the bridge to vtk files
    ranks.append(int(re.findall(r'\d+', filename)[-1]))
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
    postfix += '_roll_' + str(rolling_window)
    for rank_ind,rank in enumerate(ranks):
        positions_raw[rank_ind]['v0'] = positions_raw[rank_ind]['v0'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['v1'] = positions_raw[rank_ind]['v1'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S0'] = positions_raw[rank_ind]['S0'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S1'] = positions_raw[rank_ind]['S1'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S0full'] = positions_raw[rank_ind]['S0full'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
        positions_raw[rank_ind]['S1full'] = positions_raw[rank_ind]['S1full'].rolling(window = rolling_window, win_type = 'gaussian', closed = 'neither').mean(std=rolling_window)
else:
    postfix += ''

grid_X = np.load(input_dir + '/grid_x_' + '{:06.3f}'.format(times[0]) +'.npy')
grid_Y = np.load(input_dir + '/grid_y_' + '{:06.3f}'.format(times[0]) +'.npy')

velx = np.zeros(np.hstack([(len(times)),grid_X.shape]))
vely = np.zeros(np.hstack([(len(times)),grid_X.shape]))
S0 = np.zeros(np.hstack([(len(times)),grid_X.shape]))
S1 = np.zeros(np.hstack([(len(times)),grid_X.shape]))


for ind_t, time in enumerate(times):
    phi_rankwise = np.load(input_dir + '/phi_field' + '{:06.3f}'.format(time) +'.npy')
    time_index = int(time/0.005)-1
    for rank_ind,rank in enumerate(ranks):
        velx[ind_t][phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['v0']
        vely[ind_t][phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['v1']
        
        S0[ind_t][phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['S0full']
        S1[ind_t][phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['S1full']

# Calculating temporal mean and variances
mean_velx = np.mean(velx, axis = 0)
mean_vely = np.mean(vely, axis = 0)

var_velx = np.var(velx, axis = 0)
var_vely = np.var(vely, axis = 0)

mean_S0 = np.mean(S0, axis = 0)
mean_S1 = np.mean(S1, axis = 0)

var_S0 = np.var(S0, axis = 0)
var_S1 = np.var(S1, axis = 0)

tau_arr = np.array([0.5, 1.0])
tau_arr = np.arange(1,30,1)*0.5
vel_corr_arr = np.zeros(len(tau_arr)+1)
S_corr_arr = np.zeros(len(tau_arr)+1)
delt = 0.5 #time difference between collected data points
tau_jump = np.divide(tau_arr, delt).astype('int')

deviation_velx = velx - mean_velx
deviation_vely = vely - mean_vely

deviation_S0 = S0 - mean_S0
deviation_S1 = S1 - mean_S1

T = (times[-1] - times[0])/(times[1]-times[0]) #number of timepoints
print(tau_jump)
for ind_tau, jump in enumerate(tau_jump):
    convolution_velx = np.sum(np.multiply(deviation_velx[:-jump:jump], deviation_velx[jump::jump]), axis=0)
    convolution_vely = np.sum(np.multiply(deviation_vely[:-jump:jump], deviation_vely[jump::jump]), axis=0)
    
    point_corr_vel = (1/(2*T))*(convolution_velx/var_velx + convolution_vely/var_vely)
    vel_corr_arr[ind_tau+1] = np.mean(point_corr_vel)
    
    convolution_S0 = np.sum(np.multiply(deviation_S0[:-jump:jump], deviation_S0[jump::jump]), axis=0)
    convolution_S1 = np.sum(np.multiply(deviation_S1[:-jump:jump], deviation_S1[jump::jump]), axis=0)
    
    point_corr_S = (1/(2*T))*(convolution_S0/var_S0 + convolution_S1/var_S1)
    S_corr_arr[ind_tau+1] = np.mean(point_corr_S)

vel_corr_arr[0] = 1.0
S_corr_arr[0] = 1.0

velx_conv_zero = np.sum(np.multiply(deviation_velx, deviation_velx), axis = 0)
vely_conv_zero = np.sum(np.multiply(deviation_vely, deviation_vely), axis = 0)

vel_corr_arr[0] = np.mean((1/(2*T))*(velx_conv_zero/var_velx + vely_conv_zero/var_vely))

S0_conv_zero = np.sum(np.multiply(deviation_S0, deviation_S0), axis = 0)
S1_conv_zero = np.sum(np.multiply(deviation_S1, deviation_S1), axis = 0)

S_corr_arr[0] = np.mean((1/(2*T))*(S0_conv_zero/var_S0 + S1_conv_zero/var_S1))

tau_arr = np.hstack([0.0, tau_arr])

fig, ax = plt.subplots()
ax.plot(tau_arr, vel_corr_arr, 'b*-', label = '$C_v:$ Velocity Autocorrelation')
ax.plot(tau_arr, S_corr_arr, 'r^-', label = '$C_S:$ Nematic Autocorrelation')
ax.set_xlabel('Time')
ax.set_ylabel('Temporal Auto Correlation')
ax.legend()
plt.savefig(result_dir + '/autocorrelation.png', dpi=300)

np.save(result_dir + '/tau', tau_arr)
np.save(result_dir + '/vel_corr', vel_corr_arr)
np.save(result_dir + '/S_corr', S_corr_arr)