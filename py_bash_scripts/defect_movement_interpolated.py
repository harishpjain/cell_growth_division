import matplotlib.pyplot as plt
import sys
import numpy as np
import os
import re
import csv
import glob
import pandas as pd
#import vtk
#from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata
from matplotlib import cm

#import progress
from scipy import ndimage
import seaborn as sns; sns.set_theme()
from shapely.geometry import LineString
from matplotlib.collections import LineCollection


fig,ax = plt.subplots()

# --- Parameters
start_time = 1#0
end_time = 10
stride = 1

domain_size = np.array([100.0,100.0])

interpolation_steps = 1000
global_smooth_windows = 10


# Previous : 40steps, radius 2 and 20 rings
local_steps = 60
local_smooth_window = 10

local_radius = 6
local_rings = 60

plot_lines=True
plot_scale_interp = 0.25
plot_stride = 15

gauss_sigma = 20.0

minus = True
# --- Parameters

if len(sys.argv) > 1:
    file_pattern = sys.argv[1] + 'positions_p*.csv'
    if minus:
        out_dir = sys.argv[1] + 'minus_strain_fields'
    else:
        out_dir = sys.argv[1] + 'strain_fields'

if len(sys.argv) > 2:
    start_time = int(sys.argv[2])
    out_dir += '_start_' + str(start_time)
if len(sys.argv) > 3:
    end_time = int(sys.argv[3])
    out_dir += '_end_' + str(end_time)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
with open(out_dir+'/plus_defect_angles.csv', mode='w') as csvfile:
    csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csv_writer.writerow(['Frame','x','y','z','angle'])
    
positions_raw = []
ranks = []

for filename in glob.glob(file_pattern):
    tmp = pd.read_csv(filename)
    positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', ' vel_angle', ' nematic_angle', ' total_interaction', ' neighbours', ' confine_interaction', ' growth_rate']])
    # we also need to extract the rank to build the bridge to vtk files
    ranks.append(int(re.findall(r'\d+', filename)[-1]))
    
# Interpolation
x = np.linspace(0,domain_size[0],interpolation_steps)
y = np.linspace(0,domain_size[1],interpolation_steps)
xx,yy = np.meshgrid(x,y)
xx_vec = np.reshape(xx,interpolation_steps**2)
yy_vec = np.reshape(yy,interpolation_steps**2)

print('checkpouint 1')
count_plus = 0
count_frames = 0
for time in range(start_time,end_time):
    # Empty lists for positions, velocities and elongations
    x0 = []
    x1 = []

    S0 = []
    S1 = []

    for rank_ind,rank in enumerate(ranks):
        x0.append(positions_raw[rank_ind].iloc[time]['x0'])
        x1.append(positions_raw[rank_ind].iloc[time]['x1'])
        
        S0.append(positions_raw[rank_ind].iloc[time]['S0'])
        S1.append(positions_raw[rank_ind].iloc[time]['S1'])


    x0 = np.array(x0)
    x1 = np.array(x1)

    S0 = np.array(S0)
    S1 = np.array(S1)

    # Rewrite x0,x1 for periodic boundaries
    positions = np.stack((x0,x1),axis=1)
    positions = np.tile(positions,(5,1))
    positions[len(ranks):2*len(ranks),0] -= 100
    positions[2*len(ranks):3*len(ranks),0] += 100
    positions[3*len(ranks):4*len(ranks),1] -= 100
    positions[4*len(ranks):5*len(ranks),1] += 100

    # Copy v and S to match the positions
    
    S0 = np.tile(S0,5)
    S1 = np.tile(S1,5)

    # Construct global velocity and elongation fields by simple interpolation


    S0_glob = griddata(positions,S0,(xx_vec,yy_vec), method='linear')
    S1_glob = griddata(positions,S1,(xx_vec,yy_vec), method='linear')

    norm_S = np.sqrt(S0_glob**2 + S1_glob**2)
    S0_glob = S0_glob / norm_S
    S1_glob = S1_glob / norm_S

    S0_glob = np.reshape(S0_glob, (interpolation_steps,interpolation_steps))
    S1_glob = np.reshape(S1_glob, (interpolation_steps,interpolation_steps))

    # Smooth
    S0_glob = ndimage.uniform_filter(S0_glob, size=global_smooth_windows, mode='wrap')

    # Smooth
    S1_glob = ndimage.uniform_filter(S1_glob, size=global_smooth_windows, mode='wrap')

    print('smooth s field created')
    # Identify 0 contours of S0 and S1
    c1 = ax.contour(xx,yy,S0_glob,0,colors='red',alpha=0.5)
    c2 = ax.contour(xx,yy,S1_glob,0,colors='blue',alpha=0.5)

    defect_list = []

    c1 = c1.allsegs[0]
    c2 = c2.allsegs[0]       
    # Find intersections between zero contours
    for m in range(len(c1)):
        for n in range(len(c2)):
            line1 = LineString(c1[m])
            line2 = LineString(c2[n])

            # shapely Geometry Sequence
            if line1.intersects(line2):
                intersects = line1.intersection(line2)
                try:
                    intersects = intersects.geoms

                    for k in intersects:
                        defect_list.append(np.array([k.x,k.y]))
                except:
                    defect_list.append(np.array([intersects.x,intersects.y]))        
        

    grid = np.stack((xx_vec,yy_vec),axis=1)
    defect_array = np.array(defect_list)

    S0_dy,S0_dx = np.gradient(S0_glob)
    S0_dx = np.reshape(S0_dx,(interpolation_steps**2,1))
    S0_dy = np.reshape(S0_dy,(interpolation_steps**2,1))    

    S1_dy,S1_dx = np.gradient(S1_glob)
    S1_dx = np.reshape(S1_dx,(interpolation_steps**2,1))
    S1_dy = np.reshape(S1_dy,(interpolation_steps**2,1))       

    # Compute all defect angles
    print(S1_glob)
    omega = np.multiply(np.sign(S1_glob),(np.arctan(S0_glob/np.abs(S1_glob))/ 2.0 + np.pi/4.0))
    n1 = np.sin(omega)
    n2 = np.cos(omega)
    
    n1 = np.reshape(n1,(interpolation_steps,interpolation_steps))
    n2 = np.reshape(n2,(interpolation_steps,interpolation_steps))

    _,n1Sq_dx = np.gradient(n1**2)
    n2Sq_dy,_ = np.gradient(n2**2)
    n1n2_dy,n1n2_dx =  np.gradient(n1*n2)

    div_n_0 = n1Sq_dx + n1n2_dy
    div_n_1 = n1n2_dx + n2Sq_dy

    div_n_0 = np.reshape(div_n_0,(np.shape(div_n_0)[0]**2,1))
    div_n_1 = np.reshape(div_n_1,(np.shape(div_n_1)[0]**2,1))

    norm_div = np.sqrt(div_n_0**2+div_n_1**2)
    div_n_0 = div_n_0 / norm_div
    div_n_1 = div_n_1 / norm_div
    angles = np.arctan2(div_n_1,div_n_0)

    for coord in defect_list:

        if np.sum(np.linalg.norm(coord-defect_array,axis=1)<int(local_radius/np.sqrt(2)))>1:
            continue

        ind_local = np.unravel_index(np.linalg.norm(np.abs(grid-coord),axis=1).argmin(), np.linalg.norm(np.abs(grid-coord),axis=1).shape)
        delta = S0_dx[ind_local] * S1_dy[ind_local] - S0_dy[ind_local] * S1_dx[ind_local]    
          
        if minus:
            if delta > 0:
                continue
        else:
            if delta < 0:
                continue
        
        angle = np.mean(angles[ind_local])     
        
        with open(out_dir+'/plus_defect_angles.csv', mode='a') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            csv_writer.writerow([count_frames,coord[0],coord[1],0,angle])   
                  
    count_frames += 1       
