import os
import math
import numpy as np
import sys
import glob
import csv
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors 
from scipy.interpolate import griddata
from scipy.interpolate import interpn
from scipy.ndimage import fourier_gaussian
from scipy.ndimage import gaussian_filter
from scipy.ndimage import uniform_filter
from scipy.ndimage import spline_filter
from shapely.geometry import LineString
from matplotlib.collections import LineCollection


def defect_analyse_func(plus_def, omega, steps_x, steps_y, grid_X, grid_Y, grid_x_vec, grid_y_vec, grid, count, v0_glob, v1_glob, domain_dimension):
    #closest point to defect on a flattened grid
    close_index = np.argmin(np.linalg.norm(grid-plus_def,axis=1))
    close_point = grid[close_index]
    print('Defect is at', plus_def)
    #print('Closest grid point is', close_point)
    #print('Closest grid point index in flattened grid is', close_index)
    x_index = int(np.floor(close_index/steps_x))
    y_index = close_index % steps_y
    #print('Closest grid point index on square grid is', x_index, y_index)
    #print('Closest grid point(square) ', grid_X[x_index, y_index], grid_Y[x_index, y_index])
    
    #considering a small box around the defect
    left_ind = 1
    index_array_0 = np.arange(x_index-left_ind,x_index+left_ind+1)
    index_array_1 = np.arange(y_index-left_ind,y_index+left_ind+1)
    index_array_0[index_array_0 < 0] += steps_x
    index_array_1[index_array_1 < 0] += steps_y
    index_array_0[index_array_0 > steps_x-1] -= steps_x
    index_array_1[index_array_1 > steps_y-1] -= steps_y
    index_grid = np.meshgrid(index_array_0, index_array_1)
    
    
    grid_X_loc = grid_X[index_grid]
    grid_Y_loc = grid_Y[index_grid]
    
    omega_loc = omega[index_grid]
    v0_loc = v0_glob[index_grid]
    v1_loc = v1_glob[index_grid]
    #nematic director field plotting
    directions = np.zeros((omega_loc.shape[0]**2,2,2))
    directions[:,0,0] = grid_X_loc.ravel().squeeze()-0.05 * np.sin(omega_loc.ravel()).squeeze()
    directions[:,1,0] = grid_X_loc.ravel().squeeze()+0.05 * np.sin(omega_loc.ravel()).squeeze()
    directions[:,0,1] = grid_Y_loc.ravel().squeeze()-0.05 * np.cos(omega_loc.ravel()).squeeze()
    directions[:,1,1] = grid_Y_loc.ravel().squeeze()+0.05 * np.cos(omega_loc.ravel()).squeeze()
    
    #finding the distance to director which has the closest approach to defect
    closest_distance = np.zeros(omega_loc.shape[0]**2)
    p2_minus_p1 = directions[:,1,:] - directions[:,0,:]
    p3_minus_p1 = - directions[:,0,:] + plus_def
    closest_distance = np.abs(np.cross(p2_minus_p1, p3_minus_p1)/np.linalg.norm(p2_minus_p1, axis=1))
    close_dist_index = np.argmin(closest_distance)
    close_dist_point = np.array([grid_X_loc.ravel()[close_dist_index], grid_Y_loc.ravel()[close_dist_index]])
    angle_close = np.arctan2(close_dist_point[1]-plus_def[1],close_dist_point[0]-plus_def[0])
    defect_angle = omega_loc.ravel()[close_dist_index]
    veldef_0 = v0_loc.ravel()[close_dist_index]
    veldef_1 = v1_loc.ravel()[close_dist_index]
    a = np.array([close_dist_point[0]-plus_def[0],close_dist_point[1]-plus_def[1]])
    b = np.array([np.sin(defect_angle), np.cos(defect_angle)])
    projection_3 = np.dot(a, b)
    
    #fixing the direction of defect
    if(projection_3>0):
        defect_angle = defect_angle+np.pi
    
    # so far the defect angle is measured with respect to the vertical axis such that it is positive 
    # in clockwise direction. The range is from -pi to pi
    # Converting the angle to range 0 to 2pi such that angle is measured with respect to the horizontal
    # axis and is positive in counter clockwise direction
    defect_angle = (np.pi/2 - defect_angle)
    if(defect_angle < 0):
        defect_angle += 2*np.pi

    defectorient_x = np.cos(defect_angle)
    defectorient_y = np.sin(defect_angle)
    
    vel_angle = np.arctan2(veldef_1, veldef_0)
    if(vel_angle < 0):
        vel_angle += 2*np.pi
    angle_diff = np.abs(vel_angle-defect_angle)
    print('defect angle is ', defect_angle*360/(2*np.pi))
    print('velocity direction is ', vel_angle*360/(2*np.pi))
    print(angle_diff*360/(2*np.pi))
    dot_prod = np.dot(np.array([np.cos(defect_angle), np.sin((defect_angle))]),
                     np.array([np.cos(vel_angle), np.sin((vel_angle))]))
    if(dot_prod > 0):
        print('extensile')
    else:
        print('contractile')
    #recovering angle between defect and velocity
    diff_angle = np.arccos(dot_prod)
    """
    #rotated velocity plots
    omega_rot = omega - defect_angle
    omega_rot_loc = omega_rot[x_index-left_ind:x_index+left_ind+1, y_index-left_ind:y_index+left_ind+1]
    v0_rot_loc, v1_rot_loc = rotate_vector(v0_glob[x_index-left_ind:x_index+left_ind+1, y_index-left_ind:y_index+left_ind+1], 
                                  v1_glob[x_index-left_ind:x_index+left_ind+1, y_index-left_ind:y_index+left_ind+1],
                                  -defect_angle)
    directions_rot = np.zeros((omega_rot_loc.shape[0]**2,2,2))
    directions_rot[:,0,0] = grid_X_loc.ravel().squeeze()-0.05 * np.sin(omega_rot_loc.ravel()).squeeze()
    directions_rot[:,1,0] = grid_X_loc.ravel().squeeze()+0.05 * np.sin(omega_rot_loc.ravel()).squeeze()
    directions_rot[:,0,1] = grid_Y_loc.ravel().squeeze()-0.05 * np.cos(omega_rot_loc.ravel()).squeeze()
    directions_rot[:,1,1] = grid_Y_loc.ravel().squeeze()+0.05 * np.cos(omega_rot_loc.ravel()).squeeze()
    
    fig3, ax3 = plt.subplots()
    line_segments_rot = LineCollection(directions_rot,colors=cm.jet(0.0))
    ax3.add_collection(line_segments_rot)
    ax3.quiver(grid_X_loc, grid_Y_loc, v0_rot_loc, v1_rot_loc)
    
    """
    
    plt.close()
    
    return defect_angle, diff_angle

def roundCoordinates(center, radius, nRings, initialAngle=0, N=5):
    radii = np.linspace(0,radius,nRings)
    x = np.array([])
    y = np.array([])

    for ind_radius,r in enumerate(radii):
        nPoints = ind_radius * N + 1
        theta = np.linspace(initialAngle,initialAngle + 2*np.pi,nPoints)
        x = np.concatenate((x,center[0] + r * np.cos(theta[:-1])))
        y = np.concatenate((y,center[1] + r * np.sin(theta[:-1])))

    return x,y

def defect_angle_neo(defect_center, x_loc, y_loc, omega_loc):
    print(x_loc.shape)
    x_loc = x_loc.ravel()
    y_loc = y_loc.ravel()
    
    #Note: Omega is an angle measured with respect to the vertical line
    #because x = sin(omega) and y = cos(omega) is used to plot the director field
    omega_loc = omega_loc.ravel()
    #rescaling omega from range -pi, pi to 0, 2pi
    #omega_loc[omega_loc<0] += 2*np.pi
    
    #to indicate whether the points are above or below the defect, or left or right of defect
    y_sign = np.sign(y_loc - defect_center[1]) 
    x_sign = np.sign(x_loc - defect_center[0])
    #print(y_sign.shape)
    
    testing_number = 100
    test_angles = np.linspace(-np.pi, np.pi, testing_number)
    
    dot_product_sum = np.zeros(testing_number)
    for ind, angle in enumerate(test_angles):
        dot_product_sum[ind] = (np.sum(np.cos(omega_loc - angle)))
        
    #note the angle opposite to candidate_angle is also a candidate
    candidate_angle = test_angles[np.argmax(dot_product_sum)]
    candidate_projections = np.cos(omega_loc - angle)
    y_contribution = np.sum(np.dot(candidate_projections, y_sign))
    x_contribution = np.sum(np.dot(candidate_projections, x_sign))
    
    if(y_contribution<0):
        candidate_angle += np.pi
        if(candidate_angle > np.pi):
            candidate_angle - 2*np.pi
        
    return candidate_angle



input_dir_pos = sys.argv[1]#'i3c4a1_v1_0_neo' 
input_dir = input_dir_pos + '/global_fields_200'
result_dir = input_dir + '/nematics'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
else: 
    print('nothing')
    #sys.exit('Input directory!')
plot_scale_data = 0.7
plot_scale_interp = 0.9
plot_steps = 50
domain_dimension = np.array([100.0,100.0])
confined = False

field_smoothing = 4.0

# if multiple defects in this range, then ignore them.
# it's because these defects might try to annihilate each other
local_radius = 5 

#for interpolating to calculate 
local_steps = 50

if confined:
    mode = 'nearest'
else:
    mode = 'wrap'

times = np.load(input_dir + '/timesteps.npy')
#times = np.array([1.5, 2.5])#, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])#
jump = 10
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

count = 0 #number of defects
#storing the summed up contributions around a defect
Eyy_sum = np.zeros([20, 20])
Exy_sum = np.zeros([20, 20])
count_plus_defects = 0

plus_def_data = []
minus_def_data = []

def_angles = [] 
diff_angles = []
plusdef_dict = {}

#here I will count the nature of elongation around the defects and 
#see correlation between elongation and extensility
is_extensile_and_elongated = 0
is_extensile_and_isotropic = 0
is_contractile_and_elongated = 0
is_contractile_and_isotropic = 0

S_thres_array = np.array([7.0, 8.0, 9.0, 10.0])
is_ext_and_elon = np.zeros(len(S_thres_array))
is_ext_and_iso = np.zeros(len(S_thres_array))
is_con_and_elon = np.zeros(len(S_thres_array))
is_con_and_iso = np.zeros(len(S_thres_array))

for time in times:
    time_index = int(time/0.005)-1
    #For each timestep, we store a list of lists where the
    #lowest list stores the coordinates of the defect,
    #the orientation of the defect and difference in 
    #angle with the velocity vector
    plusdef_dict[time] = []
    
    fig = plt.figure()
    #fig.set_size_inches(6.0,6.0)
    ax = fig.add_subplot(1,1,1)
    ax.axis([0,100,0,100])

    ax.set_yticks([])
    ax.set_xticks([])

    #bounding box
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.5)
        ax.spines[axis].set_color('black')
    
    #loading data
    #S0_glob = np.load(input_dir + '/S0_glob_roll_100_' + '{:06.3f}'.format(time) +'.npy')
    #S1_glob = np.load(input_dir + '/S1_glob_roll_100_' + '{:06.3f}'.format(time) +'.npy')
    #v0_glob = np.load(input_dir + '/v0_glob_roll_100_' + '{:06.3f}'.format(time) +'.npy')
    #v1_glob = np.load(input_dir + '/v1_glob_roll_100_' + '{:06.3f}'.format(time) +'.npy')
    phi_rankwise = np.load(input_dir + '/phi_field' + '{:06.3f}'.format(time) +'.npy')
    grid_X = np.load(input_dir + '/grid_x_' + '{:06.3f}'.format(time) +'.npy')
    grid_Y = np.load(input_dir + '/grid_y_' + '{:06.3f}'.format(time) +'.npy')
    
    v0_glob = np.zeros(phi_rankwise.shape)
    v1_glob = np.zeros(phi_rankwise.shape)
    S0_glob = np.zeros(phi_rankwise.shape)
    S1_glob = np.zeros(phi_rankwise.shape)
    S0full = np.zeros(phi_rankwise.shape)
    S1full = np.zeros(phi_rankwise.shape)
    
    for rank_ind,rank in enumerate(ranks):
        v0_glob[phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['v0']
        v1_glob[phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['v1']
        S0_glob[phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['S0']
        S1_glob[phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['S1']
        S0full[phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['S0full']
        S1full[phi_rankwise==(rank)] = positions_raw[rank_ind].iloc[time_index]['S1full']
    
    #ax.contourf(grid_X,grid_Y,S0_glob,alpha=0.5)

    # Smooth field 
    S0_glob = gaussian_filter(S0_glob, field_smoothing, mode=mode)
    S1_glob = gaussian_filter(S1_glob, field_smoothing, mode=mode)
    v0_glob = gaussian_filter(v0_glob, field_smoothing, mode=mode)
    v1_glob = gaussian_filter(v1_glob, field_smoothing, mode=mode)
    
    S_mag = np.linalg.norm(np.array([S0full, S1full]), axis=0)
    ax.contourf(grid_X,grid_Y,S_mag,alpha=0.5)
    
    steps_x,steps_y = np.shape(S0_glob)
    grid_x_vec = np.reshape(grid_X, (steps_x*steps_y,1))
    grid_y_vec = np.reshape(grid_Y, (steps_x*steps_y,1))

    #coarser grid for plotting directors
    coarse_steps = int(steps_x/plot_steps)
    grid_Xc = grid_X[::coarse_steps, ::coarse_steps]
    grid_Yc = grid_Y[::coarse_steps, ::coarse_steps]
    S0_coarse = S0_glob[::coarse_steps, ::coarse_steps]
    S1_coarse = S1_glob[::coarse_steps, ::coarse_steps]
    
    # phase angle
    #omega = np.multiply(np.sign(S1_glob),(np.arctan(S0_glob/np.abs(S1_glob))/ 2.0 + np.pi/4.0))
    omega_c = np.multiply(np.sign(S1_coarse),(np.arctan(S0_coarse/np.abs(S1_coarse))/ 2.0 + np.pi/4.0))
    omega_c = omega_c.ravel()

    #nematic director field plotting
    directions = np.zeros((plot_steps**2,2,2))
    directions[:,0,0] = grid_Xc.ravel().squeeze()-0.5*plot_scale_interp * np.sin(omega_c).squeeze()
    directions[:,1,0] = grid_Xc.ravel().squeeze()+0.5*plot_scale_interp * np.sin(omega_c).squeeze()
    directions[:,0,1] = grid_Yc.ravel().squeeze()-0.5*plot_scale_interp * np.cos(omega_c).squeeze()
    directions[:,1,1] = grid_Yc.ravel().squeeze()+0.5*plot_scale_interp * np.cos(omega_c).squeeze()

    line_segments = LineCollection(directions,linewidths=plot_scale_data,colors=cm.jet(0.0))
    ax.add_collection(line_segments)
    
    # Identify 0 contours of S0 and S1
    c1 = ax.contour(grid_X,grid_Y,S0_glob,0,colors='red',alpha=0.5)
    c2 = ax.contour(grid_X,grid_Y,S1_glob,0,colors='blue',alpha=0.5)

    #fig.savefig(result_dir+'/defect_locations_'+str(count).zfill(4)+'.png', dpi=300, bbox_inches='tight')
    #print(np.asarray(c1.allsegs).shape)
    #print(len(c1.allsegs[0]))
    #if len(c1.allsegs[0]) <2 or len(c2.allsegs[0])<2:
    #    print('continuing')
    #    continue

    #print(c1.allsegs)
    defect_list = []
    print('contours identified time ' + str(time))    
    c1 = c1.allsegs[1]
    c2 = c2.allsegs[1]   

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

    grid = np.concatenate((grid_x_vec,grid_y_vec),axis=1)
    
    #grid = np.vstack((grid_x_vec,grid_y_vec))

    tmp = np.gradient(S0_glob)
    S0_dx = np.reshape(tmp[0],(steps_x*steps_y,1))
    S0_dy = np.reshape(tmp[1],(steps_x*steps_y,1))

    tmp = np.gradient(S1_glob)
    S1_dx = np.reshape(tmp[0],(steps_x*steps_y,1))
    S1_dy = np.reshape(tmp[1],(steps_x*steps_y,1))

    plus_defects = []
    minus_defects = []
    #print(len(defect_list))
    
    defect_list_isolated = []
    defect_array = np.array(defect_list)
    for coord in defect_list:
        if np.sum(np.linalg.norm(coord-defect_array,axis=1)<2*local_radius)>1:
            pass
        else:
            defect_list_isolated.append(coord)
    #print(len(defect_list_isolated))
            
    if len(defect_list_isolated) > 0:
        for coord in defect_list_isolated:
            #finding the nearest index on grid
            ind = np.linalg.norm(np.abs(grid-coord),axis=1).argmin()
            delta = S0_dy[ind] * S1_dx[ind] - S0_dx[ind] * S1_dy[ind]
            if delta > 0:
                plus_defects.append(coord)
            else:
                minus_defects.append(coord)

        plus_defects = np.array(plus_defects)
        minus_defects = np.array(minus_defects)            

        if np.shape(plus_defects)[0] > 0:
            ax.scatter(plus_defects[:,0],plus_defects[:,1],color='green',alpha=0.7)
        if np.shape(minus_defects)[0] > 0:
            ax.scatter(minus_defects[:,0],minus_defects[:,1],color='purple',alpha=0.5)  
             
    
    for plus_def in plus_defects:
            plus_def_data.append([time, plus_def[0],plus_def[1]])
    for minus_def in minus_defects:
            minus_def_data.append([time, minus_def[0], minus_def[1]])

    count += 1
    
    #computing defect orientation
    #first we take small box around the defect aligned with the domain coordinates
    defect_directions = np.zeros((len(plus_defects),2,2))
    defect_directions_2 = np.zeros((len(plus_defects),2,2))
    
    #director angle
    omega = np.multiply(np.sign(S1_glob),(np.arctan(S0_glob/np.abs(S1_glob))/ 2.0 + np.pi/4.0))
    n1 = np.sin(omega)
    n2 = np.cos(omega)

    #n1 = np.reshape(n1,(local_steps,local_steps))
    #n2 = np.reshape(n2,(local_steps,local_steps))

    _,n1Sq_dx = np.gradient(n1**2)
    n2Sq_dy,_ = np.gradient(n2**2)
    n1n2_dy,n1n2_dx =  np.gradient(n1*n2)

    div_n_0 = n1Sq_dx + n1n2_dy
    div_n_1 = n1n2_dx + n2Sq_dy
    print(div_n_0.shape)
    #div_n_0 = np.reshape(div_n_0,(np.shape(div_n_0)[0]**2,1))
    #div_n_1 = np.reshape(div_n_1,(np.shape(div_n_1)[0]**2,1))

    norm_div = np.sqrt(div_n_0**2+div_n_1**2)
    div_n_0 = div_n_0 / norm_div
    div_n_1 = div_n_1 / norm_div
    angles = np.arctan2(div_n_1,div_n_0)
    #print('max angles: ', np.max(angles))
    #print(angles.shape)
    
    v0_grad_x, v0_grad_y = np.gradient(v0_glob, domain_dimension[0]/steps_x)
    v1_grad_x, v1_grad_y = np.gradient(v1_glob, domain_dimension[0]/steps_x)
    #print(v0_grad_x.shape)
    
    E_xx = v0_grad_x
    E_yy = v1_grad_y
    E_xy = 0.5*(v0_grad_y + v1_grad_x)
    
    for pd, (plus_def) in enumerate(plus_defects):
        defect_angle, diff_angle = defect_analyse_func(plus_def, omega, steps_x, steps_y,
                                                                 grid_X, grid_Y, grid_x_vec, grid_y_vec,
                                                                 grid, count_plus_defects, v0_glob, v1_glob,
                                                                domain_dimension)
        #print(v0_glob.shape)
        #print('ok')
        #defect_analyse_velocity(plus_def, defect_angle, steps_x, steps_y, grid_X, grid_Y, grid_x_vec, grid_y_vec, grid, count, v0_glob, v1_glob)
        def_angles.append(defect_angle)
        diff_angles.append(diff_angle) 
        plusdef_dict[time].append([plus_def, defect_angle, diff_angle])
        
        ax.quiver(plus_def[0], plus_def[1], np.cos(defect_angle), np.sin(defect_angle), color='green');
        count_plus_defects += 1
    fig.savefig(result_dir+'/defect_locations_'+str(time_index+1).zfill(4)+'.png', dpi=300, bbox_inches='tight') 
    #plt.close()
    fig.clear()
    plt.close(fig)
    plusdef_info = plusdef_dict[time]
    
    cell_S0full = np.zeros(len(ranks))
    cell_S1full = np.zeros(len(ranks))
    dist_to_defect = np.zeros([len(plusdef_info), len(ranks)])

    for rank_ind,rank in enumerate(ranks):
        #print('rank is ', rank, ' and its location is', positions_raw[rank_ind].iloc[time_index]['x0'], ' & ', positions_raw[rank_ind].iloc[time_index]['x1'])
        cell_S0full[rank] = positions_raw[rank_ind].iloc[time_index]['S0full']
        cell_S1full[rank] = positions_raw[rank_ind].iloc[time_index]['S1full']

        xc = positions_raw[rank_ind].iloc[time_index]['x0']
        yc = positions_raw[rank_ind].iloc[time_index]['x1']
        for defind, defect in enumerate(plusdef_info):
            #some adjustments for periodic BC
            tempxc = positions_raw[rank_ind].iloc[time_index]['x0']
            tempyc = positions_raw[rank_ind].iloc[time_index]['x1']
            if(defect[0][0] < 50):
                if(tempxc > 75):
                    tempxc -= 100
            elif(defect[0][0] > 50):
                if(tempxc < 25):
                    tempxc += 100
            if(defect[0][1] < 50):
                if(tempyc > 75):
                    tempyc -= 100
            elif(defect[0][1] > 50):
                if(tempyc < 25):
                    tempyc += 100
            dist_to_defect[defind, rank] = np.linalg.norm(defect[0]-np.array([tempxc,tempyc]))
        
    fig3, ax3 = plt.subplots()
    #radcontour = ax.contourf(grid_X, grid_Y, radius)
    S_magcontour = ax3.contourf(grid_X, grid_Y, S_mag)#, levels = 40)
    fig3.colorbar(S_magcontour, ax=ax3)

    S_mag_neigh_avg = np.zeros(len(plusdef_info))
    cell_Smag = np.linalg.norm(np.array([cell_S0full, cell_S1full]), axis=0)
    #S_threshold = np.mean(cell_Smag) #for now it is mean, but need to be changed in the future
    S_threshold = 7.0 #8.0

    
    for defind, defect in enumerate(plusdef_info):
        #print(plusdef_info[defind][0])
        #print(np.where(dist_to_defect[defind]<12))
        #checking if mean of magnitude of S tensor of cells neighbouring the defects
        #is higher than a threshold
        mean_Smag_neighbours = np.mean(cell_Smag[np.where(dist_to_defect[defind]<12)]) #should 12 be changed?
        if(defect[2]<np.pi/2):#extensile
            my_colour = 'b'
            if(mean_Smag_neighbours < S_threshold):
                is_extensile_and_isotropic += 1
                my_marker = 'o'
            else:
                is_extensile_and_elongated += 1
                my_marker = '^'
        else:
            my_colour = 'r'
            if(mean_Smag_neighbours < S_threshold):
                is_contractile_and_isotropic += 1
                my_marker = 'o'
            else:
                is_contractile_and_elongated += 1
                my_marker = '^'
        ax3.scatter(defect[0][0], defect[0][1], c=my_colour, marker = my_marker)  

    plt.savefig(result_dir+'/defect_nature_thres7_'+str(time_index+1).zfill(4)+'.png', dpi=300, bbox_inches='tight')
    fig3.clear()
    plt.close(fig3)

    for S_ind, S_thres in enumerate(S_thres_array):
        fig5, ax5 = plt.subplots()
        S_magcontour = ax5.contourf(grid_X, grid_Y, S_mag)#, levels = 40)
        fig5.colorbar(S_magcontour, ax=ax5)
        for defind, defect in enumerate(plusdef_info):
            #print(plusdef_info[defind][0])
            #print(np.where(dist_to_defect[defind]<12))
            #checking if mean of magnitude of S tensor of cells neighbouring the defects
            #is higher than a threshold
            #if(S_ind == 0):
            mean_Smag_neighbours = np.mean(cell_Smag[np.where(dist_to_defect[defind]<12)]) #should 12 be changed?
            #    plusdef_dict[time][defind].append(mean_Smag_neighbours)
            #else:
            #    mean_Smag_neighbours = plusdef_dict[time][defind][-1]

            if(defect[2]<np.pi/2):#extensile
                my_colour = 'b'
                if(mean_Smag_neighbours < S_thres):
                    is_ext_and_iso[S_ind] += 1
                    my_marker = 'o'
                else:
                    is_ext_and_elon[S_ind] += 1
                    my_marker = '^'
            else:
                my_colour = 'r'
                if(mean_Smag_neighbours < S_thres):
                    is_con_and_iso[S_ind] += 1
                    my_marker = 'o'
                else:
                    is_con_and_elon[S_ind] += 1
                    my_marker = '^'
            ax5.scatter(defect[0][0], defect[0][1], c=my_colour, marker = my_marker)  
        plt.savefig(result_dir+'/def_nat' + postfix + '_S_thres' + str(S_thres) + '_' + 
                        str(time_index+1).zfill(4)+'.png', dpi=300, bbox_inches='tight')
        fig5.clear()
        plt.close(fig5)

    
np.save(result_dir + '/defangles' + postfix + '.npy', np.array(def_angles))
np.save(result_dir + '/diff_angles' + postfix + '.npy', np.array(diff_angles))



defect_nature = np.array([[is_extensile_and_elongated, is_extensile_and_isotropic],
                          [is_contractile_and_elongated, is_contractile_and_isotropic]])
np.save(result_dir + '/defect_nature' + postfix + '.npy', defect_nature)

print(defect_nature)

defect_nature_array = np.array([[is_ext_and_elon, is_ext_and_iso],
                          [is_con_and_elon, is_con_and_iso]])
np.save(result_dir + '/defect_nature_array' + postfix + '.npy', defect_nature_array)

print(defect_nature_array)

#def_nature_pd = pd.DataFrame(defect_nature, columns=['elongated', 'isotropic'], index = ['extensile', 'contractile'])
#def_nature_pd

#import seaborn as sn
#plt.figure(figsize = (4,3))
#sn.heatmap(def_nature_pd, annot=True)
#plt.savefig(result_dir + '/heatmap.png', dpi=300)


import pickle
with open(result_dir + '/defect_info' + postfix + '.pkl', 'wb') as f:
    pickle.dump(plusdef_dict, f) #saving defect informationm