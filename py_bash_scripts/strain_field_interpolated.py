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
from matplotlib import cm
from matplotlib import rc

#import progress

from scipy import ndimage
#import seaborn as sns; sns.set_theme()
from shapely.geometry import LineString
from matplotlib.collections import LineCollection

def roundCoordinates(center, radius, nRings, initialAngle, N=5):
    radii = np.linspace(0,radius,nRings)
    x = np.array([])
    y = np.array([])

    for ind_radius,r in enumerate(radii):
        nPoints = ind_radius * N + 1
        theta = np.linspace(initialAngle,initialAngle + 2*np.pi,nPoints)
        x = np.concatenate((x,center[0] + r * np.cos(theta[:-1])))
        y = np.concatenate((y,center[1] + r * np.sin(theta[:-1])))

    return x,y

def scale_range(x, low, high):
    """
    takes a range of numbers, and scales the values such that 
    the minimum of the numbers is equal to low and maximum is 
    equal to high"""
    x += -(np.min(x))
    x /= np.max(x) / (high - low)
    x += low
    return x

fig,ax = plt.subplots()

# --- Parameters
start_time = 100
end_time = 1000
stride = 50

save_all_files = True
save_local_images = True
plot_all = False
minus = False

# Domain and Interpolation
domain_size = np.array([100.0,100.0])
interpolation_steps = 1000
global_smooth_windows = 10

# Only consider some defects
restrict_to_extensile = False
restrict_to_contractile = False

# Round averaging. Previous : 40steps, radius 2 and 20 rings
#local_gridwidth = [10.0,10.0]
local_steps = 60
local_smooth_window = 10

local_radius = 6
local_rings = 60

#plot_box = local_radius / np.sqrt(2.0)

plot_lines=True
plot_scale_interp = 0.75
plot_stride = 15

gauss_sigma = 20.0


# --- Parameters
rescale_velocities = False
normalise_velocities = False

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
if len(sys.argv) > 4:
    if int(sys.argv[4]) == 1:
        rescale_velocities = True
    else:
        rescale_velocities = False
if len(sys.argv) > 5:
    if int(sys.argv[5]) == 1:
        normalise_velocities = True
    else:
        normalise_velocities = False

numpy_file_dir = out_dir + '/numpy_files/'
individual_images_dir = out_dir + '/individual_images/'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
if not os.path.exists(numpy_file_dir):
    os.makedirs(numpy_file_dir)
if not os.path.exists(individual_images_dir):
    os.makedirs(individual_images_dir)  

print('Start Time: ' + str(start_time))
print('End Time: '+ str(end_time))

print('Rescale: ' + str(rescale_velocities))
print('Normalise: ' + str(normalise_velocities))

if restrict_to_extensile:
    df_restriction = pd.read_csv(out_dir +'/extensile_defects.csv')
else:
    if restrict_to_contractile:
        df_restriction = pd.read_csv(out_dir +'/contractile_defects.csv')

positions_raw = []
ranks = []

for filename in glob.glob(file_pattern):
    tmp = pd.read_csv(filename)
    if((len(sys.argv) > 6) and int(sys.argv[6]) == 1):
        positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', ' vel_angle', ' nematic_angle', ' total_interaction', ' neighbours', ' confine_interaction', ' growth_rate']])
    else:
        positions_raw.append(tmp[['time','rank','x0','x1','r','S0','S1','v0','v1', 'vel_angle', 'nematic_angle', 'total_interaction', 'neighbours', 'confine_interaction', 'growth_rate']])
    # we also need to extract the rank to build the bridge to vtk files
    ranks.append(int(re.findall(r'\d+', filename)[-1]))

for rank_ind,rank in enumerate(ranks):
    positions_raw[rank_ind]= positions_raw[rank_ind].fillna(0) #fill nan values with zero
    positions_raw[rank_ind]['x0'] = ndimage.uniform_filter(positions_raw[rank_ind]['x0'], size = 10, mode = 'nearest')
    positions_raw[rank_ind]['x1'] = ndimage.uniform_filter(positions_raw[rank_ind]['x1'], size = 10, mode = 'nearest')
    positions_raw[rank_ind]['v0'] = ndimage.uniform_filter(positions_raw[rank_ind]['v0'], size = 100, mode = 'nearest')
    positions_raw[rank_ind]['v1'] = ndimage.uniform_filter(positions_raw[rank_ind]['v1'], size = 100, mode = 'nearest')
    positions_raw[rank_ind]['S0'] = ndimage.uniform_filter(positions_raw[rank_ind]['S0'], size = 100, mode = 'nearest')
    positions_raw[rank_ind]['S1'] = ndimage.uniform_filter(positions_raw[rank_ind]['S1'], size = 100, mode = 'nearest')

# Interpolation
x = np.linspace(0,domain_size[0],interpolation_steps)
y = np.linspace(0,domain_size[1],interpolation_steps)
xx,yy = np.meshgrid(x,y)
xx_vec = np.reshape(xx,interpolation_steps**2)
yy_vec = np.reshape(yy,interpolation_steps**2)

average_angle_plus = 0
average_angle_minus = 0

count_plus = 0
count_minus = 0
angles_all_times = []
count_frames = 1
count_final = 0

#times = np.load(sys.argv[1] + '/global_fields/timesteps.npy')
#for time in times:

for time in range(start_time,end_time, stride):
    # If we want to restrict our considerations, we have to check if there's anything good in the current frame
    #progress.progress(time-start_time,end_time-start_time)
    if restrict_to_extensile:
        if df_restriction.query('@count_frames==Frame').empty:
            count_frames += 1
            continue
    phi_all = []
    
    # Empty lists for positions, velocities and elongations
    x0 = []
    x1 = []
    
    v0 = []
    v1 = []

    S0 = []
    S1 = []

    for rank_ind,rank in enumerate(ranks):
        x0.append(positions_raw[rank_ind].iloc[time]['x0'])
        x1.append(positions_raw[rank_ind].iloc[time]['x1'])
        
        v0.append(positions_raw[rank_ind].iloc[time]['v0'])
        v1.append(positions_raw[rank_ind].iloc[time]['v1'])        

        S0.append(positions_raw[rank_ind].iloc[time]['S0'])
        S1.append(positions_raw[rank_ind].iloc[time]['S1'])


    x0 = np.array(x0)
    x1 = np.array(x1)
    v0 = np.array(v0)
    v1 = np.array(v1)
    S0 = np.array(S0)
    S1 = np.array(S1)


    # Rewrite x0,x1 for periodic boundaries
    positions = np.stack((x0,x1),axis=1)
    positions = np.tile(positions,(5,1))
    positions[len(ranks):2*len(ranks),0] -= 100
    positions[2*len(ranks):3*len(ranks),0] += 100
    positions[3*len(ranks):4*len(ranks),1] -= 100
    positions[4*len(ranks):5*len(ranks),1] += 100
    
    print('positions' + str(positions.shape))
    # Copy v and S to match the positions
    v0 = np.tile(v0,5)
    v1 = np.tile(v1,5)
    
    S0 = np.tile(S0,5)
    S1 = np.tile(S1,5)

    print('v0' + str(v0.shape))
    # Construct global velocity and elongation fields by simple interpolation

    v0_glob = griddata(positions,v0,(xx_vec,yy_vec), method='cubic')
    v1_glob = griddata(positions,v1,(xx_vec,yy_vec), method='cubic')

    S0_glob = griddata(positions,S0,(xx_vec,yy_vec), method='cubic')
    S1_glob = griddata(positions,S1,(xx_vec,yy_vec), method='cubic')

    norm_S = np.sqrt(S0_glob**2 + S1_glob**2)
    S0_glob = S0_glob / norm_S
    S1_glob = S1_glob / norm_S

    S0_glob = np.reshape(S0_glob, (interpolation_steps,interpolation_steps))
    S1_glob = np.reshape(S1_glob, (interpolation_steps,interpolation_steps))
    
    # Smooth
    S0_glob = ndimage.uniform_filter(S0_glob, size=global_smooth_windows, mode='wrap')

    # Smooth
    S1_glob = ndimage.uniform_filter(S1_glob, size=global_smooth_windows, mode='wrap')

    # Identify 0 contours of S0 and S1
    c1 = ax.contour(xx,yy,S0_glob,0,colors='red',alpha=0.0)
    c2 = ax.contour(xx,yy,S1_glob,0,colors='blue',alpha=0.0)

    defect_list = []

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
        

    grid = np.stack((xx_vec,yy_vec),axis=1)
    defect_array = np.array(defect_list)
    
    defect_list_distributed = []
    for coord in defect_list:
        if np.sum(np.linalg.norm(coord-defect_array,axis=1)<2*local_radius)>1:
            pass
        else:
            defect_list_distributed.append(coord)

    S0_dy,S0_dx = np.gradient(S0_glob)
    S0_dx = np.reshape(S0_dx,(interpolation_steps**2,1))
    S0_dy = np.reshape(S0_dy,(interpolation_steps**2,1))    

    S1_dy,S1_dx = np.gradient(S1_glob)
    S1_dx = np.reshape(S1_dx,(interpolation_steps**2,1))
    S1_dy = np.reshape(S1_dy,(interpolation_steps**2,1))       

    plus_defects = []
    minus_defects= []

    total_num  = len(defect_list)
    for ind_defect,coord in enumerate(defect_list_distributed):
        # Now if we do make restrictions, we can check which defect is in the extensile df
        if restrict_to_extensile or restrict_to_contractile:
            if df_restriction.query('abs(@coord[0]-x)<1 and abs(@coord[1]-y)<1 and @count_frames==Frame').empty:
                continue
        
        # Local box around the defect
        x_loc = np.linspace(coord[0]-local_radius,coord[0]+local_radius,local_steps)
        y_loc = np.linspace(coord[1]-local_radius,coord[1]+local_radius,local_steps)
        xx_loc,yy_loc = np.meshgrid(x_loc,y_loc)

        xx_loc = np.reshape(xx_loc,(local_steps)**2)
        yy_loc = np.reshape(yy_loc,(local_steps)**2)
        grid = np.stack((xx_loc,yy_loc),axis=1)
        
        # Interpolate the contributions of the global field
        v0_loc = griddata(positions,v0,(xx_loc,yy_loc), method='cubic')
        v1_loc = griddata(positions,v1,(xx_loc,yy_loc), method='cubic')
        
        S0_loc = griddata(positions,S0,(xx_loc,yy_loc), method='cubic')
        S1_loc = griddata(positions,S1,(xx_loc,yy_loc), method='cubic')
        
        # Normalize S
        norm_S = np.sqrt(S0_loc**2 + S1_loc**2)
        S0_loc = S0_loc / norm_S
        S1_loc = S1_loc / norm_S

        S0_loc = np.reshape(S0_loc,(local_steps,local_steps))    
        S1_loc = np.reshape(S1_loc,(local_steps,local_steps))
        
        v0_loc = np.reshape(v0_loc,(local_steps,local_steps))    
        v1_loc = np.reshape(v1_loc,(local_steps,local_steps))    
        
        # Save non-smoothed versions
        v0_sharp = v0_loc
        v1_sharp = v1_loc
        
        # Smooth
        S0_loc = ndimage.uniform_filter(S0_loc, size=local_smooth_window, mode='nearest')
        S1_loc = ndimage.uniform_filter(S1_loc, size=local_smooth_window, mode='nearest')

        v0_loc = ndimage.uniform_filter(v0_loc, size=local_smooth_window, mode='nearest')
        v1_loc = ndimage.uniform_filter(v1_loc, size=local_smooth_window, mode='nearest')

        # For local stress fields and visualization
        if save_all_files:
            v0_dx,v0_dy = np.gradient(v0_loc)
            v1_dx,v1_dy = np.gradient(v1_loc)

            sigma_xx = v0_dx
            sigma_yy = v1_dy
            sigma_xy = 0.5 *( v1_dx + v0_dy)

            np.save(numpy_file_dir+'/v0_loc_'+str(count_plus) + '.npy', v0_loc)
            np.save(numpy_file_dir+'/v1_loc_'+str(count_plus) + '.npy', v1_loc)
            np.save(numpy_file_dir+'/S0_loc_'+str(count_plus) + '.npy', S0_loc)
            np.save(numpy_file_dir+'/S1_loc_'+str(count_plus) + '.npy', S1_loc)

            np.save(numpy_file_dir+'/sigma_xx_'+str(count_plus) + '.npy', sigma_xx)
            np.save(numpy_file_dir+'/sigma_yy_'+str(count_plus) + '.npy', sigma_yy)
            np.save(numpy_file_dir+'/sigma_xy_'+str(count_plus) + '.npy', sigma_xy)
        
        fig2 = plt.figure()
        ax1 = fig2.add_subplot(1,2,1)
        ax2 = fig2.add_subplot(1,2,2)
        c1 = ax1.contour(np.reshape(xx_loc,(local_steps,local_steps)),np.reshape(yy_loc,(local_steps,local_steps)),S0_loc,0,colors='red',alpha=0.0)
        c2 = ax2.contour(np.reshape(xx_loc,(local_steps,local_steps)),np.reshape(yy_loc,(local_steps,local_steps)),S1_loc,0,colors='blue',alpha=0.0)
        plt.close(fig2)
        local_defect_list = []
        
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
                            local_defect_list.append(np.array([k.x,k.y]))
                    except:
                        local_defect_list.append(np.array([intersects.x,intersects.y])) 
        
        if len(local_defect_list) == 0:
            continue    
        else:
            coord = local_defect_list[0]

        # Derivatives of S for delta - determine sign of the defect
        S0_dy,S0_dx = np.gradient(S0_loc)
        S0_dx = np.reshape(S0_dx,(local_steps**2,1))
        S0_dy = np.reshape(S0_dy,(local_steps**2,1))    

        S1_dy,S1_dx = np.gradient(S1_loc)
        S1_dx = np.reshape(S1_dx,(local_steps**2,1))
        S1_dy = np.reshape(S1_dy,(local_steps**2,1))       

        ind_local = np.unravel_index(np.linalg.norm(np.abs(grid-coord),axis=1).argmin(), np.linalg.norm(np.abs(grid-coord),axis=1).shape)
        delta = S0_dx[ind_local] * S1_dy[ind_local] - S0_dy[ind_local] * S1_dx[ind_local]    
          
        
        if minus:
            if delta > 0:
                continue
        else:
            if delta < 0 :
                continue
        
        # Compute orientation of the defect
        omega = np.multiply(np.sign(S1_loc),(np.arctan(S0_loc/np.abs(S1_loc))/ 2.0 + np.pi/4.0))
        n1 = np.sin(omega)
        n2 = np.cos(omega)
        
        n1 = np.reshape(n1,(local_steps,local_steps))
        n2 = np.reshape(n2,(local_steps,local_steps))

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
        
        # Interpolation on round coordinates
        ind_angle = np.argmin(np.linalg.norm(grid-coord,axis=1))
        angle = np.mean(angles[ind_angle])
        x_round,y_round = roundCoordinates(coord, local_radius, local_rings, -(1.0)*angle)

        S0_round = griddata((xx_loc,yy_loc),np.squeeze(np.reshape(S0_loc,(local_steps**2,1))),(x_round,y_round),method='nearest')
        S1_round = griddata((xx_loc,yy_loc),np.squeeze(np.reshape(S1_loc,(local_steps**2,1))),(x_round,y_round),method='nearest')

        v0_round = griddata((xx_loc,yy_loc),np.squeeze(np.reshape(v0_sharp,(local_steps**2,1))),(x_round,y_round),method='nearest')
        v1_round = griddata((xx_loc,yy_loc),np.squeeze(np.reshape(v1_sharp,(local_steps**2,1))),(x_round,y_round),method='nearest')        

        # Maybe we need to bring all velocities to the same scale [-1,1]
        if rescale_velocities:
            max_v0 = np.max(v0_round)
            min_v0 = np.min(v0_round)
            v0_round = (2/(max_v0-min_v0))*v0_round + 1 - (2*min_v0/max_v0)

            max_v1 = np.max(v1_round)
            min_v1 = np.min(v1_round)
            v1_round = (2/(max_v1-min_v1))*v1_round + 1 - (2*min_v1/max_v1)

        if normalise_velocities:
            v_mag = np.sqrt(np.square(v0_round) + np.square(v1_round))   
            v0_round /= v_mag
            v1_round /= v_mag
        
        if save_local_images:
            fig = plt.figure()
            fig.set_size_inches(24.0,6.0)
            ax1 = fig.add_subplot(1,4,1)
            ax2 = fig.add_subplot(1,4,2)
            ax3 = fig.add_subplot(1,4,3)
            ax4 = fig.add_subplot(1,4,4)
            
            S0_vec = np.reshape(S0_loc, local_steps**2)
            S1_vec = np.reshape(S1_loc, local_steps**2)

            omega = np.multiply(np.sign(S1_vec),(np.arctan(S0_vec/np.abs(S1_vec))/ 2.0 + np.pi/4.0))
            directions = np.zeros(((np.shape(omega)[0]),2,2))
            directions[:,0,0] = xx_loc-0.25*plot_scale_interp * np.sin(omega)
            directions[:,1,0] = xx_loc+0.25*plot_scale_interp * np.sin(omega)
            directions[:,0,1] = yy_loc-0.25*plot_scale_interp * np.cos(omega)
            directions[:,1,1] = yy_loc+0.25*plot_scale_interp * np.cos(omega)

            line_segments = LineCollection(directions,linewidths=1.5,colors=cm.jet(0.0))
            ax1.tricontourf(xx_loc,yy_loc,np.reshape(sigma_xx,local_steps**2),levels=30,cmap='jet')

            ax1.add_collection(line_segments)
            ax1.set_title('xx stress for +0.5')
            ax1.set_facecolor('xkcd:white')
            ax1.axis('off')

            ax1.plot([coord[0], coord[0]+np.cos(angle)],[coord[1], coord[1]+np.sin(angle)])
            ax2.plot([coord[0], coord[0]+np.cos(angle)],[coord[1], coord[1]+np.sin(angle)])
            ax3.plot([coord[0], coord[0]+np.cos(angle)],[coord[1], coord[1]+np.sin(angle)])

            line_segments = LineCollection(directions,linewidths=1.5,colors=cm.jet(0.0))
            cs = ax2.tricontourf(xx_loc,yy_loc,np.reshape(sigma_yy,local_steps**2),levels=30,cmap='jet')

            ax2.add_collection(line_segments)
            ax2.set_title('yy stress for +0.5')
            ax2.set_facecolor('xkcd:white')
            ax2.axis('off')

            line_segments = LineCollection(directions,linewidths=1.5,colors=cm.jet(0.0))
            cs = ax3.tricontourf(xx_loc,yy_loc,np.reshape(sigma_xy,local_steps**2),levels=30,cmap='jet')

            ax3.add_collection(line_segments)
            ax3.set_title('xy stress for +0.5')
            ax3.set_facecolor('xkcd:white')
            ax3.axis('off')

            ax4.quiver(xx_loc[::5],yy_loc[::5],v0_loc[::5], v1_loc[::5])
            ax4.set_title('Velocity Field')
            ax4.set_facecolor('xkcd:white')
            ax4.axis('off')

            cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.03])#
            fig.colorbar(cs, cax = cbar_ax, orientation='horizontal')
            fig.savefig(individual_images_dir + '/stress_and_flow_' + str(count_plus) + '.png',dpi=300)
            plt.close(fig)

        
        
        if count_plus == 0:
            S0_plus = S0_round
            S1_plus = S1_round
            v0_plus = v0_round
            v1_plus = v1_round 
            count_final += 1          
        # break
        else:
            S0_plus += S0_round
            S1_plus += S1_round
            v0_plus += v0_round
            v1_plus += v1_round
            count_final += 1          

        #average_angle_plus += angle
        angles_all_times.append(angle)
        count_plus += 1        
    count_frames += 1
            
            
x_round_average,y_round_average = roundCoordinates([0,0], local_radius, local_rings,np.pi/2.0)
x_round_plot,y_round_plot = roundCoordinates([0,0], local_radius, 10, np.pi/2.0)

# For Plotting directors
S0_round_plot = griddata((x_round_average,y_round_average),S0_plus,(x_round_plot,y_round_plot),method='linear')
S1_round_plot = griddata((x_round_average,y_round_average),S1_plus,(x_round_plot,y_round_plot),method='linear')

# Go back to square coordinates for derivative computation and smoothing
#x = np.linspace(-0.5*local_radius,0.5*local_radius,0.5*local_steps)
#xx,yy = np.meshgrid(x,x)
#xx_vec = np.reshape(xx,local_steps**2)
#yy_vec = np.reshape(yy,local_steps**2)
#half_local_steps = int(local_steps)

#x_rect = np.linspace(np.min(x_round_average),np.max(x_round_average),half_local_steps)
#y_rect = np.linspace(np.min(y_round_average),np.max(y_round_average),half_local_steps)
x_rect = np.linspace(-int(local_radius/np.sqrt(2)),int(local_radius/np.sqrt(2)),local_steps)
y_rect = np.linspace(-int(local_radius/np.sqrt(2)),int(local_radius/np.sqrt(2)),local_steps)
X,Y = np.meshgrid(x_rect,y_rect)
xx_vec = np.reshape(X,local_steps**2)
yy_vec = np.reshape(Y,local_steps**2)

v0_rect = griddata((x_round_average,y_round_average),v0_plus,(X,Y),method = 'nearest')
v1_rect = griddata((x_round_average,y_round_average),v1_plus,(X,Y),method = 'nearest')
#angle_rect = griddata((x_round_average,y_round_average),velocity_angle,(X,Y),method = 'nearest')
v0_rect_sigma05 = ndimage.gaussian_filter(v0_rect, sigma=0.5, mode='nearest')
v1_rect_sigma05 = ndimage.gaussian_filter(v1_rect, sigma=0.5, mode='nearest')

v0_rect_sigma1 = ndimage.gaussian_filter(v0_rect, sigma=1.0, mode='nearest')
v1_rect_sigma1 = ndimage.gaussian_filter(v1_rect, sigma=1.0, mode='nearest')

v0_rect_sigma15 = ndimage.gaussian_filter(v0_rect, sigma=1.5, mode='nearest')
v1_rect_sigma15 = ndimage.gaussian_filter(v1_rect, sigma=1.5, mode='nearest')

#v0_rect_sigma2 = ndimage.gaussian_filter(v0_rect, sigma=2.0, mode='nearest')
#v1_rect_sigma2 = ndimage.gaussian_filter(v1_rect, sigma=2.0, mode='nearest')

v0_rect_sigma25 = ndimage.gaussian_filter(v0_rect, sigma=2.5, mode='nearest')
v1_rect_sigma25 = ndimage.gaussian_filter(v1_rect, sigma=2.5, mode='nearest')

v0_round = griddata((np.reshape(X,local_steps**2),np.reshape(Y,local_steps**2)),np.reshape(v0_rect_sigma05,local_steps**2),(x_round_average,y_round_average),method = 'nearest')
v1_round = griddata((np.reshape(X,local_steps**2),np.reshape(Y,local_steps**2)),np.reshape(v1_rect_sigma05,local_steps**2),(x_round_average,y_round_average),method = 'nearest')



v0_dy,v0_dx = np.gradient(v0_rect)
v1_dy,v1_dx = np.gradient(v1_rect)

# Restriction
v0_dy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
v0_dx[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
v1_dy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
v1_dx[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0

sigma_xx = v0_dx
sigma_yy = v1_dy
sigma_xy = 0.5 *( v0_dy + v1_dx)

#Restriction 
sigma_xx[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
sigma_yy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
sigma_xy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0

np.save(numpy_file_dir+'/average_v0.npy', v0_rect)
np.save(numpy_file_dir+'/average_v1.npy', v1_rect)
np.save(numpy_file_dir+'/average_S0.npy', S0_plus)
np.save(numpy_file_dir+'/average_S1.npy', S1_plus)

np.save(numpy_file_dir+'/average_sigma_xx.npy', sigma_xx)
np.save(numpy_file_dir+'/average_sigma_yy.npy', sigma_yy)
np.save(numpy_file_dir+'/average_sigma_xy.npy', sigma_xy)
    
    
if plot_all:    
    #### VISUALIZATION
    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(24.0,12.0)
    ax1 = fig.add_subplot(2,4,1)
    ax2 = fig.add_subplot(2,4,2)
    ax3 = fig.add_subplot(2,4,3)
    ax4 = fig.add_subplot(2,4,4)
    ax5 = fig.add_subplot(2,4,5)
    ax6 = fig.add_subplot(2,4,6)
    ax7 = fig.add_subplot(2,4,7)
    ax8 = fig.add_subplot(2,4,8)
    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax1.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx,local_steps**2),levels=30,cmap='jet')
    
    ax1.add_collection(line_segments)
    ax1.set_title('v0_dx for +0.5')
    ax1.set_facecolor('xkcd:white')
    ax1.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax2.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy,local_steps**2),levels=30,cmap='jet')
    
    ax2.add_collection(line_segments)
    ax2.set_title('v0_dy for +0.5')
    ax2.set_facecolor('xkcd:white')
    ax2.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax3.tricontourf(xx_vec,yy_vec,np.reshape(v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax3.add_collection(line_segments)
    ax3.set_title('v1_dx for +0.5')
    ax3.set_facecolor('xkcd:white')
    ax3.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax4.tricontourf(xx_vec,yy_vec,np.reshape(v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax4.add_collection(line_segments)
    ax4.set_title('v1_dy for +0.5')
    ax4.set_facecolor('xkcd:white')
    ax4.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax5.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax5.add_collection(line_segments)
    ax5.set_title('v0_dx+v1_dy for +0.5')
    ax5.set_facecolor('xkcd:white')
    ax5.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax6.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax6.add_collection(line_segments)
    ax6.set_title('v0_dx+v1_dx for +0.5')
    ax6.set_facecolor('xkcd:white')
    ax6.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax7.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax7.add_collection(line_segments)
    ax7.set_title('v0_dy+v1_dy for +0.5')
    ax7.set_facecolor('xkcd:white')
    ax7.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax8.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax8.add_collection(line_segments)
    ax8.set_title('v0_dy+v1_dx for +0.5')
    ax8.set_facecolor('xkcd:white')
    ax8.axis('off')
    #ax5.quiver(x_round_average[::5],y_round_average[::5],v0_round[::5], v1_round[::5])
    #ax5.quiver(xx_vec,yy_vec,v0_rect,v1_rect)
    #ax5.set_title('Average Velocity Field')
    #ax5.set_facecolor('xkcd:white')
    #ax5.axis('off')
    
    cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.03])#
    fig.colorbar(cs, cax = cbar_ax, orientation='horizontal')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile.png',dpi=300)
        else:
            fig.savefig(out_dir + '/average_strain_nonrestricted.png',dpi=300)
    plt.close(fig)
    
    v0_dy,v0_dx = np.gradient(v0_rect_sigma05)
    v1_dy,v1_dx = np.gradient(v1_rect_sigma05)
    
    sigma_xx = v0_dx
    sigma_yy = v1_dy
    sigma_xy = 0.5 *( v0_dy + v1_dx)
    
    # Restriction
    sigma_xx[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
    sigma_yy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
    sigma_xy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0

    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(24.0,12.0)
    ax1 = fig.add_subplot(2,4,1)
    ax2 = fig.add_subplot(2,4,2)
    ax3 = fig.add_subplot(2,4,3)
    ax4 = fig.add_subplot(2,4,4)
    ax5 = fig.add_subplot(2,4,5)
    ax6 = fig.add_subplot(2,4,6)
    ax7 = fig.add_subplot(2,4,7)
    ax8 = fig.add_subplot(2,4,8)
    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax1.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx,local_steps**2),levels=30,cmap='jet')
    
    ax1.add_collection(line_segments)
    ax1.set_title('v0_dx for +0.5')
    ax1.set_facecolor('xkcd:white')
    ax1.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax2.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy,local_steps**2),levels=30,cmap='jet')
    
    ax2.add_collection(line_segments)
    ax2.set_title('v0_dy for +0.5')
    ax2.set_facecolor('xkcd:white')
    ax2.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax3.tricontourf(xx_vec,yy_vec,np.reshape(v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax3.add_collection(line_segments)
    ax3.set_title('v1_dx for +0.5')
    ax3.set_facecolor('xkcd:white')
    ax3.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax4.tricontourf(xx_vec,yy_vec,np.reshape(v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax4.add_collection(line_segments)
    ax4.set_title('v1_dy for +0.5')
    ax4.set_facecolor('xkcd:white')
    ax4.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax5.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax5.add_collection(line_segments)
    ax5.set_title('v0_dx+v1_dy for +0.5')
    ax5.set_facecolor('xkcd:white')
    ax5.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax6.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax6.add_collection(line_segments)
    ax6.set_title('v0_dx+v1_dx for +0.5')
    ax6.set_facecolor('xkcd:white')
    ax6.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax7.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax7.add_collection(line_segments)
    ax7.set_title('v0_dy+v1_dy for +0.5')
    ax7.set_facecolor('xkcd:white')
    ax7.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax8.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax8.add_collection(line_segments)
    ax8.set_title('v0_dy+v1_dx for +0.5')
    ax8.set_facecolor('xkcd:white')
    ax8.axis('off')
    
    cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.03])#
    fig.colorbar(cs, cax = cbar_ax, orientation='horizontal')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted_sigma05.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile_sigma05.png',dpi=300)
        else:
            fig.savefig(out_dir + '/average_strain_nonrestricted_sigma05.png',dpi=300)
    plt.close(fig)
    
    v0_dy,v0_dx = np.gradient(v0_rect_sigma1)
    v1_dy,v1_dx = np.gradient(v1_rect_sigma1)
    
    sigma_xx = v0_dx
    sigma_yy = v1_dy
    sigma_xy = 0.5 *( v0_dy + v1_dx)
    
    #Restriction
    sigma_xx[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
    sigma_yy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
    sigma_xy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0

    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(24.0,12.0)
    ax1 = fig.add_subplot(2,4,1)
    ax2 = fig.add_subplot(2,4,2)
    ax3 = fig.add_subplot(2,4,3)
    ax4 = fig.add_subplot(2,4,4)
    ax5 = fig.add_subplot(2,4,5)
    ax6 = fig.add_subplot(2,4,6)
    ax7 = fig.add_subplot(2,4,7)
    ax8 = fig.add_subplot(2,4,8)
    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax1.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx,local_steps**2),levels=30,cmap='jet')
    
    ax1.add_collection(line_segments)
    ax1.set_title('v0_dx for +0.5')
    ax1.set_facecolor('xkcd:white')
    ax1.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax2.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy,local_steps**2),levels=30,cmap='jet')
    
    ax2.add_collection(line_segments)
    ax2.set_title('v0_dy for +0.5')
    ax2.set_facecolor('xkcd:white')
    ax2.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax3.tricontourf(xx_vec,yy_vec,np.reshape(v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax3.add_collection(line_segments)
    ax3.set_title('v1_dx for +0.5')
    ax3.set_facecolor('xkcd:white')
    ax3.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax4.tricontourf(xx_vec,yy_vec,np.reshape(v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax4.add_collection(line_segments)
    ax4.set_title('v1_dy for +0.5')
    ax4.set_facecolor('xkcd:white')
    ax4.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax5.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax5.add_collection(line_segments)
    ax5.set_title('v0_dx+v1_dy for +0.5')
    ax5.set_facecolor('xkcd:white')
    ax5.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax6.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax6.add_collection(line_segments)
    ax6.set_title('v0_dx+v1_dx for +0.5')
    ax6.set_facecolor('xkcd:white')
    ax6.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax7.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax7.add_collection(line_segments)
    ax7.set_title('v0_dy+v1_dy for +0.5')
    ax7.set_facecolor('xkcd:white')
    ax7.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax8.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax8.add_collection(line_segments)
    ax8.set_title('v0_dy+v1_dx for +0.5')
    ax8.set_facecolor('xkcd:white')
    ax8.axis('off')
    
    cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.03])#
    fig.colorbar(cs, cax = cbar_ax, orientation='horizontal')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted_sigma1.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile_sigma1.png',dpi=300)
        else:
            fig.savefig(out_dir + '/average_strain_nonrestricted_sigma1.png',dpi=300)
    plt.close(fig)
    
    
    
    v0_dy,v0_dx = np.gradient(v0_rect_sigma15)
    v1_dy,v1_dx = np.gradient(v1_rect_sigma15)
    
    sigma_xx = v0_dx
    sigma_yy = v1_dy
    sigma_xy = 0.5 *( v0_dy + v1_dx)

    #Restriction
    sigma_xx[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
    sigma_yy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0
    sigma_xy[np.where(np.sqrt(X**2+Y**2)>local_radius)[0],np.where(np.sqrt(X**2+Y**2)>local_radius)[1]] = 0.0   
 
    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(24.0,12.0)
    ax1 = fig.add_subplot(2,4,1)
    ax2 = fig.add_subplot(2,4,2)
    ax3 = fig.add_subplot(2,4,3)
    ax4 = fig.add_subplot(2,4,4)
    ax5 = fig.add_subplot(2,4,5)
    ax6 = fig.add_subplot(2,4,6)
    ax7 = fig.add_subplot(2,4,7)
    ax8 = fig.add_subplot(2,4,8)
    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax1.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx,local_steps**2),levels=30,cmap='jet')
    
    ax1.add_collection(line_segments)
    ax1.set_title('v0_dx for +0.5')
    ax1.set_facecolor('xkcd:white')
    ax1.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax2.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy,local_steps**2),levels=30,cmap='jet')
    
    ax2.add_collection(line_segments)
    ax2.set_title('v0_dy for +0.5')
    ax2.set_facecolor('xkcd:white')
    ax2.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax3.tricontourf(xx_vec,yy_vec,np.reshape(v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax3.add_collection(line_segments)
    ax3.set_title('v1_dx for +0.5')
    ax3.set_facecolor('xkcd:white')
    ax3.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax4.tricontourf(xx_vec,yy_vec,np.reshape(v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax4.add_collection(line_segments)
    ax4.set_title('v1_dy for +0.5')
    ax4.set_facecolor('xkcd:white')
    ax4.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax5.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax5.add_collection(line_segments)
    ax5.set_title('v0_dx+v1_dy for +0.5')
    ax5.set_facecolor('xkcd:white')
    ax5.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax6.tricontourf(xx_vec,yy_vec,np.reshape(v0_dx+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax6.add_collection(line_segments)
    ax6.set_title('v0_dx+v1_dx for +0.5')
    ax6.set_facecolor('xkcd:white')
    ax6.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax7.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dy,local_steps**2),levels=30,cmap='jet')
    
    ax7.add_collection(line_segments)
    ax7.set_title('v0_dy+v1_dy for +0.5')
    ax7.set_facecolor('xkcd:white')
    ax7.axis('off')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    cs = ax8.tricontourf(xx_vec,yy_vec,np.reshape(v0_dy+v1_dx,local_steps**2),levels=30,cmap='jet')
    
    ax8.add_collection(line_segments)
    ax8.set_title('v0_dy+v1_dx for +0.5')
    ax8.set_facecolor('xkcd:white')
    ax8.axis('off')
    
    cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.03])#
    fig.colorbar(cs, cax = cbar_ax, orientation='horizontal')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted_sigma15.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile_sigma15.png',dpi=300)
        else:
            fig.savefig(out_dir + '/average_strain_nonrestricted_sigma15.png',dpi=300)
    plt.close(fig)
    
else:
    plt.rc('axes', titlesize=20)
    
    #### VISUALIZATION
    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(6.0,12.0)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    #ax3 = fig.add_subplot(1,3,3)
    #ax1 = fig.add_axes([0.05, 0.1, 0.4, 0.8])
    #ax2 = fig.add_axes([0.50, 0.1, 0.4, 0.8])
    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    tmp = v0_dy+v1_dx
    cs1 = ax1.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')
    
    ax1.add_collection(line_segments)
    #ax1.set_title(r'Average $E_{xy}$')
    #ax1.set_facecolor('xkcd:white')
    ax1.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax1.set_ylim([np.min(yy_vec),np.max(yy_vec)])
    #ax1.axis('off')
    #ax1.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax1.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    ax1.set_xticks([])
    ax1.set_yticks([])

    
    #ax1.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax1.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax1.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax1.spines['left'].set_color('black')
    ax1.spines['right'].set_color('black')
    ax1.spines['bottom'].set_color('black')
    ax1.spines['top'].set_color('black')
    #ax1.set_xlabel('x')
    #ax1.set_ylabel('y')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    
    
    tmp = v1_dy
    cs2 = ax2.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')
    
    ax2.add_collection(line_segments)
    #ax2.set_title(r'Average $E_{yy}$')
    #ax2.set_facecolor('xkcd:white')
    ax2.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax2.set_ylim([np.min(yy_vec),np.max(yy_vec)])
    ax2.set_xticks([])
    ax2.set_yticks([])
    #ax2.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax2.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax2.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax2.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax2.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax2.spines['left'].set_color('black')
    ax2.spines['right'].set_color('black')
    ax2.spines['bottom'].set_color('black')
    ax2.spines['top'].set_color('black')
    #ax2.set_xlabel('x')
    #ax2.set_ylabel('y')
    
    #ax3.quiver(x_round_average[::5],y_round_average[::5],v0_round[::5], v1_round[::5])
    #ax3.quiver(xx_vec,yy_vec,v1_rect,v0_rect)
    #ax3.set_title('Average Velocity Field')
    #ax3.set_facecolor('xkcd:white')
    #ax3.axis('off')
    
    #cbar_ax1 = fig.add_axes([0.45, 0.05, 0.05, 0.07])#
    #cbar_ax2 = fig.add_axes([0.91, 0.1, 0.03, 0.8])
   # fig.colorbar(cs1, cax = cbar_ax1, orientation='vertical')
    #fig.colorbar(cs2, cax = cbar_ax2, orientation='vertical')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile.png',dpi=300)
        else:
            if normalise_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_normalised.png',dpi=300)
            elif rescale_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_rescaled.png',dpi=300)
            else:
                fig.savefig(out_dir + '/average_strain_nonrestricted.png',dpi=300)
    plt.close(fig)
    
    v0_dy,v0_dx = np.gradient(v0_rect_sigma05)
    v1_dy,v1_dx = np.gradient(v1_rect_sigma05)
    
    sigma_xx = v0_dx
    sigma_yy = v1_dy
    sigma_xy = 0.5 *( v0_dy + v1_dx)
    
    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(6.0,12.0)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    #ax3 = fig.add_subplot(1,3,3)
    #ax1 = fig.add_axes([0.05, 0.1, 0.4, 0.8])
    #ax2 = fig.add_axes([0.50, 0.1, 0.4, 0.8])
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))

    tmp = v0_dy+v1_dx
    cs1 = ax1.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax1.add_collection(line_segments)
    #ax1.set_title(r'Average $E_{xy}$')
    #ax1.set_facecolor('xkcd:white')
    ax1.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax1.set_ylim([np.min(yy_vec),np.max(yy_vec)])
    ax1.set_xticks([])
    ax1.set_yticks([])
    #ax1.axis('off')
    #ax1.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax1.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax1.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax1.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax1.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax1.spines['left'].set_color('black')
    ax1.spines['right'].set_color('black')
    ax1.spines['bottom'].set_color('black')
    ax1.spines['top'].set_color('black')
    #ax1.set_xlabel('x')
    #ax1.set_ylabel('y')
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))

    tmp = v1_dy

    cs2 = ax2.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax2.add_collection(line_segments)
    #ax2.set_title(r'Average $E_{yy}$')
    #ax2.set_facecolor('xkcd:white')
    ax2.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax2.set_ylim([np.min(yy_vec),np.max(yy_vec)])

    ax2.set_xticks([])
    ax2.set_yticks([])    

    #ax2.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax2.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax2.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax2.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax2.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax2.spines['left'].set_color('black')
    ax2.spines['right'].set_color('black')
    ax2.spines['bottom'].set_color('black')
    ax2.spines['top'].set_color('black')
    #ax2.set_xlabel('x')
    #ax2.set_ylabel('y')
    
    #ax3.quiver(x_round_average[::5],y_round_average[::5],v0_round[::5], v1_round[::5])
    #ax3.quiver(xx_vec,yy_vec,v1_rect_sigma05,v0_rect_sigma05)
    #ax3.set_title('Average Velocity Field')
    #ax3.set_facecolor('xkcd:white')
    #ax3.axis('off')
    
    #cbar_ax2 = fig.add_axes([0.91, 0.1, 0.03, 0.8])
   # fig.colorbar(cs1, cax = cbar_ax1, orientation='vertical')
    #fig.colorbar(cs2, cax = cbar_ax2, orientation='vertical')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted_sigma05.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile_sigma05.png',dpi=300)
        else:
            if normalise_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_normalised_sigma05.png',dpi=300)
            elif rescale_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_rescaled_sigma05.png',dpi=300)
            else:
                fig.savefig(out_dir + '/average_strain_nonrestricted_sigma05.png',dpi=300)
    plt.close(fig)
    
    v0_dy,v0_dx = np.gradient(v0_rect_sigma1)
    v1_dy,v1_dx = np.gradient(v1_rect_sigma1)
    
    sigma_xx = v0_dx
    sigma_yy = v1_dy
    sigma_xy = 0.5 *( v0_dy + v1_dx)
    
    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(6.0,12.0)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
        #ax1 = fig.add_subplot(1,2,1)
    #ax2 = fig.add_subplot(1,2,2)
    #ax3 = fig.add_subplot(1,3,3)
    #ax1 = fig.add_axes([0.05, 0.1, 0.4, 0.8])
    #ax2 = fig.add_axes([0.50, 0.1, 0.4, 0.8])
    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    tmp = v0_dy+v1_dx
    cs1 = ax1.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax1.add_collection(line_segments)
    #ax1.set_title(r'Average $E_{xy}$')
    #ax1.set_facecolor('xkcd:white')
    ax1.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax1.set_ylim([np.min(yy_vec),np.max(yy_vec)])
    #ax1.axis('off')
    #ax1.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax1.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax1.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax1.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax1.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax1.spines['left'].set_color('black')
    ax1.spines['right'].set_color('black')
    ax1.spines['bottom'].set_color('black')
    ax1.spines['top'].set_color('black')
    #ax1.set_xlabel('x')
    #ax1.set_ylabel('y')
    
    ax1.set_xticks([])
    ax1.set_yticks([])  
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))

    tmp = v1_dy

    cs2 = ax2.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax2.add_collection(line_segments)
    #ax2.set_title(r'Average $E_{yy}$')
    #ax2.set_facecolor('xkcd:white')
    ax2.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax2.set_ylim([np.min(yy_vec),np.max(yy_vec)])

    #ax2.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax2.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax2.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax2.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax2.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax2.spines['left'].set_color('black')
    ax2.spines['right'].set_color('black')
    ax2.spines['bottom'].set_color('black')
    ax2.spines['top'].set_color('black')
    #ax2.set_xlabel('x')
    #ax2.set_ylabel('y')
    
    ax2.set_xticks([])
    ax2.set_yticks([])  
    #ax3.quiver(x_round_average[::5],y_round_average[::5],v0_round[::5], v1_round[::5])
    #ax3.quiver(xx_vec,yy_vec,v1_rect_sigma05,v0_rect_sigma05)
    #ax3.set_title('Average Velocity Field')
    #ax3.set_facecolor('xkcd:white')
    #ax3.axis('off')
        
    #cbar_ax2 = fig.add_axes([0.91, 0.1, 0.03, 0.8])
   # fig.colorbar(cs1, cax = cbar_ax1, orientation='vertical')
    #fig.colorbar(cs2, cax = cbar_ax2, orientation='vertical')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted_sigma1.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile_sigma1.png',dpi=300)
        else:
            if normalise_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_normalised_sigma1.png',dpi=300)
            elif rescale_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_rescaled_sigma1.png',dpi=300)
            else:
                fig.savefig(out_dir + '/average_strain_nonrestricted_sigma1.png',dpi=300)
    plt.close(fig)
    
    
    
    v0_dy,v0_dx = np.gradient(v0_rect_sigma15)
    v1_dy,v1_dx = np.gradient(v1_rect_sigma15)
    
    sigma_xx = v0_dx
    sigma_yy = v1_dy
    sigma_xy = 0.5 *( v0_dy + v1_dx)
    
    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(6.0,12.0)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)    #ax1 = fig.add_subplot(1,2,1)
    #ax2 = fig.add_subplot(1,2,2)
    #ax3 = fig.add_subplot(1,3,3)
    #ax1 = fig.add_axes([0.05, 0.1, 0.4, 0.8])
    #ax2 = fig.add_axes([0.49, 0.1, 0.4, 0.8])

    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    tmp = v0_dy+v1_dx
    cs1 = ax1.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax1.add_collection(line_segments)
    #ax1.set_title(r'Average $E_{xy}$')
    #ax1.set_facecolor('xkcd:white')
    ax1.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax1.set_ylim([np.min(yy_vec),np.max(yy_vec)])
    #ax1.axis('off')
    #ax1.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax1.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax1.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax1.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax1.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax1.spines['left'].set_color('black')
    ax1.spines['right'].set_color('black')
    ax1.spines['bottom'].set_color('black')
    ax1.spines['top'].set_color('black')
    #ax1.set_xlabel('x')
    #ax1.set_ylabel('y')
    ax1.set_xticks([])
    ax1.set_yticks([])  
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))

    tmp = v1_dy

    cs2 = ax2.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax2.add_collection(line_segments)
    #ax2.set_title(r'Average $E_{yy}$')
    #ax2.set_facecolor('xkcd:white')
    ax2.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax2.set_ylim([np.min(yy_vec),np.max(yy_vec)])

    #ax2.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax2.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax2.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax2.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax2.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax2.spines['left'].set_color('black')
    ax2.spines['right'].set_color('black')
    ax2.spines['bottom'].set_color('black')
    ax2.spines['top'].set_color('black')
    #ax2.set_xlabel('x')
    #ax2.set_ylabel('y')
    ax2.set_xticks([])
    ax2.set_yticks([])  
    #ax3.quiver(x_round_average[::5],y_round_average[::5],v0_round[::5], v1_round[::5])
    #ax3.quiver(xx_vec,yy_vec,v1_rect_sigma05,v0_rect_sigma05)
    #ax3.set_title('Average Velocity Field')
    #ax3.set_facecolor('xkcd:white')
    #ax3.axis('off')
    
    #cbar_ax2 = fig.add_axes([0.91, 0.1, 0.03, 0.8])
   # fig.colorbar(cs1, cax = cbar_ax1, orientation='vertical')
    #fig.colorbar(cs2, cax = cbar_ax2, orientation='vertical')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted_sigma15.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile_sigma15.png',dpi=300)
        else:
            if normalise_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_normalised_sigma15.png',dpi=300)
            elif rescale_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_rescaled_sigma15.png',dpi=300)
            else:
                fig.savefig(out_dir + '/average_strain_nonrestricted_sigma15.png',dpi=300)
    plt.close(fig)

    v0_dy,v0_dx = np.gradient(v0_rect_sigma25)
    v1_dy,v1_dx = np.gradient(v1_rect_sigma25)
    
    sigma_xx = v0_dx
    sigma_yy = v1_dy
    sigma_xy = 0.5 *( v0_dy + v1_dx)
    
    ## Start with +0.5
    fig = plt.figure()
    fig.set_size_inches(6.0,12.0)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)    #ax1 = fig.add_subplot(1,2,1)
    #ax2 = fig.add_subplot(1,2,2)
    #ax3 = fig.add_subplot(1,3,3)
    #ax1 = fig.add_axes([0.05, 0.1, 0.4, 0.8])
    #ax2 = fig.add_axes([0.49, 0.1, 0.4, 0.8])

    
    
    omega = np.multiply(np.sign(S1_round_plot),(np.arctan(S0_round_plot/np.abs(S1_round_plot))/ 2.0 + np.pi/4.0))
    directions = np.zeros(((np.shape(omega)[0]),2,2))
    directions[:,0,0] = x_round_plot-0.25*plot_scale_interp * np.sin(omega)
    directions[:,1,0] = x_round_plot+0.25*plot_scale_interp * np.sin(omega)
    directions[:,0,1] = y_round_plot-0.25*plot_scale_interp * np.cos(omega)
    directions[:,1,1] = y_round_plot+0.25*plot_scale_interp * np.cos(omega)
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))
    tmp = v0_dy+v1_dx
    cs1 = ax1.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax1.add_collection(line_segments)
    #ax1.set_title(r'Average $E_{xy}$')
    #ax1.set_facecolor('xkcd:white')
    ax1.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax1.set_ylim([np.min(yy_vec),np.max(yy_vec)])
    #ax1.axis('off')
    #ax1.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax1.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax1.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax1.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax1.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax1.spines['left'].set_color('black')
    ax1.spines['right'].set_color('black')
    ax1.spines['bottom'].set_color('black')
    ax1.spines['top'].set_color('black')
    #ax1.set_xlabel('x')
    #ax1.set_ylabel('y')
    ax1.set_xticks([])
    ax1.set_yticks([])  
    
    line_segments = LineCollection(directions,linewidths=3,colors=cm.jet(0.0))

    tmp = v1_dy

    cs2 = ax2.tricontourf(xx_vec,yy_vec,np.reshape(tmp,local_steps**2),levels=30,cmap='jet')    
    ax2.add_collection(line_segments)
    #ax2.set_title(r'Average $E_{yy}$')
    #ax2.set_facecolor('xkcd:white')
    ax2.set_xlim([np.min(xx_vec),np.max(xx_vec)])
    ax2.set_ylim([np.min(yy_vec),np.max(yy_vec)])

    #ax2.set_xticks([np.min(xx_vec),0.0,np.max(xx_vec)])
    #ax2.set_xticklabels([str(int(np.min(xx_vec))), str(0), str(int(np.min(xx_vec)))])
    #ax2.set_yticks([np.min(yy_vec),0.0,np.max(yy_vec)])
    #ax2.set_yticklabels([str(int(np.min(yy_vec))), str(0), str(int(np.min(yy_vec)))])    
    #ax2.tick_params(direction='out', top=False, right=False, bottom=True, left=True) 
    ax2.spines['left'].set_color('black')
    ax2.spines['right'].set_color('black')
    ax2.spines['bottom'].set_color('black')
    ax2.spines['top'].set_color('black')
    #ax2.set_xlabel('x')
    #ax2.set_ylabel('y')
    ax2.set_xticks([])
    ax2.set_yticks([])  
    #ax3.quiver(x_round_average[::5],y_round_average[::5],v0_round[::5], v1_round[::5])
    #ax3.quiver(xx_vec,yy_vec,v1_rect_sigma05,v0_rect_sigma05)
    #ax3.set_title('Average Velocity Field')
    #ax3.set_facecolor('xkcd:white')
    #ax3.axis('off')
    
    #cbar_ax2 = fig.add_axes([0.91, 0.1, 0.03, 0.8])
   # fig.colorbar(cs1, cax = cbar_ax1, orientation='vertical')
    #fig.colorbar(cs2, cax = cbar_ax2, orientation='vertical')
    if restrict_to_extensile:
        fig.savefig(out_dir + '/average_strain_restricted_sigma25.png',dpi=300)
    else:
        if restrict_to_contractile:
            fig.savefig(out_dir + '/average_strain_restricted_contractile_sigma25.png',dpi=300)
        else:
            if normalise_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_normalised_sigma25.png',dpi=300)
            elif rescale_velocities:
                fig.savefig(out_dir + '/average_strain_nonrestricted_rescaled_sigma25.png',dpi=300)
            else:
                fig.savefig(out_dir + '/average_strain_nonrestricted_sigma25.png',dpi=300)
    plt.close(fig)
print('Total number: ' + str(count_final))
