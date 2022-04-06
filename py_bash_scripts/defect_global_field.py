
import os
import math
import numpy as np
import sys
import glob
import csv
import re
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors 
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from shapely.geometry import LineString
from matplotlib.collections import LineCollection

from joblib import Parallel, delayed
import multiprocessing as mp

#----------PARAMETERS----------

plot_scale_data = 0.7
plot_scale_interp = 0.9
plot_steps = 50
domain_dimension = np.array([100.0,100.0])
confined = False

if confined:
    mode = 'nearest'
else:
    mode = 'wrap'

plot_window = 10

def compute(folder):
    input_dir = folder 
    result_dir = folder + '/nematics'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    else: 
        print('nothing')
        #sys.exit('Input directory!')
        
    if os.path.exists(result_dir+'/plus_defects.csv'):
        os.remove(result_dir+'/plus_defects.csv')
    if os.path.exists(result_dir+'/minus_defects.csv'):
        os.remove(result_dir+'/minus_defects.csv')

    with open(result_dir+'/plus_defects.csv', mode='w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['Frame','x','y','z'])
    with open(result_dir+'/minus_defects.csv', mode='w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['Frame','x','y','z'])
        
    times = np.load(input_dir + '/timesteps.npy')
    print(times)
    file_pattern = input_dir + '/S0_glob_*.npy'
    
    # For plotting: 
    x = np.linspace(0,domain_dimension[0],plot_steps)
    y = np.linspace(0,domain_dimension[1],plot_steps)
    X,Y = np.meshgrid(x,y)
    xx = np.reshape(X,(plot_steps**2,1))
    yy = np.reshape(Y,(plot_steps**2,1))
    
    
    # Each filename is also one timestep
    count = 0
    for time in times[1:]:
        fig = plt.figure()
        fig.set_size_inches(6.0,6.0)
        ax = fig.add_subplot(1,1,1)
        ax.axis([0,100,0,100])
            
        ax.set_yticks([])
        ax.set_xticks([])
            
            #bounding box
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(0.5)
            ax.spines[axis].set_color('black')
        
        S0_glob = np.load(input_dir + '/S0_glob_' + '{:06.3f}'.format(time) +'.npy')
        S1_glob = np.load(input_dir + '/S1_glob_' + '{:06.3f}'.format(time) +'.npy')
        grid_X = np.load(input_dir + '/grid_x_' + '{:06.3f}'.format(time) +'.npy')
        grid_Y = np.load(input_dir + '/grid_y_' + '{:06.3f}'.format(time) +'.npy')
    
        ax.contourf(grid_X,grid_Y,(S0_glob),alpha=0.5)
        # Smooth
        S0_glob = gaussian_filter(S0_glob, 10.0, mode=mode)
    
        # Smooth
        S1_glob = gaussian_filter(S1_glob, 10.0, mode=mode)
    
        print('filtered time ' + str(time))
        steps_x,steps_y = np.shape(S0_glob)
        
        grid_x_vec = np.reshape(grid_X, (steps_x*steps_y,1))
        grid_y_vec = np.reshape(grid_Y, (steps_x*steps_y,1))
        
        S0_interp = griddata((np.squeeze(grid_x_vec),np.squeeze(grid_y_vec)),np.reshape(S0_glob,(steps_x*steps_y,1)),(xx,yy))
        S1_interp = griddata((np.squeeze(grid_x_vec),np.squeeze(grid_y_vec)),np.reshape(S1_glob,(steps_x*steps_y,1)),(xx,yy))
    
        # phase angle
        omega = np.multiply(np.sign(S1_interp),(np.arctan(S0_interp/np.abs(S1_interp))/ 2.0 + np.pi/4.0))
        omega = np.reshape(omega,(plot_steps**2,1))
    
        directions = np.zeros((plot_steps**2,2,2))
        directions[:,0,0] = xx.squeeze()-0.5*plot_scale_interp * np.sin(omega).squeeze()
        directions[:,1,0] = xx.squeeze()+0.5*plot_scale_interp * np.sin(omega).squeeze()
        directions[:,0,1] = yy.squeeze()-0.5*plot_scale_interp * np.cos(omega).squeeze()
        directions[:,1,1] = yy.squeeze()+0.5*plot_scale_interp * np.cos(omega).squeeze()
    
        line_segments = LineCollection(directions,linewidths=plot_scale_data,colors=cm.jet(0.0))
        ax.add_collection(line_segments)
        #print('saving')
        
        # Identify 0 contours of S0 and S1
        c1 = ax.contour(grid_X,grid_Y,S0_glob,0,colors='red',alpha=0.5)
        c2 = ax.contour(grid_X,grid_Y,S1_glob,0,colors='blue',alpha=0.5)
        
        #fig.savefig(result_dir+'/defect_locations_'+str(count).zfill(4)+'.png', dpi=300, bbox_inches='tight')
        print(np.asarray(c1.allsegs).shape)
        print(len(c1.allsegs[0]))
        if len(c1.allsegs[0]) <2 or len(c2.allsegs[0])<2:
            print('continuing')
            continue
        
        defect_list = []
        print('contours identified time ' + str(time))    
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
        
        grid = np.concatenate((grid_x_vec,grid_y_vec),axis=1)
        
        tmp = np.gradient(S0_glob)
        S0_dx = np.reshape(tmp[0],(steps_x*steps_y,1))
        S0_dy = np.reshape(tmp[1],(steps_x*steps_y,1))
    
        tmp = np.gradient(S1_glob)
        S1_dx = np.reshape(tmp[0],(steps_x*steps_y,1))
        S1_dy = np.reshape(tmp[1],(steps_x*steps_y,1))
        
        plus_defects = []
        minus_defects = []
        print(len(defect_list))
        if len(defect_list) > 0:
            for coord in defect_list:
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
    
        fig.savefig(result_dir+'/defect_locations_'+str(count).zfill(4)+'.png', dpi=300, bbox_inches='tight')        
        with open(result_dir+'/plus_defects.csv', mode='a') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for plus_def in plus_defects:
                csv_writer.writerow([count,plus_def[0],plus_def[1],0])
        
        with open(result_dir+'/minus_defects.csv', mode='a') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            for minus_def in minus_defects:
                csv_writer.writerow([count,minus_def[0],minus_def[1],0])        
        
        plt.close(fig)
        count += 1

noCores = mp.cpu_count()
Parallel(n_jobs=len(sys.argv)-1)(delayed(compute)(folder) for folder in sys.argv[1:])
