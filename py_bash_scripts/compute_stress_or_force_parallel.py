import matplotlib.pyplot as plt
import sys
import numpy as np
import os
import re
import glob
import pandas as pd
import multiprocessing as mp
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from scipy.interpolate import griddata
from scipy import ndimage
from scipy.ndimage import gaussian_filter
import load_files as load_files
import time
import queue

import dask
from dask.distributed import Client
from dask_jobqueue import SLURMCluster
#print("Number of processors: ", mp.cpu_count())

"""
This code is valid only for square domains. 

Returns:
    _type_: _description_
"""

# --- Parameters
## Change the vbalues below based on how you run the simulations
endTimestep = np.inf
stride = 100#25#100
a_ = 1.0#1.5#1.0 #1.5
Ca_ = 20e-2
In_ = 10e-2
a_rep_ = 1.0
a_adh_ = 1.5
potential_type_ = 2

if len(sys.argv) > 3:
    #print('using customised Ca')
    Ca_ = float(sys.argv[3])

if len(sys.argv) > 4:
    
    if (sys.argv[4] == 'old'):
        potential_type_ = 1
        #print('using new potential')
    elif (sys.argv[4] == 'new'):
        potential_type_ = 2
        #print('using new potential')
        

if len(sys.argv) > 5:
    #print('using customised a or a_rep')
    a_ = float(sys.argv[5])
    a_rep_ = float(sys.argv[5])
    #print(a_rep_)

if len(sys.argv) > 5:
    #print('using customised a_adh')
    a_adh_ = float(sys.argv[6])
    #print(a_adh_)

pbc=True
domain_size = np.array([100.0,100.0])
confined = False
new_interaction = False #False - old B*w interaction and True - new interaction

quantity = 'force'
if int(sys.argv[2]) == 1: 
    quantity = 'stress'  #also saves energy
    #print('checkpoint 1: stress')
if int(sys.argv[2]) == 2:
    #print('checkpoint 1: force from stored mu')
    quantity = 'force'
if int(sys.argv[2]) == 3:
    #print('checkpoint 1: free energy density')
    quantity = 'free_energy_density'
    


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
        

"""
if len(sys.argv) > 3:
    endTimestep = int(sys.argv[3])
else:
    endTimestep = np.inf

if len(sys.argv) > 4:
    stride = int(sys.argv[4])
else:
    stride = 100
"""

dx=0.5
# Interpolation
x_ = np.linspace(0,domain_size[0],interpolation_steps)
y_ = np.linspace(0,domain_size[1],interpolation_steps)
xx_,yy_ = np.meshgrid(x_,y_)
xx_ = np.reshape(xx_,(interpolation_steps**2,1))
yy_ = np.reshape(yy_,(interpolation_steps**2,1))

dx = domain_size[0]/interpolation_steps #assuming dx = dy
np.save(out_dir + '/grid_x',np.reshape(xx_,(interpolation_steps,interpolation_steps)))
np.save(out_dir + '/grid_y',np.reshape(yy_,(interpolation_steps,interpolation_steps)))


def calculate_parallel_stress(parameters):
    indtlist = parameters['indtlist']
    folder = parameters['folder']
    stride = parameters['stride']
    interpolation_steps = parameters['interpolation_steps']
    Ca_ = parameters['Ca_']
    In_ = parameters['In_']
    a_adh_ = parameters['a_adh_']
    a_rep_ = parameters['a_rep_']
    a_ = parameters['a_']
    dx = parameters['dx']
    potential_type_ = parameters['potential_type_']
    domain_size = parameters['domain_size']
    
    # Interpolation
    x = np.linspace(0,domain_size[0],interpolation_steps)
    y = np.linspace(0,domain_size[1],interpolation_steps)
    xx,yy = np.meshgrid(x,y)
    xx = np.reshape(xx,(interpolation_steps**2,1))
    yy = np.reshape(yy,(interpolation_steps**2,1))

    print('checkpoint 3')
    
    file_pattern = sys.argv[1] + 'neo_positions_p*.csv'
    out_dir = sys.argv[1] + 'stress_fields_500'
    
    times_ = np.load(folder + '/global_fields_200/timesteps.npy')
    
    
    positions_raw = []
    ranks = []
    count = 0
    for filename in glob.glob(file_pattern):
        count += 1
        tmp = pd.read_csv(filename)#, engine='python')#-fwf')#, engine='python')
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
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    #row_indices = np.arange(0,positions_raw[0].index.shape[0],stride,dtype=int) uncomment this if 0 problem is fixed
    row_indices = np.arange(stride-1,positions_raw[0].index.shape[0],stride,dtype=int)

    # we now have one particular timepoint
    #print((times_==time))
    #print(np.where(times_==time))
    #print(neighbours_arr[np.where(times_==time)].shape)
    #print(neighbours_arr[np.where(times_==time)])
    
    neighbours_arr = load_files.get_neighbour_relations(folder, 0, 90000, stride) #reading out neighbour relations
    
    for indt in indtlist:
        time = times_[indt]
        neighbours = neighbours_arr[np.where(times_==time)][0]
        
        phi_cell = np.zeros([len(ranks), interpolation_steps, interpolation_steps])
        #mu_cell = np.zeros([len(ranks), interpolation_steps, interpolation_steps])
        print(time)
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
            #mu = vtk_to_numpy(data.GetPointData().GetArray(1))

            phi_interp = griddata(points[:,0:2],phi,(xx,yy),method='nearest')
            phi_interp = np.reshape(phi_interp,(interpolation_steps,interpolation_steps))
            
            #mu_interp = griddata(points[:,0:2],mu,(xx,yy),method='nearest')
            #mu_interp = np.reshape(mu_interp,(interpolation_steps,interpolation_steps))
            
            phi_cell[rank] = phi_interp
            #mu_cell[rank] = mu_interp
        
        sigma_00 = np.zeros(phi_cell[0].shape)
        sigma_11 = np.zeros(phi_cell[0].shape)
        sigma_01 = np.zeros(phi_cell[0].shape)
        free_energy = np.zeros(phi_cell[0].shape)

        #curvature_cell0 = np.zeros(phi_cell[0].shape)
        #curvature_cell1 = np.zeros(phi_cell[0].shape)
        print('checkpoint 5')
        for cell_n in ranks:
            grad_phi_ = grad_phi(phi_cell[cell_n], dx=dx, pbc=True)
            
            isotropic = np.zeros(phi_cell[0].shape)
            
            f, mu = f_mu_CH(phi_cell[cell_n], grad_phi_, epsilon=0.15, Ca=Ca_)
            f_ = np.copy(f)
            isotropic = f
            isotropic -= 2*mu*(phi_cell[cell_n]+1)/2  #note these 2's are because of conversion of range from [-1,1] to [0,1]
            
            
            #for cell_k in [cell_k for cell_k in rank if cell_k!=cell_n]:#considers all cells 
            for cell_k in np.where(neighbours[cell_n]==1)[0]: #considers only neighbours, dependent on method of neighbour determination
                if potential_type_ == 2:
                    fi, mui = f_mu_IN_nk_neo(phi_cell[cell_n], phi_cell[cell_k], In=In_, a_rep=a_rep_, a_adh = a_adh_)
                elif potential_type_ == 1:
                    fi, mui = f_mu_IN_nk(phi_cell[cell_n], phi_cell[cell_k], In=In_, a=a_)
                f_ += np.copy(fi)
                isotropic += fi
                isotropic -= 2*mui*(phi_cell[cell_n]+1)/2
            
            #if method == 2: #doesn't calculate mu but uses mu_cell
            #    isotropic -= 2*mu_cell[cell_n]*(phi_cell[cell_n]+1)/2
            epsilon=0.15
            sigma_00 += isotropic - epsilon*grad_phi_[0]**2
            sigma_11 += isotropic - epsilon*grad_phi_[1]**2
            sigma_01 += -epsilon*grad_phi_[0]*grad_phi_[1]
            free_energy += f_

            #if cell_n == 0:
            #    curvature_cell0 = curvature(phi_cell[cell_n], dx=dx, pbc=True)
            #if cell_n == 1:
            #    curvature_cell1 = curvature(phi_cell[cell_n], dx=dx, pbc=True)
        print('checkpoint 6')        
        np.save(out_dir + '/free_energy' + '_{:06.3f}'.format(time), free_energy)
            
        #np.save(out_dir + '/sigma_00' + '_{:06.3f}'.format(time), sigma_00)
        #np.save(out_dir + '/sigma_11' + '_{:06.3f}'.format(time), sigma_11)
        #np.save(out_dir + '/sigma_01' + '_{:06.3f}'.format(time), sigma_01)

        #np.save(out_dir + '/curvature_cell0' + '_{:06.3f}'.format(time), curvature_cell0)
        #np.save(out_dir + '/curvature_cell1' + '_{:06.3f}'.format(time), curvature_cell1)
            
    return 0
print('checkpoint 2d')

processes=20

times_ = np.load(sys.argv[1] + '/global_fields_200/timesteps.npy')
nt = len(times_)
ntlist = np.arange(nt)
parameter_set = []
for i in range(processes):
    parameters = {}
    parameters['indtlist'] = ntlist[ntlist%processes==i]
    print(ntlist[ntlist%processes==i])
    parameters['folder'] = sys.argv[1]
    parameters['stride'] = stride
    parameters['interpolation_steps'] = interpolation_steps
    parameters['Ca_'] = Ca_
    parameters['In_'] = In_
    parameters['a_adh_'] = a_adh_
    parameters['a_rep_'] = a_rep_
    parameters['a_'] = a_
    parameters['dx'] = dx
    parameters['potential_type_'] = potential_type_
    parameters['domain_size'] = domain_size
    parameter_set.append(parameters)

pool = mp.Pool(processes=20)#mp.cpu_count())
result = pool.map(calculate_parallel_stress, parameter_set)


print('checkpoint 2b')
if quantity == 'strdess':
    ncores = 31
    pool = mp.Pool(processes=ncores)
    chunksize=int(len(times_)/(ncores+1))+1
    print('checkpoint 2c')
    result = pool.map(calculate_parallel_stress, times_, chunksize=5)
    print('checkpoint 2d')

"""
pool = mp.Pool(21)
pool.map(calculate_parallel_stress, range(len(times_)))
"""
"""
#https://www.digitalocean.com/community/tutorials/python-multiprocessing-example
def parallel_queue(todo, done):
    while True:
        print('checkpoint parallel_queue')
        try:
            '''
                try to get task from the queue. get_nowait() function will 
                raise queue.Empty exception if the queue is empty. 
                queue(False) function would do the same task also.
            '''
            indt = todo.get_nowait()
            print(f"{indt} is what we will do now")
            calculate_parallel_stress(indt)
        except queue.Empty:
            break
        else:
            '''
                if no exception has been raised, 
            '''
            print(f"{indt} performed by {mp.current_process.name}")
            done.put(indt)
            time.sleep(0.5)
    return True

number_of_tasks = len(times_)
number_of_processes = 30
todo=mp.Queue()
done=mp.Queue()

for i in range(number_of_tasks):
    todo.put(i)

#creating processes
processes=[]
for indp in range(number_of_processes):
    p = mp.Process(target=parallel_queue, args=(todo, done))
    processes.append(p)
    p.start()
    
#completing processes
for p in processes:
    p.join()
    
while not done.empty():
    print(done.get())   
"""