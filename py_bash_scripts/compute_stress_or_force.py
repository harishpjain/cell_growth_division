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
import load_files as load_files

"""
This code is valid only for square domains. 

Returns:
    _type_: _description_
"""

print(sys.argv)

# --- Parameters
## Change the vbalues below based on how you run the simulations
endTimestep = np.inf
stride = 100#25#100
a_ = 1.0#1.5#1.0 #1.5
Ca_ = 20e-2
In_ = 10e-2

if len(sys.argv) > 3:
    print('using customised Ca')
    Ca_ = float(sys.argv[3])
pbc=True
domain_size = np.array([100.0,100.0])
confined = False

quantity = 'force'
if int(sys.argv[2]) == 1:
    quantity = 'stress'
    print('checkpoint 1: stress')
if int(sys.argv[2]) == 2:
    print('checkpoint 1: force from stored mu')
    quantity = 'force'
if int(sys.argv[2]) == 3:
    print('checkpoint 1: free energy density')
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
        
        


print(file_pattern)

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


positions_raw = []
ranks = []
count = 0
for filename in glob.glob(file_pattern):
    print(count)
    print(filename)
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

# Interpolation
x = np.linspace(0,domain_size[0],interpolation_steps)
y = np.linspace(0,domain_size[1],interpolation_steps)
xx,yy = np.meshgrid(x,y)
xx = np.reshape(xx,(interpolation_steps**2,1))
yy = np.reshape(yy,(interpolation_steps**2,1))
np.save(out_dir + '/grid_x',np.reshape(xx,(interpolation_steps,interpolation_steps)))
np.save(out_dir + '/grid_y',np.reshape(yy,(interpolation_steps,interpolation_steps)))

dx = domain_size[0]/interpolation_steps #assuming dx = dy
    
print(Ca_)
print(a_)
print(In_)

reader = vtk.vtkXMLUnstructuredGridReader()
print(stride)
print(positions_raw[0].index.shape)
print(positions_raw[0].index)
#row_indices = np.arange(0,positions_raw[0].index.shape[0],stride,dtype=int) uncomment this if 0 problem is fixed
row_indices = np.arange(stride-1,positions_raw[0].index.shape[0],stride,dtype=int)
times = []
times_ = np.load(sys.argv[1] + '/global_fields_200/timesteps.npy')
count = 0

neighbours_arr = load_files.get_neighbour_relations(sys.argv[1], 0, 30000, stride) #reading out neighbour relations


print('checkpoint 2')
def grad_phi(phi_n, dx=dx, pbc=True):
    if pbc:
        grad_phix = np.gradient(np.pad(phi_n, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        grad_phiy = np.gradient(np.pad(phi_n, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
    else:
        grad_phix = np.gradient(phi_n, dx, axis=0)
        grad_phiy = np.gradient(phi_n, dx, axis=1)
    
    return np.array([grad_phix, grad_phiy])

def g(phi_n):
    #return 3/(2*(1-phi_n)*(1+phi_n))
    return 1

def W(phi_n):
    return ((phi_n**2 - 1)**2)/4 

def grad_W(phi_n):
    return phi_n**3 - phi_n

def w_int(phi_k, a=1.5):
    return 1 - ((a+1)*(phi_k-1)**2/4) + (a*(phi_k-1)**4/16)

def w_int_der(phi_k, a=1.5):
    return -(2*(a+1)*(phi_k-1)/4) + (4*a*(phi_k-1)**3/16)

def laplace_phi2(gradphi_n, dx=dx, pbc = True):
    if pbc:
        grad_phixx = np.gradient(np.pad(gradphi_n[0], 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
        grad_phiyy = np.gradient(np.pad(gradphi_n[1], 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
    else:
        grad_phixx = np.gradient(gradphi_n[0], dx, axis=0)
        grad_phiyy = np.gradient(gradphi_n[1], dx, axis=1)
    return grad_phixx + grad_phiyy

def f_mu_IN_nk(phi_n, phi_k, In=10e-2, a=1.5):
    #returns free energy and chemical potential of interaction energy
    w = w_int(phi_k, a)
    w_der = w_int_der(phi_n, a)
    f = (1/In)*((phi_n+1)/2)*(w)
    mu = (1/In)*((0.5*w) + w_der*(phi_k+1)/2)
    return f, mu

def f_mu_CH(phi_n, gradphi_n, epsilon=0.15, Ca=20e-2, dx=dx):
    #returns free energy and chemical potential of Cahn Hilliard energy
    f = (1/Ca)*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))
    mu = (1/Ca)*((grad_W(phi_n)/epsilon) - (epsilon*(laplace_phi2(gradphi_n, dx=dx))))
    return f, mu

def f_CH(phi_n, gradphi_n, epsilon=0.15, Ca=20e-2):
    #return (1/Ca)*(g(phi_n))*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))
    return (1/Ca)*((W(phi_n)/epsilon) + (epsilon*(gradphi_n[0]**2+gradphi_n[1]**2)/2))

def f_IN_nk(phi_n, phi_k, In=10e-2, a=1.5):
    #returns free energy and chemical potential of interaction energy
    w = w_int(phi_k, a)
    f = (1/In)*((phi_n+1)/2)*(w)
    return f

if quantity == 'stress':
    #####for ind in row_indices:
    for time in times_:
        # we now have one particular timepoint
        #####time = positions_raw[0].iloc[ind]['time']
        print(time)
        #if time < 5:
        #    continue
        #if not os.path.exists(sys.argv[1] + 'data/phase_p0_' + '{:06.3f}'.format(time) + '.vtu'):
        #    continue
        times.append(time)
        neighbours = neighbours_arr[np.where(times_==time)][0]
        
        phi_cell = np.zeros([len(ranks), interpolation_steps, interpolation_steps])
        #mu_cell = np.zeros([len(ranks), interpolation_steps, interpolation_steps])

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
        
        for cell_n in ranks:
            grad_phi_ = grad_phi(phi_cell[cell_n], dx=dx, pbc=True)
            
            isotropic = np.zeros(phi_cell[0].shape)
            
            f, mu = f_mu_CH(phi_cell[cell_n], grad_phi_, epsilon=0.15, Ca=Ca_)
            f_ = np.copy(f)
            isotropic = f
            isotropic -= 2*mu*(phi_cell[cell_n]+1)/2  #note these 2's are because of conversion of range from [-1,1] to [0,1]
            
            
            #for cell_k in [cell_k for cell_k in rank if cell_k!=cell_n]:#considers all cells 
            for cell_k in np.where(neighbours[cell_n]==1)[0]: #considers only neighbours, dependent on method of neighbour determination
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
            
        np.save(out_dir + '/free_energy' + '_{:06.3f}'.format(time), free_energy)
            
        np.save(out_dir + '/sigma_00' + '_{:06.3f}'.format(time), sigma_00)
        np.save(out_dir + '/sigma_11' + '_{:06.3f}'.format(time), sigma_11)
        np.save(out_dir + '/sigma_01' + '_{:06.3f}'.format(time), sigma_01)
            
    np.save(out_dir + '/timesteps',np.array(times))
    np.save(out_dir + '/ranks',np.array(ranks))
    
print('checkpoint 3')
if quantity == 'force':
    print('checkpoint 4')
    for ind in row_indices:
        # we now have one particular timepoint
        time = positions_raw[0].iloc[ind]['time']
        print(time)
        #if time < 5:
        #    continue
        #if not os.path.exists(sys.argv[1] + 'data/phase_p0_' + '{:06.3f}'.format(time) + '.vtu'):
        #    continue
        times.append(time)
        neighbours = neighbours_arr[np.where(times_==time)][0]
        
        force_0 = np.zeros([interpolation_steps, interpolation_steps])
        force_1 = np.zeros([interpolation_steps, interpolation_steps])

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
            mu = vtk_to_numpy(data.GetPointData().GetArray(1))

            phi_interp = griddata(points[:,0:2],phi,(xx,yy),method='nearest')
            phi_interp = np.reshape(phi_interp,(interpolation_steps,interpolation_steps))
            
            mu_interp = griddata(points[:,0:2],mu,(xx,yy),method='nearest')
            mu_interp = np.reshape(mu_interp,(interpolation_steps,interpolation_steps))
            
            smoothen = True #shift htis later
            phi_interp = gaussian_filter(phi_interp, 3, mode='wrap')
            mu_interp = gaussian_filter(mu_interp, 3, mode='wrap')
            if pbc:
                grad_mu_interpx = np.gradient(np.pad(mu_interp, 2, mode="wrap"), dx, axis=0)[2:-2, 2:-2]
                grad_mu_interpy = np.gradient(np.pad(mu_interp, 2, mode="wrap"), dx, axis=1)[2:-2, 2:-2]
                
            #Force is divergence of stress = \partial_j \sigma_ij = F_i
            #divergence of stress is also equal to \sum_n -(\phi_i * \partial_i \mu_n) = F_i
            
            
            
            
            force_0 += phi_interp*grad_mu_interpx
            force_1 += phi_interp*grad_mu_interpy
        
        np.save(out_dir + '/force_0' + '_{:06.3f}'.format(time), force_0)
        np.save(out_dir + '/force_1' + '_{:06.3f}'.format(time), force_1)

            
    np.save(out_dir + '/timesteps',np.array(times))
    np.save(out_dir + '/ranks',np.array(ranks))
    
if quantity == 'free_energy_density':
    #####for ind in row_indices:
    for time in times_:
        # we now have one particular timepoint
        #####time = positions_raw[0].iloc[ind]['time']
        print(time)
        #if time < 5:
        #    continue
        #if not os.path.exists(sys.argv[1] + 'data/phase_p0_' + '{:06.3f}'.format(time) + '.vtu'):
        #    continue
        times.append(time)
        neighbours = neighbours_arr[np.where(times_==time)][0]
        
        phi_cell = np.zeros([len(ranks), interpolation_steps, interpolation_steps])
        #mu_cell = np.zeros([len(ranks), interpolation_steps, interpolation_steps])

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
        
        free_energy = np.zeros(phi_cell[0].shape)
        
        for cell_n in ranks:
            grad_phi_ = grad_phi(phi_cell[cell_n], dx=dx, pbc=True)
            isotropic = np.zeros(phi_cell[0].shape)
            
            f = f_CH(phi_cell[cell_n], grad_phi_, epsilon=0.15, Ca=Ca_)

            #for cell_k in [cell_k for cell_k in rank if cell_k!=cell_n]:#considers all cells 
            for cell_k in np.where(neighbours[cell_n]==1)[0]: #considers only neighbours
                f += f_IN_nk(phi_cell[cell_n], phi_cell[cell_k], In=In_, a=a_)
            
            #if method == 2: #doesn't calculate mu but uses mu_cell
            #    isotropic -= 2*mu_cell[cell_n]*(phi_cell[cell_n]+1)/2
            epsilon=0.15
            free_energy += f
            
        np.save(out_dir + '/free_energy' + '_{:06.3f}'.format(time), free_energy)
            
    np.save(out_dir + '/timesteps',np.array(times))
    np.save(out_dir + '/ranks',np.array(ranks))