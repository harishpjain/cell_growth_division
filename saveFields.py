import numpy as np
from numpy.lib.npyio import save
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk

import os
import sys
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import rank_filter
import matplotlib.pyplot as plt
from shapely.geometry import LineString
import math
import matplotlib.ticker as ticker

plt.style.use('seaborn-bright')
#plt.style.use('dark_background')
positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "total_int": 10,\
     "neighbours": 11, "growth_rate": 12}

outdir = '/scratch/ws/1/haja565a-workspace2/master_thesis/output'
savedir = '/scratch/ws/1/haja565a-workspace2/quant/'
expName = sys.argv[1]
totalCells = int(sys.argv[2])
ranks = list(range(0,totalCells))

dt = 0.005

if not os.path.exists(savedir+expName):
    os.makedirs(savedir+expName)

T = np.genfromtxt(outdir + expName + '/positions_p0.csv', delimiter=',',skip_header=1)[:,0] #time

save_interval = 1000 #40 #all data is saved every 40 timesteps

timesteps = len(T)
savesteps = int(timesteps/save_interval)
interpolation_steps = 100
domain_size = 100
x = np.linspace(0, domain_size,interpolation_steps)
xx, yy = np.meshgrid(x,x)

cx, cy = xx-50, yy-50
angpos = np.arctan2(cy,cx) #angular posiiton of the pixel

positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "total_int": 10,\
     "neighbours": 11, "confine_int": 12, "growth_rate": 13}

Vxtot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #Vxtot
Vytot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #Vytot
S0tot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #S0tot
S1tot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #S1tot
inttot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #interaction potential
phitot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) 
dowtot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #dow for finding signs of defects
Vradtot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #Radial Velocity
Vorttot = np.zeros([savesteps, interpolation_steps, interpolation_steps]) #Orthoradial velocity

for t in range(savesteps):
    Vx = []
    Vy = []
    S0 = []
    S1 = []
    interaction_pot = []
    phi_all = []
    reader = vtk.vtkXMLUnstructuredGridReader()
    for rank_ind,rank in enumerate(ranks):
        time = T[t*save_interval]
        ind = int(time/dt)
        filename = outdir + sys.argv[1] + '/data/phase_p' + str(rank) + '_' + '{:06.3f}'.format(time) + '.vtu'
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

        positionfile = outdir + sys.argv[1] + '/positions_p' + str(rank) + '.csv'
        positiondata = np.genfromtxt(positionfile, delimiter=',',skip_header=1)

        #S0.append(positiondata[ind][5])
        #S1.append(positiondata[ind][6])
        #interaction_pot.append(positiondata[ind][10])
        Vx.append(positiondata[ind][7])
        Vy.append(positiondata[ind][8])

    # after this we have a list IN THE SAME ORDER AS ranks with all phi
    # now the axis 0 here is the rank axis which we want to remove
    phi_all = np.array(phi_all)

    #S0 = np.array(S0)
    #S1 = np.array(S1)
    #interaction_pot = np.array(interaction_pot)
    Vx = np.array(Vx)
    Vy = np.array(Vy)
    # global phasefield, given by a lot of 1s and something in between

    phi_glob = np.max(phi_all,axis=0)
    # this is more interesting -> this locally gives the rank the contributed to the max, ie the cell
    rank_max = np.argmax(phi_all,axis=0)

    Vxglob = Vx[rank_max] * np.sign(phi_glob)
    Vyglob = Vy[rank_max] * np.sign(phi_glob)
    """
    normV = np.sqrt(Vxglob**2 + Vyglob**2) #velocity magnitude
    Uxglob = np.divide(Vxglob, normV, out=np.zeros_like(Vxglob), where=normV!=0) #unit x vel
    Uyglob = np.divide(Vyglob, normV, out=np.zeros_like(Vyglob), where=normV!=0) #unit y vel
    #Vxglobf = gaussian_filter(Vxglob, sigma=40, mode = 'constant')#, mode='nearest') #filter velx
    #Vyglobf = gaussian_filter(Vyglob, sigma=40, mode = 'constant')#, mode='nearest') #filter vely
    """
    Vxglobf = gaussian_filter(Vxglob, sigma=80, mode = 'constant')#, mode='nearest') #filter velx
    Vyglobf = gaussian_filter(Vyglob, sigma=80, mode = 'constant')#, mode='nearest') #filter vely
    """
    vorticity = np.gradient(Vyglobf, axis = 1) - np.gradient(Vxglobf, axis = 0) 
    vorticity = np.multiply(vorticity, np.sign(phi_glob)) #vorticity

    S0_glob = S0[rank_max] * phi_glob
    S1_glob = S1[rank_max] * phi_glob
    interaction_pot_glob = interaction_pot[rank_max] * phi_glob

    S0_glob = gaussian_filter(S0_glob, sigma=80, mode='nearest')
    S1_glob = gaussian_filter(S1_glob, sigma=80, mode='nearest')

    a = np.gradient(S0_glob, axis = 0)
    b = np.gradient(S0_glob, axis = 1)
    c = np.gradient(S1_glob, axis = 0)
    d = np.gradient(S1_glob, axis = 1)
    dow = a*d-b*c

    Vxtot[t] = Vxglobf
    Vytot[t] = Vyglobf
    S0tot[t] = S0_glob
    S1tot[t] = S1_glob
    inttot[t] = interaction_pot_glob
    """
    Vradtot[t] = Vxglobf*np.cos(angpos)+Vyglobf*np.sin(angpos)
    Vorttot[t] = Vxglobf*np.sin(angpos)-Vyglobf*np.cos(angpos)

    del reader

dist_cent = np.hypot(cx, cy)

print("calculation done. Now saving")

#np.save(savedir + expName + "/Vx.npy", Vxtot)
#np.save(savedir + expName + "/Vy.npy", Vytot)
#np.save(savedir + expName + "/S0.npy", S0tot)
#np.save(savedir + expName + "/S1.npy", S1tot)
#np.save(savedir + expName + "/interaction.npy", inttot)
np.save(savedir + expName + "/Vrad.npy", Vradtot)
np.save(savedir + expName + "/Vort.npy", Vorttot)
np.save(savedir + expName + "/dist_cent.npy", dist_cent)