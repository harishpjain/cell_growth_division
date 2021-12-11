import numpy as np
from numpy.lib.npyio import save
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
import matplotlib.tri as mtri
import os
import sys
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import rank_filter
import matplotlib.pyplot as plt
from shapely.geometry import LineString
import math
import matplotlib.ticker as ticker
import matplotlib.colors as colors

plt.style.use('seaborn-bright')
#plt.style.use('dark_background')
positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "total_int": 10,\
     "neighbours": 11, "growth_rate": 12}

outdir = '/scratch/ws/1/haja565a-workspace2/master_thesis/output'
savedir = "/scratch/ws//1/haja565a-workspace2/DeformationField/"
expName = sys.argv[1]
totalCells = int(sys.argv[2])
time = float(sys.argv[3])
ranks = list(range(0,totalCells))

print("saving" + expName)

saveVTK = False

dt = 0.005

if not os.path.exists(savedir+expName):
    os.makedirs(savedir+expName)


interpolation_steps = 1000
domain_size = 100
x = np.linspace(0, domain_size,interpolation_steps)
xx, yy = np.meshgrid(x,x)

cx, cy = xx-50, yy-50
angpos = np.arctan2(cy,cx) #angular posiiton of the pixel

positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "total_int": 10,\
     "neighbours": 11, "confine_int": 12, "growth_rate": 13}

#ageplot = False
#if(ageplot):
age = np.load("/scratch/ws//1/haja565a-workspace2/quant/" + expName + "/age.npy")
agelist = []

Vx = []
Vy = []
S0 = []
S1 = []
interaction_pot = []
phi_all = []
reader = vtk.vtkXMLUnstructuredGridReader()
for rank_ind,rank in enumerate(ranks):
    #time = T[t*save_interval]
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

    #if(ageplot):
    agelist.append(age[rank, ind])
    
    positionfile = outdir + sys.argv[1] + '/positions_p' + str(rank) + '.csv'
    positiondata = np.genfromtxt(positionfile, delimiter=',',skip_header=1)

    S0.append(positiondata[ind][5])
    S1.append(positiondata[ind][6])
    interaction_pot.append(positiondata[ind][10])
    Vx.append(positiondata[ind][7])
    Vy.append(positiondata[ind][8])
    

# after this we have a list IN THE SAME ORDER AS ranks with all phi
# now the axis 0 here is the rank axis which we want to remove
phi_all = np.array(phi_all)

S0 = np.array(S0)
S1 = np.array(S1)
interaction_pot = np.array(interaction_pot)
Vx = np.array(Vx)
Vy = np.array(Vy)
# global phasefield, given by a lot of 1s and something in between

phi_glob = np.max(phi_all,axis=0)
# this is more interesting -> this locally gives the rank the contributed to the max, ie the cell
rank_max = np.argmax(phi_all,axis=0)
print(rank_max)




agelist = np.array(agelist)
#phi_temp = 2*phi_glob-1.
#agefield = agelist[rank_max]*(np.tanh(10*(phi_glob-0.8))+1.0)*0.5
agefield = agelist[rank_max]*np.sign(phi_glob-0.9)

"""
if(ageplot):
    agelist = np.array(agelist)
    phi_temp = 2*phi_glob-1.0
    agefield = agelist[rank_max]*(phi_temp)

    figa, axa = plt.subplots()
    Ages = axa.pcolormesh(xx, yy, agefield, norm =colors.LogNorm(vmin=1, vmax=30000), cmap = "gist_rainbow")
    #Ages = axa.pcolormesh(xx, yy, agefield*dt, norm =colors.LogNorm(vmin=0.005, vmax=10000*dt), cmap = "gist_rainbow")
    axa.set_xlabel("X")
    axa.set_ylabel("Y")
    axa.axis('square')
    #Ages = axa.pcolormesh(xx, yy, agefield, norm =colors.TwoSlopeNorm(vmin=0, vcenter=1, vmax=20000), cmap = "gist_rainbow")
    figa.colorbar(Ages, ax=axa, label = "Age(in timesteps)")
    #figa.colorbar(Ages, ax=axa, label = "Age")
    plt.savefig(savedir + expName + "/age_log"+ '{:05d}'.format(ind) + '.png')

    figa, axa = plt.subplots()
    #Ages = axa.pcolormesh(xx, yy, agefield, norm =colors.LogNorm(vmin=1, vmax=100), cmap = "gist_rainbow")
    Ages = axa.pcolormesh(xx, yy, agefield, norm =colors.Normalize(vmin=1, vmax=30000), cmap = "gist_rainbow")
    axa.set_xlabel("X")
    axa.set_ylabel("Y")
    axa.axis('square')
    figa.colorbar(Ages, ax=axa, label = "Age(in timesteps)")
    plt.savefig(savedir + expName + "/age"+ '{:05d}'.format(ind) + '.png')
"""

Vxglob = Vx[rank_max] * np.sign(phi_glob)
Vyglob = Vy[rank_max] * np.sign(phi_glob)

#normV = np.sqrt(Vxglob**2 + Vyglob**2) #velocity magnitude
#Uxglob = np.divide(Vxglob, normV, out=np.zeros_like(Vxglob), where=normV!=0) #unit x vel
#Uyglob = np.divide(Vyglob, normV, out=np.zeros_like(Vyglob), where=normV!=0) #unit y vel
#Vxglobf = gaussian_filter(Vxglob, sigma=40, mode = 'constant')#, mode='nearest') #filter velx
#Vyglobf = gaussian_filter(Vyglob, sigma=40, mode = 'constant')#, mode='nearest') #filter vely
#Vxglobf = gaussian_filter(Vxglob, sigma=80, mode = 'constant')#, mode='nearest') #filter velx
#Vyglobf = gaussian_filter(Vyglob, sigma=80, mode = 'constant')#, mode='nearest') #filter vely
Vxglobf = gaussian_filter(Vxglob, sigma=40, mode = 'constant')#, mode='nearest') #filter velx
Vyglobf = gaussian_filter(Vyglob, sigma=40, mode = 'constant')#, mode='nearest') #filter vely
#vorticity = np.gradient(Vyglobf, axis = 1) - np.gradient(Vxglobf, axis = 0) 
#vorticity = np.multiply(vorticity, np.sign(phi_glob)) #vorticity


if(saveVTK):
    np.savetxt(savedir+expName+'/Vx'+ '{:06.3f}'.format(time) +'.csv', Vxglobf, delimiter =",")
    np.savetxt(savedir+expName+'/Vy'+ '{:06.3f}'.format(time) +'.csv', Vyglobf, delimiter =",")


S0_glob = S0[rank_max] * phi_glob
S1_glob = S1[rank_max] * phi_glob
interaction_pot_glob = interaction_pot[rank_max] *np.sign(phi_glob-0.9)

S0_glob = gaussian_filter(S0_glob, sigma=80, mode='nearest')
S1_glob = gaussian_filter(S1_glob, sigma=80, mode='nearest')

#a = np.gradient(S0_glob, axis = 0)
#b = np.gradient(S0_glob, axis = 1)
#c = np.gradient(S1_glob, axis = 0)
#d = np.gradient(S1_glob, axis = 1)
#dow = a*d-b*c

Vradtot = Vxglobf*np.cos(angpos)+Vyglobf*np.sin(angpos)
Vorttot = Vxglobf*np.sin(angpos)-Vyglobf*np.cos(angpos)

dist_cent = np.hypot(cx, cy)

print("calculation done. Now saving")

np.save(savedir + expName + "/age"+ '{:06.3f}'.format(time) + ".npy", agefield)
np.save(savedir + expName + "/Vx"+ '{:06.3f}'.format(time) + ".npy", Vxglobf)
np.save(savedir + expName + "/Vx"+ '{:06.3f}'.format(time) + ".npy", Vxglobf)
np.save(savedir + expName + "/Vy"+ '{:06.3f}'.format(time) + ".npy", Vyglobf)
np.save(savedir + expName + "/S0"+ '{:06.3f}'.format(time) + ".npy", S0_glob)
np.save(savedir + expName + "/S1"+ '{:06.3f}'.format(time) + ".npy", S1_glob)
np.save(savedir + expName + "/interaction"+ '{:06.3f}'.format(time) + ".npy", interaction_pot_glob)
np.save(savedir + expName + "/phi"+ '{:06.3f}'.format(time) + ".npy", phi_glob)
#np.save(savedir + expName + "/Vrad.npy", Vradtot)
#np.save(savedir + expName + "/Vort.npy", Vorttot)
np.save(savedir + expName + "/dist_cent"+ '{:06.3f}'.format(time) + ".npy", dist_cent)

print("saving" + expName + "successful")
