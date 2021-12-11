import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import sys
from matplotlib.ticker import MaxNLocator

from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import griddata


plt.style.use('seaborn-bright')
#plt.style.use('dark_background')
positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "total_int": 10,\
     "neighbours": 11, "confine_int": 12, "growth_rate": 13}

outdir = '/scratch/ws/1/haja565a-workspace2/master_thesis/output'
savedir = '/scratch/ws/1/haja565a-workspace2/quant/'
expName = sys.argv[1]
totalCells = int(sys.argv[2])
dt = 0.005

if not os.path.exists(savedir+expName):
    os.makedirs(savedir+expName)

T = np.genfromtxt(outdir + expName + '/positions_p0.csv', delimiter=',',skip_header=1)[:,0] #time

timesteps = len(T)
Time = np.zeros([totalCells, timesteps])
print(timesteps)

X = np.zeros([totalCells, timesteps]) #posX
Y = np.zeros([totalCells, timesteps]) #posy
S0 = np.zeros([totalCells, timesteps]) #S0
S1 = np.zeros([totalCells, timesteps]) #S1
angle = np.zeros([totalCells, timesteps]) #angle of elongation axis
Vx = np.zeros([totalCells, timesteps]) #velX
Vy = np.zeros([totalCells, timesteps]) #velY

radius = np.zeros([totalCells, timesteps]) #radius
total_int = np.zeros([totalCells, timesteps]) #total_interaction
neighbours = np.zeros([totalCells, timesteps], dtype=int) #neighbours
growth_rate = np.zeros([totalCells, timesteps]) #growth_rate
confine_int = np.zeros([totalCells, timesteps])

print("Reading Data:" + expName)
for cell in range(totalCells):
    X[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,2] 
    Y[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,3] 
    S0[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,5] 
    S1[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,6]
    Vx[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,7]/dt 
    Vy[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,8]/dt 
    angle[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,9] 
    radius[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,4] 
    total_int[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,10] 
    confine_int[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,12]
    neighbours[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,11]
    growth_rate[cell] = np.genfromtxt(outdir + expName + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,13] 
    #growth_rate[cell][0:-1] = (np.square(radius[cell][1:])/(np.square(radius[cell][0:-1]))-1)/dt
    Time[cell] = T

np.save(savedir + expName + "/S0.npy", S0)
np.save(savedir + expName + "/S1.npy", S1)

np.save(savedir + expName + "/Vx.npy", Vx)
np.save(savedir + expName + "/Vy.npy", Vy)
np.save(savedir + expName + "/angle.npy", angle)
np.save(savedir + expName + "/radius.npy", radius)
np.save(savedir + expName + "/total_int.npy", total_int)
np.save(savedir + expName + "/confine_int.npy", confine_int)
np.save(savedir + expName + "/neighbours.npy", neighbours)
np.save(savedir + expName + "/X.npy", X)
np.save(savedir + expName + "/Y.npy", Y)
np.save(savedir + expName + "/growth_rate.npy", growth_rate)
np.save(savedir + expName + "/T.npy", T)

print("Finished Saving Data:" + expName)