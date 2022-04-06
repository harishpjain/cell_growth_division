import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import sys
from matplotlib.ticker import MaxNLocator
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import griddata


outdir = '/beegfs/ws/1/haja565a-workspacebeegfs/phd18122201/output'
#savedir = '/beegfs/ws/1/haja565a-workspacebeegfs/phd18122201/'
#expName = sys.argv[1]
outdir = sys.argv[1]
totalCells = int(sys.argv[2])
dt = 0.005

if not os.path.exists(outdir + '/quant'):
    os.makedirs(outdir + '/quant')

T = np.genfromtxt(outdir + '/positions_p0.csv', delimiter=',',skip_header=1)[:,0] #time

timesteps = len(T)
Time = np.zeros([totalCells, timesteps])
print('The simulation ran for ' + str(timesteps) + ' timesteps')
print('saving data')

savePosition = True
saveRadius = True
saveSTensor = True
saveVelocity = True
saveVelocityOrientation = True
saveNematicOrientation = True
saveTotalInteraction = True
saveConfinementInteraction = True
saveNeighbourNumber = True
saveGrowthRate = True

if(savePosition):
    print('saving positions')
    X = np.zeros([totalCells, timesteps])
    Y = np.zeros([totalCells, timesteps])
    for cell in range(totalCells):
        X[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,2]
        Y[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,3]
    np.save(outdir + '/quant' + "/X.npy", X)
    np.save(outdir + '/quant' + "/Y.npy", Y)

if(saveVelocity):
    print('saving velocity')
    Vx = np.zeros([totalCells, timesteps]) #velX
    Vy = np.zeros([totalCells, timesteps]) #velY
    for cell in range(totalCells):
        Vx[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,7]/dt 
        Vy[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,8]/dt 
    np.save(outdir + '/quant' + "/velXSim.npy", Vx)
    np.save(outdir + '/quant' + "/velYSim.npy", Vy)

if(saveSTensor):
    print('saving S Tensor')
    S0 = np.zeros([totalCells, timesteps])
    S1 = np.zeros([totalCells, timesteps]) 
    for cell in range(totalCells):
        S0[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,5] 
        S1[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,6] 
    np.save(outdir + '/quant' + "/S0.npy", S0)
    np.save(outdir + '/quant' + "/S1.npy", S1)

if(saveRadius):
    print('saving radius')
    R = np.zeros([totalCells, timesteps])
    for cell in range(totalCells):
        R[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,4]
    np.save(outdir + '/quant' + "/radius.npy", R)

if(saveVelocityOrientation):
    print('saving velocity orientations')
    theta = np.zeros([totalCells, timesteps])
    for cell in range(totalCells):
        theta[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,9]
    np.save(outdir + '/quant' + "/velocityAngle.npy", theta)

if(saveVelocityOrientation):
    print('saving nematic orientations')
    alpha = np.zeros([totalCells, timesteps])
    for cell in range(totalCells):
        alpha[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,10]
    np.save(outdir + '/quant' + "/nematicAngle.npy", alpha)

if(saveTotalInteraction or saveConfinementInteraction):
    print('saving interactions')
    totalInt = np.zeros([totalCells, timesteps])
    confineInt = np.zeros([totalCells, timesteps])
    for cell in range(totalCells):
        totalInt[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,11]
        confineInt[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,13]
    np.save(outdir + '/quant' + "/totalInteraction.npy", totalInt)
    np.save(outdir + '/quant' + "/confinementInteraction.npy", confineInt)

if(saveNeighbourNumber):
    print('saving number of neighbours')
    neighbours = np.zeros([totalCells, timesteps], dtype = int)
    for cell in range(totalCells):
        neighbours[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,12]
    np.save(outdir + '/quant' + "/neighbours.npy", neighbours)


if(saveGrowthRate):
    print('saving growth rate')
    growthRate = np.zeros([totalCells, timesteps])
    for cell in range(totalCells):
        growthRate[cell] = np.genfromtxt(outdir + '/positions_p' + str(cell) +'.csv', delimiter=',',skip_header=1)[:,14]
    np.save(outdir + '/quant' + "/growthRate.npy", growthRate)