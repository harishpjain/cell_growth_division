import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import sys
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy import stats
from scipy.optimize import least_squares
import matplotlib as mpl

plt.style.use('seaborn-bright')

outdir = '/scratch/ws/1/haja565a-workspace2/master_thesis/output'
savedir = '/scratch/ws/1/haja565a-workspace2/quant/'
expName = sys.argv[1]
totalCells = int(sys.argv[2])
dt = 0.005

if not os.path.exists(savedir+expName):
    os.makedirs(savedir+expName)

print("Reading Data:" + expName)
Vx = np.load(savedir + expName + "/Vx.npy")
Vy = np.load(savedir + expName + "/Vy.npy")
angle = np.load(savedir + expName + "/angle.npy")
radius = np.load(savedir + expName + "/radius.npy")
total_int = np.load(savedir + expName + "/total_int.npy")
confine_int = np.load(savedir + expName + "/confine_int.npy")
neighbours = np.load(savedir + expName + "/neighbours.npy")
X = np.load(savedir + expName + "/X.npy")
Y = np.load(savedir + expName + "/Y.npy")
growth_rate = np.load(savedir + expName + "/growth_rate.npy")
T = np.load(savedir + expName + "/T.npy")
timesteps = len(T)

Vx[radius<0.1]  = 0
Vy[radius<0.1]  = 0
X[radius<0.1]  = 0
Y[radius<0.1]  = 0

numCells = np.count_nonzero(radius>0.01, axis = 0) #change if cells are very small

age = np.load(savedir + expName  + "/age.npy")
age = age*dt



colony_volume_1 = np.sum(np.pi*(radius**2), axis=0)
fillIndex = np.where(colony_volume_1>=0.99*max(colony_volume_1))[0][0]

#making list of all ages and total int before confinement is filled
age_list = age[:,:fillIndex].flatten()[radius[:,:fillIndex].flatten()>0.1]
total_int_list = total_int[:,:fillIndex].flatten()[radius[:,:fillIndex].flatten()>0.1]
neighbours_list = neighbours[:,:fillIndex].flatten()[radius[:,:fillIndex].flatten()>0.1]
volume_list = np.pi*radius[::fillIndex].flatten()[radius[::fillIndex].flatten()>0.01]**2

total_int_list2 = total_int.flatten()[radius.flatten()>0.1]
neighbours_list2 = neighbours.flatten()[radius.flatten()>0.1]

print(np.count_nonzero(np.isnan(total_int_list)))
#age_list = age.flatten()
#total_int_list = total_int.flatten()


#age_list_sort = age_list[np.argsort(age_list)]
#total_int_list_sort = total_int_list[np.argsort(age_list)]

#ages = np.arange(0.005, np.max(age_list), 0.005)
ages = np.arange(0.005, 50.0, 0.005)
print(ages)
total_int_mean_age = np.zeros(len(ages))
total_int_std_age = np.zeros(len(ages))

#age_list = np.array((age_list/0.005), dtype = int)
#age_list = age_list*0.005

for t in range(len(ages)):
    #total_int_mean_age[t] = np.mean((total_int_list[age_list == ages[t]])[~np.isnan(total_int_list[age_list == ages[t]])])
    #total_int_std_age[t] = np.std((total_int_list[age_list == ages[t]])[~np.isnan(total_int_list[age_list == ages[t]])])
    total_int_mean_age[t] = np.mean(total_int_list[(age_list < ages[t] + 0.0025)  & (age_list > ages[t] - 0.0025)])
    total_int_std_age[t] = np.std(total_int_list[(age_list < ages[t] + 0.0025)  & (age_list > ages[t] - 0.0025)])
#for t in range(len(ages)):
#    total_int_mean_age[t] = np.mean(total_int_list[(age_list == ages[t]) &(~np.isnan(total_int_list))])
#    total_int_std_age[t] = np.std(total_int_list[(age_list == ages[t]) &(~np.isnan(total_int_list))])

print(age_list)

neighs = np.arange(0, 9, 1, dtype = int)

total_int_mean_neighs = np.zeros(len(neighs))
total_int_std_neighs = np.zeros(len(neighs))
total_neighs_list = [[] for i in range(len(neighs))]




for t in range(len(neighs)):
    total_int_mean_neighs[t] = np.mean(total_int_list[neighbours_list == neighs[t]])
    total_int_std_neighs[t] = np.std(total_int_list[neighbours_list == neighs[t]])
    total_neighs_list[neighs[t]].extend(total_int_list[neighbours_list == neighs[t]])
#total_int_mean_neighs = np.mean(total_int_list[neighbours_list == neighs])
#total_int_std_neighs = np.std(total_int_list[neighbours_list == neighs])

fig, ax = plt.subplots()
ax.errorbar(ages[::200], total_int_mean_age[::200], yerr = total_int_std_age[::200], capsize=4)
ax.set_xlabel('Age')
ax.set_ylabel('Total Interactions')
plt.savefig(savedir+expName+"/neoagestress.png")

print(ages[::100])
print(total_int_mean_age[::100])

fig, ax = plt.subplots()
ax.errorbar(neighs, total_int_mean_neighs, yerr = total_int_std_neighs, capsize=4)
ax.set_xlabel('Number of neighbours')
ax.set_ylabel('Total Interactions')
plt.savefig(savedir+expName+"/neoneighstress.png")

"""
fig, ax = plt.subplots()
ax.violinplot(total_neighs_list, showmeans=False, showmedians=False,
        showextrema=False)
ax.scatter(neighs, total_int_mean_neighs)
plt.savefig(savedir+expName+"/neoneighstressviolin.png")
"""