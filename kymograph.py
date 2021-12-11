import numpy as np
import matplotlib.pyplot as plt
from os import sys

readdir = '/scratch/ws/1/haja565a-workspace2/quant/'
expName = sys.argv[1]

Vrad = np.load(readdir + expName + "/Vrad.npy")
Vort = np.load(readdir + expName + "/Vort.npy")
dist_cent = np.load(readdir + expName + "/dist_cent.npy")

totaltimes = Vrad.shape[0]
timegap = 1000 #100
time = np.arange(0, totaltimes*timegap, timegap)

#tt, dd = np.meshgrid(time, dist_cent)

def radial_profile(data, center): #averages a quantity radially
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

center = (50, 50)
dist_1d = radial_profile(dist_cent, center) #List of unique dist_from_center
Vrad_arr = np.zeros([totaltimes, len(dist_1d)])
Vort_arr = np.zeros([totaltimes, len(dist_1d)])
for i in range(totaltimes):
    Vrad_arr[i] = radial_profile(Vrad[i], center)
    Vort_arr[i] = radial_profile(Vort[i], center)

tt, dd = np.meshgrid(time, dist_1d)
print(dist_1d)

fig, ax = plt.subplots(1,2)
ax[0].pcolormesh(dd, tt, Vrad_arr.T)
ax[1].pcolormesh(dd, tt, Vort_arr.T)
fig.tight_layout()
plt.savefig(readdir+expName+"/kymo.png")
