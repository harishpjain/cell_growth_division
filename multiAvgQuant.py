import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import sys
import matplotlib as mpl
plt.style.use('seaborn-bright')
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter1d

savedir = '/scratch/ws/1/haja565a-workspace2/quant/'

Simuls =  ['700v11', '700v11b', '700v11c', '700v12', '700v12b', '700v12c', '700v13', '700v13b', '700v13c', '700v14', '700v14b', '700v14c', '700v15', \
    '700v15b', '700v15c', '700h13', '700h13vb', '700h13vc', '700a13b']
Vels = {'v0':0.0, 'v0_25':0.25, 'v0_5':0.5, 'v0_75':0.75, 'v1':1.0, 'v1_25':1.25, 'v1_5':1.5}
VelsList = ['v0', 'v0_25', 'v0_5', 'v0_75', 'v1', 'v1_25', 'v1_5']
VelsExp = {'v0':['700a13b'], 'v0_25':['700v11', '700v11b', '700v11c'], 'v0_5':['700v12', '700v12b', '700v12c'], \
    'v0_75':['700v13', '700v13b', '700v13c'], 'v1':['700h13', '700h13vb', '700h13vc'], 'v1_25':[ '700v14', '700v14b', '700v14c'], \
    'v1_5':['700v15', '700v15b', '700v15c']}

Dictrad = {}
DictT = {}

for key in Simuls:
    allrad = np.load(savedir + key + "/radius.npy")
    Dictrad[key] = np.sum(allrad**2, axis=0)**0.5
    DictT[key] = np.load(savedir + key + "/T.npy")
    print(key)
    print(len(DictT[key]))

timeidx = 15000
timearray = DictT["700v11"]
VelsAvg = {}
for key in VelsList:
    VelsAvg[key] = np.zeros(timeidx)
    for exp in VelsExp[key]:
        VelsAvg[key] += Dictrad[exp][ :15000]
    VelsAvg[key] = VelsAvg[key]/len(VelsExp[key])

my_colors = ['blue', 'green', 'red', 'magenta', 'orange', 'purple', 'grey']
fig, ax = plt.subplots()#figsize = (5,10))
ind = 0
for key in VelsList:
    print(VelsAvg[key].shape)
    print(timearray[:15000].shape)
    ax.plot(timearray[:15000], VelsAvg[key][:15000], color=my_colors[ind],linewidth=0.5,  label = Vels[key])
    ind+=1
ax.legend()
ax.set_ylabel("mean colony radius")
ax.set_xlabel("time")
plt.savefig(savedir + 'velavg2.png', dpi=300)

