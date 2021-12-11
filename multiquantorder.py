import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter1d

plt.style.use('seaborn-bright')
savedir = '/scratch/ws/1/haja565a-workspace2/quant/'

expNames = [  '700g12','700g13', '700g14','700g15', '700g16','700g17']#, ]#'700g15',, '700g18', '700g19']#'700j17', '700j18', '700j19','700g11',
expConRad = {'700j17': 25, '700j18': 20, '700j19': 15, '700g13': 25, '700g14': 20, '700g16': 15, '700g11': 35, '700g12': 30, '700g15': 17.5,\
    '700g17': 12.5, '700g18': 10, '700g19': 7.5}
expLi = {'700j17': 7500, '700j18': 7500, '700j19': 7500, '700g13': 10000, '700g14': 10000, '700g16': 10000, '700g11': 10000, '700g12': 10000, '700g15': 10000,\
    '700g17': 10000, '700g18': 10000, '700g19': 10000}

Dictrotorder = {}
Dicttransorder = {}
Dicttime = {}

for exp in expNames:
    Dictrotorder[exp] = np.load(savedir + exp + "/rotOrder.npy")
    Dicttransorder[exp] = np.load(savedir + exp + "/transOrder.npy")
    Dicttime[exp] = np.load(savedir + exp + "/timearray.npy")

fig, ax = plt.subplots()
for exp in expNames:
    ax.plot(Dicttime[exp], gaussian_filter(Dictrotorder[exp], sigma = 100), label = r"$r_c = %s$"%expConRad[exp] )#$L_i = %s$"%expLi[exp] + ", " + r"
ax.legend()
ax.set_xlabel('time')
ax.set_ylabel(r'$O_R$')
ax.set_xlim(0,200)
ax.set_ylim(-1,1)
ax.grid()
plt.savefig(savedir + 'rotOrdermix.png')