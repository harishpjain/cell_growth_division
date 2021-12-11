import numpy as np
from numpy.lib.npyio import save
import os
import sys
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import rank_filter
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors

savedir = "/scratch/ws//1/haja565a-workspace2/DeformationField/"
expName = sys.argv[1]
time = float(sys.argv[2])
dt=0.005
ind = int(time/dt)
interpolation_steps = 1000
domain_size = 100
x = np.linspace(0, domain_size,interpolation_steps)
xx, yy = np.meshgrid(x,x)

agefield = np.load(savedir + expName + "/age"+ '{:06.3f}'.format(time) + ".npy")
agefield *= 0.005
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.2, 1]
        return np.ma.masked_array(np.interp(value, x, y))

figa, axa = plt.subplots()
#divnorm=mpl.colors.TwoSlopeNorm(vmin=0., vcenter=1., vmax=10000)
midnorm = MidpointNormalize(vmin=0., vcenter=0.005, vmax=40)
Ages = axa.pcolormesh(xx, yy, agefield, norm =midnorm, cmap = mpl.cm.get_cmap('CMRmap_r'))
#Ages = axa.pcolormesh(xx, yy, agefield, norm =colors.logNorm(vmin=1, vmax=30000), cmap = "gist_rainbow")
#axa.set_xlabel("X")
#axa.set_ylabel("Y")
angle = np.linspace(0,2*np.pi,300)
radius = 45
x = np.cos(angle)*radius+50
y = np.sin(angle)*radius+50
axa.plot(x,y, color='k')
axa.axis('square')
axa.set_xticks([])
axa.set_yticks([])
plt.tight_layout()
plt.savefig(savedir + expName + "/ageneonil"+ '{:05d}'.format(ind) + '.png', dpi=250)
figa.colorbar(Ages, ax=axa, label = "Age")#,orientation='horizontal')
plt.tight_layout()
plt.savefig(savedir + expName + "/ageneo"+ '{:05d}'.format(ind) + '.png', dpi=250)