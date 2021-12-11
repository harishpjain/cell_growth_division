import numpy as np
from numpy.lib.npyio import save
import os
import sys
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import rank_filter
import matplotlib.pyplot as plt
from shapely.geometry import LineString
#from shapely.geometry import GeometryCollection
#from shapely.geometry import Point
import math
import matplotlib.ticker as ticker
import matplotlib as mpl
from matplotlib import colors


savedir = "/scratch/ws//1/haja565a-workspace2/DeformationField/"
expName = sys.argv[1]
time = float(sys.argv[2])

delt = 0.005
ind = int(time/delt)

interpolation_steps = 1000
domain_size = 100
x = np.linspace(0, domain_size,interpolation_steps)
xx, yy = np.meshgrid(x,x)

Vxglobf = np.load(savedir + expName + "/Vx"+ '{:06.3f}'.format(time) + ".npy")
Vyglobf = np.load(savedir + expName + "/Vy"+ '{:06.3f}'.format(time) + ".npy")
S0_glob = np.load(savedir + expName + "/S0"+ '{:06.3f}'.format(time) + ".npy")
S1_glob = np.load(savedir + expName + "/S1"+ '{:06.3f}'.format(time) + ".npy")
int_glob = np.load(savedir + expName + "/interaction"+ '{:06.3f}'.format(time) + ".npy")
phi_glob = np.load(savedir + expName + "/phi"+ '{:06.3f}'.format(time) + ".npy")

vorticity = np.gradient(Vyglobf, axis = 1) - np.gradient(Vxglobf, axis = 0) 
vorticity = np.multiply(vorticity, np.sign(phi_glob)) #vorticity

normV = np.sqrt(Vxglobf**2 + Vyglobf**2) #velocity magnitude
Uxglob = np.divide(Vxglobf, normV, out=np.zeros_like(Vxglobf), where=normV!=0) #unit x vel
Uyglob = np.divide(Vyglobf, normV, out=np.zeros_like(Vyglobf), where=normV!=0) #unit y vel

a = np.gradient(S0_glob, axis = 0)
b = np.gradient(S0_glob, axis = 1)
c = np.gradient(S1_glob, axis = 0)
d = np.gradient(S1_glob, axis = 1)
dow = a*d-b*c

norm_S = np.sqrt(S0_glob**2 + S1_glob**2) #magnitude of Deformation?? Is that right?

S0_globn = S0_glob / norm_S
S1_globn = S1_glob / norm_S

omega = np.multiply(np.sign(S1_globn),(np.arctan(S0_glob/np.abs(S1_glob))/ 2.0 + np.pi/4.0))
eig_v0 = np.sin(omega)
eig_v1 = np.cos(omega)


print("stage 2 over") 


#lines of trivial S Tensor
fig, ax = plt.subplots()
c1 = ax.contour(xx,yy,S0_glob,0,colors='red',alpha=1.0)
c2 = ax.contour(xx,yy,S1_glob,0,colors='blue',alpha=1.0)    
plt.savefig(savedir+expName+'/conti'+ '{:06.3f}'.format(time) +'.png', dpi = 200, format='png')
defect_list = []

"""
c1 = c1.allsegs[1]
c2 = c2.allsegs[1]
# Find intersections between zero contours
for m in range(len(c1)):
    for n in range(len(c2)):
        line1 = LineString(c1[m])
        line2 = LineString(c2[n])
        # shapely Geometry Sequence
        if line1.intersects(line2):
            intersects = line1.intersection(line2)
            try:
                for k in intersects:
                    x,y = k.xy
                    if(len(x)==1):
                        defect_list.append(np.array([x,y]))
            except:
                x,y = intersects.xy
                if(len(x)==1):
                    defect_list.append(np.array([x,y]))

print("stage 3 over") 
"""
fig13, ax13 = plt.subplots()
ConVor = ax13.contourf(xx, yy, vorticity, cmap = 'bwr')
VelF = ax13.quiver(xx[::20, ::20], yy[::20, ::20], Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig13.colorbar(ConVor, ax=ax13, label = "vorticity")
ax13.set_title("Vorticity field with \n velocity direction field")
plt.tight_layout()
plt.savefig(savedir + expName + '/vortex_'+ '{:05d}'.format(ind) + '.png', dpi=150, format='png')

fig14, ax14 = plt.subplots()
inhibition_limit = 10000.0
normint_glob = int_glob/inhibition_limit
class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.2, 1]
        return np.ma.masked_array(np.interp(value, x, y))
midnorm = MidpointNormalize(vmin=0, vcenter=0.0001, vmax=2)
##I1 = ax14.contourf(xx, yy, int_glob, cmap = "Reds", alpha = 0.8, norm =colors.Normalize(vmin=1, vmax=16000))#, vmin = 0, vmax=16000)
I1 = ax14.pcolormesh(xx, yy, normint_glob, norm =midnorm, cmap = 'CMRmap_r')
#I1 = ax14.pcolormesh(xx, yy, int_glob, norm =colors.Normalize(vmin=0, vmax=30000), cmap = 'nipy_spectral_r')
#Qvel = ax12[0,1].quiver(xx[::20, ::20], yy[::20, ::20],  Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
#C = fig14.colorbar(I1, ax=ax14, label="Total interactions")#, ticks = [0, 4000, 8000, 12000, 16000])# vmin = 0, vmax = 16000)
#C.ax.set_yticklabels(['0', '4000', '8000', '12000', '16000']) 
#C.mappable.set_clim(vmin =1, vmax = 16000) 
#C.ax.set_ylim(0,16000)
#mpl.set_clim(C, vmin = 0, vmax = 16000)
#C.set_clim(vmin = 0, vmax = 16000)
#ax14.set_title("Total Interactions")
#ax14.set_xlabel('X')
#ax14.set_ylabel('Y')
ax14.axis('square')
angle = np.linspace(0,2*np.pi,300)
radius = 45
x = np.cos(angle)*radius+50
y = np.sin(angle)*radius+50
ax14.plot(x,y, color='k')
ax14.set_xticks([])
ax14.set_yticks([])
plt.savefig(savedir + expName + '/interactionsneonil_'+ '{:05d}'.format(ind) + '.png', dpi=250, format='png')
plt.tight_layout()
C = fig14.colorbar(I1, ax=ax14, label=r'$\frac{T_i}{L}$')#,orientation='horizontal')#, ticks = [0, 4000, 8000, 12000, 16000])# vmin = 0, vmax = 16000)
plt.tight_layout()
plt.savefig(savedir + expName + '/interactionsneo_'+ '{:05d}'.format(ind) + '.png', dpi=250, format='png')


"""
fig12, ax12 = plt.subplots(2,2, figsize = (12,10))
ConVor = ax12[0,0].contourf(xx, yy, vorticity, cmap = 'bwr')
VelF = ax12[0,0].quiver(xx[::20, ::20], yy[::20, ::20], Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig12.colorbar(ConVor, ax=ax12[0,0], label = "vorticity")
ax12[0,0].set_title("Vorticity field with \n velocity direction field")

I1 = ax12[0, 1].contourf(xx, yy, int_glob, cmap = "Reds", alpha = 0.8)
#Qvel = ax12[0,1].quiver(xx[::20, ::20], yy[::20, ::20],  Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig12.colorbar(I1, ax=ax12[0, 1], label="Total interactions")
ax12[0,1].set_title("Total Interactions")

#interaction, Q director, defects
I = ax12[1, 0].contourf(xx, yy, int_glob, cmap = "Reds", alpha = 0.3)
Q = ax12[1, 0].quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Cont2 = ax10[1].contourf(xx, yy, norm_S, cmap = 'bwr', alpha=0.3)
for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(np.mean(dow[indx-5:indx+6, indy-5:indy+6])>0):
        ax12[1, 0].scatter(point[0], point[1], marker="o", alpha=0.5, color='g')
    else:
        ax12[1, 0].scatter(point[0], point[1], marker="^", alpha=0.5, color='k')
fig12.colorbar(I, ax=ax12[1, 0], label="Total interactions")
ax12[1,0].set_title("Nematic director field \n with defects and total interactions")

Cont = ax12[1,1].contourf(xx, yy, vorticity, cmap = 'bwr', alpha=0.6)
Q = ax12[1,1].streamplot(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], color='k', linewidth = 0.1, density = 5, arrowsize=0)
for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(np.mean(dow[indx-5:indx+6, indy-5:indy+6])>0):
        ax12[1,1].scatter(point[0], point[1], marker="o", alpha=0.5, color='g')
    else:
        ax12[1,1].scatter(point[0], point[1], marker="^", alpha=0.5, color='k')
fig12.colorbar(Cont, ax=ax12[1,1], label = "Vorticity")
ax12[1,1].set_title("Nematic director streamfield \n with defects and vorticity field")

plt.tight_layout()
plt.savefig(savedir + expName + '/total_'+ '{:05d}'.format(ind) + '.png', dpi=150, format='png')




fig15, ax15 = plt.subplots()
I = ax15.contourf(xx, yy, int_glob, cmap = "Reds", alpha = 0.3)
Q = ax15.quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Cont2 = ax10[1].contourf(xx, yy, norm_S, cmap = 'bwr', alpha=0.3)
for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(np.mean(dow[indx-5:indx+6, indy-5:indy+6])>0):
        ax15.scatter(point[0], point[1], marker="o", alpha=0.5, color='g')
    else:
        ax15.scatter(point[0], point[1], marker="^", alpha=0.5, color='k')
fig15.colorbar(I, ax=ax15, label="Total interactions")
plt.tight_layout()
plt.savefig(savedir + expName + '/defects_'+ '{:05d}'.format(ind) + '.png', dpi=150, format='png')


"""