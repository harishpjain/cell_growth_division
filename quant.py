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
from scipy.signal import lfilter
from scipy.signal import savgol_filter

import scipy
from scipy.optimize import least_squares
import matplotlib as mpl

def fitfunc1(x, a, b):
    return -(np.exp(-a*x+b)-a)
    #return a * np.exp(b * x) + c
def fitfunc2(x, a, b):
    return -a/(b*x) + b
def fitfunc3(x, a):
    return -a/(x) +1.0
def fitfunc4(x, b):
    return -(np.exp(-x+b)-1.0)
def fit_gauss(x, a, b):
     return ((1)/(np.sqrt(2*np.pi*b*b)))*(np.exp(-0.5*((x-a)**2/(b*b))))


from datetime import datetime
tic1 = datetime.now()

readpos = True
readvel = True
readrad = True
readint = True
readangle = True
readneighbours = True
readgrowth = True
plottrajec = True
plotneighbours = True
plotradius= True
plotquant = True #int, rad, vel vs distance from center
plotorder = True #
multitime_avg = True#find number neighbours averaged over a window
sparsetime_avg = True #find number of neighbours vs time smoothly
plotFill = True # a marker that colony has filled confinement

import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.use('PDF')


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
for cell in range(totalCells):#redundant!!
     Time[cell] = T
print(type(total_int))

X[radius<0.1]  = 0
Y[radius<0.1]  = 0
Vx[radius<0.1]  = 0
Vy[radius<0.1]  = 0

print("Finished Reading Data:" + expName)
volume = np.pi*radius**2
tic2 = datetime.now()
print(tic2-tic1)
#alive_mask = radius>0.01 #true when cells have non zero radius
dist_center = np.sqrt(np.square(X-50)+np.square(Y-50)) #distance of cell center from domain center, center is 50,50 from 0,0 of simulation domain

#numCells = np.zeros(timesteps, dtype=int) #stores number of cells at a given time
colony_volume_1 = np.zeros(timesteps) #sum of volume of all cells

trans_order = np.zeros(timesteps)
rot_order = np.zeros(timesteps)
V_mag = np.sqrt(Vx**2+ Vy**2) #magnitude of velocity  

#Ux = Vx/V_mag #unit x velocity
Ux = np.divide(Vx, V_mag, out=np.zeros_like(Vx), where=V_mag!=0)

#Uy = Vy/V_mag #unit y velocity
Uy = np.divide(Vy, V_mag, out=np.zeros_like(Vy), where=V_mag!=0)

Cx = X-50 #shifted coordinates
Cy = Y-50

angular_disp = np.arctan2(Cy, Cx)+np.pi/2 #angle of point measured from center of domain
ang_posX = np.cos(angular_disp)
ang_posY = np.sin(angular_disp)


numCells = np.count_nonzero(radius>0.01, axis = 0) #change if cells are very small
divEvents = np.gradient(numCells)
divEventsbool = divEvents>0
divIndex = np.where(divEventsbool)[0]

colony_volume_1 = np.sum(volume, axis=0)
fillIndex = np.where(colony_volume_1>=0.99*max(colony_volume_1))[0][0] #index where confinement is filled!!
print(fillIndex)


trans_order = np.sqrt(np.sum(Ux, axis=0)**2+np.sum(Uy, axis=0)**2)/numCells

rot_order = np.sum(Ux*ang_posX+Uy*ang_posY, axis=0)/numCells

colony_radius_1 = np.sqrt(colony_volume_1/np.pi)

nematic_order = np.zeros(timesteps)
aliveCells = [[] for _ in range(timesteps)]
age = np.ones([totalCells, timesteps], dtype = int)*0
#for t in range(timesteps):
#     [aliveCells[t]] = list(np.nonzero(radius[:,t]>0.01))
#     nematic_order[t] = np.sqrt(np.sum(np.cos(2*angle[aliveCells[t],t]))**2+np.sum(np.sin(2*angle[aliveCells[t],t]))**2)/numCells[t]

[aliveCells[0]] = list(np.nonzero(radius[:,0]>0.01))
nematic_order[0] = np.sqrt(np.sum(np.cos(2*angle[aliveCells[0],0]))**2+np.sum(np.sin(2*angle[aliveCells[0],0]))**2)/numCells[0]
for t in range(timesteps-1):
     [aliveCells[t+1]] = list(np.nonzero(radius[:,t+1]>0.01))
     nematic_order[t+1] = np.sqrt(np.sum(np.cos(2*angle[aliveCells[t+1],t+1]))**2+np.sum(np.sin(2*angle[aliveCells[t+1],t+1]))**2)/numCells[t+1]
     for cell in range(totalCells):
          if(radius[cell,t+1]>0.1):
               if((radius[cell, t+1]**2)<(0.5*radius[cell, t]**2)):
                    age[cell, t+1] = 0
               else:
                    age[cell, t+1] = age[cell, t] +1

np.save(savedir+expName+"/age.npy", age)
vel_avg = np.sum(V_mag, axis = 0)/numCells

age = age*dt
# Finding quanties for cells that are alive, and averaging then over the averaging_period
averaging_period = int(timesteps/40)

select_interactions = []
select_dist_center = []
select_velocity = []
select_velx = []
select_vely = []
select_radius = []
select_ang_disp = []
time = []


for t in range(int(timesteps/averaging_period)):
     numAliveCells = numCells[t*averaging_period] #number of cells alive at start of averaging period
     for cell in aliveCells[t*averaging_period]:
          time.append(t*averaging_period*dt)
     select_dist_center.extend(np.mean(dist_center[aliveCells[t*averaging_period], t*averaging_period: (t+1)*averaging_period], axis=1))
     select_interactions.extend(np.mean(total_int[aliveCells[t*averaging_period], t*averaging_period: (t+1)*averaging_period], axis=1))
     select_velocity.extend(np.mean(V_mag[aliveCells[t*averaging_period], t*averaging_period: (t+1)*averaging_period], axis=1))
     select_radius.extend(np.mean(radius[aliveCells[t*averaging_period], t*averaging_period: (t+1)*averaging_period], axis=1))
     select_velx.extend(np.mean(Vx[aliveCells[t*averaging_period], t*averaging_period: (t+1)*averaging_period], axis=1))
     select_vely.extend(np.mean(Vy[aliveCells[t*averaging_period], t*averaging_period: (t+1)*averaging_period], axis=1))
     select_ang_disp.extend(np.mean(Vy[aliveCells[t*averaging_period], t*averaging_period: (t+1)*averaging_period], axis=1))

#interpolating data for Kymographs
select_velx = np.array(select_velx)
select_vely = np.array(select_vely)
select_ang_disp = np.array(select_ang_disp)
vel_rad = select_velx*np.cos(select_ang_disp) + select_vely*np.sin(select_ang_disp)
vel_ortho = select_velx*np.sin(select_ang_disp) - select_vely*np.cos(select_ang_disp)
timearray = T[::averaging_period]
distarray = np.linspace(0, max(select_dist_center), 100)
tt, dd = np.meshgrid(timearray, distarray)
vel_rad_interp = griddata((time, select_dist_center), vel_rad, (tt,dd), method='nearest')
vel_ortho_interp = griddata((time, select_dist_center), vel_ortho, (tt,dd), method='nearest')


saveOrder = True
if(saveOrder):
     np.save(savedir+expName+'/rotOrder.npy', rot_order)
     np.save(savedir+expName+'/transOrder.npy', trans_order)
     np.save(savedir+expName+'/timearray.npy', T)


tic3 = datetime.now()
print(tic3-tic2)

#****************************************************#

print("plotting ...")

if(plotradius == True):

     #bar plot of number of cells vs time
     fig2, ax2 = plt.subplots()
     ax2.scatter(T, numCells, marker='.')
     ax2.set_xlabel('time')
     ax2.set_ylabel('Cell Number')
     ax2.set_title('Number of Cells')
     plt.savefig(savedir+expName+'/numCells.png', dpi=150, format='png')
     plt.close("all")

     #colony radius 1 vs time 
     fig3, ax3 = plt.subplots()
     ax3.plot(T, colony_radius_1)
     ax3.set_ylabel('colony radius')
     ax3.set_xlabel('time')
     ax3.set_title('Colony Radius')
     ax3.grid()
     plt.savefig(savedir+expName+'/colonyRadius.png', dpi=150, format='png')
     plt.close("all")

          #colony radius 1 vs time 
     fig3a, ax3a = plt.subplots()
     n = 200
     b = [1.0/n]*n
     a = 1
     #smooth_bound_vel = np.gradient(lfilter(b, a, (colony_radius_1/dt)))
     #smooth_bound_vel = savgol_filter(np.gradient(colony_radius_1/dt, 101, 2)
     #smooth_bound_vel = np.gradient(savgol_filter(colony_radius_1/dt, 101, 2))
     #kernel_size = 10
     #kernel = np.ones(kernel_size) / kernel_size
     #convolved_radius = np.convolve(colony_radius_1, kernel, mode='same')
     colony_radius_gap = colony_radius_1[:7500:30]
     avg_smooth = np.zeros(len(colony_radius_gap))
     avg_smooth[0:3] = colony_radius_gap[0:3]
     avg_smooth[-3:] = colony_radius_gap[-3:]
     avg_smooth[3:-3] = (colony_radius_gap[0:-6]+colony_radius_gap[1:-5]+colony_radius_gap[2:-4]+colony_radius_gap[3:-3]+colony_radius_gap[4:-2]+colony_radius_gap[5:-1]+colony_radius_gap[6:])/7
     ax3a.plot(T[:7500:30], np.gradient(avg_smooth/dt))

     #ax3a.plot(T, smooth_bound_vel)
     #ax3a.plot(T, gaussian_filter(np.gradient(colony_radius_1), sigma=100)/dt)
     #ax3a.plot(T, np.gradient(gaussian_filter1d(colony_radius_1, sigma=3))/dt)
     #ax3a.plot(T, gaussian_filter(np.gradient(gaussian_filter((colony_radius_1),sigma=10)/dt), sigma = 10))
     ax3a.set_ylabel('colony radius rate')
     ax3a.set_xlabel('time')
     ax3a.set_title('Colony Radius rate')
     ax3a.grid()
     plt.savefig(savedir+expName+'/colonyRadiusRate.png', dpi=150, format='png')
     plt.close("all")

     #colony volume 1 vs time
     fig4, ax4 = plt.subplots()
     ax4.plot(T, colony_volume_1)
     ax4.set_xlabel('time')
     ax4.set_ylabel('colony volume')
     ax4.set_title('Colony Volume')
     ax4.grid()
     plt.savefig(savedir+expName+'/colonyVolume.png', dpi=150, format='png')
     plt.close("all")

          #colony volume 1 vs time
     fig4, ax4 = plt.subplots()
     ax4.plot(T, np.log(colony_volume_1))
     ax4.set_xlabel('time')
     ax4.set_ylabel(r'$\log{(V_c)}$')
     ax4.set_xticks(T[divIndex], minor=True)
     ax4.tick_params(which='minor', length=4, color='r')
     #ax4.set_title('Colony Volume')
     ax4.grid()
     ax4.set_xlim(0,18)
     plt.savefig(savedir+expName+'/colonyVolumelog.png', dpi=150, format='png')
     plt.close("all")

     #volume stack plots vs time
     fig5, ax5 = plt.subplots()
     ax5.stackplot(T, volume, baseline = 'zero')
     ax5.set_xlabel('time')
     ax5.set_ylabel('volume')
     ax5.set_title('Volume of Cells')
     plt.savefig(savedir+expName+'/volumeStack.png', dpi=150, format='png')
     plt.close("all")

     #radius stack plots vs time
     fig6, ax6 = plt.subplots()
     ax6.stackplot(T, radius, baseline = 'zero')
     ax6.set_xlabel('time')
     ax6.set_ylabel('radius')
     ax6.set_title('Radius of Cells')
     plt.savefig(savedir+expName+'/radiusStack.png', dpi=150, format='png')
     plt.close("all")

if(plotorder == True):
     #translational order parameter vs time
     fig8, ax8 = plt.subplots()
     #ax8.plot(T[50::], trans_order[50::], 'r')
     ax8.plot(T[50::], gaussian_filter(trans_order[50::], sigma=5), 'b')
     ax8.set_xlabel('time')
     ax8.set_ylabel(r'$o_T$')
     ax8.plot(T[fillIndex]*np.ones(100), np.linspace(-0.2,1.2, 100), c='r', linewidth = 3.0, alpha = 0.5)
     ax8.set_xticks(T[divIndex], minor=True)
     ax8.tick_params(which='minor', length=4, color='r')
     ax8.grid()
     ax8.set_ylim([-0.1, 1.1])
     #ax8.set_title('Translational Order Parameter(gaussian filter of size 3)')
     plt.savefig(savedir+expName+'/transOrder.png', dpi=150, format='png')
     plt.close("all")

     #rotational order parameter vs time
     fig9, ax9 = plt.subplots()
     #ax9.plot(T[50::], rot_order[50::], 'r')
     ax9.plot(T[50::], gaussian_filter(rot_order[50::], sigma=5), 'b')
     ax9.set_xlabel('time')
     ax9.set_ylabel(r'$o_R$')
     ax9.set_ylim([-1.0, 1.0])
     ax9.set_xticks(T[divIndex], minor=True)
     ax9.tick_params(which='minor', length=4, color='r')
     ax9.plot(T[fillIndex]*np.ones(100), np.linspace(-1,1, 100), c='r', linewidth = 3.0, alpha = 0.5)
     #ax9.set_title('Rotational Order Parameter(gaussian filter of size 3)')
     ax9.grid()
     plt.savefig(savedir+expName+'/rotOrder.png', dpi=150, format='png')
     plt.close("all")

     #nematic order parameter vs time
     fig12, ax12 = plt.subplots()
     ax12.plot(T[50::], gaussian_filter(nematic_order[50::], sigma=3))
     ax12.set_xlabel('time')
     ax12.set_ylabel('Nematic Order Parameter')
     ax12.set_xticks(T[divIndex], minor=True)
     ax12.grid()
     ax12.tick_params(which='minor', length=4, color='r')
     ax12.plot(T[fillIndex]*np.ones(100), np.linspace(-0.2,1.2, 100), c='r', linewidth = 3.0, alpha = 0.5)
     ax12.set_title('Nematic Order(Gaussian filter of size 3)')
     plt.savefig(savedir+expName+'/nematicOrder.png', dpi=150, format='png')
     plt.close("all")

if(plottrajec== True):
     #trajectory of 4 cells
     fig10, ax10 = plt.subplots()
     startTime = 2000
     k = np.random.choice(aliveCells[startTime], size=4, replace=False)
     gap = 1000
     num = len(X[k[0], startTime::gap])
     for i in range(num):
          ax10.scatter(X[k[0], startTime::gap][i], Y[k[0], startTime::gap][i], color = 'r', alpha = i*1.0/num, marker=".")
          ax10.scatter(X[k[1], startTime::gap][i], Y[k[1], startTime::gap][i], color = 'g', alpha = i*1.0/num, marker=".")
          ax10.scatter(X[k[2], startTime::gap][i], Y[k[2], startTime::gap][i], color = 'b', alpha = i*1.0/num, marker=".")
          ax10.scatter(X[k[3], startTime::gap][i], Y[k[3], startTime::gap][i], color = 'k', alpha = i*1.0/num, marker=".")
     ax10.set_xlabel('X')
     ax10.axis('square')
     ax10.set_ylabel('Y')
     ax10.set_xlim([0, 100])
     ax10.set_ylim([0, 100])
     ax10.set_title("Trajectories")
     plt.savefig(savedir+expName+'/trajec.png', dpi=150, format='png')
     plt.close("all")

     fig10b, ax10b = plt.subplots()
     ax10b.plot(X[k[0], startTime::], Y[k[0], startTime::], 'r')
     ax10b.plot(X[k[1], startTime::], Y[k[1], startTime::], 'g')
     ax10b.plot(X[k[2], startTime::], Y[k[2], startTime::], 'b')
     ax10b.plot(X[k[3], startTime::], Y[k[3], startTime::], 'k')
     plt.axis('square');
     ax10b.set_xlabel('X')
     ax10b.set_ylabel('Y')
     ax10b.set_xlim([0, 100])
     ax10b.set_ylim([0, 100])
     ax10b.set_title("Trajectories")
     plt.savefig(savedir+expName+'/trajec2.png', dpi=150, format='png')
     plt.close("all")

if(plotquant==True):
     #average velocity of all cells
     fig11, ax11 = plt.subplots()
     #ax11.plot(T, gaussian_filter(vel_avg, sigma=1), 'r')
     #ax11.plot(T, gaussian_filter(vel_avg, sigma=2), 'g')
     ax11.plot(T, gaussian_filter(vel_avg, sigma=3), 'b')
     #ax11.plot(T, vel_avg)
     ax11.set_xlabel(r'$t$')
     ax11.set_ylabel(r'$|v|_{avg}$')
     ax11.set_xticks(T[divIndex], minor=True)
     ax11.tick_params(which='minor', length=4, color='r')
     ax11.grid()
     #ax11.set_title('Average velocity')
     ax11.plot(T[fillIndex]*np.ones(100), np.linspace(np.min(vel_avg),np.max(vel_avg), 100), c='r', linewidth = 3.0, alpha = 0.5)
     plt.savefig(savedir+expName + '/avgVel.png', dpi=150, format='png')
     plt.close("all")

     ##average growth rate over a time window say 100 steps and scatter 
     # plot it with respect to time and distance from center
     fig13, ax13 = plt.subplots()
     A13 = ax13.scatter(time, select_dist_center, c=select_interactions, cmap='plasma', marker='_')
     ax13.set_xlabel('time')
     ax13.set_ylabel('distance from center')
     ax13.set_title('average total interactions')
     fig13.colorbar(A13, label='total interactions')
     ax13.set_xticks(T[divIndex], minor=True)
     ax13.tick_params(which='minor', length=4, color='r')
     ax13.set_title('Total interactions')
     ax13.plot(T[fillIndex]*np.ones(100), np.linspace(np.min(select_dist_center),np.max(select_dist_center), 100), c='r', linewidth = 3.0, alpha = 0.5)
     plt.savefig(savedir+expName+'/int_dist.png', dpi=150, format='png')
     plt.close("all")

     fig14, ax14 = plt.subplots()
     A14 = ax14.scatter(time, select_dist_center, c=select_velocity, cmap='plasma', marker='_')
     ax14.set_xlabel('time')
     ax14.set_ylabel('distance from center')
     ax14.set_title('average velocity')
     ax14.set_xticks(T[divIndex], minor=True)
     ax14.tick_params(which='minor', length=4, color='r')
     fig14.colorbar(A14, label=r'$v_{avg}$')
     ax14.plot(T[fillIndex]*np.ones(100), np.linspace(np.min(select_dist_center),np.max(select_dist_center), 100), c='r', linewidth = 3.0, alpha = 0.5)
     ax14.set_title('Velocity')
     plt.savefig(savedir+expName+'/vel_dist.png', dpi=150, format='png')
     plt.close("all")

     fig15, ax15 = plt.subplots()
     A15 = ax15.scatter(time, select_dist_center, c=select_radius, cmap='plasma', marker='_')
     ax15.set_xlabel('time')
     ax15.set_ylabel('distance from center')
     ax15.set_title('average radius')
     ax15.set_xticks(T[divIndex], minor=True)
     ax15.tick_params(which='minor', length=4, color='r')
     ax15.plot(T[fillIndex]*np.ones(100), np.linspace(np.min(select_dist_center),np.max(select_dist_center), 100), c='r', linewidth = 3.0, alpha = 0.5)
     fig15.colorbar(A15, label=r'$r_{avg}$')
     ax15.set_title('Radius')
     plt.savefig(savedir+expName+'/rad_dist.png', dpi=150, format='png')
     plt.close("all")

     fig13l, ax13l = plt.subplots()
     A13l = ax13l.scatter(time, select_dist_center, c=np.log(select_interactions), cmap='plasma', marker='_')
     ax13l.set_xlabel('time')
     ax13l.set_ylabel('distance from center')
     ax13l.set_title('log average total interactions')
     ax13l.plot(T[fillIndex]*np.ones(100), np.linspace(np.min(select_dist_center),np.max(select_dist_center), 100), c='r', linewidth = 3.0, alpha = 0.5)
     fig13l.colorbar(A13l, label='log(total interactions)')
     ax13l.set_title('Logarithm of Total interactions')
     plt.savefig(savedir+expName+'/int_distlog.png', dpi=150, format='png')
     plt.close("all")

     fig14l, ax14l = plt.subplots()
     A14l = ax14l.scatter(time, select_dist_center, c=np.log(select_velocity), cmap='plasma', marker='_')
     ax14l.set_xlabel('time')
     ax14l.set_ylabel('distance from center')
     ax14l.set_title('log average velocity')
     ax14l.set_xticks(T[divIndex], minor=True)
     ax14l.tick_params(which='minor', length=4, color='r')
     ax14l.plot(T[fillIndex]*np.ones(100), np.linspace(np.min(select_dist_center),np.max(select_dist_center), 100), c='r', linewidth = 3.0, alpha = 0.5)
     fig14l.colorbar(A14l, label=r'$\log(v_{avg})$')
     ax14l.set_title('Logarithm of Velocity')
     plt.savefig(savedir+expName+'/vel_distlog.png', dpi=150, format='png')
     plt.close("all")

     fig15l, ax15l = plt.subplots()
     A15l = ax15l.scatter(time, select_dist_center, c=np.log(select_radius), cmap='plasma', marker='_')
     ax15l.set_xlabel('time')
     ax15l.set_ylabel('distance from center')
     ax15l.set_title('log average radius')
     fig15l.colorbar(A15l, label=r'$\log(r_{avg})$')
     ax15l.plot(T[fillIndex]*np.ones(100), np.linspace(np.min(select_dist_center),np.max(select_dist_center), 100), c='r', linewidth = 3.0, alpha = 0.5)
     ax15l.set_title('Logarithm of radius')
     ax15l.set_xticks(T[divIndex], minor=True)
     ax15l.tick_params(which='minor', length=4, color='r')
     plt.savefig(savedir+expName+'/rad_distlog.png', dpi=150, format='png')
     plt.close("all")



if(plotneighbours==True):
     if(sparsetime_avg): #for a few disparate timesteps
          #number of neighbours
          fig16, ax16 = plt.subplots()
          timen = [3000, 4000, 5000, 6000, 7000]
          #timen = [10000, 12000, 14000, 16000, 18000]
          neighs = [1, 2, 3, 4, 5, 6, 7, 8, 9]
          neighbours_num = np.zeros([len(timen), len(neighs)], dtype=int)
          prob = np.zeros([len(timen), len(neighs)])

          for t in range(len(timen)):
               for n in range(len(neighs)):
                    neighbours_num[t, n]= np.count_nonzero(neighbours[aliveCells[timen[t]], timen[t]]==neighs[n])
               for n in range(len(neighs)):
                    prob[t,n] = float(neighbours_num[t, n])/np.sum(neighbours_num[t,:])
               print(np.sum(neighbours_num[t,:]))
          xbaraxis = np.arange(1, len(neighs)+1, 1)
          neigh_sum = np.sum(neighbours_num, axis=1)*1.0
          width = 0.8/len(timen)
          mycolors = ['r', 'g', 'b', 'k', 'm', 'c', 'y', 'darkorange', 'lime', 'slategray']
          for i in range(len(timen)):
               ax16.bar(xbaraxis+((i-(len(timen)/2))*width), prob[i, :], width, label = "t = " + str(timen[i]*dt), color = mycolors[i])
          ax16.set_xlabel("Number of Neighbours")
          ax16.set_ylabel("Probability")
          ax16.legend(fontsize=12)
          ax16.set_xticks(xbaraxis, neighs)
          ax16.xaxis.set_major_locator(MaxNLocator(integer=True))
          #ax16.yaxis.set_major_locator(MaxNLocator(integer=True))
          plt.savefig(savedir+expName+'/neighbours.png', dpi=150, format='png')
          plt.close("all")


          #number of neighbours inside confinement
          fig16b, ax16b = plt.subplots()
          #timen = [4000, 7000, 10000, 13000, 16000]
          #timen = [10000, 12000, 14000, 16000, 18000]
          neighs = [1, 2, 3, 4, 5, 6, 7, 8, 9]
          neighbours_num = np.zeros([len(timen), len(neighs)], dtype=int)
          prob = np.zeros([len(timen), len(neighs)])
          
          for t in range(len(timen)):
               for n in range(len(neighs)):
                    neighbours_num[t, n]= np.count_nonzero((neighbours[aliveCells[timen[t]], timen[t]]==neighs[n]) & (confine_int[aliveCells[timen[t]], timen[t]]<=0.1))
                    
               for n in range(len(neighs)):
                    prob[t,n] = np.divide(float(neighbours_num[t, n]),np.sum(neighbours_num[t,:]))    
          xbaraxis = np.arange(1, len(neighs)+1, 1)

          width = 0.8/len(timen)
          for i in range(len(timen)):
               ax16b.bar(xbaraxis+((i-(len(timen)/2))*width), prob[i, :], width, label = "t = " + str(timen[i]*dt), color = mycolors[i])
          ax16b.set_xlabel("Number of Neighbours")
          ax16b.set_ylabel("Probability")
          ax16b.legend(fontsize=12)
          ax16b.set_xticks(xbaraxis, neighs)
          ax16b.xaxis.set_major_locator(MaxNLocator(integer=True))
          #ax16b.yaxis.set_major_locator(MaxNLocator(integer=True))
          plt.savefig(savedir+expName+'/neighbours2.png', dpi=150, format='png')
          plt.close("all")

          fig16p, ax16p = plt.subplots()
          for i in range(len(timen)):
               popt1, pcov1 = curve_fit(fit_gauss, neighs,  prob[i, :])
               ax16p.plot(np.linspace(1,9,100), fit_gauss(np.linspace(1,9,100), *popt1), label = "t = " + str(timen[i]*dt))
               ax16p.scatter(neighs, prob[i, :], label = "t = " + str(timen[i]*dt))
          ax16p.set_xlabel("Number of Neighbours")
          ax16p.set_ylabel("Probability")
          ax16p.legend(fontsize=12)
          plt.savefig(savedir+expName+'/neighbours3.png', dpi=300, format='png')
          plt.close("all")

          #fit gamma functions
          poptsarr = np.zeros(len(timen))
          pconvarr = np.zeros(len(timen))

          #prob1 = prob[5,:]
          nei = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])

          """
          fig, ax = plt.subplots()
          for i in range(len(timen)):
               poptsarr, pconvarr = curve_fit(fitfunc1, nei, np.array(prob[i, :], dtype=float))
               ax.plot(np.linspace(0,10), fitfunc1(np.linspace(0,10), *poptsarr), label = str(timen[i]))
          #popt, pconv = curve_fit(fitfunc1, nei, prob1)
          #ax.plot(np.linspace(0,10), fitfunc1(np.linspace(0,10), *popt))
          ax.legend()
          ax.grid()
          plt.savefig(savedir+expName+'/neighboursdist.png', dpi=200, format='png')
          """

     #multitime_average
     if(multitime_avg):
          neighs = [1, 2, 3, 4, 5, 6, 7, 8]
          colors = ['r', 'b', 'g', 'y', 'orange', 'purple', 'k', 'm']
          neighbours_num = np.zeros([timesteps, len(neighs)], dtype=float)
          prob = np.zeros([timesteps, len(neighs)])
          for j in neighs:
               neighbours_num[:, j-1] = np.count_nonzero(neighbours==j, axis = 0)
          neighbours_sum = np.sum(neighbours_num, axis = 1)#, keepdims=True)
          prob = np.divide(neighbours_num.T, neighbours_sum).T
          #prob = gaussian_filter(prob, sigma = 10)

          
          fig17, ax17 = plt.subplots()
          for j in neighs:
               ax17.plot(T, prob[:,j-1], label = str(j), c = colors[j-1])
          ax17.legend()
          ax17.set_xlabel("time")
          ax17.set_ylabel("probability")
          ax17.set_title("Probabiliity")
          ax17.set_xticks(T[divIndex], minor=True)
          ax17.tick_params(which='minor', length=4, color='r')
          ax17.plot(T[fillIndex]*np.ones(100), np.linspace(0,1, 100), c='r', linewidth = 3.0, alpha = 0.5)
          plt.savefig(savedir+expName+"/neighbours_smooth1.png", dpi=150, format="png")

          #smoothened
          prob = gaussian_filter1d(prob, sigma = 10, axis = 0)
          fig17b, ax17b = plt.subplots()
          for j in neighs:
               ax17b.plot(T, prob[:,j-1], label = str(j), c = colors[j-1])
          ax17b.legend()
          ax17b.set_xlabel("time")
          ax17b.set_ylabel("probability")
          ax17b.set_title("Probabiliity with filter")
          ax17b.set_xticks(T[divIndex], minor=True)
          ax17b.tick_params(which='minor', length=4, color='r')
          ax17b.plot(T[fillIndex]*np.ones(100), np.linspace(0,1, 100), c='r', linewidth = 3.0, alpha = 0.5)
          plt.savefig(savedir+expName+"/neighbours_smooth1b.png", dpi=150, format="png")
          
          #finding average number of neighbours
          #neigh_avg = np.zeros(timesteps)
          neigh_avg_all = np.sum(prob*np.array(neighs), axis=1)
          np.save(savedir+expName+"/neighavgall.npy", neigh_avg_all)


          #ignoring cells that touch confinmenet
          neighbours_num = np.zeros([timesteps, len(neighs)], dtype=float)
          prob = np.zeros([timesteps, len(neighs)])
          for j in neighs:
               neighbours_num[:, j-1] = np.count_nonzero((neighbours==j) & (confine_int < 0.1), axis = 0)
          neighbours_sum = np.sum(neighbours_num, axis = 1)#, keepdims=True)
          prob = np.divide(neighbours_num.T, neighbours_sum).T

          fig17c, ax17c = plt.subplots()
          for j in neighs:
               ax17c.plot(T, prob[:,j-1], label = str(j), c = colors[j-1])
          ax17c.legend()
          ax17c.set_xlabel("time")
          ax17c.set_ylabel("probability")
          ax17c.set_title("Probabiliity")
          ax17c.set_xticks(T[divIndex], minor=True)
          ax17c.tick_params(which='minor', length=4, color='r')
          ax17c.plot(T[fillIndex]*np.ones(100), np.linspace(0,1, 100), c='r', linewidth = 3.0, alpha = 0.5)
          plt.savefig(savedir+expName+"/neighbours_smooth2.png", dpi=150, format="png")

          prob = gaussian_filter1d(prob, sigma = 10, axis = 0)
          fig17d, ax17d = plt.subplots()
          for j in neighs:
               ax17d.plot(T, prob[:,j-1], label = str(j), c = colors[j-1])
          ax17d.legend()
          ax17d.set_xlabel("time")
          ax17d.set_ylabel("probability")
          ax17d.set_title("Probabiliity with filter")
          ax17d.set_xticks(T[divIndex], minor=True)
          ax17d.tick_params(which='minor', length=4, color='r')
          ax17d.plot(T[fillIndex]*np.ones(100), np.linspace(0,1, 100), c='r', linewidth = 3.0, alpha = 0.5)
          plt.savefig(savedir+expName+"/neighbours_smooth2b.png", dpi=150, format="png")

          neigh_avg_inner = np.sum(prob*np.array(neighs), axis=1)
          np.save(savedir+expName+"/neighavgin.npy", neigh_avg_inner)
          fig, ax = plt.subplots()
          ax.plot(T, neigh_avg_all, label = "All cells")
          ax.plot(T, neigh_avg_inner, label = "Interior cells")
          ax.grid()
          ax.legend()
          ax.set_xlabel('time')
          ax.set_ylabel(r'$<neighbours>$')
          plt.savefig(savedir+expName+"/neighavg.png", dpi=150, format="png")

          """
          avgRadius = np.sum(radius, axis=0)/numCells
          rad_neigh = np.zeros([timesteps, len(neighs)], dtype=float)
          for j in neighs:
               rad_neigh[:, j-1] = np.sum(radius[(neighbours==j) & (confine_int < 0.1)], axis = 0)
          rad_neigh = rad_neigh/neighbours_num
          An_by_A = np.zeros([timesteps, len(neighs)], dtype=float)
          for j in neighs:
               An_by_A[:,j-1] = (rad_neigh[:,j-1]/avgRadius)
          fig, ax = plt.subplots()
          ax.plot(T[10000:], An_by_A[10000:])
          plt.savefig(savedir+expName+"/An_by_A.png", dpi=150, format="png")
          """

fillIndex = 7000
print("Fill Index is set to")
print(fillIndex)

fig20, ax20 = plt.subplots()
AS = ax20.hist2d(age.flatten()[age.flatten()>0], total_int.flatten()[age.flatten()>0], bins = 20, norm=LogNorm(), cmap=plt.cm.jet)
ax20.set_xlim(left = 0)
ax20.set_xlabel("Age")
ax20.set_ylabel("Total Interactions")
fig20.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stress.png", dpi=150, format="png")

fig20b, ax20b = plt.subplots()
AS = ax20b.hist2d(age.flatten()[age.flatten()>0], total_int.flatten()[age.flatten()>0], bins = 20, cmap=plt.cm.jet)
ax20b.set_xlim(left = 0)
ax20b.set_xlabel("Age")
ax20b.set_ylabel("Total Interactions")
fig20b.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stressnl.png", dpi=150, format="png")
fig21, ax21 = plt.subplots()

age_xx = age[:,:fillIndex].flatten()[age[:,:fillIndex].flatten()>0]
int_yy = np.abs(total_int[:,:fillIndex].flatten()[age[:,:fillIndex].flatten()>0])
neigh_xx = neighbours[:,:fillIndex].flatten()[age[:,:fillIndex].flatten()>0]
inhibition_limit = 10000.0
normint_yy = int_yy/inhibition_limit


fig20c, ax20c = plt.subplots()
AS = ax20c.hist2d(age_xx, normint_yy, bins = 31, cmin = 250, vmin = 1, vmax =3500, range = [[0, 25], [0, 1.5000]], norm=LogNorm(), cmap=plt.cm.jet)
ax20c.set_xlim(left = 0, right=25)
ax20c.set_ylim(0, 1.5000)
ax20c.set_xlabel("Age")
ax20c.set_ylabel(r"$\frac{T_i}{L}$")
fig20c.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stress_fill.png", dpi=150, format="png")


fig20d, ax20d = plt.subplots()
AS = ax20d.hist2d(age_xx, normint_yy, bins = 31, cmin = 250, vmin = 1, vmax = 3500, range = [[0, 25], [0, 1.5000]],cmap='CMRmap_r')# cmap=plt.cm.jet)
ax20d.set_xlim(left = 0, right=25)
ax20d.set_ylim(0, 1.5000)
ax20d.set_xlabel("Age")
ax20d.set_ylabel(r"$\frac{T_i}{L}$")
fig20d.colorbar(AS[3], label = "Count", ticks=[0, 500, 1000, 1500, 2000, 2500, 3000, 3500])
plt.savefig(savedir+expName+"/age_stressnl_fill.png", dpi=150, format="png")

'''
#return a-(np.exp(-a*x+b))
popt1, pcov1 = curve_fit(fitfunc1, age_xx[age_xx>1.0], normint_yy[age_xx>1.0])
print(popt1)
fig20dx, ax20dx = plt.subplots()
AS = ax20dx.hist2d(age_xx, normint_yy, bins = 31, cmin = 250, vmin = 1, vmax = 3500, range = [[0, 25], [0, 1.5000]], cmap=plt.cm.jet)
ax20dx.plot(np.linspace(0,25, 100), fitfunc1(np.linspace(0,25,100), *popt1), color='k', label = r'$a-e^{-ax+b}$' + ", a = " + '%.3f' %popt1[0] + ", b = " + '%.3f' %popt1[1])
ax20dx.set_xlim(left = 0, right=25)
ax20dx.set_ylim(0, 1.5000)
ax20dx.set_xlabel("Age")
ax20dx.set_ylabel("$T_i/L_i$")
ax20dx.legend()
fig20dx.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stressnl_fill_fit1.png", dpi=150, format="png")

#return -a/(b*x) + b used for fit
popt1, pcov1 = curve_fit(fitfunc2, age_xx[age_xx>1.0], normint_yy[age_xx>1.0])
fig20dx, ax20dx = plt.subplots()
AS = ax20dx.hist2d(age_xx, normint_yy, bins = 31, cmin = 250, vmin = 1, vmax = 3500, range = [[0, 25], [0, 1.5000]], cmap=plt.cm.jet)
ax20dx.plot(np.linspace(0,25, 100), fitfunc2(np.linspace(0,25,100), *popt1), color='k', label = r'$b-\frac{a}{bx}$' + ", a = " + '%.3f' %popt1[0] + ", b = " + '%.3f' %popt1[1])
ax20dx.set_xlim(left = 0, right=25)
ax20dx.set_ylim(0, 1.5000)
ax20dx.set_xlabel("Age")
ax20dx.set_ylabel("$T_i/L_i$")
ax20dx.legend()
fig20dx.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stressnl_fill_fit2.png", dpi=150, format="png")

#return -a/(x)+1.0  used for fit
popt1, pcov1 = curve_fit(fitfunc3, age_xx[age_xx>1.0], normint_yy[age_xx>1.0])
fig20dx, ax20dx = plt.subplots()
AS = ax20dx.hist2d(age_xx, normint_yy, bins = 31, cmin = 250, vmin = 1, vmax = 3500, range = [[0, 25], [0, 1.5000]], cmap=plt.cm.jet)
ax20dx.plot(np.linspace(0,25, 100), fitfunc3(np.linspace(0,25,100), *popt1), color='k', label = r'$1-\frac{a}{x}$' + ", a = " + '%.3f' %popt1[0])
ax20dx.set_xlim(left = 0, right=25)
ax20dx.set_ylim(0, 1.5000)
ax20dx.set_xlabel("Age")
ax20dx.set_ylabel("$T_i/L_i$")
ax20dx.legend()
fig20dx.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stressnl_fill_fit3.png", dpi=150, format="png")

#return -(np.exp(-x+b)-1.0)
popt1, pcov1 = curve_fit(fitfunc4, age_xx[age_xx>1.0], normint_yy[age_xx>1.0])
fig20dx, ax20dx = plt.subplots()
AS = ax20dx.hist2d(age_xx, normint_yy, bins = 31, cmin = 250, vmin = 1, vmax = 3500, range = [[0, 25], [0, 1.5000]], cmap=plt.cm.jet)
ax20dx.plot(np.linspace(0,25, 100), fitfunc4(np.linspace(0,25,100), *popt1), color='k', label =r'$1-e^{-x+b}$' + ", b = " + '%.3f' %popt1[0])
ax20dx.set_xlim(left = 0, right=25)
ax20dx.set_ylim(0, 1.5000)
ax20dx.set_xlabel("Age")
ax20dx.set_ylabel("$T_i/L_i$")
ax20dx.legend()
fig20dx.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stressnl_fill_fit4.png", dpi=150, format="png")
'''

popt1, pcov1 = curve_fit(fitfunc1, age_xx[age_xx>1.0], normint_yy[age_xx>1.0])
fig20dx, ax20dx = plt.subplots()
A20dx = ax20dx.plot(np.linspace(0,25, 100), fitfunc1(np.linspace(0,25,100), *popt1))
plt.savefig(savedir+expName+"/age_stress_fit_exp.png", dpi=150, format="png")

popt1, pcov1 = curve_fit(fitfunc2, age_xx[age_xx>1.0], int_yy[age_xx>1.0])
fig20dx, ax20dx = plt.subplots()
A20dx = ax20dx.plot(np.linspace(0,25, 100), fitfunc2(np.linspace(0,25,100), *popt1))
plt.savefig(savedir+expName+"/age_stress_fit_rec.png", dpi=150, format="png")


p=500
#popt1, pcov1 = curve_fit(fitfunc1, age[:,::p].flatten()[age[:,::p].flatten()>0], total_int[:,::p].flatten()[age[:,::p].flatten()>0])

fig20e, ax20e = plt.subplots()
ASs = ax20e.scatter(age.flatten()[age.flatten()>0], total_int.flatten()[age.flatten()>0], label='.', s = 0.025)
#ax20e.plot(np.linspace(5,timesteps, 10000), fitfunc1(np.linspace(5,timesteps, 10000), *popt1), 'r-',
#         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt1))
ax20e.set_xlabel("Age")
ax20e.set_ylabel("Total Interactions")
plt.savefig(savedir+expName+"/age_stress_scat.png", dpi=150, format="png")


#popt2, pcov2 = curve_fit(fitfunc1, age[:,:fillIndex:p].flatten()[age[:,:fillIndex:p].flatten()>0], total_int[:,:fillIndex:p].flatten()[age[:,:fillIndex:p].flatten()>0])
fig20f, ax20f = plt.subplots()
ASs = ax20f.scatter(age[:,:fillIndex].flatten()[age[:,:fillIndex].flatten()>0], total_int[:,:fillIndex].flatten()[age[:,:fillIndex].flatten()>0], label='.', s = 0.025)
#ax20f.plot(np.linspace(5,timesteps, 10000), fitfunc1(np.linspace(5,timesteps, 10000), *popt1), 'r-',
#         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt1))
ax20f.set_xlabel("Age")
ax20f.set_ylabel("Total Interactions")
plt.savefig(savedir+expName+"/age_stress_scatfill.png", dpi=150, format="png")


#plotting age vs stress plots only for interior cells
nonconfineIndex = confine_int < 0.1
fig20g, ax20g = plt.subplots()
AS = ax20g.hist2d(age[nonconfineIndex].flatten()[age[nonconfineIndex].flatten()>0], total_int[nonconfineIndex].flatten()[age[nonconfineIndex].flatten()>0], bins = 20, norm=LogNorm(), cmap=plt.cm.jet)
ax20g.set_xlim(left = 0)
ax20g.set_xlabel("Age")
ax20g.set_ylabel("Total Interactions")
fig20g.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stress_con.png", dpi=150, format="png")



fig21, ax21 = plt.subplots()
RS = ax21.hist2d(radius.flatten()[radius.flatten()>0.1], total_int.flatten()[radius.flatten()>0.1], bins = 20, norm=LogNorm(), cmap=plt.cm.jet)
ax21.set_xlim(4, 7)
ax21.set_xlabel("Radius")
ax21.set_ylabel("Total Interactions")
fig21.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/rad_stress.png", dpi=150, format="png")

fig21b, ax21b = plt.subplots()
RS = ax21b.hist2d(radius.flatten()[radius.flatten()>0.1], total_int.flatten()[radius.flatten()>0.1], bins = 20, cmap=plt.cm.jet)
ax21b.set_xlim(4, 7)
ax21b.set_xlabel("Radius")
ax21b.set_ylabel("Total Interactions")
fig21b.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/rad_stressnl.png", dpi=150, format="png")
#Heatmap of growth rate as function of number of cells and distance from center

fig21c, ax21c = plt.subplots()
RS = ax21c.hist2d(radius[:,:fillIndex].flatten()[radius[:,:fillIndex].flatten()>0.1], total_int[:,:fillIndex].flatten()[radius[:,:fillIndex].flatten()>0.1], bins = 20, norm=LogNorm(), cmap=plt.cm.jet)
ax21c.set_xlim(4, 7)
ax21c.set_xlabel("Radius")
ax21c.set_ylabel("Total Interactions")
fig21c.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/rad_stress_fill.png", dpi=150, format="png")

fig21d, ax21d = plt.subplots()
RS = ax21d.hist2d(radius[:,:fillIndex].flatten()[radius[:,:fillIndex].flatten()>0.1], total_int[:,:fillIndex].flatten()[radius[:,:fillIndex].flatten()>0.1], bins = 20, cmap=plt.cm.jet)
ax21d.set_xlim(4, 7)
ax21d.set_xlabel("Radius")
ax21d.set_ylabel("Total Interactions")
fig21d.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/rad_stressnl_fill.png", dpi=150, format="png")
#Heatmap of growth rate as function of number of cells and distance from center

fig21e, ax21e = plt.subplots()
RS = ax21e.hist2d(radius[nonconfineIndex].flatten()[radius[nonconfineIndex].flatten()>0.1], total_int[nonconfineIndex].flatten()[radius[nonconfineIndex].flatten()>0.1], bins = 20, norm=LogNorm(), cmap=plt.cm.jet)
ax21e.set_xlim(4, 7)
ax21e.set_xlabel("Radius")
ax21e.set_ylabel("Total Interactions")
fig21e.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/rad_stress_con.png", dpi=150, format="png")

fig22, ax22 = plt.subplots()
RS = ax22.hist2d(neighbours.flatten()[radius.flatten()>0.1], total_int.flatten()[radius.flatten()>0.1], bins = [9,20], norm=LogNorm(), cmap=plt.cm.jet)
ax22.set_xlim(left = 0, right = 8)
#ax22.set_xlabel("neighbours")
ax22.set_ylabel("Total Interactions")
xticks = [i for i in range(9)]
xtick_labels = ['%d' % (f) for f in range(9)]
ax22.set_xticks([])
ax22.set_xticklabels([])
fig22.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/nei_stress.png", dpi=150, format="png")

fig22b, ax22b = plt.subplots()
RS = ax22b.hist2d(neighbours.flatten()[radius.flatten()>0.1], total_int.flatten()[radius.flatten()>0.1], bins = [9,20], cmap=plt.cm.jet)
ax22b.set_xlim(left = 0, right = 8)
#ax22b.set_xlabel("neighbours")
ax22b.set_ylabel("Total Interactions")
ax22b.set_xticks([])
ax22b.set_xticklabels([])
#for label in ax22b.xaxis.get_xticklabels():
#    label.set_horizontalalignment('right')
fig22b.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/nei_stressnl.png", dpi=150, format="png")
#Heatmap of growth rate as function of number of cells and distance from center


fig22c, ax22c = plt.subplots()
RS = ax22c.hist2d(neigh_xx, normint_yy, bins = [9,20], cmin = 500, vmin = 1, vmax = 40000, norm=LogNorm(), cmap=plt.cm.jet)
ax22c.set_xlim(left = 0, right = 8)
#ax22c.set_xlabel("neighbours")
ax22c.set_ylabel(r'$\frac{T_i}{L}$')
ax22c.set_xticks([])
ax22c.set_xticklabels([])
#for label in ax22c.xaxis.get_xticklabels():
#    label.set_horizontalalignment('center')
fig22c.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/nei_stress_fill.png", dpi=150, format="png")

fig22d, ax22d = plt.subplots()
RS = ax22d.hist2d(neigh_xx, normint_yy, bins = [9,20], cmin = 500, vmin = 1, vmax = 40000, cmap=plt.cm.jet)
ax22d.set_xlim(left = 0, right = 8)
#ax22d.set_xlabel("neighbours")
ax22d.set_ylabel(r'$\frac{T_i}{L}$')
ax22d.set_xticks([])
ax22d.set_xticklabels([])
#for label in ax22d.xaxis.get_xticklabels():
#    label.set_horizontalalignment('center')
fig22d.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/nei_stressnl_fill.png", dpi=150, format="png")

"""
fig20d, ax20d = plt.subplots()
AS = ax20d.hist2d(age_xx, normint_yy, bins = 31, cmin = 250, vmin = 1, vmax = 3500, range = [[0, 25], [0, 1.5000]], cmap=plt.cm.jet)
ax20d.set_xlim(left = 0, right=25)
ax20d.set_ylim(0, 1.5000)
ax20d.set_xlabel("Age")
ax20d.set_ylabel("$T_i/L_i$")
fig20d.colorbar(AS[3], label = "Count")
plt.savefig(savedir+expName+"/age_stressnl_fill.png", dpi=150, format="png")
"""

fig22e, ax22e = plt.subplots()
RS = ax22e.hist2d(neighbours[nonconfineIndex].flatten()[radius[nonconfineIndex].flatten()>0.1], total_int[nonconfineIndex].flatten()[radius[nonconfineIndex].flatten()>0.1], bins = [9,20], norm=LogNorm(), cmap=plt.cm.jet)
ax22e.set_xlim(left = 0, right = 8)
#ax22e.set_xlabel("neighbours")
ax22e.set_ylabel("Total Interactions")
ax22e.set_xticks([])
ax22e.set_xticklabels([])
#for tick in ax22.xaxis.get_major_ticks():
#    tick.label.set_horizontalalignment('left')
fig22e.colorbar(RS[3], label = "Count")
plt.savefig(savedir+expName+"/nei_stress_con.png", dpi=150, format="png")

tic4 = datetime.now()
print(tic4-tic3)


#https://matplotlib.org/stable/gallery/statistics/time_series_histogram.html#sphx-glr-gallery-statistics-time-series-histogram-py

#plot cell trajectories

#MSD

#vorticity field: w = dv_y/dx - dv_x/dy

#nematic orientational order <cos 2*theta>

#aspect ratio




#kymographs
figk, axk = plt.subplots(1,2)
pc1 = axk[0].pcolormesh(tt, dd, vel_rad_interp, cmap = 'Spectral', norm =mpl.colors.Normalize(vmin=-0.2, vmax=0.2))
pc2 = axk[1].pcolormesh(tt, dd, vel_ortho_interp, cmap = 'Spectral', norm =mpl.colors.Normalize(vmin=-0.2, vmax=0.2))
axk[0].set_xlabel("time")
axk[1].set_xlabel("time")
axk[0].set_ylabel("distance from center")
axk[1].set_ylabel("distance from center")
axk[0].set_title("Radial Vel")
axk[1].set_title("Orthoradial Vel")
figk.colorbar(pc1, ax = axk[0],norm =mpl.colors.Normalize(vmin=-0.2, vmax=0.2))
figk.colorbar(pc2, ax = axk[1],norm =mpl.colors.Normalize(vmin=-0.2, vmax=0.2))
figk.tight_layout()
plt.savefig(savedir+expName+'/kymograph.png', dpi=150, format='png')

figk, axk = plt.subplots(1,2)
pc1 = axk[0].pcolormesh(tt, dd, gaussian_filter(vel_rad_interp, sigma = 3), cmap = 'Spectral', norm = mpl.colors.Normalize(vmin = -0.2, vmax = 0.2))
pc2 = axk[1].pcolormesh(tt, dd, gaussian_filter(vel_ortho_interp, sigma = 3), cmap = 'Spectral', norm = mpl.colors.Normalize(vmin = -0.2, vmax = 0.2))
axk[0].set_xlabel("time")
axk[1].set_xlabel("time")
axk[0].set_ylabel("distance from center")
axk[1].set_ylabel("distance from center")
axk[0].set_title("Radial Vel")
axk[1].set_title("Orthoradial Vel")
figk.colorbar(pc1, ax = axk[0], norm = mpl.colors.Normalize(vmin = -0.2, vmax = 0.2))
figk.colorbar(pc2, ax = axk[1], norm = mpl.colors.Normalize(vmin = -0.2, vmax = 0.2))
figk.tight_layout()
plt.savefig(savedir+expName+'/kymographfilt.png', dpi=150, format='png')

figkl, axkl = plt.subplots(1,2)
pc1l = axkl[0].pcolormesh(tt, dd, np.log(vel_rad_interp), cmap = 'Spectral')
pc2l = axkl[1].pcolormesh(tt, dd, np.log(vel_ortho_interp), cmap = 'Spectral')
axkl[0].set_xlabel("time")
axkl[1].set_xlabel("time")
axkl[0].set_ylabel("distance from center")
axkl[1].set_ylabel("distance from center")
axkl[0].set_title("Log Radial Vel")
axkl[1].set_title("Log Orthoradial Vel")
figkl.colorbar(pc1l, ax = axkl[0])
figkl.colorbar(pc2l, ax = axkl[1])
figkl.tight_layout()
plt.savefig(savedir+expName+'/kymographlog.png', dpi=150, format='png')