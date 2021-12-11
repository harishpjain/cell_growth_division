import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib as mpl
import csv
plt.style.use('seaborn-bright')
#mpl.rcParams['text.usetex'] = True
#mpl.use('PDF')
positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "Total_int_potential": 10}

expNames = ['r701', 'r702', 'r703', 'r704','r705','r706','r707','r708','r709','r710',\
            'r711', 'r712', 'r713', 'r714','r715','r716','r717','r718','r719','r720',\
            'r700', 'r70a', 'r70b', 'r70c', 'r70d', '700k14', '700k15', '700k16', \
                '700k17', '700k18', '700k19', '700k11', '700k12', '700k13']
expRad0 = {}
expTime0 = {}
expTI = {}

directory = "/scratch/ws/1/haja565a-workspace2/master_thesis/output"

for exp in expNames:
    
    expRad0[exp] = np.genfromtxt(directory + exp + '/positions_p0.csv', delimiter=',',skip_header=1)[:,4]
    expTime0[exp] = np.genfromtxt(directory + exp + '/positions_p0.csv', delimiter=',',skip_header=1)[:,0]
    expTI[exp] = np.genfromtxt(directory + exp + '/positions_p0.csv', delimiter=',',skip_header=1)[:,10]

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k14'], (expRad0['700k14']), label = r"$L_i = 2000$", color = 'g')
axes.plot(expTime0['700k15'], (expRad0['700k15']), label = r"$L_i = 6000$", color = 'b')
axes.plot(expTime0['700k16'], (expRad0['700k16']), label = r"$L_i = 10000$", color = 'k')
axes.plot(expTime0['700k17'], (expRad0['700k17']), label = r"$L_i = 15000$", color = 'purple')
axes.plot(expTime0['700k18'], (expRad0['700k18']), label = r"$L_i = 20000$", color = 'm')
axes.plot(expTime0['700k19'], (expRad0['700k19']), label = r"$L_i = 25000$", color = 'r')
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadius_700kx.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k11'], (expRad0['700k11']), label = r"$I_n = 0.005$", color = 'purple')
axes.plot(expTime0['700k12'], (expRad0['700k12']), label = r"$I_n = 0.01$", color = 'm')
axes.plot(expTime0['700k13'], (expRad0['700k13']), label = r"$I_n = 0.025$", color = 'r')
axes.plot(expTime0['700k16'], (expRad0['700k16']), label = r"$I_n = 0.05$", color = 'k')
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadius_700kxI.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k14'], np.pi*(expRad0['700k14']**2), label = r"$L_i = 2000$", color = 'g')
axes.plot(expTime0['700k15'], np.pi*(expRad0['700k15']**2), label = r"$L_i = 6000$", color = 'b')
axes.plot(expTime0['700k16'], np.pi*(expRad0['700k16']**2), label = r"$L_i = 10000$", color = 'k')
axes.plot(expTime0['700k17'], np.pi*(expRad0['700k17']**2), label = r"$L_i = 15000$", color = 'purple')
axes.plot(expTime0['700k18'], np.pi*(expRad0['700k18']**2), label = r"$L_i = 20000$", color = 'm')
axes.plot(expTime0['700k19'], np.pi*(expRad0['700k19']**2), label = r"$L_i = 25000$", color = 'r')
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume_700kx.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k11'], np.pi*(expRad0['700k11']**2), label = r"$I_n = 0.005$", color = 'purple')
axes.plot(expTime0['700k12'], np.pi*(expRad0['700k12']**2), label = r"$I_n = 0.01$", color = 'm')
axes.plot(expTime0['700k13'], np.pi*(expRad0['700k13']**2), label = r"$I_n = 0.025$", color = 'r')
axes.plot(expTime0['700k16'], np.pi*(expRad0['700k16']**2), label = r"$I_n = 0.05$", color = 'k')
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume_700kxI.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k14'], np.log(np.pi*(expRad0['700k14']**2)), label = r"$L_i = 2000$", color = 'g')
axes.plot(expTime0['700k15'], np.log(np.pi*(expRad0['700k15']**2)), label = r"$L_i = 6000$", color = 'b')
axes.plot(expTime0['700k16'], np.log(np.pi*(expRad0['700k16']**2)), label = r"$L_i = 10000$", color = 'k')
axes.plot(expTime0['700k17'], np.log(np.pi*(expRad0['700k17']**2)), label = r"$L_i = 15000$", color = 'purple')
axes.plot(expTime0['700k18'], np.log(np.pi*(expRad0['700k18']**2)), label = r"$L_i = 20000$", color = 'm')
axes.plot(expTime0['700k19'], np.log(np.pi*(expRad0['700k19']**2)), label = r"$L_i = 25000$", color = 'r')
axes.set_xlabel('time')
axes.set_ylabel('log volume')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolumeLog_700kx.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k11'], np.log(np.pi*(expRad0['700k11']**2)), label = r"$I_n = 0.005$", color = 'purple')
axes.plot(expTime0['700k12'], np.log(np.pi*(expRad0['700k12']**2)), label = r"$I_n = 0.01$", color = 'm')
axes.plot(expTime0['700k13'], np.log(np.pi*(expRad0['700k13']**2)), label = r"$I_n = 0.025$", color = 'r')
axes.plot(expTime0['700k16'], np.log(np.pi*(expRad0['700k16']**2)), label = r"$I_n = 0.05$", color = 'k')
axes.set_xlabel('time')
axes.set_ylabel('log volume')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolumeLog_700kxI.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k14'][10:-100], np.gradient(gaussian_filter(expRad0['700k14'][10:-100], sigma=20),edge_order=2)/0.005, label = r"$L_i = 2000$", color = 'g')
axes.plot(expTime0['700k15'][10:-100], np.gradient(gaussian_filter(expRad0['700k15'][10:-100], sigma=20),edge_order=2)/0.005, label = r"$L_i = 6000$", color = 'b')
axes.plot(expTime0['700k16'][10:-100], np.gradient(gaussian_filter(expRad0['700k16'][10:-100], sigma=20),edge_order=2)/0.005, label = r"$L_i = 10000$", color = 'k')
axes.plot(expTime0['700k17'][10:-100], np.gradient(gaussian_filter(expRad0['700k17'][10:-100], sigma=20),edge_order=2)/0.005, label = r"$L_i = 15000$", color = 'purple')
axes.plot(expTime0['700k18'][10:-100], np.gradient(gaussian_filter(expRad0['700k18'][10:-100], sigma=20),edge_order=2)/0.005, label = r"$L_i = 20000$", color = 'm')
axes.plot(expTime0['700k19'][10:-100], np.gradient(gaussian_filter(expRad0['700k19'][10:-100], sigma=20),edge_order=2)/0.005, label = r"$L_i = 25000$", color = 'r')
axes.set_xlabel('time')
axes.set_ylabel(r'$\dot r$')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadiusRate_700kx.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k14'], (expTI['700k14']), label = r"$L_i = 2000$", color = 'g')
axes.plot(expTime0['700k15'], (expTI['700k15']), label = r"$L_i = 6000$", color = 'b')
axes.plot(expTime0['700k16'], (expTI['700k16']), label = r"$L_i = 10000$", color = 'k')
axes.plot(expTime0['700k17'], (expTI['700k17']), label = r"$L_i = 15000$", color = 'purple')
axes.plot(expTime0['700k18'], (expTI['700k18']), label = r"$L_i = 20000$", color = 'm')
axes.plot(expTime0['700k19'], (expTI['700k19']), label = r"$L_i = 25000$", color = 'r')
axes.set_xlabel('time')
axes.set_ylabel('Total Interactions')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingInteractions_700kx.png', dpi=300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['700k11'], (expTI['700k11']), label = r"$I_n = 0.005$", color = 'purple')
axes.plot(expTime0['700k12'], (expTI['700k12']), label = r"$I_n = 0.01$", color = 'm')
axes.plot(expTime0['700k13'], (expTI['700k13']), label = r"$I_n = 0.025$", color = 'r')
axes.plot(expTime0['700k16'], (expTI['700k16']), label = r"$I_n = 0.05$", color = 'k')
axes.set_xlabel('time')
axes.set_ylabel('Total Interactions')
axes.legend()
axes.grid()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingInteractions_700kxI.png', dpi=300, format='png')

"""
fig, axes = plt.subplots(1,1, figsize = ((8,8)))

#axes.plot(expTime0['701'], np.log(np.pi*expRad0['701']**2), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['702'], np.log(np.pi*expRad0['702']**2), label = r"$L_i = 2000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['703'], np.log(np.pi*expRad0['703']**2), label = r"$L_i = 4000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['704'], np.log(np.pi*expRad0['704']**2), label = r"$L_i = 6000$", color = 'k', linestyle = 'dotted')
axes.plot(expTime0['705'], np.log(np.pi*expRad0['705']**2), label = r"$L_i = 8000$", color = 'purple', linestyle = 'dotted')
axes.plot(expTime0['706'], np.log(np.pi*expRad0['706']**2), label = r"$L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['707'], np.log(np.pi*expRad0['707']**2), label = r"$L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['708'], np.log(np.pi*expRad0['708']**2), label = r"$L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['709'], np.log(np.pi*expRad0['709']**2), label = r"$L_i = 25000$", color = 'k', linestyle = 'dashed')
axes.plot(expTime0['710'], np.log(np.pi*expRad0['710']**2), label = r"$L_i = 30000$", color = 'purple', linestyle = 'dashed')
axes.plot(expTime0['70a'], np.log(np.pi*expRad0['70a']**2), label = r"$L_i = 50000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['70b'], np.log(np.pi*expRad0['70b']**2), label = r"$L_i = 100000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['70c'], np.log(np.pi*expRad0['70c']**2), label = r"$L_i = 250000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['70d'], np.log(np.pi*expRad0['70d']**2), label = r"$L_i = 500000$", color = 'k', linestyle = 'dashdot')
axes.plot(expTime0['700'], np.log(np.pi*expRad0['700']**2), label = r"$L_i = 1000000$", color = 'purple', linestyle = 'dashdot')
#r"\bf{phase field} $\phi$"
#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel(r'$\log$(volume)')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingLogVolume_70x.pdf', dpi = 300, format='pdf')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))

#axes.plot(expTime0['701'], (np.pi*expRad0['701']**2), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['702'], (np.pi*expRad0['702']**2), label = r"$L_i = 2000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['703'], (np.pi*expRad0['703']**2), label = r"$L_i = 4000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['704'], (np.pi*expRad0['704']**2), label = r"$L_i = 6000$", color = 'k', linestyle = 'dotted')
axes.plot(expTime0['705'], (np.pi*expRad0['705']**2), label = r"$L_i = 8000$", color = 'purple', linestyle = 'dotted')
axes.plot(expTime0['706'], (np.pi*expRad0['706']**2), label = r"$L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['707'], (np.pi*expRad0['707']**2), label = r"$L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['708'], (np.pi*expRad0['708']**2), label = r"$L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['709'], (np.pi*expRad0['709']**2), label = r"$L_i = 25000$", color = 'k', linestyle = 'dashed')
axes.plot(expTime0['710'], (np.pi*expRad0['710']**2), label = r"$L_i = 30000$", color = 'purple', linestyle = 'dashed')
axes.plot(expTime0['70a'], (np.pi*expRad0['70a']**2), label = r"$L_i = 50000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['70b'], (np.pi*expRad0['70b']**2), label = r"$L_i = 100000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['70c'], (np.pi*expRad0['70c']**2), label = r"$L_i = 250000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['70d'], (np.pi*expRad0['70d']**2), label = r"$L_i = 500000$", color = 'k', linestyle = 'dashdot')
axes.plot(expTime0['700'], (np.pi*expRad0['700']**2), label = r"$L_i = 1000000$", color = 'purple', linestyle = 'dashdot')


#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume_70x.pdf', dpi=300, format='pdf')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))

axes.plot(expTime0['711'], np.log(np.pi*expRad0['711']**2), label = r"$I_n = 2.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['712'], np.log(np.pi*expRad0['712']**2), label = r"$I_n = 2.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['713'], np.log(np.pi*expRad0['713']**2), label = r"$I_n = 2.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['714'], np.log(np.pi*expRad0['714']**2), label = r"$I_n =5.0e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['715'], np.log(np.pi*expRad0['715']**2), label = r"$I_n =5.0e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['716'], np.log(np.pi*expRad0['716']**2), label = r"$I_n =5.0e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['717'], np.log(np.pi*expRad0['717']**2), label = r"$I_n =7.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['718'], np.log(np.pi*expRad0['718']**2), label = r"$I_n =7.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['719'], np.log(np.pi*expRad0['719']**2), label = r"$I_n =7.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['720'], np.log(np.pi*expRad0['720']**2), label = r"$I_n =2.5e^{-2}, L_i = 30000$", color = 'k', linestyle = 'dotted')

#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('$\log$(volume)')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingLogVolume_71x.pdf', dpi = 300, format='pdf')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))

axes.plot(expTime0['711'], (np.pi*expRad0['711']**2), label = r"$I_n =2.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['712'], (np.pi*expRad0['712']**2), label = r"$I_n =2.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['713'], (np.pi*expRad0['713']**2), label = r"$I_n =2.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['714'], (np.pi*expRad0['714']**2), label = r"$I_n =5.0e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['715'], (np.pi*expRad0['715']**2), label = r"$I_n =5.0e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['716'], (np.pi*expRad0['716']**2), label = r"$I_n =5.0e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['717'], (np.pi*expRad0['717']**2), label = r"$I_n =7.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['718'], (np.pi*expRad0['718']**2), label = r"$I_n =7.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['719'], (np.pi*expRad0['719']**2), label = r"$I_n =7.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['720'], (np.pi*expRad0['720']**2), label = r"$I_n =2.5e^{-2}, L_i = 30000$", color = 'k', linestyle = 'dotted')

#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume_71x.pdf', dpi = 300, format='pdf')


#radius
fig, axes = plt.subplots(1,1, figsize = ((8,8)))

#axes.plot(expTime0['701'], np.log(expRad0['701']**2), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['702'], np.log(expRad0['702']), label = r"$L_i = 2000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['703'], np.log(expRad0['703']), label = r"$L_i = 4000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['704'], np.log(expRad0['704']), label = r"$L_i = 6000$", color = 'k', linestyle = 'dotted')
axes.plot(expTime0['705'], np.log(expRad0['705']), label = r"$L_i = 8000$", color = 'purple', linestyle = 'dotted')
axes.plot(expTime0['706'], np.log(expRad0['706']), label = r"$L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['707'], np.log(expRad0['707']), label = r"$L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['708'], np.log(expRad0['708']), label = r"$L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['709'], np.log(expRad0['709']), label = r"$L_i = 25000$", color = 'k', linestyle = 'dashed')
axes.plot(expTime0['710'], np.log(expRad0['710']), label = r"$L_i = 30000$", color = 'purple', linestyle = 'dashed')
axes.plot(expTime0['70a'], np.log(expRad0['70a']), label = r"$L_i = 50000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['70b'], np.log(expRad0['70b']), label = r"$L_i = 100000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['70c'], np.log(expRad0['70c']), label = r"$L_i = 250000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['70d'], np.log(expRad0['70d']), label = r"$L_i = 500000$", color = 'k', linestyle = 'dashdot')
axes.plot(expTime0['700'], np.log(expRad0['700']), label = r"$L_i = 1000000$", color = 'purple', linestyle = 'dashdot')
#r"\bf{phase field} $\phi$"
#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel(r'$\log (radius)$')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingLogRadius_70x.pdf', dpi = 300, format='pdf')


fig, axes = plt.subplots(1,1, figsize = ((8,8)))

#axes.plot(expTime0['701'], (expRad0['701']), label = r"$L_i = 10", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['702'], (expRad0['702']), label = r"$L_i = 2000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['703'], (expRad0['703']), label = r"$L_i = 4000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['704'], (expRad0['704']), label = r"$L_i = 6000$", color = 'k', linestyle = 'dotted')
axes.plot(expTime0['705'], (expRad0['705']), label = r"$L_i = 8000$", color = 'purple', linestyle = 'dotted')
axes.plot(expTime0['706'], (expRad0['706']), label = r"$L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['707'], (expRad0['707']), label = r"$L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['708'], (expRad0['708']), label = r"$L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['709'], (expRad0['709']), label = r"$L_i = 25000$", color = 'k', linestyle = 'dashed')
axes.plot(expTime0['710'], (expRad0['710']), label = r"$L_i = 30000$", color = 'purple', linestyle = 'dashed')
axes.plot(expTime0['70a'], (expRad0['70a']), label = r"$L_i = 50000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['70b'], (expRad0['70b']), label = r"$L_i = 100000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['70c'], (expRad0['70c']), label = r"$L_i = 250000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['70d'], (expRad0['70d']), label = r"$L_i = 500000$", color = 'k', linestyle = 'dashdot')
axes.plot(expTime0['700'], (expRad0['700']), label = r"$L_i = 1000000$", color = 'purple', linestyle = 'dashdot')


#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadius_70x.pdf', dpi=300, format='pdf')



fig, axes = plt.subplots(1,1, figsize = ((8,8)))

axes.plot(expTime0['711'], np.log(expRad0['711']), label = r"$I_n = 2.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['712'], np.log(expRad0['712']), label = r"$I_n = 2.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['713'], np.log(expRad0['713']), label = r"$I_n = 2.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['714'], np.log(expRad0['714']), label = r"$I_n =5.0e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['715'], np.log(expRad0['715']), label = r"$I_n =5.0e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['716'], np.log(expRad0['716']), label = r"$I_n =5.0e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['717'], np.log(expRad0['717']), label = r"$I_n =7.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['718'], np.log(expRad0['718']), label = r"$I_n =7.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['719'], np.log(expRad0['719']), label = r"$I_n =7.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['720'], np.log(expRad0['720']), label = r"$I_n =2.5e^{-2}, L_i = 30000$", color = 'k', linestyle = 'dotted')

#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('$\log$(radius)')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingLogRadius_71x.pdf', dpi = 300, format='pdf')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))

axes.plot(expTime0['711'], (expRad0['711']), label = r"$I_n =2.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['712'], (expRad0['712']), label = r"$I_n =2.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['713'], (expRad0['713']), label = r"$I_n =2.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['714'], (expRad0['714']), label = r"$I_n =5.0e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['715'], (expRad0['715']), label = r"$I_n =5.0e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['716'], (expRad0['716']), label = r"$I_n =5.0e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['717'], (expRad0['717']), label = r"$I_n =7.5e^{-2}, L_i = 10000$", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['718'], (expRad0['718']), label = r"$I_n =7.5e^{-2}, L_i = 15000$", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['719'], (expRad0['719']), label = r"$I_n =7.5e^{-2}, L_i = 20000$", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['720'], (expRad0['720']), label = r"$I_n =2.5e^{-2}, L_i = 30000$", color = 'k', linestyle = 'dotted')

#axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadius_71x.pdf', dpi = 300, format='pdf')

"""