import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import csv
plt.style.use('seaborn')

positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "Total_int_potential": 10}

expNames = ['641', '642', '643', '644', '645','646', '647', '648', '649', '651','652', '653', '654',\
            '621', '622', '623', '624', '625','626', '627', '628', '629',\
            '631', '632', '633', '634', '635','636', '637', '638', '639']
expRad0 = {}
expTime0 = {}
expMuInMax = {'641':2000, '642':3000, '643':4000, '644':5000, '645':6000,'646':7000, '647':8000, '648':9000, '649':10000\
                , '651':400,'652':800, '653':1200, '654':1600, \
                '621':3750, '622':5000, '623':6250, '624':1875, '625':2500,'626':3125, '627':1250, '628':1666, '629':2083}

directory = "/scratch/ws/1/haja565a-workspace2/master_thesis/outputring"
for exp in expNames:
    
    expRad0[exp] = np.genfromtxt(directory + exp + '/positions_p0.csv', delimiter=',',skip_header=1)[:,4]
    expTime0[exp] = np.genfromtxt(directory + exp + '/positions_p0.csv', delimiter=',',skip_header=1)[:,0]

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
axes.plot(expTime0['651'], expRad0['651'], label = "400", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['652'], expRad0['652'], label = "800", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['653'], expRad0['653'], label = "1200", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['654'], expRad0['654'], label = "1600", color = 'k', linestyle = 'dashdot')
axes.plot(expTime0['641'], expRad0['641'], label = "2000", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['642'], expRad0['642'], label = "3000", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['643'], expRad0['643'], label = "4000", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['644'], expRad0['644'], label = "5000", color = 'k', linestyle = 'dotted')
axes.plot(expTime0['645'], expRad0['645'], label = "6000", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['646'], expRad0['646'], label = "7000", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['647'], expRad0['647'], label = "8000", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['648'], expRad0['648'], label = "9000", color = 'k', linestyle = 'dashed')
axes.plot(expTime0['649'], expRad0['649'], label = "10000", color = 'm', linestyle = 'dashdot')

axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadius_64x65x.png', dpi = 300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
axes.plot(expTime0['651'], np.pi*expRad0['651']**2, label = "400", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['652'], np.pi*expRad0['652']**2, label = "800", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['653'], np.pi*expRad0['653']**2, label = "1200", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['654'], np.pi*expRad0['654']**2, label = "1600", color = 'k', linestyle = 'dashdot')
axes.plot(expTime0['641'], np.pi*expRad0['641']**2, label = "2000", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['642'], np.pi*expRad0['642']**2, label = "3000", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['643'], np.pi*expRad0['643']**2, label = "4000", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['644'], np.pi*expRad0['644']**2, label = "5000", color = 'k', linestyle = 'dotted')
axes.plot(expTime0['645'], np.pi*expRad0['645']**2, label = "6000", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['646'], np.pi*expRad0['646']**2, label = "7000", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['647'], np.pi*expRad0['647']**2, label = "8000", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['648'], np.pi*expRad0['648']**2, label = "9000", color = 'k', linestyle = 'dashed')
axes.plot(expTime0['649'], np.pi*expRad0['649']**2, label = "10000", color = 'm')

axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume_64x65x.png', dpi = 300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))

axes.plot(expTime0['621'], expRad0['621'], label = "3750 and 0.4", color='b')
axes.plot(expTime0['622'], expRad0['622'], label = "5000 and 0.4", color='g')
axes.plot(expTime0['623'], expRad0['623'], label = "6250 and 0.4", color='k')
axes.plot(expTime0['624'], expRad0['624'], label = "1875 and 0.8", color='m')
axes.plot(expTime0['625'], expRad0['625'], label = "2500 and 0.8", color='r')
axes.plot(expTime0['626'], expRad0['626'], label = "3125 and 0.8", color='purple')
axes.plot(expTime0['627'], expRad0['627'], label = "1250 and 1.2", color='orange')
axes.plot(expTime0['628'], expRad0['628'], label = "1666.66 and 1.2", color='fuchsia')
axes.plot(expTime0['629'], expRad0['629'], label = "2083.33 and 1.2", color='brown')

axes.set_title("Effect of inhibition and interaction")
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadius2_62x.png', dpi = 300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
axes.plot(expTime0['621'], np.pi*expRad0['621']**2, label = "3750 and 0.4", color='b')
axes.plot(expTime0['622'], np.pi*expRad0['622']**2, label = "5000 and 0.4", color='g')
axes.plot(expTime0['623'], np.pi*expRad0['623']**2, label = "6250 and 0.4", color='r')
axes.plot(expTime0['624'], np.pi*expRad0['624']**2, label = "1875 and 0.8", color='fuchsia')
axes.plot(expTime0['625'], np.pi*expRad0['625']**2, label = "2500 and 0.8", color='brown')
axes.plot(expTime0['626'], np.pi*expRad0['626']**2, label = "3125 and 0.8", color='m')
axes.plot(expTime0['627'], np.pi*expRad0['627']**2, label = "1250 and 1.2", color='k')
axes.plot(expTime0['628'], np.pi*expRad0['628']**2, label = "1666.66 and 1.2", color='orange')
axes.plot(expTime0['629'], np.pi*expRad0['629']**2, label = "2083.33 and 1.2", color='purple')

axes.set_title("Effect of inhibition and interaction")
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume2_62x.png', dpi = 300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))

axes.plot(expTime0['631'], expRad0['631'], label = "3750 and 0.4", color='b')
axes.plot(expTime0['632'], expRad0['632'], label = "5000 and 0.4", color='g')
axes.plot(expTime0['633'], expRad0['633'], label = "6250 and 0.4", color='r')
axes.plot(expTime0['634'], expRad0['634'], label = "1875 and 0.8", color='fuchsia')
axes.plot(expTime0['635'], expRad0['635'], label = "2500 and 0.8", color='orange')
axes.plot(expTime0['636'], expRad0['636'], label = "3125 and 0.8", color='purple')
axes.plot(expTime0['637'], expRad0['637'], label = "1250 and 1.2", color='brown')
axes.plot(expTime0['638'], expRad0['638'], label = "1666.66 and 1.2", color='k')
axes.plot(expTime0['639'], expRad0['639'], label = "2083.33 and 1.2", color='m')

axes.set_title("Effect of inhibition and interaction")
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingRadius3_63x.png', dpi = 300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))

axes.plot(expTime0['631'], np.pi*expRad0['631']**2, label = "3750 and 0.4", color='b')
axes.plot(expTime0['632'], np.pi*expRad0['632']**2, label = "5000 and 0.4", color='g')
axes.plot(expTime0['633'], np.pi*expRad0['633']**2, label = "6250 and 0.4", color='r')
axes.plot(expTime0['634'], np.pi*expRad0['634']**2, label = "1875 and 0.8", color='fuchsia')
axes.plot(expTime0['635'], np.pi*expRad0['635']**2, label = "2500 and 0.8", color='m')
axes.plot(expTime0['636'], np.pi*expRad0['636']**2, label = "3125 and 0.8", color='k')
axes.plot(expTime0['637'], np.pi*expRad0['637']**2, label = "1250 and 1.2", color='brown')
axes.plot(expTime0['638'], np.pi*expRad0['638']**2, label = "1666.66 and 1.2", color='purple')
axes.plot(expTime0['639'], np.pi*expRad0['639']**2, label = "2083.33 and 1.2", color='orange')

axes.set_title("Effect of inhibition and interaction")
axes.set_xlabel('time')
axes.set_ylabel('radius')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume3_63x.png', dpi = 300, format='png')



#changed order of plots legends
fig, axes = plt.subplots(1,1, figsize = ((8,8)))
axes.plot(expTime0['627'], np.pi*expRad0['627']**2, label = "1250, 1.2 and 1500", color='b')
axes.plot(expTime0['628'], np.pi*expRad0['628']**2, label = "1666.66, 1.2 and 2000", color='g')
axes.plot(expTime0['624'], np.pi*expRad0['624']**2, label = "1875, 0.8 and 1500", color='r')
axes.plot(expTime0['629'], np.pi*expRad0['629']**2, label = "2083.33, 1.2 and 2500", color='fuchsia')
axes.plot(expTime0['625'], np.pi*expRad0['625']**2, label = "2500, 0.8 and 2000", color='m')
axes.plot(expTime0['626'], np.pi*expRad0['626']**2, label = "3125, 0.8 and 2500", color='orange')
axes.plot(expTime0['621'], np.pi*expRad0['621']**2, label = "3750, 0.4 and 1500", color='k')
axes.plot(expTime0['622'], np.pi*expRad0['622']**2, label = "5000, 0.4 and 2000", color='purple')
axes.plot(expTime0['623'], np.pi*expRad0['623']**2, label = "6250, 0.4 and 2500", color='brown')

axes.set_title("Effect of inhibition and interaction (Legend: New Inhibition Limit, In and Old Inhibition Limit)")
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume4_62x.png', dpi = 300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
axes.plot(expTime0['637'], np.pi*expRad0['637']**2, label = "1250, 1.2 and 1500", color='b')
axes.plot(expTime0['638'], np.pi*expRad0['638']**2, label = "1666.66, 1.2 and 2000", color='g')
axes.plot(expTime0['634'], np.pi*expRad0['634']**2, label = "1875, 0.8 and 1500", color='r')
axes.plot(expTime0['639'], np.pi*expRad0['639']**2, label = "2083.33, 1.2 and 2500", color='fuchsia')
axes.plot(expTime0['635'], np.pi*expRad0['635']**2, label = "2500, 0.8 and 2000", color='m')
axes.plot(expTime0['636'], np.pi*expRad0['636']**2, label = "3125, 0.8 and 2500", color='orange')
axes.plot(expTime0['631'], np.pi*expRad0['631']**2, label = "3750, 0.4 and 1500", color='k')
axes.plot(expTime0['632'], np.pi*expRad0['632']**2, label = "5000, 0.4 and 2000", color='purple')
axes.plot(expTime0['633'], np.pi*expRad0['633']**2, label = "6250, 0.4 and 2500", color='brown')

axes.set_title("Effect of inhibition and interaction (Legend: New Inhibition Limit, In and Old Inhibition Limit)")
axes.set_xlabel('time')
axes.set_ylabel('volume')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingVolume5_63x.png', dpi = 300, format='png')


#changed order of plots legends logarithm
fig, axes = plt.subplots(1,1, figsize = ((8,8)))
axes.plot(expTime0['627'], np.log(np.pi*expRad0['627']**2), label = "1250, 1.2 and 1500", color='b')
axes.plot(expTime0['628'], np.log(np.pi*expRad0['628']**2), label = "1666.66, 1.2 and 2000", color='g')
axes.plot(expTime0['624'], np.log(np.pi*expRad0['624']**2), label = "1875, 0.8 and 1500", color='r')
axes.plot(expTime0['629'], np.log(np.pi*expRad0['629']**2), label = "2083.33, 1.2 and 2500", color='fuchsia')
axes.plot(expTime0['625'], np.log(np.pi*expRad0['625']**2), label = "2500, 0.8 and 2000", color='m')
axes.plot(expTime0['626'], np.log(np.pi*expRad0['626']**2), label = "3125, 0.8 and 2500", color='orange')
axes.plot(expTime0['621'], np.log(np.pi*expRad0['621']**2), label = "3750, 0.4 and 1500", color='k')
axes.plot(expTime0['622'], np.log(np.pi*expRad0['622']**2), label = "5000, 0.4 and 2000", color='purple')
axes.plot(expTime0['623'], np.log(np.pi*expRad0['623']**2), label = "6250, 0.4 and 2500", color='brown')

axes.set_title("Effect of inhibition and interaction (Legend: New Inhibition Limit, In and Old Inhibition Limit)")
axes.set_xlabel('time')
axes.set_ylabel('log(volume)')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingLogVolume4_62x.png', dpi = 300, format='png')

fig, axes = plt.subplots(1,1, figsize = ((8,8)))
axes.plot(expTime0['651'], np.log(np.pi*expRad0['651']**2), label = "400", color = 'r', linestyle = 'dashdot')
axes.plot(expTime0['652'], np.log(np.pi*expRad0['652']**2), label = "800", color = 'g', linestyle = 'dashdot')
axes.plot(expTime0['653'], np.log(np.pi*expRad0['653']**2), label = "1200", color = 'b', linestyle = 'dashdot')
axes.plot(expTime0['654'], np.log(np.pi*expRad0['654']**2), label = "1600", color = 'k', linestyle = 'dashdot')
axes.plot(expTime0['641'], np.log(np.pi*expRad0['641']**2), label = "2000", color = 'r', linestyle = 'dotted')
axes.plot(expTime0['642'], np.log(np.pi*expRad0['642']**2), label = "3000", color = 'g', linestyle = 'dotted')
axes.plot(expTime0['643'], np.log(np.pi*expRad0['643']**2), label = "4000", color = 'b', linestyle = 'dotted')
axes.plot(expTime0['644'], np.log(np.pi*expRad0['644']**2), label = "5000", color = 'k', linestyle = 'dotted')
axes.plot(expTime0['645'], np.log(np.pi*expRad0['645']**2), label = "6000", color = 'r', linestyle = 'dashed')
axes.plot(expTime0['646'], np.log(np.pi*expRad0['646']**2), label = "7000", color = 'g', linestyle = 'dashed')
axes.plot(expTime0['647'], np.log(np.pi*expRad0['647']**2), label = "8000", color = 'b', linestyle = 'dashed')
axes.plot(expTime0['648'], np.log(np.pi*expRad0['648']**2), label = "9000", color = 'k', linestyle = 'dashed')
axes.plot(expTime0['649'], np.log(np.pi*expRad0['649']**2), label = "10000", color = 'm')

axes.set_title("Effect of inhibition limit for constant In")
axes.set_xlabel('time')
axes.set_ylabel('log(volume)')
axes.legend()
plt.savefig('/scratch/ws/1/haja565a-workspace2/CellInRingLogVolume_64x65x.png', dpi = 300, format='png')