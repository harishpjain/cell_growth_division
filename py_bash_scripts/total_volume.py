import numpy as np
import matplotlib.pyplot as plt

out_dir = '/beegfs/ws/1/haja565a-workspacebeegfs/out10b/'
quant_dir = out_dir + 'quant/'

cells = 150
radius = np.load(quant_dir + 'radius.npy')
#time = np.load(quant_dir + 'time.npy')
time = np.linspace(0.0, len(radius[0])*0.005, len(radius[0]))

volume = np.sum(np.pi*(radius**2), axis=0)
num_alive_cells = np.count_nonzero(radius>0.01, axis = 0)
print(time)
print(volume)
print(volume.shape)
print(num_alive_cells.shape)
print(num_alive_cells)

fig, ax = plt.subplots()
ax.plot(time, volume)
ax.plot(time, np.ones(len(time))*10000, linestyle='dotted')
ax.set_xlabel("time")
ax.set_ylabel("volume")
plt.savefig(out_dir + 'quant/total_volume.png', dpi=300)

fig, ax = plt.subplots()
ax.plot(time, num_alive_cells)
ax.set_xlabel("time")
ax.set_ylabel("total number of cells")
plt.savefig(out_dir + 'quant/total_cells.png', dpi=300)