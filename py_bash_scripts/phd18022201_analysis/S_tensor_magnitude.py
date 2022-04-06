import numpy as np
import matplotlib.pyplot as plt

out_dir = '/beegfs/ws/1/haja565a-workspacebeegfs/phd18022201/'

a_prefix = 'bout_3000t_Ca_10_In_2_5_a_'
a = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9])

Ca_prefix = 'bout_3000t_a_1_5_In_2_5_Ca_'
Ca = np.array([20.0e-2, 15.0e-2, 10.0e-2, 5.0e-2])

In_prefix = 'bout_3000t_a_1_5_Ca_10_In_'
In = np.array([10e-2, 7.5e-2, 5.0e-2, 2.5e-2, 1.0e-2, 0.75e-2, 0.5e-2])

a_lambda = np.zeros(len(a))
Ca_lambda = np.zeros(len(Ca))
In_lambda = np.zeros(len(In))

a_int= np.zeros(len(a)) #store total interactions of first cell for last timestep
Ca_int = np.zeros(len(Ca))
In_int = np.zeros(len(In))

for i, In_value in enumerate(In):
    quant_dir = out_dir + In_prefix + str(i+1) + '/quant/'
    S0 = np.load(quant_dir + 'S0.npy')[0,-1] #first cell last time step
    S1 = np.load(quant_dir + 'S1.npy')[0,-1]
    In_lambda[i] = np.sqrt(S0**2 + S1**2)
    In_int[i] = np.load(quant_dir + 'totalInteraction.npy')[1,-1]
    print(np.max(np.load(quant_dir + 'totalInteraction.npy')))

for i, Ca_value in enumerate(Ca):
    quant_dir = out_dir + Ca_prefix + str(i+1) + '/quant/'
    S0 = np.load(quant_dir + 'S0.npy')[0,-1] #first cell last time step
    S1 = np.load(quant_dir + 'S1.npy')[0,-1]
    Ca_lambda[i] = np.sqrt(S0**2 + S1**2)
    Ca_int[i] = np.load(quant_dir + 'totalInteraction.npy')[1,-1]
    print(np.max(np.load(quant_dir + 'totalInteraction.npy')))


for i, a_value in enumerate(a):
    quant_dir = out_dir + a_prefix + str(i+1) + '/quant/'
    S0 = np.load(quant_dir + 'S0.npy')[0,-1] #first cell last time step
    S1 = np.load(quant_dir + 'S1.npy')[0,-1]
    a_lambda[i] = np.sqrt(S0**2 + S1**2)
    a_int[i] = np.load(quant_dir + 'totalInteraction.npy')[1,-1]
    print(np.max(np.load(quant_dir + 'totalInteraction.npy')))

fig, ax = plt.subplots()
ax.plot(In, In_int)
ax.set_xlabel("In")
ax.set_ylabel("total interaction")
ax.set_title('Ca = 0.1, a = 1.5')
plt.savefig(out_dir + 'interaction_In.png', dpi=300)

fig, ax = plt.subplots()
ax.plot(a, a_int)
ax.set_xlabel("a")
ax.set_title('Ca = 0.1, In = 0.025')
ax.set_ylabel("total interaction")
plt.savefig(out_dir + 'interaction_a.png', dpi=300)

fig, ax = plt.subplots()
ax.plot(Ca, Ca_int)
ax.set_xlabel("Ca")
ax.set_title('In = 0.025, a = 1.5')
ax.set_ylabel("total interaction")
plt.savefig(out_dir + 'interaction_Ca.png', dpi=300)
