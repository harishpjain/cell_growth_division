import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import os
import sys
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

times = np.linspace(15.0,28.0,14, endpoint=True)
#delt = 0.005

#ranks = list(range(0,128))

#interpolation_steps = 1000
#x = np.linspace(0,100,interpolation_steps)
#xx, yy = np.meshgrid(x,x)
#positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "Total_int_potential": 10}
#write_idx = 0
for time in times:
    print(time)
    reader = vtk.vtkXMLUnstructuredGridReader()

    delt = 0.005
    ind = int(time/delt)
    ranks = list(range(0,180))
    S0 = []
    S1 = []
    interaction_pot = []
    interpolation_steps = 1000
    x = np.linspace(0,100,interpolation_steps)
    phi_all = []
    xx, yy = np.meshgrid(x,x)
    positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "Total_int_potential": 10}


    for rank_ind,rank in enumerate(ranks):
        filename = sys.argv[1] + '/data/phase_p' + str(rank) + '_' + '{:06.3f}'.format(time) + '.vtu'
        reader.SetFileName(filename)
        reader.Update()
        data = reader.GetOutput()

        # grid points
        points = data.GetPoints()
        points = vtk_to_numpy(points.GetData())
        # values 
        phi = vtk_to_numpy(data.GetPointData().GetArray(0))


        phi_interp = griddata(points[:,0:2],phi,(xx,yy),method='nearest')

        phi_interp = np.reshape(phi_interp,(interpolation_steps,interpolation_steps))
        phi_all.append(0.5 * phi_interp + 0.5)

        positionfile = sys.argv[1] + '/positions_p' + str(rank) + '.csv'
        positiondata = np.genfromtxt(positionfile, delimiter=',',skip_header=1)

        S0.append(positiondata[ind][5])
        S1.append(positiondata[ind][6])
        interaction_pot.append(positiondata[ind][10])
        #S0.append(positions_raw[rank_ind].iloc[ind]['S0'])
        #S1.append(positions_raw[rank_ind].iloc[ind]['S1'])

    # after this we have a list IN THE SAME ORDER AS ranks with all phi
    # now the axis 0 here is the rank axis which we want to remove
    phi_all = np.array(phi_all)

    S0 = np.array(S0)
    S1 = np.array(S1)
    interaction_pot = np.array(interaction_pot)
    # global phasefield, given by a lot of 1s and something in between

    phi_glob = np.max(phi_all,axis=0)
    # this is more interesting -> this locally gives the rank the contributed to the max, ie the cell
    rank_max = np.argmax(phi_all,axis=0)



    S0_glob = S0[rank_max] * phi_glob
    S1_glob = S1[rank_max] * phi_glob
    interaction_pot_glob = interaction_pot[rank_max] * phi_glob
    norm_S = np.sqrt(S0_glob**2 + S1_glob**2)
    S0_globn = S0_glob / norm_S
    S1_globn = S1_glob / norm_S
    S0_globn_rot = -S1_globn #rotation by pi/2
    S1_globn_rot = S0_globn

    omega = np.multiply(np.sign(S1_globn),(np.arctan(S0_globn/np.abs(S1_globn))/ 2.0 + np.pi/4.0))
    S0_shift = np.sin(omega)
    S1_shift = np.cos(omega)

    S0_globr = S0_glob/np.amax(S0_glob)
    S1_globr = S0_glob/np.amax(S1_glob)
    #S_mag = np.hypot(S0_glob, S1_glob)
    #S_mag = np.sqrt(S0_glob**2+S1_glob**2)

    fig1, ax1 = plt.subplots()
    I = ax1.contourf(xx, yy, interaction_pot_glob, cmap = "Reds")
    Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], S0_shift[::20, ::20], S1_shift[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
    #Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], S0_globn_rot[::20, ::20], S1_globn_rot[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
    #Q = ax1.quiver(xx[::10, ::10], yy[::10, ::10], S0_globn_rot[::10, ::10], S1_globn_rot[::10, ::10], norm_S[::10, ::10], units='width', pivot='mid', headaxislength=0, headlength=0)
    plt.colorbar(I)
    #ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
    #             label='Quiver key, length = 10', labelpos='E')
    plt.savefig('/scratch/ws/1/haja565a-workspace2/DeformationField/ring/400d56/' + '{:06.3f}'.format(time) + '_' + sys.argv[2] + '.png', dpi = 300, format='png')
    fig2, ax2 = plt.subplots()
    plt.close()
