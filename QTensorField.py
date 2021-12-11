import numpy as np
from numpy.lib.npyio import save
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk



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

time = float(sys.argv[2])
outdir = "/scratch/ws//1/haja565a-workspace2/master_thesis/output"
reader = vtk.vtkXMLUnstructuredGridReader()
savedir = "/scratch/ws//1/haja565a-workspace2/DeformationField/"
expName = sys.argv[1]

saveVTK = True
"""
if(abs(time-81.0) < 1.0):
    saveVTK = True
else:
    saveVTK = False
if not os.path.exists(savedir+expName):
    os.makedirs(savedir+expName)
"""
gaussian_blur = True
plotkymograph = False

delt = 0.005
ind = int(time/delt)
ranks = list(range(0,int(sys.argv[3])))
S0 = []
S1 = []
interaction_pot = []
Vx = []
Vy = []
interpolation_steps = 1000
domain_size = 100
x = np.linspace(0, domain_size,interpolation_steps)
phi_all = []
xx, yy = np.meshgrid(x,x)
positions_columns = {'time': 0, 'rank': 1, "posx": 2, "posy": 3, "radius": 4, "S0": 5, "S1": 6, "velx": 7, "vely": 8, "angle": 9, "Total_int_potential": 10}


for rank_ind,rank in enumerate(ranks):
    filename = outdir + sys.argv[1] + '/data/phase_p' + str(rank) + '_' + '{:06.3f}'.format(time) + '.vtu'
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

    positionfile = outdir + sys.argv[1] + '/positions_p' + str(rank) + '.csv'
    positiondata = np.genfromtxt(positionfile, delimiter=',',skip_header=1)

    S0.append(positiondata[ind][5])
    S1.append(positiondata[ind][6])
    interaction_pot.append(positiondata[ind][10])
    Vx.append(positiondata[ind][7])
    Vy.append(positiondata[ind][8])


    #S0.append(positions_raw[rank_ind].iloc[ind]['S0'])
    #S1.append(positions_raw[rank_ind].iloc[ind]['S1'])
print("stage 1 over")

# after this we have a list IN THE SAME ORDER AS ranks with all phi
# now the axis 0 here is the rank axis which we want to remove
phi_all = np.array(phi_all)

S0 = np.array(S0)
S1 = np.array(S1)
interaction_pot = np.array(interaction_pot)
Vx = np.array(Vx)
Vy = np.array(Vy)
# global phasefield, given by a lot of 1s and something in between

phi_glob = np.max(phi_all,axis=0)
# this is more interesting -> this locally gives the rank the contributed to the max, ie the cell
rank_max = np.argmax(phi_all,axis=0)


Vxglob = Vx[rank_max] * np.sign(phi_glob)
Vyglob = Vy[rank_max] * np.sign(phi_glob)

normV = np.sqrt(Vxglob**2 + Vyglob**2) #velocity magnitude
Uxglob = np.divide(Vxglob, normV, out=np.zeros_like(Vxglob), where=normV!=0) #unit x vel
Uyglob = np.divide(Vyglob, normV, out=np.zeros_like(Vyglob), where=normV!=0) #unit y vel
#Vxglobf = gaussian_filter(Vxglob, sigma=40, mode = 'constant')#, mode='nearest') #filter velx
#Vyglobf = gaussian_filter(Vyglob, sigma=40, mode = 'constant')#, mode='nearest') #filter vely
Vxglobf = gaussian_filter(Vxglob, sigma=80, mode = 'constant')#, mode='nearest') #filter velx
Vyglobf = gaussian_filter(Vyglob, sigma=80, mode = 'constant')#, mode='nearest') #filter vely
vorticity = np.gradient(Vyglobf, axis = 1) - np.gradient(Vxglobf, axis = 0) 
vorticity = np.multiply(vorticity, np.sign(phi_glob)) #vorticity

## saving velocity to vtk
if(saveVTK):
    np.savetxt(savedir+expName+'/Vx'+ '{:06.3f}'.format(time) +'.csv', Vxglobf, delimiter =",")
    np.savetxt(savedir+expName+'/Vy'+ '{:06.3f}'.format(time) +'.csv', Vyglobf, delimiter =",")
    """
    Vxg = Vxglobf.reshape(interpolation_steps, interpolation_steps, 1)
    Vyg = Vyglobf.reshape(interpolation_steps, interpolation_steps, 1)
    Vx_vtk = numpy_to_vtk(Vxglobf, deep=True, array_type=vtk.VTK_FLOAT)
    Vy_vtk = numpy_to_vtk(Vyglobf, deep=True, array_type=vtk.VTK_FLOAT)

    imageVTK = vtk.vtkImageData()
    imageVTK.SetSpacing(100/999, 100/999, 1)
    imageVTK.SetOrigin([0.0, 0.0, 0.0])
    imageVTK.SetDimensions(1000, 1000, 1)
    imageVTK.GetPointData().SetVectors([Vx_vtk, Vy_vtk])
    """

if(plotkymograph):
    rad_vel = Vxglobf*np.cos(angpos)+Vyglobf*np.sin(angpos)
    ortho_vel = Vxglobf*np.sin(angpos)-Vyglobf*np.cos(angpos)
    dist_cent = np.hypot(cx, cy)


S0_glob = S0[rank_max] * phi_glob
S1_glob = S1[rank_max] * phi_glob
interaction_pot_glob = interaction_pot[rank_max] * phi_glob

S0_glob = gaussian_filter(S0_glob, sigma=80, mode='nearest')
S1_glob = gaussian_filter(S1_glob, sigma=80, mode='nearest')

a = np.gradient(S0_glob, axis = 0)
b = np.gradient(S0_glob, axis = 1)
c = np.gradient(S1_glob, axis = 0)
d = np.gradient(S1_glob, axis = 1)
dow = a*d-b*c

norm_S = np.sqrt(S0_glob**2 + S1_glob**2) #magnitude of Deformation?? Is that right?

S0_globn = S0_glob / norm_S
S1_globn = S1_glob / norm_S
#S0_globn_rot = -S1_globn #rotation by pi/2
#S1_globn_rot = S0_globn

omega = np.multiply(np.sign(S1_globn),(np.arctan(S0_glob/np.abs(S1_glob))/ 2.0 + np.pi/4.0))
eig_v0 = np.sin(omega)
eig_v1 = np.cos(omega)

#S0_globr = S0_glob/np.amax(S0_glob)
#S1_globr = S0_glob/np.amax(S1_glob)

print("stage 2 over") 


#lines of trivial S Tensor
fig, ax = plt.subplots()
c1 = ax.contour(xx,yy,S0_glob,0,colors='red',alpha=1.0)
c2 = ax.contour(xx,yy,S1_glob,0,colors='blue',alpha=1.0)    
plt.savefig(savedir+expName+'/conti'+ '{:06.3f}'.format(time) +'.png', dpi = 200, format='png')
defect_list = []


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
#vorticity contour plot + velocity field
figVor, axVor = plt.subplots(1,1)
ConVor = axVor.contourf(xx, yy, vorticity, cmap = 'bwr')
VelF = axVor.quiver(xx[::20, ::20], yy[::20, ::20], Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
figVor.colorbar(ConVor, ax=axVor)
plt.savefig(savedir+expName+'/vorticity_' + '{:05d}'.format(ind) + '.png', dpi=300, format='png')

#velocity field 
figVel, axVel = plt.subplots(1,1)
Qvel = axVel.quiver(xx[::20, ::20], yy[::20, ::20],  Vxglob[::20, ::20], Vyglob[::20, ::20], units='x', scale=0.5)
plt.savefig(savedir+expName+'/velocity_' + '{:05d}'.format(ind) + '.png', dpi=300, format='png')


fig10, ax10 = plt.subplots(1,2, figsize=(10,4))
Cont = ax10[0].contourf(xx, yy, vorticity, cmap = 'bwr')
VelF = ax10[0].quiver(xx[::20, ::20], yy[::20, ::20], Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig10.colorbar(Cont, ax=ax10[0])
I = ax10[1].contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha = 0.3)
Q = ax10[1].quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Cont2 = ax10[1].contourf(xx, yy, norm_S, cmap = 'bwr', alpha=0.3)
for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(dow[indx, indy]>0):
        ax10[1].scatter(point[0], point[1], marker="o", alpha=0.5, color='g')
    else:
        ax10[1].scatter(point[0], point[1], marker="^", alpha=0.5, color='k')

fig10.colorbar(I, ax=ax10[1])
plt.savefig(savedir + expName + '/defect_' + '{:05d}'.format(ind) + '.png', dpi=300, format='png')

fig11, ax11 = plt.subplots(1,1)
Cont = ax11.contourf(xx, yy, vorticity, cmap = 'bwr', alpha=0.3)
Q = ax11.streamplot(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], linewidth = 0.1, density = 10, arrowsize=0)
for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(dow[indx, indy]>0):
        ax11.scatter(point[0], point[1], marker="o", alpha=0.5, color='g')
    else:
        ax11.scatter(point[0], point[1], marker="^", alpha=0.5, color='k')

fig11.colorbar(Cont, ax=ax11)
plt.savefig(savedir + expName + '/stream_'+ '{:05d}'.format(ind) + '.png', dpi=300, format='png')
"""

fig12, ax12 = plt.subplots(2,2, figsize = (12,10))
ConVor = ax12[0,0].contourf(xx, yy, vorticity, cmap = 'bwr')
VelF = ax12[0,0].quiver(xx[::20, ::20], yy[::20, ::20], Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig12.colorbar(ConVor, ax=ax12[0,0], label = "vorticity")
ax12[0,0].set_title("Vorticity field with \n velocity direction field")

I1 = ax12[0, 1].contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha = 0.8)
#Qvel = ax12[0,1].quiver(xx[::20, ::20], yy[::20, ::20],  Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig12.colorbar(I1, ax=ax12[0, 1], label="Total interactions")
ax12[0,1].set_title("Total Interactions")

#interaction, Q director, defects
I = ax12[1, 0].contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha = 0.3)
Q = ax12[1, 0].quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Cont2 = ax10[1].contourf(xx, yy, norm_S, cmap = 'bwr', alpha=0.3)
for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(dow[indx, indy]>0):
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
    if(dow[indx, indy]>0):
        ax12[1,1].scatter(point[0], point[1], marker="o", alpha=0.5, color='g')
    else:
        ax12[1,1].scatter(point[0], point[1], marker="^", alpha=0.5, color='k')
fig12.colorbar(Cont, ax=ax12[1,1], label = "Vorticity")
ax12[1,1].set_title("Nematic director streamfield \n with defects and vorticity field")

plt.tight_layout()
plt.savefig(savedir + expName + '/total_'+ '{:05d}'.format(ind) + '.png', dpi=300, format='png')

fig13, ax13 = plt.subplots()
ConVor = ax13.contourf(xx, yy, vorticity, cmap = 'bwr')
VelF = ax13.quiver(xx[::20, ::20], yy[::20, ::20], Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig13.colorbar(ConVor, ax=ax13, label = "vorticity")
ax13.set_title("Vorticity field with \n velocity direction field")
plt.tight_layout()
plt.savefig(savedir + expName + '/vortex_'+ '{:05d}'.format(ind) + '.png', dpi=300, format='png')

fig14, ax14 = plt.subplots()
I1 = ax14.contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha = 0.8)
#Qvel = ax12[0,1].quiver(xx[::20, ::20], yy[::20, ::20],  Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
fig14.colorbar(I1, ax=ax14, label="Total interactions")
ax14.set_title("Total Interactions")
plt.tight_layout()
plt.savefig(savedir + expName + '/interactions_'+ '{:05d}'.format(ind) + '.png', dpi=300, format='png')

fig15, ax15 = plt.subplots()
I = ax15.contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha = 0.3)
Q = ax15.quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Cont2 = ax10[1].contourf(xx, yy, norm_S, cmap = 'bwr', alpha=0.3)
for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(dow[indx, indy]>0):
        ax15.scatter(point[0], point[1], marker="o", alpha=0.5, color='g')
    else:
        ax15.scatter(point[0], point[1], marker="^", alpha=0.5, color='k')
fig15.colorbar(I, ax=ax15, label="Total interactions")
plt.tight_layout()
plt.savefig(savedir + expName + '/defects_'+ '{:05d}'.format(ind) + '.png', dpi=300, format='png')
"""
fig1, ax1 = plt.subplots()
I = ax1.contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha = 0.3)
Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], S0_globn_rot[::20, ::20], S1_globn_rot[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Q = ax1.quiver(xx[::10, ::10], yy[::10, ::10], S0_globn_rot[::10, ::10], S1_globn_rot[::10, ::10], norm_S[::10, ::10], units='width', pivot='mid', headaxislength=0, headlength=0)
plt.colorbar(I)

ax1.xaxis.set_major_locator(ticker.NullLocator())
ax1.yaxis.set_major_locator(ticker.NullLocator())
#ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
#             label='Quiver key, length = 10', labelpos='E')
plt.savefig(savedir + expName + '/Qfield' + '{:06.3f}'.format(time) + '.png', dpi = 300, format='png')

for point in defect_list:
    indx = math.floor(point[0]*(interpolation_steps/domain_size))
    indy = math.floor(point[1]*(interpolation_steps/domain_size))
    if(dow[indx, indy]>0):
        ax1.scatter(point[0], point[1], marker="^", alpha=0.5, color='g')
    else:
        ax1.scatter(point[0], point[1], marker="o", alpha=0.5, color='k')
plt.savefig(savedir + expName + '/TopDef_' + '{:06.3f}'.format(time) +'.png', dpi = 300, format='png')
"""

"""
fig9, ax9 = plt.subplots(1,2, figsize=(10,4))
I = ax9[0].contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha=0.3)
VelF = ax9[0].quiver(xx[::20, ::20], yy[::20, ::20], Uxglob6[::20, ::20], Uyglob6[::20, ::20], units='x', scale=0.5)
fig9.colorbar(I, ax=ax9[0])
Cont = ax9[1].contourf(xx, yy, vorticity6, cmap = 'bwr')
#Vels = ax9[1].streamplot(xx[::, ::], yy[::, ::], Vxglob6[::, ::], Vyglob6[::, ::], density=10, linewidth=0.3, arrowsize=0, color = 'k')
#Uxglob6 = gaussian_filter(Uxglob6, sigma=40, mode='nearest')
#Uyglob6 = gaussian_filter(Uyglob6, sigma=40, mode='nearest')
VelF = ax9[1].quiver(xx[::40, ::40], yy[::40, ::40], Uxglob6[::40, ::40], Uyglob6[::40, ::40], units='x', scale=0.5)
fig9.colorbar(Cont, ax=ax9[1])
plt.savefig(savedir + expName + '/stream' + '{:06.3f}'.format(time)  + '.png', dpi = 300, format='png')
"""


"""
fig1, ax1 = plt.subplots()
I = ax1.contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha=0.3)
#Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
VelF = ax1.quiver(xx[::20, ::20], yy[::20, ::20], Uxglob[::20, ::20], Uyglob[::20, ::20], units='x', scale=0.5)
#Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], S0_globn_rot[::20, ::20], S1_globn_rot[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Q = ax1.quiver(xx[::10, ::10], yy[::10, ::10], S0_globn_rot[::10, ::10], S1_globn_rot[::10, ::10], norm_S[::10, ::10], units='width', pivot='mid', headaxislength=0, headlength=0)
plt.colorbar(I)
#ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
#             label='Quiver key, length = 10', labelpos='E')
plt.savefig('/scratch/ws/1/haja565a-workspace2/DeformationField/ring/vel' + '{:06.3f}'.format(time) + '_' + sys.argv[1] + '.png', dpi = 300, format='png')


figst, axst = plt.subplots()
#I = axst.contourf(xx, yy, interaction_pot_glob, cmap = "Reds", alpha=0.3)
Cont = axst.contourf(xx, yy, vorticity, cmap = 'bwr')
Vels = axst.streamplot(xx[::, ::], yy[::, ::], Uxglob[::, ::], Uyglob[::, ::], density=5, linewidth=0.1, arrowsize=0)
VelF = axst.quiver(xx[::40, ::40], yy[::40, ::40], Uxglob[::40, ::40], Uyglob[::40, ::40], units='x', scale=0.5)

plt.colorbar(Cont)
plt.savefig('/scratch/ws/1/haja565a-workspace2/DeformationField/ring/stream' + '{:06.3f}'.format(time) + '_' + sys.argv[1] + '.png', dpi = 300, format='png')
"""

#plt.savefig(savedir + expName + '/TopDef_' + '{:06.3f}'.format(time) +'.png', dpi = 300, format='png')
#Vels = ax9[1].streamplot(xx[::, ::], yy[::, ::], Vxglob6[::, ::], Vyglob6[::, ::], density=10, linewidth=0.3, arrowsize=0, color = 'k')
#Uxglob6 = gaussian_filter(Uxglob6, sigma=40, mode='nearest')
#Uyglob6 = gaussian_filter(Uyglob6, sigma=40, mode='nearest')
#VelF = ax10[1].quiver(xx[::40, ::40], yy[::40, ::40], Uxglob6[::40, ::40], Uyglob6[::40, ::40], units='x', scale=0.5)


"""

#find indices of degenerate points in S_Tensor_Field -> where eig_v1 = 0
S0_zero = np.ma.masked_where(abs(S0_glob)<0.001, S0_glob)  #ignoring all zero values(masked) of S0 i.e where we assume there is no cell
S1_zero = np.ma.masked_where(abs(S1_glob)>0.001, S1_glob)     #condition for degeneracy -> S1 = zero(not masked)
S0_zero_mask = abs(S0_glob)>0.01
S1_zero_mask = abs(S1_glob)<0.01
degen_mask = S0_zero_mask & S1_zero_mask



a = np.gradient(S0_glob, axis = 0)
b = np.gradient(S0_glob, axis = 1)
c = np.gradient(S1_glob, axis = 0)
d = np.gradient(S1_glob, axis = 1)
dow = a*d-b*c


plus_dow = dow>0.0
minus_dow = dow<0.0

plus_def = S0_zero_mask & S1_zero_mask & plus_dow
minus_def = S1_zero_mask & S1_zero_mask & minus_dow


fig1, ax1 = plt.subplots()
I = ax1.contourf(xx, yy, interaction_pot_glob, cmap = "Reds")
Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
S = ax1.scatter(xx[plus_def], yy[plus_def])
#Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], S0_globn_rot[::20, ::20], S1_globn_rot[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Q = ax1.quiver(xx[::10, ::10], yy[::10, ::10], S0_globn_rot[::10, ::10], S1_globn_rot[::10, ::10], norm_S[::10, ::10], units='width', pivot='mid', headaxislength=0, headlength=0)
plt.colorbar(I)
#ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
#             label='Quiver key, length = 10', labelpos='E')
plt.savefig('/scratch/ws/1/haja565a-workspace2/DeformationField/ring/' + 'plushalf2{:06.3f}'.format(time) + '_' + sys.argv[2] + '.png', dpi = 300, format='png')

fig1, ax1 = plt.subplots()
I = ax1.contourf(xx, yy, interaction_pot_glob, cmap = "Reds")
Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], eig_v0[::20, ::20], eig_v1[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
S = ax1.scatter(xx[minus_def], yy[minus_def])
#Q = ax1.quiver(xx[::20, ::20], yy[::20, ::20], S0_globn_rot[::20, ::20], S1_globn_rot[::20, ::20], norm_S[::20, ::20], units='width', pivot='mid', headaxislength=0, headlength=0)
#Q = ax1.quiver(xx[::10, ::10], yy[::10, ::10], S0_globn_rot[::10, ::10], S1_globn_rot[::10, ::10], norm_S[::10, ::10], units='width', pivot='mid', headaxislength=0, headlength=0)
plt.colorbar(I)
#ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
#             label='Quiver key, length = 10', labelpos='E')
plt.savefig('/scratch/ws/1/haja565a-workspace2/DeformationField/ring/' + 'minushalf2{:06.3f}'.format(time) + '_' + sys.argv[2] + '.png', dpi = 300, format='png')







degen_dow = dow[S0_zero_mask & S1_zero_mask]

print(degen_dow)


"""

#fig2, ax2 = plt.subplots()
#I = ax2.contourf(xx, yy, interaction_pot_glob)
#ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
#             label='Quiver key, length = 10', labelpos='E')
#plt.savefig('/scratch/ws/1/haja565a-workspace2/InteractionPotentialfield', dpi = 300)
#fig1, ax1 = plt.subplots()
#Q = ax1.quiver(xx, yy, S0_globr, S1_globr, norm_S, units='width')#, units='x', pivot='tip')
#ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
#             label='Quiver key, length = 10', labelpos='E')
#plt.savefig('Deformation field.jpg', dpi = 150)
"""
fig, axs = plt.subplots(2,2, figsize=(10, 10))
#plt.figsize(10,10)
Q1 = axs[0,0].quiver(xx, yy, S0_globn, S1_globn, norm_S, units='width')#, units='x', pivot='tip')
Q2 = axs[0,1].quiver(xx, yy, S0_globn, S1_globn, norm_S, units='x')#, units='x', pivot='tip')
Q3 = axs[1,0].quiver(xx, yy, S0_glob, S1_glob, norm_S, units='width')#, units='x', pivot='tip')
Q4 = axs[1,1].quiver(xx, yy, S0_glob, S1_glob, norm_S, units='x')#, units='x', pivot='tip')
#ax1.quiverkey(Q, X=0.3, Y=1.1, U=10,
#             label='Quiver key, length = 10', labelpos='E')
plt.savefig('/scratch/ws/1/haja565a-workspace2/Deformation field2.jpg', dpi = 150)
"""

print("stage 4 over") 
