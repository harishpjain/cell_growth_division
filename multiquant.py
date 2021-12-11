import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import sys
import matplotlib as mpl
plt.style.use('seaborn-bright')
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter1d
#mpl.rcParams['text.usetex'] = True
#mpl.use('PDF')

savedir = '/scratch/ws/1/haja565a-workspace2/quant/'


SimulsL = ["700a11b", '700a12', '700a13b', '700a14b', "700a15b"]
SimulsIn = ["700a13b",  '700a17b','700a16b', '700a18']
SimulsCa = ["700a13b", '700a19b']

#Simuls = ["700a13b", "700h13"]


Simuls = ["700a11b", '700a12', '700a13b', '700a14b', "700a15b", '700a16b', '700a17b', '700a18']

Simuls = ["700h11", "700h12c", "700h13", "700h14c" ,"700h15", "700h16c", "700h17", "700h18"]
DictL = {"700a11b":5000, '700a12':7500, '700a13b':10000, '700a14b':12500, "700a15b":15000, '700a16b':10000, '700a17b':10000, '700a18':10000}#, '700a19b':10000}
DictIn = {"700a11b":5e-2, '700a12':5e-2, '700a13b':5e-2, '700a14b':5e-2, "700a15b":5e-2,  '700a17b':2.5e-2, '700a16b':1e-2,'700a18':0.5e-2, '700a19b':5e-2}
DictCa = {"700a11b":10e-2, '700a12':10e-2, '700a13b':10e-2, '700a14b':10e-2, "700a15b":10e-2, '700a16b':10e-2, '700a17b':10e-2, '700a18':10e-2, '700a19b':15e-2}


SimulsL = ["700h11", "700h12c", "700h13", "700h14c" ,"700h15"]
SimulsIn = ["700h13", "700h17", "700h16c", "700h18"]#"700h13"
DictL = {"700h11":5000, "700h12c":7500, '700h13':10000, "700h14c":12500,  "700h15":15000, '700h16c':10000, '700h17':10000, '700h18':10000, '700h19':10000}
DictIn = {"700h11":5e-2, '700h13':5e-2, "700h15":5e-2, '700h17':2.5e-2, '700h16c':1e-2, '700h18':0.5e-2, '700h19':5e-3}

Simuls = ["700a13b", "700v11","700v12","700v13","700h13","700v14","700v15"]
SimulsV = ["700a13b", "700v11","700v12","700v13","700h13","700v14","700v15"]
DictV = {"700a13b": 0, "700v11":0.25, "700v12": 0.5 ,"700v13": 0.75 ,"700h13":1, "700v14": 1.25 ,"700v15": 1.5}
Dictrad = {}
DictT = {}
DictVx = {}
DictVy = {}
DictVmag = {}
DictNumCells = {}
DictX = {}
DictY = {}
DictS0 = {}
DictS1 = {}
DictSmag = {}
Dictneigh_avg = {}

#for key in DictL.keys():
for key in Simuls:
    Dictrad[key] = np.load(savedir + key + "/radius.npy")
    DictT[key] = np.load(savedir + key + "/T.npy")
    DictX[key] = np.load(savedir + key + "/X.npy")
    DictY[key] = np.load(savedir + key + "/Y.npy")
    DictVx[key] = np.load(savedir + key + "/Vx.npy")
    DictVy[key] = np.load(savedir + key + "/Vy.npy")
    #DictVx[key] = np.gradient(DictX[key], axis = 1)/0.005
    #DictVy[key] = np.gradient(DictY[key], axis = 1)/0.005
    DictVmag[key] = np.linalg.norm([DictVx[key], DictVy[key]], axis = 0)
    DictNumCells[key] = np.count_nonzero(Dictrad[key]>0.01, axis = 0)
    DictS0[key] = np.load(savedir + key + "/S0.npy")
    DictS1[key] = np.load(savedir + key + "/S1.npy")
    DictSmag[key] = np.linalg.norm([DictS0[key], DictS1[key]], axis = 0)
    #Dictneigh_avg[key] = np.load(savedir + key + "/neighavgin.npy")
"""
fig, ax = plt.subplots()
ax.plot(DictT["700a13b"], Dictneigh_avg["700a13b"], label = "non migrating")
ax.plot(DictT["700h13"], Dictneigh_avg["700h13"], label = "migrating")
ax.legend(fontsize=14)
ax.set_xlabel("time")
ax.set_ylabel(r"$<neighbours>$")
ax.set_xlim(0,50)
ax.grid()
plt.savefig(savedir +'neigh_avg_700ac_noac_growth.png', dpi=300)
"""
"""
fig,ax = plt.subplots()
for key in SimulsL:
    ax.plot(DictT[key], Dictneigh_avg[key], label = r"$L_i = %s$"%DictL[key])
ax.legend()
ax.set_xlabel("time")
ax.set_ylabel(r"$<neighbours>$")
ax.set_xlim(0,100)
ax.grid()
plt.savefig(savedir +'neigh_avg_700aL.png')

fig,ax = plt.subplots()
for key in SimulsIn:
    ax.plot(DictT[key], Dictneigh_avg[key], label = r"$I_n = %s$"%DictIn[key])
ax.legend()
ax.set_xlabel("time")
ax.set_ylabel(r"$<neighbours>$")
ax.set_xlim(0,100)
ax.grid()
plt.savefig(savedir +'neigh_avg_700aI.png')
"""
"""
fig, ax = plt.subplots()
ax.plot(DictT['700h13'], gaussian_filter(DictVmag['700h13'][0], sigma = 10), 'r')
ax.plot(DictT['700h13'], gaussian_filter(DictVmag['700h13'][1], sigma = 10), 'g')
ax.plot(DictT['700h13'], gaussian_filter(DictVmag['700h13'][2], sigma = 10), 'b')
ax.plot(DictT['700h13'], gaussian_filter(DictVmag['700h13'][3], sigma = 10), 'k')
ax.plot(DictT['700h13'], gaussian_filter(DictVmag['700h13'][4], sigma = 10), 'm')
plt.savefig(savedir +'temp.png')
"""
"""
figvm, axvm = plt.subplots()
for key in SimulsL:
    axvm.plot(DictT[key], gaussian_filter(np.sum(DictVmag[key], axis = 0)/DictNumCells[key], sigma=10), label = r"$L_i = %s$"%DictL[key] + ", " + r"$I_n = %s$"%DictIn[key] )
for key in SimulsIn:
    axvm.plot(DictT[key], gaussian_filter(np.sum(DictVmag[key], axis = 0)/DictNumCells[key], sigma=10), label = r"$L_i = %s$"%DictL[key] + ", " + r"$I_n = %s$"%DictIn[key])
axvm.set_xlabel(r"$t$")
axvm.set_ylabel(r"$|v|$")
axvm.grid()
axvm.legend()
axvm.set_xlim(0, 150)
plt.savefig(savedir + "Speed700hxx.png", dpi=150)

figvm, axvm = plt.subplots()
for key in SimulsL:
    axvm.plot(DictT[key], gaussian_filter(np.sum(DictSmag[key], axis = 0)/DictNumCells[key], sigma=10), label = r"$L_i = %s$"%DictL[key] + ", " + r"$I_n = %s$"%DictIn[key] )
for key in SimulsIn:
    axvm.plot(DictT[key], gaussian_filter(np.sum(DictSmag[key], axis = 0)/DictNumCells[key], sigma=10), label = r"$L_i = %s$"%DictL[key] + ", " + r"$I_n = %s$"%DictIn[key])
axvm.set_xlabel(r"$t$")
axvm.set_ylabel(r"$\lambda$")
axvm.set_ylim(0,5)
axvm.grid()
axvm.legend()
axvm.set_xlim(0, 150)
plt.savefig(savedir + "Deformation700hxx.png", dpi=150)
"""

"""
figLv, axLv = plt.subplots()
for key in SimulsL:
    axLv.plot(DictT[key], np.sum(np.pi*Dictrad[key]**2, axis=0), label = r"$L = %s$"%DictL[key])
axLv.set_xlabel(r"$t$")
axLv.set_ylabel(r"$V_c$")
axLv.grid()
axLv.legend()
axLv.set_xlim(0, 100)
#plt.savefig(savedir + "ColonyVolume700axxL.png", dpi=150)
plt.savefig(savedir + "ColonyVolume700hxxL.png", dpi=150)

figIv, axIv = plt.subplots()
for key in SimulsIn:
    axIv.plot(DictT[key], np.sum(np.pi*Dictrad[key]**2, axis=0), label = r"$In = %s$"%DictIn[key])
axIv.set_xlabel(r"$t$")
axIv.set_ylabel(r"$V_c$")
axIv.grid()
axIv.legend()
axIv.set_xlim(0, 100)
#plt.savefig(savedir + "ColonyVolume700axxI.png", dpi=150)
plt.savefig(savedir + "ColonyVolume700hxxI.png", dpi=150)
"""
"""
figCv, axCv = plt.subplots()
for key in SimulsCa:
    axCv.plot(DictT[key], np.sum(np.pi*Dictrad[key]**2, axis=0), label = DictCa[key])
axCv.set_xlabel(r"$t$")
axCv.set_ylabel(r"$V_c$")
axCv.grid()
axCv.legend()
axCv.set_xlim(0, 100)
#plt.savefig(savedir + "ColonyVolume700axxC.png", dpi=150)
plt.savefig(savedir + "ColonyVolume700hxxC.png", dpi=150)
"""

my_colors = ['blue', 'green', 'red', 'magenta', 'orange', 'purple', 'grey']
figLr, axLr = plt.subplots()
indexc = 0
for key in SimulsV:
    axLr.plot(DictT[key], np.sum(Dictrad[key]**2, axis=0)**0.5, label = r"$v = %s$"%DictV[key], color = my_colors[indexc])
    indexc+=1
axLr.set_xlabel(r"$t$", fontsize=14)
axLr.set_ylabel("colony radius", fontsize=14)
axLr.grid()
axLr.legend(fontsize=14)
axLr.set_xlim(0, 50)
axLr.set_ylim(bottom=0.0)
axLr.plot(np.linspace(0,100, 20), np.ones(20)*45, color = 'k', linewidth=1, linestyle= "dotted")
plt.savefig(savedir + "ColonyRadius700vel.png", dpi=200)

"""
figLr, axLr = plt.subplots()
for key in SimulsL:
    axLr.plot(DictT[key], np.sum(Dictrad[key]**2, axis=0)**0.5, label = r"$L = %s$"%DictL[key])
axLr.set_xlabel(r"$t$", fontsize=14)
axLr.set_ylabel("colony radius", fontsize=14)
axLr.grid()
axLr.legend(fontsize=14)
axLr.set_xlim(0, 100)
axLr.set_ylim(bottom=0.0)
axLr.plot(np.linspace(0,100, 20), np.ones(20)*45, color = 'k', linewidth=1, linestyle= "dotted")
#plt.savefig(savedir + "ColonyRadius700axxL.png", dpi=200)
plt.savefig(savedir + "ColonyRadius700hxxL.png", dpi=200)

figIr, axIr = plt.subplots()
for key in SimulsIn:
    axIr.plot(DictT[key], np.sum(Dictrad[key]**2, axis=0)**0.5, label = r"$In = %s$"%DictIn[key])
axIr.set_xlabel(r"$t$", fontsize=14)
axIr.set_ylabel("colony radius", fontsize=14)
axIr.grid()
axIr.legend(fontsize=14)
axIr.set_xlim(0, 100)
axIr.set_ylim(bottom=0.0)
axIr.plot(np.linspace(0,100, 20), np.ones(20)*45, color = 'k', linewidth=1, linestyle= "dotted")
#plt.savefig(savedir + "ColonyRadius700axxI.png", dpi=200)
plt.savefig(savedir + "ColonyRadius700hxxI.png", dpi=200)
"""
"""
figCr, axCr = plt.subplots()
for key in SimulsCa:
    axCr.plot(DictT[key], np.sum(Dictrad[key]**2, axis=0)**0.5, label = DictCa[key])
axCr.set_xlabel("time")
axCr.set_ylabel("Colony Radius")
axCr.grid()
axCr.legend()
axCr.set_xlim(0, 100)
axCr.set_ylim(bottom=0.0)
#plt.savefig(savedir + "ColonyRadius700axxC.png", dpi=150)
plt.savefig(savedir + "ColonyRadius700hxxI.png", dpi=150)

#printing fill fractions
print("printing fill fractions")
print("For Li")
figLf, axLf = plt.subplots()
for key in SimulsL:
    print(np.sum(np.pi*Dictrad[key]**2, axis=0)[-1]/np.sum(np.pi*Dictrad["700a15b"]**2, axis=0)[-1])
    axLf.scatter(DictL[key], np.sum(np.pi*Dictrad[key]**2, axis=0)[-1])
axLf.set_xlabel("Inhibition Limit")
axLf.set_ylabel("Colony Volume")
plt.savefig(savedir+ "ColonyVolume700axxFracL.png", dpi=150)


print("For In")
figIf, axIf = plt.subplots()
for key in SimulsIn:
    print(np.sum(np.pi*Dictrad[key]**2, axis=0)[-1]/np.sum(np.pi*Dictrad["700a13b"]**2, axis=0)[-1])
    axIf.scatter(DictIn[key], np.sum(np.pi*Dictrad[key]**2, axis=0)[-1])
axIf.set_xlabel("I_n")
axIf.set_ylabel("Colony Volume")
plt.savefig(savedir+ "ColonyVolume700axxFracI.png", dpi=150)


#printing slopes of colony radius vs time graph in the linear region
print("printing slopes of linear region of colony radius vs time")
print("for Li")
#for time from 400 to 4400
figLs, axLs = plt.subplots()
for key in SimulsL:
    print(str(DictL[key]) + ":" + str(((np.sum(Dictrad[key]**2, axis=0)**0.5)[4400] - (np.sum(Dictrad[key]**2, axis=0)**0.5)[400])/20))
    axLs.scatter(DictL[key], ((np.sum(Dictrad[key]**2, axis=0)**0.5)[4400] - (np.sum(Dictrad[key]**2, axis=0)**0.5)[400])/20)
axLs.set_xlabel("Li")
axLs.set_ylabel("Slope of Radius vs Time in Linear Region")
plt.savefig(savedir+ "ColonyRadius700axxSlopeL.png", dpi=150)

print("for In")
#for time from 400 to 4400
figIs, axIs = plt.subplots()
for key in SimulsIn:
    print(str(DictIn[key]) + ":" + str(((np.sum(Dictrad[key]**2, axis=0)**0.5)[4400] - (np.sum(Dictrad[key]**2, axis=0)**0.5)[400])/20))
    axIs.scatter(DictIn[key], ((np.sum(Dictrad[key]**2, axis=0)**0.5)[4400] - (np.sum(Dictrad[key]**2, axis=0)**0.5)[400])/20)
axIs.set_xlabel("I_n")
axIs.set_ylabel("Slope of Radius vs Time in Linear Region")
plt.savefig(savedir+ "ColonyRadius700axxSlopeI.png", dpi=150)
"""