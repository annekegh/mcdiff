#!/bin/python

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from outreading import read_F_D_edges,read_many_profiles

def plotsettings():
    plt.rc(('xtick','ytick','axes'), labelsize=24.0)
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['figure.subplot.left']=0.15
    plt.rcParams['figure.subplot.bottom']=0.14
    plt.rcParams['legend.fontsize'] = 18
def plotsettingsax(ax):
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        #line.set_color('green')
        line.set_markersize(10)
        line.set_markeredgewidth(2)
plotsettings()


import sys
list_filenames = sys.argv[1:]
print list_filenames
list_filenames.sort()

# guess name
figname="fig_"
if len(list_filenames) > 0:
   name = list_filenames[0]
   if "mito" in name: figname = "fig_mito_"
   elif "HexdWat" in name: figname = "fig_HexdWat_"
   elif "popc" in name: figname = "fig_popc_"
   if "outrad" in name:
       figname += "outrad_"

# lagtimes?
lts = []
dtc = 1.
nbins = 100
for filename in list_filenames:
    words = filename.split("/")
    parts = words[-1].split(".")
    if "outrad" in filename:
        shift = int(parts[4])
    else:
        shift = int(parts[2])
    #lts.extend([shift*dtc for i in range(nbins)])
    lts.extend([shift*dtc])
lts = np.array(lts)
# sort
#arg = np.argsort(lts)
#lts = lts[arg]
#list_filenames = list_filenames[arg]

F,D,edges = read_many_profiles(list_filenames)
#print F
#print F.ravel()  # this is row by row

print "F",F.shape
print "lts",lts
# F has shape nbins x nprofiles

#######
nlt = len(lts)
A = np.zeros((nbins*nlt,2*nbins))
res = np.zeros((nbins*nlt),float)
for i,lt in enumerate(lts):
    for j in range(nbins):
        res[i*nbins+j] = D[j,i]
        A[i*nbins+j,j] = 1.
        A[i*nbins+j,j+nbins] = 1./lt
#print A

print A[:6,:6]
print A[-6:,-6:]
print A[:6,-6:]
print "res",res

print A.shape, res.shape
x,residues,rank,s = np.linalg.lstsq(A,res)
print x
print x.shape
#print "residues",residues

xD = edges[:len(D)]
Deff = x[:nbins]
t0 = x[nbins:]/x[:nbins]

plt.figure()
plt.plot(xD,D)
plt.plot(xD,Deff,lw=5,color='black')
plt.savefig(figname+"Dreal.png")

plt.figure()
plt.plot(xD,Deff,lw=5,color='black')
plt.plot(xD,t0,color='red')
plt.legend(["Dreal","t0"])
plt.savefig(figname+"timezero.png") 


#######
print done

Fr = D.ravel()

p = np.polyfit(lts,Fr,1)
print p
plt.figure()
plt.plot(lts,Fr,"o")
plt.savefig("fig_lts.png")

p = np.polyfit(1./lts,Fr,1)
print p
plt.figure()
plt.plot(1./lts,Fr,"o")
plt.savefig("fig_lts1.png")

