"""Script to plot histogram of xyz coords
AG, August 13, 2012"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from functionshistogram import read_coor
from .functionsdiffusion import fit_sqrt_vs_time


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


#=====================
# ANALYZE
#=====================
def analyze_dist_multi_1D(list_x,nstep,outdir,dtc):
    "Take the average over the D estimates in each of these files"""
    assert len(list_x) > 0
    nfiles = len(list_x)
    t = np.arange(0,nstep*dtc,dtc)

    alldist = np.zeros((nstep,nfiles),float)
    allD = np.zeros((nfiles),float)

    print("="*5)
    print("Results")
    print("making %i fits" %nfiles)
    print("fit from %i steps, i.e. time %f ps, actually time %f ps" %(nstep,nstep*dtc,(nstep-1)*dtc))
    for i in range(nfiles):
        dist = calc_dist_1D(np.array(read_coor(list_x[i])))
        d = dist[:nstep,:]   # only fit the first nstep time steps (use all atoms)
        average = np.mean(d,1)
        alldist[:,i] = average

        p = np.polyfit(t,average**2,1)
        allD[i] = p[0]   # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                                # = 1e-8 meter**2/second = 1e-4 cm**2/s
        # if I want many figures:
        #fit_sqrt_vs_time(np.mean(dist_xy,1),dtc,outdir+"/fig_dist.xy.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_z,1), dtc,outdir+"/fig_dist.z.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_r,1), dtc,outdir+"/fig_dist.r.%i.png"%i,title=str(i))

    m = np.mean(alldist,axis=1)
    #s = np.std(alldist,axis=1)

    fit_sqrt_vs_time(m,dtc,outdir+"/fig_dist.average.png",title="average of %i"%nfiles)
    print_output_D_ave_1D(allD)

def analyze_dist_multi(list_x,list_y,list_z,nstep,outdir,dtc):
    "Take the average over the D estimates in each of these files"""
    assert len(list_x) > 0
    assert len(list_x) == len(list_y)
    assert len(list_x) == len(list_z)
    nfiles = len(list_x)
    t = np.arange(0,nstep*dtc,dtc)

    alldist = np.zeros((nstep,nfiles,3),float)
    allD = np.zeros((nfiles,3),float)

    print("="*5)
    print("Results")
    print("making %i fits" %nfiles)
    print("fit from %i steps, i.e. time %f ps, actually time %f ps" %(nstep,nstep*dtc,(nstep-1)*dtc))
    for i in range(nfiles):
        dist_xy,dist_z,dist_r = calc_dist_from3files(list_x[i],list_y[i],list_z[i])
        for it,dist in enumerate([dist_xy,dist_z,dist_r]):
            d = dist[:nstep,:]   # only fit the first nstep time steps (use all atoms)
            average = np.mean(d,1)
            alldist[:,i,it] = average
            #print i,list_x[i],it,average[:10]

            p = np.polyfit(t,average**2,1)
            allD[i,it] = p[0]   # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                                # = 1e-8 meter**2/second = 1e-4 cm**2/s
        # if I want many figures:
        #fit_sqrt_vs_time(np.mean(dist_xy,1),dtc,outdir+"/fig_dist.xy.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_z,1), dtc,outdir+"/fig_dist.z.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_r,1), dtc,outdir+"/fig_dist.r.%i.png"%i,title=str(i))

    m_xy = np.mean(alldist[:,:,0],1)
    #s_xy = np.std(tot_xy,1)
    m_z  = np.mean(alldist[:,:,1],1)
    #s_z  = np.std(tot_z,1)
    m_r  = np.mean(alldist[:,:,2],1)
    #s_r  = np.std(tot_r,1)

    fit_sqrt_vs_time(m_xy,dtc,outdir+"/fig_dist.xy.average.png",title="average of %i"%nfiles)
    fit_sqrt_vs_time(m_z, dtc,outdir+"/fig_dist.z.average.png",title="average of %i"%nfiles)
    fit_sqrt_vs_time(m_r, dtc,outdir+"/fig_dist.r.average.png",title="average of %i"%nfiles)

    print_output_D_ave(allD)

def print_output_D_ave(allD):
    print("="*20)
    print("Diffusion estimates")
    print(allD.shape)
    #print "xy,  z,  r"
    #print allD

    print("="*5)
    print("xy")
    print(allD[:,0]/4.)  # 2D
    print("z")
    print(allD[:,1]/2.)  # 1D
    print("r")
    print(allD[:,2]/6.)  # 3D

    print("="*5)
    print("Diffusion constant")
    print("   %20s   %20s" %("Dmean[e-4cm^2/s]","Dstd"))
    print("xy %20.10f   %20.10f" %(np.mean(allD[:,0])/4., np.std(allD[:,0])/4.))
    print("z  %20.10f   %20.10f" %(np.mean(allD[:,1])/2., np.std(allD[:,1])/2.))
    print("r  %20.10f   %20.10f" %(np.mean(allD[:,2])/6., np.std(allD[:,2])/6.))
    print("="*5)

def print_output_D_ave_1D(allD):
    print("="*20)
    print("Diffusion estimates")
    print(allD.shape)
    #print allD

    print("="*5)
    print(allD[:]/2.)  # 1D

    print("="*5)
    print("Diffusion constant")
    print("   %20s   %20s" %("Dmean[e-4cm^2/s]","Dstd"))
    print("r  %20.10f   %20.10f" %(np.mean(allD[:])/2., np.std(allD[:])/2.))
    print("="*5)


#=========================================

def calc_dist_fromdir(dirname):
    x = np.array(read_coor(dirname+"/smoothcrd-x/out"))
    y = np.array(read_coor(dirname+"/smoothcrd-y/out"))
    z = np.array(read_coor(dirname+"/smoothcrd-z/out"))
    dist_xy,dist_z,dist_r = calc_dist(x,y,z)
    return dist_xy,dist_z,dist_r

def calc_dist_from3files(xfile,yfile,zfile):
    x = np.array(read_coor(xfile))
    y = np.array(read_coor(yfile))
    z = np.array(read_coor(zfile))
    dist_xy,dist_z,dist_r = calc_dist(x,y,z)   # use ALL data in the file!
    return dist_xy,dist_z,dist_r

def calc_dist_1D(x,):   # kind of oversampling!
    # coor format: x[timestemp,atom]
    nstep = len(x)  # number of time steps
    natom = x.shape[1]  # number of atoms
    alldist  = np.zeros((nstep,natom),float)

    for lt in range(1,nstep):  # lt is lagtime
        # fill in alldist[i]
        diff  = abs(x[:-lt,:] - x[lt:,:])
        alldist[lt,:]  = np.mean(diff,axis=0)   #z/(self.nstep-lt)

    return alldist

def calc_dist(x,y,z):   # kind of oversampling!
    # coor format: x[timestemp,atom]
    nstep = len(x)  # number of time steps
    natom = x.shape[1]  # number of atoms
    alldist_xy = np.zeros((nstep,natom),float)
    alldist_z  = np.zeros((nstep,natom),float)

    for lt in range(1,nstep):  # lt is lagtime
        # fill in alldist[i]
        diff_xy = (x[:-lt,:] - x[lt:,:])**2 + (y[:-lt,:] - y[lt:,:])**2
        diff_z  = abs(z[:-lt,:] - z[lt:,:])
        alldist_xy[lt,:] = np.mean(diff_xy,axis=0)  #xy2/(self.nstep-lt)
        alldist_z[lt,:]  = np.mean(diff_z,axis=0)   #z/(self.nstep-lt)

    dist_xy = np.sqrt(alldist_xy)
    dist_z  = alldist_z
    dist_r  = np.sqrt(alldist_xy+alldist_z**2)
    return dist_xy,dist_z,dist_r

#=========================================

