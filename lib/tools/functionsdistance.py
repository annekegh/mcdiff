"""Script to plot histogram of xyz coords
AG, August 13, 2012"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from functionshistogram import read_coor
from functionsdiffusion import fit_sqrt_vs_time, plot_vs_time


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
# notations
# dtc  --  time interval between saved coordinates
#=====================


#=====================
# ANALYZE
#=====================

def analyze_dist_multi_1D(list_x,nstep,outdir,dtc):
    "Take the average over the D estimates in each of these files"""
    assert len(list_x) > 0
    nfiles = len(list_x)
    t = np.arange(0,nstep*dtc,dtc)

    alldist2 = np.zeros((nstep,nfiles),float)
    allD = np.zeros((nfiles),float)

    print "="*5
    print "Results"
    print "making %i fits" %nfiles
    print "fit from %i steps, i.e. time %f ps, actually time %f ps" %(nstep,nstep*dtc,(nstep-1)*dtc)
    for i in range(nfiles):
        dist2 = calc_dist_1D(np.array(read_coor(list_x[i])))
        d2 = dist2[:nstep,:]   # only fit the first nstep time steps (use all atoms)
        average = np.mean(d2,1)
        alldist2[:,i] = average

        p = np.polyfit(t,average,1)
        allD[i] = p[0]   # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                                # = 1e-8 meter**2/second = 1e-4 cm**2/s
        # if I want many figures:
        #fit_sqrt_vs_time(np.mean(dist_xy,1),dtc,outdir+"/fig_dist.xy.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_z,1), dtc,outdir+"/fig_dist.z.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_r,1), dtc,outdir+"/fig_dist.r.%i.png"%i,title=str(i))

    m2 = np.mean(alldist2,axis=1)
    #s2 = np.std(alldist2,axis=1)

    fit_sqrt_vs_time(np.sqrt(m2),dtc,outdir+"/fig_dist.average.png",title="average of %i"%nfiles)
    print_output_D_ave_1D(allD)



def analyze_dist_multi(list_x,list_y,list_z,nstep,outdir,dtc):
    "Take the average over the D estimates in each of these files"""
    assert len(list_x) > 0
    assert len(list_x) == len(list_y)
    assert len(list_x) == len(list_z)
    nfiles = len(list_x)
    t = np.arange(0,nstep*dtc,dtc)

    alldist2 = np.zeros((nstep,nfiles,3),float)
    allD = np.zeros((nfiles,3),float)

    print "="*5
    print "Results"
    print "making %i fits" %nfiles
    print "fit from %i steps, i.e. time %f ps, actually time %f ps" %(nstep,nstep*dtc,(nstep-1)*dtc)
    for i in range(nfiles):
        dist2_xy,dist2_z,dist2_r,weight = calc_dist_from3files(list_x[i],list_y[i],list_z[i])
        for it,dist2 in enumerate([dist2_xy,dist2_z,dist2_r]):
            d2 = dist2[:nstep,:]   # only fit the first nstep time steps (use all atoms)
            average = np.mean(d2,1)
            alldist2[:,i,it] = average
            #print i,list_x[i],it,average[:10]

            p = np.polyfit(t,average,1)
            allD[i,it] = p[0]   # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                                # = 1e-8 meter**2/second = 1e-4 cm**2/s
        # if I want many figures:
        #fit_sqrt_vs_time(np.mean(dist_xy,1),dtc,outdir+"/fig_dist.xy.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_z,1), dtc,outdir+"/fig_dist.z.%i.png"%i,title=str(i))
        #fit_sqrt_vs_time(np.mean(dist_r,1), dtc,outdir+"/fig_dist.r.%i.png"%i,title=str(i))

    m2_xy = np.mean(alldist2[:,:,0],1)
    #s_xy = np.std(tot_xy,1)
    m2_z  = np.mean(alldist2[:,:,1],1)
    #s_z  = np.std(tot_z,1)
    m2_r  = np.mean(alldist2[:,:,2],1)
    #s_r  = np.std(tot_r,1)

    fit_sqrt_vs_time(np.sqrt(m2_xy),dtc,outdir+"/fig_dist.xy.average.png",title="average of %i"%nfiles)
    fit_sqrt_vs_time(np.sqrt(m2_z), dtc,outdir+"/fig_dist.z.average.png",title="average of %i"%nfiles)
    fit_sqrt_vs_time(np.sqrt(m2_r), dtc,outdir+"/fig_dist.r.average.png",title="average of %i"%nfiles)

    print_output_D_ave(allD)

def print_output_D_ave(allD):
    print "="*20
    print "Diffusion estimates"
    print allD.shape
    #print "xy,  z,  r"
    #print allD

    print "="*5
    print "xy"
    print allD[:,0]/4.  # 2D
    print "z"
    print allD[:,1]/2.  # 1D
    print "r"
    print allD[:,2]/6.  # 3D

    print "="*5
    print "Diffusion constant"
    print "   %20s   %20s" %("Dmean[e-4cm^2/s]","Dstd")
    print "xy %20.10f   %20.10f" %(np.mean(allD[:,0])/4., np.std(allD[:,0])/4.)
    print "z  %20.10f   %20.10f" %(np.mean(allD[:,1])/2., np.std(allD[:,1])/2.)
    print "r  %20.10f   %20.10f" %(np.mean(allD[:,2])/6., np.std(allD[:,2])/6.)
    print "="*5

def print_output_D_ave_1D(allD):
    print "="*20
    print "Diffusion estimates"
    print allD.shape
    #print allD

    print "="*5
    print allD[:]/2.  # 1D

    print "="*5
    print "Diffusion constant"
    print "   %20s   %20s" %("Dmean[e-4cm^2/s]","Dstd")
    print "r  %20.10f   %20.10f" %(np.mean(allD[:])/2., np.std(allD[:])/2.)
    print "="*5


#=========================================


def analyze_dist_multi_BIS1(list_x,list_y,list_z,dtc):
    "Collect distances in all these files"""
    assert len(list_x) > 0
    assert len(list_x) == len(list_y)
    assert len(list_x) == len(list_z)
    nfiles = len(list_x)
    #nfiles = 4

    print "="*5
    print "Results"
    #alllagtimes = []
    #alldist_xy = []
    #alldist_z = []
    #alldist_r = []
    allD = np.zeros((nfiles,3),float)

    for i in range(nfiles):
        lagtimes,dist2_xy,dist2_z,dist2_r = collect_dist_from3files(list_x[i],list_y[i],list_z[i],dtc)
        print "file",i,"shape",dist2_xy.shape
        for it,dist2 in enumerate([dist2_xy,dist2_z,dist2_r]):

            p = np.polyfit(lagtimes,dist2,1)
            allD[i,it] = p[0]   # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                                # = 1e-8 meter**2/second = 1e-4 cm**2/s
            #alllagtimes.append(lagtimes)
            #alldist_xy.append(dist_xy)
            #alldist_z.append(dist_z)
            #alldist_r.append(dist_r)

    #lagtimes = np.array(alllagtimes).ravel()
    #dist_xy = np.array(alldist_xy).ravel()
    #dist_z = np.array(alldist_z).ravel()
    #dist_r = np.array(alldist_r).ravel()

    print_output_D_ave(allD)

    #import matplotlib.cm as cm
    #plt.figure()
    #plt.hexbin(lagtimes,dist_r,bins='log', cmap=cm.jet)
    #cb = plt.colorbar()
    #cb.set_label('log10(N)')
    #plt.savefig("lag_vs_r.png")



def collect_dist_from3files(xfile,yfile,zfile,dtc):
    x = np.array(read_coor(xfile))
    y = np.array(read_coor(yfile))
    z = np.array(read_coor(zfile))
    lagtimes,dist2_xy,dist2_z,dist2_r = collect_dist(x,y,z,dtc)
    return lagtimes,dist2_xy,dist2_z,dist2_r

def collect_dist(x,y,z,dtc):
    lagtimes = []
    dist2_xy = []
    dist2_z = []
    dist2_r = []
    natom = x.shape[1]
    for lt in range(1,len(x)):  # lt is lagtime   #TODO I adapted this not
        diff_xy = (x[:-lt,:] - x[lt:,:])**2 + (y[:-lt,:] - y[lt:,:])**2
        diff_z  = (z[:-lt,:] - z[lt:,:])**2
        lagtimes.extend(lt*np.ones(len(diff_xy)*natom)*dtc)
        dist2_xy.extend(diff_xy.ravel().tolist())
        dist2_z.extend(diff_z.ravel().tolist())
        #print lt, len(diff_xy), len(dist_xy)
    lagtimes = np.array(lagtimes)
    dist2_xy = np.array(dist2_xy)
    dist2_z = np.array(dist2_z)
    dist2_r = dist2_xy+dist2_z
    return lagtimes,dist2_xy,dist2_z,dist2_r


#=========================================

def analyze_dist_timedependent(list_x,list_y,list_z,dtc,nstep,figname,zoom=50,shift=1):
    "Collect distances in all these files"""
    assert len(list_x) > 0
    assert len(list_x) == len(list_y)
    assert len(list_x) == len(list_z)
    nfiles = len(list_x)

    nstep = 1000
    print "="*5
    print "Results"

    allddr = np.zeros((nstep,nfiles,6),float)  # six of those    
    for i in range(nfiles):
        ddr,weight = calc_timedependent_D_from3files(list_x[i],list_y[i],list_z[i],dtc,shift=shift)
        for it,ddist in enumerate(ddr):
            allddr[:,i,it] = np.mean(ddr[it],axis=1)
            print "file",i,"type",it,"ddr",ddr[it].shape

    for it,label in enumerate(["xy1","xy2","z1","z2","r1","r2"]):
        dim = [2,2,1,1,3,3,][it]
        plot_vs_time(np.mean(allddr[:,:,it],axis=1) / (2.*dim) ,dtc,figname+".type-%s.shift"%label +str(shift)+".png")

def calc_timedependent_D_fromdir(dirname,dtc,figname,zoom=100,shift=1):
    x = np.array(read_coor(dirname+"/smoothcrd-x/out"))
    y = np.array(read_coor(dirname+"/smoothcrd-y/out"))
    z = np.array(read_coor(dirname+"/smoothcrd-z/out")) 
    ddr,weight = calc_timedependent_D(x,y,z,dtc,shift=shift)
    meanddr = np.mean(ddr,axis=1)  # average over atoms
    plot_vs_time(meanddr,dtc,figname+".png")
    plot_vs_time(meanddr[:zoom],dtc,figname+".zoom.png")
    return ddr,weight

def calc_timedependent_D_from3files(xfile,yfile,zfile,dtc,shift=1):
    x = np.array(read_coor(xfile))
    y = np.array(read_coor(yfile))
    z = np.array(read_coor(zfile))
    ddr,weight = calc_timedependent_D(x,y,z,dtc,shift=shift)
    return ddr,weight

def calc_timedependent_D(x,y,z,dtc,shift=1):
    """
    First approach:
        Calculate [dr**2(t+shift) - dr**2(t)] / shift as a function of shift,
        averaged over available times
    Second approach:
        Calculate [dr**2(t+shift) - dr**2(t)] / shift as a function of shift,
        averaged over available shifts"""
    nstep = x.shape[0]  # number of time steps
    natom = x.shape[1]  # number of atoms

    dist2_xy,dist2_z,dist2_r,weight = calc_dist(x,y,z)   # includes averaging

    ddr = []

    for it,dist2 in enumerate([dist2_xy,dist2_z,dist2_r]):
        # dist2 is Delta r**2(t), dimension nstep x natom
        #label = ["xy","z","r"][it]
        #plot_vs_time(dr2,dtc,"dr2.%s.png"%label)

        # first approach
        ddr1 = np.zeros((nstep,natom),float)
        ddr1[shift:,:] = (dist2[shift:,:]-dist2[:-shift,:]) / (shift*dtc)
        ddr1[0,:] = ddr1[1,:]  # do trick to fix first timestep

        # second approach
        ddr2 = np.zeros((nstep,natom),float)
        for sh in range(1,nstep):
            ddr2[sh,:] = np.mean(dist2[sh:,:]-dist2[:-sh,:],axis=0) / (sh*dtc)
        ddr2[0,:] = ddr2[1,:]  # do trick to fix first timestep

        # third approach
        #ddr3 = np.zeros((nstep,natom),float)
        #for shift in range(1,nstep):
        #    lt = np.arange(dtc,(nstep-shift+1)*dtc,dtc)
        #    ddr3[shift,:] = np.mean((dr2[shift:,:]-dr2[:-shift,:]).transpose() / (6.*lt),axis=1) 
        #ddr3[0,:] = ddr3[1,:]  # do trick to fix first timestep
        ddr.extend([ddr1,ddr2])

    return ddr,weight  # this is [ddxy1,ddxy2,ddz1,ddz2,ddr1,ddr2]


#=========================================

def calc_dist_lt_from3files(xfile,yfile,zfile,lt,rv=False):
    x = np.array(read_coor(xfile,rv=rv,axis=0))
    y = np.array(read_coor(yfile,rv=rv,axis=1))
    z = np.array(read_coor(zfile,rv=rv,axis=2))
    d2_xy,d2_z,d2_r = calc_dist_lt(x,y,z,lt)
    return d2_xy,d2_z,d2_r

def calc_dist_lt(x,y,z,lt):
    natom = x.shape[1]
    dist2_xy = (x[:-lt,:] - x[lt:,:])**2 + (y[:-lt,:] - y[lt:,:])**2
    dist2_z  = (z[:-lt,:] - z[lt:,:])**2
    dist2_r = dist2_xy + dist2_z
    return dist2_xy,dist2_z,dist2_r

def analyze_dist_conditionalz(list_x,list_y,list_z,dtc,zpbc,lt,figname,rv=False):
    "Collect distances in all these files"""
    assert len(list_x) > 0
    assert len(list_x) == len(list_y)
    assert len(list_x) == len(list_z)
    nfiles = len(list_x)

    nstep = 1000
    print "="*5
    print "Results"

    alldist2_xy = []
    alldist2_z = []
    alldist2_r = []
    allz_init = []
    allz_final = []
    for i in range(nfiles):
        dist2_xy,dist2_z,dist2_r = calc_dist_lt_from3files(list_x[i],list_y[i],list_z[i],lt,rv=rv)
        z = np.array(read_coor(list_z[i],rv=rv))

        alldist2_xy.append(dist2_xy)
        alldist2_z.append(dist2_z)
        alldist2_r.append(dist2_r)

        zinit = z[:-lt]    # initial z
        zinit -= zpbc*np.floor(zinit/zpbc+0.5)
        allz_init.append(zinit)

        zfinal = z[lt:]    # final z
        zfinal -= zpbc*np.floor(zfinal/zpbc+0.5)
        allz_final.append(zfinal)
    analyze_dist_conditionalz_essence(alldist2_xy,alldist2_z,alldist2_r,allz_init,allz_final,lt,dtc,figname)

def analyze_dist_conditionalz_essence(alldist2_xy,alldist2_z,alldist2_r,allz_init,allz_final,lt,dtc,figname):

    print "XXXXXX lt*dtc",lt*dtc
    if lt*dtc < 0.6:
        distmax = 2.
    elif lt*dtc < 6.:  #ps
        distmax = 5.
    elif lt*dtc < 60.:
        distmax = 10.
    else:
        distmax = 20.

    for c,allz in enumerate([allz_init,allz_final]):
      print "="*10
      labelz = ["zinit","zfinal"][c]
      print labelz
      for i,dist2 in enumerate([alldist2_xy,alldist2_z,alldist2_r]):
        label = labelz+"-"+["xy","z","r"][i]
        factor = [np.sqrt(2),1.,np.sqrt(3)][i]
        k,xedges,yedges = np.histogram2d(np.array(np.sqrt(dist2)).ravel(),np.array(allz).ravel(),[80,40])
        #print dist2_r.shape,z.shape
        print "hist,xedges,yedges",k.shape, xedges.shape, yedges.shape
        for n,donorm in enumerate([False,True]): #enumerate([1.,np.sum(k,axis=0)]):
            Label = label+"-%s" %(["nonorm","norm"][n])
            if i == 2:  # r
                kk = (xedges[1:]+xedges[:-1]).reshape((-1,1))/2.*k
            else:
                import copy
                kk = copy.deepcopy(k)
            if donorm:  # normalized
                norm = np.sum(kk,axis=0)
                kk = kk / norm
            import matplotlib.pyplot as plt
            plt.figure()
            # factor rescaling xedges/factor
            plt.contourf(xedges[:-1],yedges[:-1],kk.transpose()) #locator=ticker.LogLocator()

            # expectation value
            E = np.sum((xedges[1:]+xedges[:-1]).reshape((-1,1))/2.*kk,axis=0)/np.sum(kk,axis=0)
            plt.plot(E,yedges[:-1],color="k")
            #D = E**2 / 2. /factor**2/lt*10 ########################################################3
            #plt.plot(D,yedges[:-1],color="grey",lw=2)

            plt.xlabel("distance [A]")
            plt.xlim(xmax=distmax)
            plt.ylabel("z [A]")
            plt.title("root mean square distance: "+Label)
            plt.colorbar()
            plt.savefig(figname+"_hist2d.%s.lt%i.png"%(Label,lt))


#=========================================


def calc_dist_fromdir(dirname):
    x = np.array(read_coor(dirname+"/smoothcrd-x/out"))
    y = np.array(read_coor(dirname+"/smoothcrd-y/out"))
    z = np.array(read_coor(dirname+"/smoothcrd-z/out"))
    dist2_xy,dist2_z,dist2_r,weight = calc_dist(x,y,z)
    return dist2_xy,dist2_z,dist2_r,weight

def calc_dist_from3files(xfile,yfile,zfile,rv=False):
    x = np.array(read_coor(xfile,rv=rv,axis=0))
    y = np.array(read_coor(yfile,rv=rv,axis=1))
    z = np.array(read_coor(zfile,rv=rv,axis=2))
    dist2_xy,dist2_z,dist2_r,weight = calc_dist(x,y,z)
    return dist2_xy,dist2_z,dist2_r,weight

def calc_dist_1D(x,):   # kind of oversampling!
    # coor format: x[timestemp,atom]
    nstep = len(x)  # number of time steps
    natom = x.shape[1]  # number of atoms
    dist2  = np.zeros((nstep,natom),float)

    for lt in range(1,nstep):  # lt is lagtime
        # fill in alldist[i]
        diff  = (x[:-lt,:] - x[lt:,:])**2
        dist2[lt,:]  = np.mean(diff,axis=0)   #z/(self.nstep-lt)

    return dist2

def calc_dist(x,y,z):
    # kind of oversampling!
    # Use ALL data in the files !!!
    # coor format: x[timestemp,atom]
    nstep = x.shape[0]  # number of time steps
    natom = x.shape[1]  # number of atoms
    dist2_xy = np.zeros((nstep,natom,),float)
    dist2_z  = np.zeros((nstep,natom,),float)
    weight  = np.zeros((nstep),float)

    for lt in range(1,nstep):  # lt is lagtime
        diff2_xy = (x[:-lt,:] - x[lt:,:])**2 + (y[:-lt,:] - y[lt:,:])**2
        diff2_z  = (z[:-lt,:] - z[lt:,:])**2
        dist2_xy[lt,:] = np.mean(diff2_xy,axis=0)  # xy2/(self.nstep-lt)
        dist2_z[lt,:]  = np.mean(diff2_z,axis=0)   #z/(self.nstep-lt)
        weight[lt] = len(diff2_xy)

    dist2_r = dist2_xy+dist2_z
    return dist2_xy,dist2_z,dist2_r,weight

