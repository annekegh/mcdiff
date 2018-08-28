#!/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from outreading import read_F_D_edges,read_many_profiles,read_many_profiles_Drad

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
#plotsettings()


import sys
assert len(sys.argv)>=3
system = sys.argv[1]
assert "/" not in system
list_filenames = sys.argv[2:]
assert len(list_filenames)>0
print list_filenames

figname = "figs-rescale/fig_"+system+"_"
# guess name
#figname="fig_"
#if len(list_filenames) > 0:
#    name = list_filenames[0]
#    if "HexdWat" in name: figname = "fig_HexdWat_"
#    elif "Hexd" in name:  figname = "fig_Hexd_"
#    elif "Wat" in name:   figname = "fig_Wat_"
#    elif "popc" in name:  figname = "fig_popc_"
#    elif "mito" in name:  figname = "fig_mito_"
#    elif "uwat" in name:  figname = "fig_uwat_"
#if "rad" in name:
    #figname += "rad_"
    #doradial = True
#else: doradial = False
if "rad" in system:
    doradial = True
else: doradial = False


def extract_lagtimes_from_filenames(list_filenames):
    # lagtimes?
    lts = []
    for filename in list_filenames:
        words = filename.split("/")
        parts = words[-1].split(".")
        #if "outrad" in filename:
        #    shift = int(parts[4])
        #else:
        shift = int(parts[2])
        #lts.extend([shift*dtc for i in range(nbins)])
        lts.extend([shift*dtc])
    lts = np.array(lts)
    return lts

dtc = 1.
lts = extract_lagtimes_from_filenames(list_filenames)

# read profiles
F,D,edges,Fst,Dst = read_many_profiles(list_filenames)   # Fst=None, Dst=None
if doradial:
    Drad,RE,Dradst = read_many_profiles_Drad(list_filenames)    #Dradst=None

F = np.array(F).T
D = np.array(D).T
if doradial:
    D = np.array(Drad).T    # !!!!!!! I am overwriting !!!!
edges = np.array(edges).T

print "Inputdata:"
print "F  ",F.shape
print "D  ",D.shape
print "lts",lts
# F has shape nbins x nprofiles
nbins = len(F)

#######


def fit_line_inverse(lt,D):
    """ FIRST
    a+2*Dreal*t = 2*D*t
    D = a/2/t + Dreal
    Dreal = true D, long lt
    D = the one I have, erronously

    Fitting
    D[i](lt) = A[i] + B[i]/lt
    D[i](lt) = A[i] * (1.+t0[i]/lt)
    D[i](lt) = (A[i]*lt + B[i]) / lt
    with two coefficients (profiles) A and B that are determined with least square fit
    A = Dreal
    B = A*t0, so t0 = B/Dreal
    """
    nlt = len(lts)    # number of lag times
    nbins = D.shape[0]
    A = np.zeros((nbins*nlt,2*nbins))   # coeffs
    res = np.zeros((nbins*nlt),float)   # results
    for i,lt in enumerate(lts):
        for j in range(nbins):
            res[i*nbins+j] = D[j,i]
            A[i*nbins+j,j] = 1.
            A[i*nbins+j,j+nbins] = 1./lt
    if False:
        print "A"
        #print "A",A
        print A[:6,:6]
        print A[-6:,-6:]
        print A[:6,-6:]
        #print "res",res

    print "--"*10
    print "Doing least square fit"
    print "A",A.shape, "res",res.shape
    x,residues,rank,s = np.linalg.lstsq(A,res)

    print "x",x.shape
    print "residues",residues.shape
    print "rank",rank
    #print "s",s
    print "--"*10
    # extract profiles from fit
    Dreal = x[:nbins]
    t0 = x[nbins:]/x[:nbins]

    # error analysis
    ## Calculate vector of residuals
    #res = as.matrix(women$weight-bh[1]-bh[2]*women$height)
    #res = Y-(bh[0]+X[:,1]*bh[1])
    # residues
    ## Define n and k parameters
    n = len(D)
    n = len(A)   # number of data points ???
    n = 4
    k = 2   # number of fitted parameters


    #data = np.array(list(df1))[1:,3:5].astype('float')
    #data = data with two columns
 #   data = np.concatenate(????,res,1)
 #   nrow = len(nbins)   #data.shape[0]
 #
 #   intercept = np.ones( (nrow,1) )
 #   b2 = data[:,0].reshape(-1, 1)
 #
 #   X = np.concatenate((intercept, b2), axis=1)
 #   Y = data[:,1].T
 #
 #   ## X and Y arrays must have the same number of columns for the matrix multiplication to work:
 #   print(X.shape)
 #   print(Y.shape)

    ## Calculate Variance-Covariance Matrix
 #   VCV = np.true_divide(1,n-k)*np.dot(np.dot(residues.T,residues),np.linalg.inv(np.dot(X.T,X)))
    VCV = np.true_divide(1,n-k)*np.dot(np.dot(residues.T,residues),np.linalg.inv(np.dot(A.T,A)))

    ## Standard errors of the estimated coefficients
    stderr = np.sqrt(np.diagonal(VCV))

    #print VCV.shape
    #print stderr

    return Dreal, t0

"""
def fit_line_gerhard(lt,D):
=> No, this suggestion is not okay

    Fitting SECOND
    D[i](lt) = 1 / (A[i] + B[i]*lt)
             = c[i] / (1.+d[i]*lt)   = Dreal / (1 + lt/t0)
    1/D[i](lt) = A[i] + B[i]*lt = 1./c[i] + lt * d[i]/c[i]
               = (1+lt/t0) / Dreal = 1/Dreal + lt / (t0*Dreal)
    with two coefficients (profiles) A and B that are determined with least square fit
    Dreal = 1/A
    B = 1/(t0*Dreal), so t0 = 1/B/Dreal = A/B


    nlt = len(lts)
    nbins = D.shape[0]
    A = np.zeros((nbins*nlt,2*nbins))
    res = np.zeros((nbins*nlt),float)
    for i,lt in enumerate(lts):
        for j in range(nbins):
            res[i*nbins+j] = 1./D[j,i]
            A[i*nbins+j,j] = 1.
            A[i*nbins+j,j+nbins] = lt

    print "--"*10
    print "Doing least square fit"
    print "A",A.shape, "res",res.shape
    x,residues,rank,s = np.linalg.lstsq(A,res)

    print "x",x.shape
    print "residues",residues.shape
    print "rank",rank
    #print "s",s
    print "--"*10
    # extract profiles from fit
    Dreal = 1./x[:nbins]
    t0 = x[:nbins]/x[nbins:]
    return Dreal, t0

=> No, this suggestion is not okay
 """


Dreal,t0 = fit_line_inverse(lts,D)



print "Dreal",Dreal.shape
print "Dreal[:10]",Dreal[:10]

##################################
# Do some plotting
##################################
xD = edges[:nbins,0]
print "--"*20
plt.figure()
plt.plot(xD,D,label="various lt")
plt.plot(xD,Dreal,lw=5,color='black',label="lt to infty (Dreal)")
plt.xlabel("z [A]")
plt.ylabel("D [A^2/ps]")
plt.legend()
plt.savefig(figname+"Dreal.png")

plt.figure()
plt.plot(xD,Dreal,lw=5,color='black',label="Dreal [A^2/ps]")
plt.plot(xD,t0,color='red',label="t0 [ps]")
plt.xlabel("z [A]")
plt.legend()
plt.savefig(figname+"timezero.png") 

##################################

plt.figure()
plt.plot(lts,D.T,"o-")
plt.xlabel("lt [ps]")
plt.ylabel("D [A^2/ps]")
plt.savefig(figname+"lts.png")


plt.figure()
plt.plot(1./lts,D.T,"o-")
plt.xlabel("1/lt [1/ps]")
plt.ylabel("D [A^2/ps]")
plt.savefig(figname+"lts1.png")

# do not plot all

from plot import plotsettings
plotsettings()
plotx = np.zeros(len(lts)+1)
plotx[-1] = 0.
plotx[:-1] = 1./lts
plotD = np.zeros((len(D),len(lts)+1))
plotD[:,-1] = Dreal
plotD[:,:-1] = D
plotD = (plotD.T)[:,:len(D)/2:5]  # has a transpose already # rows-> different lts, cols-> different bin
words = figname.split("_")
title = words[1].upper()


plt.figure()
plt.plot(plotx,plotD,"o-")
plt.plot(np.zeros(plotD.shape[1]),plotD[-1,:],"o",color='k',markersize=8)
#plt.plot(1./lts,(D.T)[:,:len(D)/2:5],"o-")   # subsample the first half of the profile
plt.rc(('xtick','ytick','axes'), labelsize=18.0)
plt.locator_params(axis='y',nbins=6) #to specify number of ticks on both or any single axes
#plt.locator_params(axis='x',nbins=10)
plt.xlabel(r"$1/t$ (1/ps)")
plt.ylabel(r"$D$ ($\AA^2$/ps)")
plt.xlim(xmin=0)
plt.title(title,fontsize=24)
plt.tight_layout()
plt.savefig(figname+"lts1.somelines.png")


##################################

msd = 2*D*lts
plt.figure()
plt.subplot(121)
plt.plot(lts,msd.T,"o-")
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.xlabel("lt [ps]")
plt.ylabel("msd [A^2]")
plt.subplot(122)
plt.loglog(lts,msd.T,"o-")
plt.xlabel("lt [ps]")
plt.ylabel("msd [A^2]")
plt.savefig(figname+"msd.png")

##################################


# other fitting

p = np.polyfit(lts,D.T,1)
#print p
p = np.polyfit(1./lts,D.T,1)
#print p

##################################
# Writing profiles
##################################

def print_profile(filename,v,edges,startline="None"):
    f = file(filename,"w+")
    if startline is not None:
        print >> f,startline
    for i in range(len(v)):
        print >> f, "%8d %13.5e %13.5e  %13.5f"%( i,edges[i],edges[i+1],v[i])
    print >> f, "="*10
    print "file written...",filename
    f.close()

def print_profile_2vecs(filename,v,d,edges,startline=None):
    assert len(v)==len(d)   # doing PBC
    f = file(filename,"w+")
    if startline is not None:
        print >> f,startline
    for i in range(len(v)):
        print >> f, "%8d %13.5e %13.5e  %13.5f  %13.5f"%( i,edges[i],edges[i+1],v[i],d[i])
    print >> f, "="*10
    print "file written...",filename
    f.close()

##################################

Dreal_file = figname+"Dreal.dat"
t0_file = figname+"t0.dat"

if not doradial:
    startline = "   index  bin-str  bin-end"
    print_profile_2vecs(Dreal_file,F[:,0],Dreal,edges[:,0],startline=startline)
else:
    startline = "   index  bin-str  bin-end  diffusion-coefficient-at[i]"
    print_profile(Dreal_file,Dreal,edges[:,0],startline=startline)

print_profile(t0_file,t0,edges[:,0])


# check if writing profiles worked well
#from outreading import read_F_D_edges, read_Drad
#F,D,edges = read_F_D_edges(Dreal_file)
#if doradial:
#    read_Drad(Dreal_file)


