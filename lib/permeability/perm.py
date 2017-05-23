"""Script to compare different methods to extract profiles:
calc flux
AG, April 9, 2013
AG, April 25, 2013
AG, Jan 11, 2016
AG, Sept 2016"""

#lt = 10

import numpy as np
import matplotlib.pyplot as plt
import mcdiff
from mcdiff.outreading import read_F_D_edges, read_Drad
from mcdiff.permeability.deff import calc_Dave_midF, calc_Dave_notmidF
from analyzeprofiles import construct_rate_matrix_from_F_D

##### UNITS #####
# F -- in kBT
# D -- in angstrom**2/ps
# edges -- in angstrom

# st, end -> using bins [st,st+1,...,end-1]
# this means that number of used bins = end-1-st+1 = end-st

# TODO
# rate*lagtime
# lagtime in ps
# rate in 1/dt
# so do *dt or /dt somewhere 

#################################
##### PERMEABILITIES        #####
#################################

### 1) NORMAL PERMEABILITY

def calc_permeability(F,D,dx,st,end,edges=None,ref=0,doprint=False):
    # F in units kBT, D in units A^2/ps, dx in units A
    # choose reference:
    # ref = 0 is default #### ASSUME BULK value at bin 0
    # ref = st would be another plausible assumption

    """ from bin st until bin end, using bins [st,st+1,...,end-1]
    number of used bins = end-1 - st + 1 = end-st (typical Python way of slicing)
    edges[st] is between bin [st-1] and bin [st]
    edges[end] is between bin [end-1] and bin [end]
    h = edges[end]-edges[st] = length of covered bins = height
    e.g. all bins: st=0, end=len(F)
    """
    if edges is None:
        edges = np.arange(len(F)+1)*dx   # is accurate if dx is accurate...
    h = edges[end]-edges[st]  # angstrom, is accurate
    #print "st,end",st,end

    aveD = [(D[0]+D[-1])/2.] + ((D[1:]+D[:-1])/2.).tolist()  # assume PBC
    aveD = np.array(aveD)    # angstrom**2/ps

    Fref = F[ref]
    part = np.exp(-(F-Fref))      # no units, partition function, not normalized
    dRdx = 1./(part*aveD)         # permeation resistance per unit length, unit ps/A^2
    R = np.sum(dRdx[st:end]) *dx  # resistance to permeability, unit ps/A
                                  # integrate from x=-25.5=edges[st] to x=25.5=edges[end]
    P = 1./R                      # permeability, unit A/ps
    Deff = h*P                    # effective D, in A**2/ps

    # effective length
    #heff = D[ref]/P  # in units A
    heff = aveD[ref]/P  # in units A   TODO
    diff_h = heff-h  # in units A     # this looks weird on graph??
    ratio_h = heff/h
    if doprint:
        ##print "st,end %3i %3i"%(st,end),
        print "st,end,h %7.2f %7.2f %7.2f"%(edges[st],edges[end],h),
        print "P",P, "Deff",Deff, "heff",heff,"R",R, "dRdx",dRdx[st]
    return P,Deff,heff,diff_h

### 2) RADIAL PERMEABILITY

def calc_permeability_radial_h(F,Drad,dx,radius,st,end,edges=None,doprint=True):
    # F in units kBT, D in units A^2/ps, dx in units A
    """ from bin st until bin end
    see function calc_permeability
    """
    if edges is None:
        edges = np.arange(len(F)+1)*dx
    h = edges[end]-edges[st]  # angstrom

    # naieve:   (I might not have the correct end bin)
    #P,Deff,R = calc_permeability_radial(F[st:end],Drad[st:end],dx,radius)

    #Deff = calc_Dave_midF(F,Drad,st=st,end=end)   # average Drad, unit A^2/ps  # TODO
    Deff = calc_Dave_notmidF(F,Drad,st=st,end=end)   # average Drad, unit A^2/ps
    P = Deff/radius                               # permeability, unit A/ps
    R = radius/Deff                               # permeability resistance, unit ps/A

    if doprint:
        print "st,end,h %7.2f %7.2f %7.2f"%(edges[st],edges[end],h),
        print "P",P, "Deff",Deff, "R",R

    return P,Deff,R

def calc_permeability_radial(F,Drad,dx,radius):
    #Deff = calc_Dave_midF(F,Drad)   # average Drad, unit A^2/ps    # TODO
    Deff = calc_Dave_notmidF(F,Drad)   # average Drad, unit A^2/ps
    P = Deff/radius                 # permeability, unit A/ps
    R = radius/Deff                 # permeability resistance, unit ps/A
    return P,Deff,R

########
### 3) OXYGEN TRANSPORT PARAMETER

def calc_oxygentransportparameter(F,D,edges):
    """compute W ~ D*exp(-F)"""
    # F in units kBT, D in units A^2/ps
    aveD = [(D[0]+D[-1])/2.] + ((D[1:]+D[:-1])/2.).tolist()
    aveD = np.array(aveD)
    part = np.exp(-(F-min(F)))
    edges_mid = (edges[:-1]+edges[1:])/2.
    product = aveD*part     # in unit A^2/ps
    return product,edges_mid


### 4) EQUILIBRIUM DISTRIBUTION PROFILE WHEN MEASURING PERMEABILITY

def calc_permeability_distribution(F,D,dx,dt,st,end,figname=None):
    # F in units kBT, D in units A^2/ps, dx in units A
    rate = construct_rate_matrix_from_F_D(F,D,dx,dt)   # in 1/dt, PBC

    rate = rate[st:end,st:end]

    # solve subrate*p = vec
    subrate = rate[1:-1,1:-1]
    # rhs of equation 
    vec = np.zeros(len(subrate),float)
    vec[0] = -rate[1,0]
    # p = subrate^-1 * p
    r1 = np.linalg.inv(subrate)
    p = np.dot(r1,vec)
    # plug in vector of original size, with boundaries
    P = np.zeros(len(rate),float)
    P[1:-1] = p
    P[0] = 1.
    #print "prob",P
    if figname is not None:
        plt.figure()
        plt.plot(P)
        plt.xlabel("bins")
        plt.ylabel("prob")
        plt.savefig(figname)
        print "file written...",figname
    return P


