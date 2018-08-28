"""Script to compare different methods to extract profiles:
calc flux
AG, April 9, 2013
AG, April 25, 2013
AG, Jan 11, 2016"""

import numpy as np
import mcdiff
from mcdiff.outreading import read_F_D_edges, read_Drad
from mcdiff.utils import construct_rate_matrix_from_F_D
import matplotlib.pyplot as plt

##### UNITS #####
# F -- in kBT
# D -- in angstrom**2/ps
# edges -- in angstrom

# TODO
# rate*lagtime
# lagtime in ps
# rate in 1/dt
# so do *dt or /dt somewhere 

#################################
# average D, effective D
#################################

def calc_Dave_notmidF(F,D,st=None,end=None,edges=None):
    """Compute average diffusion constant, in parallel layers, e.g. Drad
         Dave = sum_i  D_i  exp(-F_i) ,

    where D_i is actually D_(i+1/2)
    not shifting F nor D to the middle of the bins

    using bins [st:end] = [st,st+1,...,end-1]
    number of used bins = end-st
    when st=None,end=None, then this means [:]"""
    # F in units kBT, D in units A^2/ps
    assert len(F) == len(D)
    v = F[st:end]     # when st=None,end=None, then this means [:]
    d = D[st:end]

    # average Drad:
    part = np.exp(-(v-min(v)))
    Dave = np.sum(d*part)/np.sum(part)   # a weighted sum, normalized

    # for some printing:  # I should make it nicer
    if edges is not None:
        assert len(edges)>0
        x_st = edges[st]
        x_end = edges[end]
        h = edges[end]-edges[st]
        # print for postprocessing
        print "st,end %3i %3i"%(st,end),
        print "st,end,h %7.2f %7.2f %7.2f"%(x_st,x_end,h),
    #print "Dave",Dave
    return Dave

def calc_Dave_midF(F,D,st=None,end=None):
    """Compute average D, parallel routes, take mid-values of F
    See calc_Dave_notmidF
         Dave = sum_i  D_i  exp(-F_i) ,

    where D_i is actually D_(i+1/2)
    where F_i is now also F_(i+1/2)"""
    # F in units kBT, D in units A^2/ps
    assert len(F) == len(D)
    assert len(F) > 1 # I cannot take average F if only one value
    # value F at edge of bin
    Fmid = (0.5*(F[:-1]+F[1:])).tolist() + [0.5*(F[0]+F[-1])]   # I assume here periodic boundary conditions!
    Fmid = np.array(Fmid)
    v = Fmid[st:end]     # when st=None,end=None, then this means [:]
    d = D[st:end]
    return calc_Dave_notmidF(v,d)

def calc_Deff_1(F,D,st=None,end=None,ave=0,dz=None):
    """Compute effective D, serial, e.g. layered medium

       Deff**(-1) = sum_i ( D_i exp(-F_i) )**(-1)

  TODO and norm

    This is a formula from literature: Berendsen, Attila Szabo, ... """
    # F in units kBT, D in units A^2/ps
    assert len(F) == len(D)
    assert len(F) > 1 # I cannot take average F if only one value
    # value F at edge of bin
    Fmid = (0.5*(F[:-1]+F[1:])).tolist() + [0.5*(F[0]+F[-1])]
    Fmid = np.array(Fmid)
    #v = Fmid[st:end]     # when st=None,end=None, then this means [:]
    d = D[st:end]

    # independent of reference F
    v  = Fmid[st:end]
    v -= min(v)    # for numerical stability
    part = np.exp(-v)
    if ave == 1:
        Deff = np.sum(part)/np.sum(part/d)   # TODO which???    # TODO I used to compute this like this one !!!! TODO TODO ALERT
    elif ave == 2:
        Deff = 1./np.sum(part)/np.sum(1./(d*part))
    elif ave == 0:
        #h = len(part)*dz
        N = len(part)
        Deff = N**2/ np.sum(part)/ np.sum(1./(d*part))   # this would be most correct???

    elif ave == -1:
        h = len(part)*dz
        N = len(part)
        Deff = N/ np.sum(1./(d*part))

    # use reference Fmid[0]
    v  = Fmid[st:end]-F[0]   # reference is water bin
    part = np.exp(-v)
    if ave == 3:
        Deff = len(part)/np.sum(part/d)     # this one makes most sense, I guess,
                                            # but it is somewhat strange to have the reference OUTSIDE, to need a reference even.
    elif ave == 4:
        Deff = 1./len(part)/np.sum(1./(d*part))

    # use reference min(Fmid), in practice this is F-middle-membrane
    v = Fmid[st:end]
    v = v-min(v)
    #print min(v)  apparantly min(v) is zero the whole time,
    # so this is the same as just Fmid[st:end]
    part = np.exp(-v)
    if ave == 5:
        Deff = len(part)/np.sum(part/d)
    elif ave == 6:
        Deff = 1./len(part)/np.sum(1./(d*part))

    if dz is not None:
        h = dz*len(part)
        P = Deff/h
        print "P ave",ave,P,

    # effective D
    #print "Deff_1",Deff
    return Deff

def calc_Deff_2(F,D,dx,dt):
    """Compute effective D through series of layers, formula Ghysels 2016"""
    # F in units kBT, D in units A^2/ps
    assert len(F) == len(D)
    assert len(F) > 1  # can't take Fmidst otherwise
    # effective D
    # 1) part average D
    Dave = calc_Dave_midF(F,D)
    # 2) part with rate matrix
    Dinh = calc_Dinh(F,D,dx,dt)
    return Dave + Dinh


def calc_Dinh(F,D,dx,dt):
    """Compute inhomogeneous contribution to Deff"""
    # compute the rate matrix
    rate = construct_rate_matrix_from_F_D(F,D,dx,dt)  # in 1/dt
    # compute pseudo-inverse
    # symmetrize and pseudo-inverse
    a = np.exp(-(F-min(F))/2.)
    M = np.diag(a)   #/ np.sum(a**2)
    M1 = np.diag(1./a) #* np.sum(a**2)
    rateS = np.dot(np.dot(M1,rate),M)      # symmetrized rate, in 1/dt
    rateS1 = np.linalg.pinv(rateS,1.e-7)   # pseudo-inverse, in dt
    #print "rateS",rateS[:4,:4]

    # gradient-like vector
    gradD = np.zeros(D.shape)
    for m in range(-1,len(gradD)-1):
        gradD[m] = D[m]*a[m+1] - D[m-1]*a[m-1]   # in units A^2/ps
    #print "gradD",gradD

    # combine
    # double product is in ((A^2)/ps)**2 * dt = A^4/ps^2*dt
    Dinh = np.dot(np.dot(gradD,rateS1),gradD)/np.sum(a**2)
    Dinh *= dt    # put dt in units A
    Dinh /= dx**2 # finally, should be divided by dx**2 to arrive at MSD/(2t) interpretation
    #print "Dinh",Dinh
    return Dinh

#====================================
# TRYING
# Jan 2018
#====================================

# TODO clean up implementation try-out formula of Hummer

def calc_Dpar_logs(F,D,Drad,dz,dt,st=None,end=None,edges=None):
    """Compute average diffusion constant, in parallel layers, e.g. Drad
         Dave = sum_i  D_i  exp(-F_i) ,

    where D_i is actually D_(i+1/2)
    not shifting F nor D to the middle of the bins

    using bins [st:end] = [st,st+1,...,end-1]
    number of used bins = end-st
    when st=None,end=None, then this means [:]"""
    # F in units kBT, D in units A^2/ps
    assert len(F) == len(D)
    #v = F[st:end]     # when st=None,end=None, then this means [:]
    #d = D[st:end]
    #drad = Drad[st:end]

    # average Drad:
    #part = np.exp(-(v-min(v)))
    #Dave = np.sum(d*part)/np.sum(part)   # a weighted sum, normalized
    drad = Drad[st+1:end]
    # absorbing boundaries
    rate = construct_rate_matrix_from_F_D(F,D,dz,dt,pbc=False,st=st,end=end,side="both")  # in 1/dt
    #rate = construct_rate_matrix_from_F_D(F,D,dz,dt,pbc=False,st=st,end=end,side=None)  # in 1/dt
    print "rate",rate.shape
    vals,U = np.linalg.eig(rate)  # U contains eigenvectors as columns
    print vals
    U1 = np.linalg.inv(U)
    Ddiag = np.diag(drad)
    UDdiagU1 = np.dot(U,np.dot(Ddiag,U1))

    "o-ow, I have three positive eigenvals?????"

    Dpar = 0.
    # first term
    Dpar_1 = 0.
    for i in range(len(vals)):
        Dpar_1 += 1./vals[i]*U1[-1,i]*UDdiagU1[i,i]*U[i,0]
    print "part1",Dpar_1
    # second term
    Dpar_2 = 0.
    for i in range(len(vals)):
        for j in range(len(vals)):
              if i != j:
                  term = np.log(abs(vals[i]/vals[j]))/(vals[i]-vals[j])
                  term *= U1[-1,i]*UDdiagU1[i,j]*U[j,0]
                  Dpar_2 += term
    print "part2",Dpar_2
    Dpar = Dpar_1 + Dpar_2
    print "sum",Dpar
    rate1 = np.linalg.inv(rate)
    Dpar /= rate1[-1,0]
    print "final",Dpar
    return Dpar


def calc_Dpar_ratios(F,D,Drad,t,dz,dt,st=None,end=None,edges=None):
    """Compute average diffusion constant, in parallel layers, e.g. Drad
         Dave = sum_i  D_i  percentage-time-spent-in-bin-i

    where D_i is actually D_(i+1/2)

    using bins [st:end] = [st,st+1,...,end-1]
    number of used bins = end-st
    when st=None,end=None, then this means [:]"""
    # F in units kBT, D in units A^2/ps
    assert len(F) == len(D)
    v = F[st:end]     # when st=None,end=None, then this means [:]
    d = D[st:end]
    drad = Drad[st:end]

    # average Drad:
    #part = np.exp(-(v-min(v)))
    #Dave = np.sum(d*part)/np.sum(part)   # a weighted sum, normalized

    # ow.... watch out !!!
    drad = Drad[st+1:end]

    # absorbing boundaries
    rate = construct_rate_matrix_from_F_D(F,D,dz,dt,pbc=False,st=st,end=end,side="both")  # in 1/dt
    # reflective boundaries
    #rate = construct_rate_matrix_from_F_D(F,D,dz,dt,pbc=False,st=st,end=end,side=None)  # in 1/dt

    #print "rate",rate.shape
    n = len(rate)
    vals,U = np.linalg.eig(rate)  # U contains eigenvectors as columns
    #print "vals",vals
    U1 = np.linalg.inv(U)
    Ddiag = np.diag(drad)
    #UDdiagU1 = np.dot(U,np.dot(Ddiag,U1))
    U1DdiagU = np.dot(U1,np.dot(Ddiag,U))

    maxval = max(vals)
    exp_valt = np.exp((vals-maxval)*t)
    Deff_t = np.zeros((n,n),)
    for i in range(n):
        for j in range(n):
          #  for k in range(n):
          #      Deff_t[i,j] += drad[k] * np.sum(     # summation over l
          #           U[k,:]*U1[:,j]*U[i,:]*U1[:,k] * t * exp_valt )
            Deff_t[i,j] = np.sum(U1[:,j]*U[i,:]* np.diag(U1DdiagU) * t * exp_valt )
#    print Deff_t
    for i in range(n):
        for j in range(n):
                for l in range(n):
                  for m in range(n):
                    if m != l:
                       Deff_t[i,j] += U1[l,j]*U[i,m]*U1DdiagU[m,l] * (exp_valt[l]-exp_valt[m])/(vals[l]-vals[m])

#    print Deff_t
    import scipy.linalg
    Deff_t /= scipy.linalg.expm2((rate-np.diag(np.ones(n))*maxval)*t) * t
    print "t",t
    print Deff_t
    return Deff_t

