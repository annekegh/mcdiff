#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# An Ghysels (August 2012)
#

import numpy as np
import scipy
import scipy.linalg
import numpy.linalg


#------------------------
# UNITS
# F -- free energy, in kBT
# D -- diffusion constant, in A^2/ps
# v -- equal to F, in kBT
# w -- log( D/(dx^2/dt) ), has no unit
# exp(w) -- D/(dx^2/dt), has no unit, amounts to D in units dx^2/dt
# rate -- in 1/dt
# lagtime -- in ps or in dt   # TODO CONFIRM  ook elders nakijken
#------------------------

#------------------------
# EXTRA FUNCTIONS
#------------------------

def init_rate_matrix(n,v,w,pbc,st=None,end=None,side=None):
    # st -- absorbing or reflective bin to the left
    # end -- absorbing or reflective bin to the right
    # side -- which side is absorbing (both, left, right)
    if pbc:
        if st is not None or end is not None or side is not None:
            print "you are asking too much:"
            print "asking for rate matrix with pbc=True AND with absorption/reflection in bins=",st,end
            raise NotImplemented
        return init_rate_matrix_pbc(n,v,w,)  # PBC

    else:
        if st is None and end is None and side is None:
            rate = init_rate_matrix_nopbc(n,v,w)  # NOPBC, left=right=reflective
        else:
            rate = init_rate_matrix_pbc(n,v,w)  # PBC   # because I want correct corner elements
            rate[-1,0] = 0.     #remove PBC
            rate[0,-1] = 0.     #remove PBC

            if st is None:
                st = -1
            if end is None:
                end = n
            assert st >= -1 and st <n-1
            assert end <=n  and end >1
            assert st+1<end

            rate = rate[st+1:end,st+1:end]   # taking a cut
            # sides are automatically absorbing when taking a cut
            if side not in ["left","both"]:
                rate[0,0] = -rate[1,0]      # make left reflective
            if side not in ["right","both"]:
                rate[-1,-1] = -rate[-2,-1]  # make right reflective
        return rate


def init_rate_matrix_pbc(n,v,w):
    """initialize rate matrix from potential vector v and diffusion
    vector w = log(D(i)/delta^2)"""
    assert len(v) == n  # number of bins
    assert len(w) == n
    rate = np.float64(np.zeros((n,n)))  # high precision

    # off-diagonal elements
    diffv = v[1:]-v[:-1] #length n-1  # diffv[i] = v[i+1]-v[i]
    exp1 = w[:n-1]-0.5*diffv
    exp2 = w[:n-1]+0.5*diffv
    rate.ravel()[n::n+1] = np.exp(exp1)[:n-1]
    rate.ravel()[1::n+1] = np.exp(exp2)[:n-1]
    #this amounts to doing:
    #for i in range(n-1):
    #    rate[i+1,i] = np.exp(w[i]-0.5*(v[i+1]-v[i]))
    #    rate[i,i+1] = np.exp(w[i]-0.5*(v[i]-v[i+1]))

    # corners    # periodic boundary conditions
    rate[0,-1]  = np.exp(w[-1]-0.5*(v[0]-v[-1]))
    rate[-1,0]  = np.exp(w[-1]-0.5*(v[-1]-v[0]))
    rate[0,0]   = - rate[1,0] - rate[-1,0]
    rate[-1,-1] = - rate[-2,-1] - rate[0,-1]

    # diagonal elements
    for i in range(1,n-1):
        rate[i,i] = - rate[i-1,i] - rate[i+1,i]
    return rate

def init_rate_matrix_nopbc(n,v,w):
    """initialize rate matrix from potential vector v and diffusion
    vector w = log(D(i)/delta^2)"""
    assert len(v) == n  # number of bins
    assert len(w)+1 == n
    rate = np.float64(np.zeros((n,n)))  # high precision

    # off-diagonal elements
    diffv = v[1:]-v[:-1] #length n-1  # diffv[i] = v[i+1]-v[i]
    exp1 = w[:n-1]-0.5*diffv
    exp2 = w[:n-1]+0.5*diffv
    rate.ravel()[n::n+1] = np.exp(exp1)[:n-1]
    rate.ravel()[1::n+1] = np.exp(exp2)[:n-1]
    #this amounts to doing:
    #for i in range(n-1):
    #    rate[i+1,i] = np.exp(w[i]-0.5*(v[i+1]-v[i]))
    #    rate[i,i+1] = np.exp(w[i]-0.5*(v[i]-v[i+1]))

    # corners   # reflecting boundaries (is equal to a hard wall)
    rate[0,0]   = - rate[1,0]
    rate[-1,-1] = - rate[-2,-1]

    # diagonal elements
    for i in range(1,n-1):
        rate[i,i] = - rate[i-1,i] - rate[i+1,i]
    return rate

##### CONSTRUCT #####
# in the following:
##### UNITS #####
# F -- in kBT
# D -- in angstrom**2/ps
# edges -- in angstrom

# TODO
# rate*lagtime
# lagtime in ps
# rate in 1/dt
# so do *dt or /dt somewhere 

# input: F and D profiles
# output: rate matrix, propagator, extra matrices

def construct_rate_matrix_from_F_D(F,D,dx,dt,pbc=True,st=None,end=None,side=None):
    """compute the rate matrix"""
    # dx in angstrom
    # dt in ps
    # rate in 1/dt
    # F in units kBT, D in units A^2/ps
    from mcdiff.utils import init_rate_matrix
    n = len(F)
    wunit = np.log(dx**2/dt)
    w = np.log(D)-wunit     # w has no real unit
    v = F                   # v in unit kBT

    if pbc or st is not None or end is not None or side is not None:
        rate = init_rate_matrix(n,v,w,pbc,st=st,end=end,side=side)  # PBC or ABSORBING
    else:
        rate = init_rate_matrix(n,v,w[:-1],False,st=st,end=end,side=side)  # NOPBC

    # check: verify normalization
    #print "norm prop: np.sum(norm**2)", np.sum(np.sum(rate,0)**2)

    return rate   # in 1/dt

def construct_propagator_from_F_D(F,D,Drad,dz,dr,dt,lagtime,lmax,dim_rad,pbc=True,
    st=None,end=None,side=None,power=None,radialintegral=None):
    # F -- in kBT
    # D -- in unit=angstrom**2/ps, D/unit has no dimension
    #      convert to w = ln(D/unit)
    # dz, dr -- angstrom
    # dt -- in ps
    # dim_rad -- len(redges), where redges is [0.,dr,2dr,..]
    # propagator  --  no unit, is per r-bin per z-bin
    # normal
    # TODO radialintegral=2 for <r**2> but nothing is done with this option
    rate = construct_rate_matrix_from_F_D(F,D,dz,dt,pbc=pbc,st=st,end=end,side=side)   # in 1/dt
    n = len(rate)

    # radial
    from mcdiff.twod import setup_bessel_functions, propagator_radial_diffusion
    wradunit = np.log(dr**2/dt)
    wrad = np.log(Drad)-wradunit
    if st is None: st=-1
    if end is None: end=len(F)
    wrad = wrad[st+1:end]

    bessel0_zeros,bessels = setup_bessel_functions(lmax,dim_rad)
    if radialintegral is not None:
        from mcdiff.permeability.msd import setup_bessel_functions_integral
        besselsintegral = setup_bessel_functions_integral(lmax,dim_rad)

    if power is None:
      if radialintegral is None:
        propagator = propagator_radial_diffusion(n,dim_rad,rate,wrad,lagtime,
           lmax,bessel0_zeros,bessels,)   # no unit, is per r-bin per z-bin
        return propagator
      else: raise NotImplementedError

    elif power in [-1,-2]:
      if radialintegral is None:
        from mcdiff.permeability.msd import timeintegral_propagator_radial_diffusion
        mat = timeintegral_propagator_radial_diffusion(n,dim_rad,rate,wrad,lagtime,
           lmax,bessel0_zeros,bessels,power=power)   # in 1/time**power per r-bin per z-bin
        return mat
      else:
        from mcdiff.permeability.msd import radialintegral_propagator_radial_diffusion
        mat = radialintegral_propagator_radial_diffusion(len(rate),dim_rad,rate,wrad,lmax,bessel0_zeros,
                   besselsintegral,power)
        return mat

#=================================

# this is now obsolete!!!
# this differs from the other definition, where bins -1 and N are absorbing, by taking a cut
def make_rate_absorbing_old(rate,side):
    if side not in ["left","right"]:
        raise ValueError("side is not known:",side)
    if side == "left":
        rate[1,0] = 0.    # absorbing in bin 0, nothing goes from 0 to 1
    elif side == "right":
        rate[-2,-1] = 0.  # absorbing in bin N-1, nothing goes from N-1 to N-2
    return rate

def string_energy(vec,k,pbc):
    if False:
        v = np.exp(vec)   # TODO this is as in paper GH
    else:
        v = vec
    diff = v[1:]-v[:-1]
    energy = k/2.*np.sum(diff**2)
    if pbc:
        energy += k/2.*(v[0]-v[-1])**2
    return energy

def string_vecs(n,pbc):
    M = np.diag(np.ones((n),float)*2)
    M.ravel()[1::n+1] = -1.
    M.ravel()[n::n+1] = -1.
    if pbc:
        M[0,-1] = -1.
        M[-1,0] = -1.
    else:
        M[0,0] = 1.
        M[-1,-1] = 1.
    vals,vecs = np.linalg.eigh(M)
    return vecs

def log_likelihood(n,ilag,transition,lagtime,rate):
    """calculate log-likelihood from rate matrix and transition matrix
    assuming time step lagtime"""
    # lagtime -- in dt  # TODO confirm elders
    # rate -- rate matrix, in 1/dt
    # calc propagator as matrix exponential
    propagator = scipy.linalg.expm2(lagtime*rate)
    # sum over all transitions
    # GH: in case of numerical issues try: np.log(np.abs(propagator[i,j]))
    log_like = np.float64(0.0)  # use high precision
    tiny = 1e-10

    # PUT CUT-OFF
    b = transition[ilag,:,:]*np.log(propagator.clip(tiny))
    val = np.sum(b)
    #print propagator
    #count = np.sum(propagator<tiny)
    #vals,vecs = np.linalg.eig(propagator)
    #line = ""
    #for v in vals:
    #    if v.imag < 1e-10: VAL=v.real
    #    line += str(VAL)+" "
    #print "val",val, count, line

    # PUT NAN
    #b = np.where(transition[ilag,:,:]>0,
    #        np.where(propagator>-tiny, transition[ilag,:,:] * np.log(abs(propagator)), float("Nan")),
    #        0)
    #val = np.sum(b)
    # only sum over non-zero transition values (kind of unuseful if using numpy arrays
    #a = np.where(transition[ilag,:,:]>0, transition[ilag,:,:] * np.log(propagator), 0)
    # this should work
    #a = np.where(propagator>tiny, transition[ilag,:,:] * np.log(abs(propagator)), 0)
    # this is not giving any warnings
    #a = np.sum(np.ma.masked_invalid(transition[ilag,:,:] * np.ma.log(abs(propagator))))

    # INDUCE CRASH whenever propagator gets negative
    #np.seterr(invalid='raise')
    #try:
    #    b = np.log(propagator)
    #except:
    #    return None
    #val = np.sum(b)

    log_like += val
    return log_like


def log_like_lag(num_bin,num_lag, v,w,lagtimes,transition, pbc):
    """calculate log-likelihood summed over all umbrella windows"""
    log_like = np.float64(0.0)
    rate = init_rate_matrix(num_bin,v,w,pbc)

    for ilag in range(num_lag):
        # add several contributions
        ll = log_likelihood(num_bin,ilag,transition,lagtimes[ilag],rate)
        if ll is None:
            return None
        else:
           log_like+=ll

    return log_like


# TODO
def calc_overlap_basis(p_basis):
    L = len(p_basis)
    overlap = np.zeros((L,L),float)
    for i in xrange(L):
        for j in np.arange(j,L):
            a = np.sum(p_basis[:,i]*p_basis[:,j])
            overlap[i,j] = a
            overlap[j,i] = a
    return overlap

def calc_rhs_diffusionequation_basisfunctions(p_basis,rate):
    L = len(p_basis)
    rhs = np.zeros((L,L),float)
    for i in xrange(L):
        for j in xrange(L):
            a = np.sum(p_basis[:,i]*p_basis[:,j]*rate[:,j])  # I can do this better  TODO might be wrong
            rhs[i,j] = a
    return rhs


def log_likelihood_basis(transition,p_basis,rate,lagtime):
    overlap = calc_overlap_basis(p_basis)
    rhs = calc_rhs_diffusionequation_basisfunctions(p_basis,rate)
    vals,vecs = np.linalg.eig(np.dot(np.inv(overlap),rhs))

    # lagtime:
    prop = np.sum(vecs*np.exp(lagtime*vals),axis=X)

    pass
    # TODO


