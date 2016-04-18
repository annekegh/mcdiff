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
# rate -- in 1/dt   # TODO CONFIRM  ja, dit is zo.
# lagtime -- in ps or in dt   # TODO CONFIRM  ook elders nakijken
#------------------------

#------------------------
# EXTRA FUNCTIONS
#------------------------

def init_rate_matrix(n,v,w,pbc):
    if pbc:
        return init_rate_matrix_pbc(n,v,w)
    else:
        return init_rate_matrix_nopbc(n,v,w)

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


