#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# An Ghysels (August 2012)
#

import numpy as np
import scipy
import scipy.linalg


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
    import numpy.linalg
    vals,vecs = np.linalg.eigh(M)
    return vecs

def log_likelihood(n,ilag,transition,lagtime,rate):
    """calculate log-likelihood from rate matrix and transition matrix
    assuming time step lagtime"""
    # calc propagator as matrix exponential
    propagator = scipy.linalg.expm2(lagtime*rate)
    # sum over all transitions
    # in case of numerical issues try: np.log(np.abs(propagator[i,j]))
    log_like = np.float64(0.0)
    a = np.where(transition[ilag,:,:]>0, transition[ilag,:,:] * np.log(abs(propagator)), 0)
    log_like += a.sum()  #use high precision
    return log_like
# TODO use tiny


def log_like_lag(num_bin,num_lag, v,w,lagtimes,transition, pbc):
    """calculate log-likelihood summed over all umbrella windows"""
    log_like = np.float64(0.0)
    rate = init_rate_matrix(num_bin,v,w,pbc)

    for ilag in range(num_lag):
        # add several distributions
        log_like += log_likelihood(num_bin,ilag,transition,lagtimes[ilag],rate)

    return log_like


