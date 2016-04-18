#!/usr/bin/env python
#

import numpy as np
import sys
import scipy
from scipy import linalg, special

#=============================
# some testing at bottom of file


def rad_log_like_lag(dim_trans,dim_rad, num_lag, rate, wrad, lagtimes, transition,
            rad,lmax,bessel0_zeros,bessels, epsilon ):
    """calculate log-likelihood summed over different lag times"""
    tiny = np.empty((dim_rad,dim_trans,dim_trans))
    tiny.fill(1.e-32) # lower bound of propagator (to avoid NaN's)
    log_like = np.float64(0.0)
    # add contributions of different lag times
    for ilag in range(num_lag):
        # use elementwise maximum with tiny to avoid NaN errors
        lnpropagator =  np.log(np.maximum(propagator_radial_diffusion(dim_trans,dim_rad,
                          rate,wrad,lagtimes[ilag],lmax,bessel0_zeros,bessels),tiny))
        # sum up log likelihood
        log_like += np.sum( transition[ilag,:,:,:] * lnpropagator[:,:,:] )

    # smoothness prior for log(D)
    if (epsilon > 0.0):
        for i in range(num_bin-1):
            log_like += - ( ( w[i+1] - w[i] ) / epsilon ) ** 2 

    return log_like

def setup_bessel_functions(lmax,redges,):
    """set up Bessel functions first type zero-th order J_0(b_l x)
    redges  --  bin edges for radial bins, starts with 0.
    dim_w  --  number of transitions between bins (size diffusion vector)   TODO
    rmax  --  radius where all functions are zero
    lmax  --  number of Bessel functions taken into account, each labeled by the zero by which the argument is rescaled
    Output
    bessels  --  in units 1/dr"""

    # first lmax zeros of 0th order Bessel first type
    bessel0_zeros = scipy.special.jn_zeros(0,lmax)   # no unit
    # 1st order Bessel first type in those zeros
    bessel1_inzeros = scipy.special.j1(bessel0_zeros)   # no unit

    dim_rad = len(redges)  # transitions between radial bins
    rmax = np.float64(dim_rad)  # in units [dr]   # TODO dimensions!!!!!!!!
    bessels = np.zeros((lmax,dim_rad),dtype=np.float64)

    r = np.arange(dim_rad,dtype=np.float64)+0.5   # in units [dr]
    #print "r",r
    #print "rmax",rmax

    # set up Bessel functions
    for l in range(lmax):
        bessels[l,:] = 2*r*scipy.special.j0(r/rmax*bessel0_zeros[l]) / bessel1_inzeros[l]**2 /rmax**2
        # in units r [dr] / rmax**2 [dr**2], so in units [1/dr]
    if False:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.plot(bessels.transpose())
        plt.savefig("bessels.%i.png"%lmax)

    return bessel0_zeros,bessels

def propagator_radial_diffusion(n,dim_rad,rate,wrad,lagtime,
           lmax,bessel0_zeros,bessels,):
    """calculate propagator for radial diffusion as matrix exponential
    rate  --  rate matrix for 1-D diffusion in z-direction, in [1/dt]
    wrad  --  ln Drad, radial diffusion coefficient, dimension n
              Drad = exp(wrad), in [dr**2/dt]
    n  --  dim_trans, dimension transition matrix, usually number of bins in z-direction
    dim_rad  --  dimension transition matrix, always equal to len(redges)
    bessels0_zeros  --  first lmax zeros, no unit
    bessels  --  dimension lmax x dim_rad, in units [1/dr]

    rate_l  --  rate matrix including sink equation, in [1/dt]
    propagator  --  in 1/dr (I think)"""  # TODO unit double check?

    rmax = np.float64(dim_rad)  # in units [dr]

    # initialize arrays
    rate_l = np.zeros((n,n),dtype=np.float64)   # N x N
    propagator = np.zeros((dim_rad,n,n),dtype=np.float64) # dim_rad x N x N
    # set up sink
    sink = np.zeros((n),dtype=np.float64)   # N
    # loop over l (index of Bessel function)
    for l in range(lmax):
        sink = np.exp(wrad)*bessel0_zeros[l]**2/rmax**2 # sink term D_par(i) * (b_l)**2
        # in units np.exp(wrad) [dr**2/dt] / rmax**2 [dr**2], so in units [1/dt]

        rate_l[:,:] = rate[:,:]                 # take rate matrix for 1-D diffusion
        rate_l.ravel()[::n+1] -= sink           # and add sink term
        mat_exp = scipy.linalg.expm2(lagtime*rate_l) # matrix exponential

        # increment propagator by solution of sink equation for each l
        # propagator to arrive in radial bin k, contribution from Bessel function l
        # bessels in [1/dr], mat_exp is "kind of" [1/dz], so propagator in [1/dr/dz]
        # so propagator in [1/dr] (the 1/dz is more of a probability/bin, so gives no unit)
        for k in range(dim_rad):
            propagator[k,:,:] += bessels[l,k] * mat_exp[:,:]
    # TODO normalize?
    #propagator /= np.sum(np.sum(propagator,axis=0),axis=0)
    return propagator

#=============================
# TESTING
#=============================
if False:
    from mcdiff.utils import init_rate_matrix
    n = 30
    D = 1.   # in angstrom**2/ps

    dt = 1.  # ps
    dz = 1.  # angstrom
    wunit = np.log(dz**2/dt)
    redges = np.arange(0,50.1,0.5)
    dim_rad = len(redges)
    #rmax = max(redges)
    dr = redges[1]-redges[0]
    wradunit = np.log(dr**2/dt)

    w = np.zeros((n),float)
    w[0] = np.log(D)-wunit
    v = np.zeros((n),float)
    rate = init_rate_matrix(n,v,w,True)

    wrad = np.zeros((n),float)
    wrad[0] = np.log(D)-wradunit

    lagtime = 10.
    lmax = 50

    bessel0_zeros,bessels = setup_bessel_functions(lmax,redges,) #rmax=rmax)
    propagator = propagator_radial_diffusion(n,dim_rad,rate,wrad,lagtime,
           lmax,bessel0_zeros,bessels,) #rmax=rmax)
    print propagator

if False:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for i in range(n):
        plt.figure()
        plt.plot(propagator[:,:,i])
        plt.title("in start z-bin %i" %i)
        plt.savefig("propagator.3-%i.png"%i)
    for i in range(n):
        plt.figure()
        plt.plot(propagator[:,i,:].transpose())
        plt.title("in end z-bin %i" %i)
        plt.savefig("propagator.2-%i.png"%i)
    for i in range(dim_rad):
        plt.figure()
        plt.plot(propagator[i,:,:].transpose())
        plt.title("in end r-bin %i" %i)
        plt.savefig("propagator.1-%i.png"%i)
    for i in range(dim_rad):
        plt.figure()
        plt.contourf(propagator[i,:,:].transpose())
        plt.title("in end r-bin %i" %i)
        plt.savefig("propagator.1s-%i.png"%i)
    for i in range(n):
        plt.figure()
        plt.contourf(propagator[:,i,:].transpose())
        plt.title("in end z-bin %i" %i)
        plt.savefig("propagator.2s-%i.png"%i)
    for i in range(n):
        plt.figure()
        plt.contourf(propagator[:,:,i].transpose())
        plt.title("in start z-bin %i" %i)
        plt.savefig("propagator.3s-%i.png"%i)


# check normalization
if False:
    print "Normalization"
    a = np.sum(propagator,axis=0)  #*dr # already in units 1/dr/dz
    print "sum-axis-0",a
    b = np.sum(a,axis=0)   # *dz # already in units dz
    print "sum-axis-0-0",b
    print np.sum(b)
