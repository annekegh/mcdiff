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

def setup_bessel_functions(lmax,dim_rad):
    """set up Bessel functions first type zero-th order J_0(b_l x)
    Input
      lmax  --  number of Bessel functions taken into account, each labeled
                by the zero by which the argument is rescaled
      dim_rad  --  len(redges), number of radial bins
    Output
      bessel0_zeros  --  no unit
      bessels  --  no unit: in units per-r-bin"""

    # first lmax zeros of 0th order Bessel first type
    bessel0_zeros = scipy.special.jn_zeros(0,lmax)   # no unit
    # 1st order Bessel first type in those zeros
    bessel1_inzeros = scipy.special.j1(bessel0_zeros)   # no unit

    # rmax = dim_rad in high precision ~ radius where all functions are zero
    rmax = np.float64(dim_rad)
    r = np.arange(dim_rad,dtype=np.float64)+0.5
    #print "rmax",rmax   # this is float(dim_rad)
    #print "r",r   # this is [0,1,2,...dim_rad-1]
    #actually, this is [0.5,1.5,2.5,...,dim_rad-1+0.5] so always the middle of the dim_rad bins

    # set up Bessel functions
    bessels = np.zeros((lmax,dim_rad),dtype=np.float64)
    for l in range(lmax):
        bessels[l,:] = 2*r*scipy.special.j0(r/rmax*bessel0_zeros[l]) / bessel1_inzeros[l]**2 /rmax**2
        # in units r [dr] / rmax**2 [dr**2], so in units [1/dr]
        # so in units 'per r-bin' !

    return bessel0_zeros,bessels


def plot_bessels(lmax,redges,figname,color=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    #from mcdiff.plot import plotsettings
    #plotsettings()
    # construct
    dr = redges[1]-redges[0]
    dim_rad = len(redges)
    bessel0_zeros,bessels = setup_bessel_functions(lmax,dim_rad)
    # plot
    plt.plot(redges,np.zeros(len(redges)),color='grey') #,linewidth=2)
    if color is not None:
        plt.plot(redges,bessels.transpose(),color=color)
    else:
        plt.plot(redges,bessels.transpose())
    plt.xlabel("r [A]")
    plt.ylabel("per r-bin bin, with dr = %5.2fA"%dr)
    plt.title("bessel functions")
    plt.savefig("%s.lmax%i.png"%(figname,lmax))
    print("picture printed:...", figname)


def propagator_radial_diffusion(n,dim_rad,rate,wrad,lagtime,
           lmax,bessel0_zeros,bessels,):
    """calculate propagator for radial diffusion as matrix exponential
    n  --  dim_trans, dimension transition matrix, usually number of bins in z-direction
    dim_rad  --  dimension transition matrix, always equal to len(redges)
    rate  --  rate matrix for 1-D diffusion in z-direction, in [1/dt]
    wrad  --  ln Drad, radial diffusion coefficient, dimension n
              Drad = exp(wrad), in [dr**2/dt]
    lagtime  --  should be in units [dt]
    bessels0_zeros  --  first lmax zeros, no unit
    bessels  --  dimension lmax x dim_rad, no unit, in unit 'per r-bin'

    rate_l  --  rate matrix including sink equation, in [1/dt]
    propagator  --  no unit, is per r-bin per z-bin"""

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
        mat_exp = scipy.linalg.expm2(lagtime*rate_l) # matrix exponential, no unit

        # increment propagator by solution of sink equation for each l
        # propagator to arrive in radial bin k, contribution from Bessel function l
        # bessels is 'per r-bin', no unit
        # mat_exp is 'per z-bin', no unit
        # so propagator is 'per r-bin per z-bin', no unit
        for k in range(dim_rad):
            propagator[k,:,:] += bessels[l,k] * mat_exp[:,:]   # no unit

    # TODO normalize? some probability might flow away after long times
    #propagator /= np.sum(np.sum(propagator,axis=0),axis=0)
    return propagator

#=============================
# TESTING
#=============================
def test_propagator_twod():
    from mcdiff.utils import init_rate_matrix
    n = 30
    D = 1.   # in angstrom**2/ps

    dt = 1.  # ps
    dz = 1.  # angstrom
    dr = 0.5 # angstrom
    wunit = np.log(dz**2/dt)
    redges = np.arange(0,50.1,dr)  # in angstrom
    dim_rad = len(redges)
    #rmax = max(redges)
    wradunit = np.log(dr**2/dt)

    w = np.zeros((n),float)
    w[0] = np.log(D)-wunit
    v = np.zeros((n),float)
    rate = init_rate_matrix(n,v,w,True)

    wrad = np.zeros((n),float)
    wrad[0] = np.log(D)-wradunit

    lagtime = 10.
    lmax = 50

    bessel0_zeros,bessels = setup_bessel_functions(lmax,dim_rad)
    propagator = propagator_radial_diffusion(n,dim_rad,rate,wrad,lagtime,
           lmax,bessel0_zeros,bessels,)  # no unit, is 'per r-bin per z-bin'
    print("propagator-2d",propagator)

    # make plots
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
    # propagator is 'per end-r-bin per end-z-bin', no unit
    print("Normalization")
    a = np.sum(propagator,axis=0)  # is per end-z-bin
    print("sum-axis-0  ",a)
    b = np.sum(a,axis=0)   # is total probability (for each initial z-bin) 
    print("sum-axis-0-0",b)
    print("sum         ",np.sum(b))

