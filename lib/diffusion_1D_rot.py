#!/usr/bin/env python
#
# $Id: diffusion_1D_rot.py,v 1.3 2012/08/30 15:29:40 hummer Exp $
# $Log: diffusion_1D_rot.py,v $
# Revision 1.3  2012/08/30 15:29:40  hummer
# enhanced commenting
#
# Revision 1.2  2012/08/24 13:22:14  hummer
# made code more "pythonic" which resulted in significant speedup
#
# Revision 1.1  2012/08/23 20:58:49  hummer
# Initial revision
#
# python program to determine angular diffusion coefficient from transition counts
#
# It is assumed that the underlying dynamics is diffusive.  With this assumption,
# a discrete 1D diffusion model is fit to the data in each window using a Bayesian or
# maximum likelihood approach (depending on choice of temp).  In calculating the
# likelihood function, it is assumed that the dynamics (1) in the radial direction
# is diffusive with an input potential of mean force Vperp_i and diffusion coefficient
# Dperp_i, both of which are read in at start time.  It is further assumed that
# the dynamics in the angular direction is diffusive with a diffusion coefficient
# Dpar_i. The likelihood is determined based on given transition
# counts T(k,j,i) = (number of transitions from radial bin i to
# [radial bin j,angular bin k]). The output diffusion
# coefficient is in units of (deltax**2/timestep).
# Angular bins are uniform in [-1,1] corresponding to the cosine of the polar angle
# theta
# Prior enforcing smoothness is added if epsilon>0
# log_like_smooth = - Sum_i [ ( log(D_i)- log(D_{i-1}) ) / epsilon ] ** 2
#
# Input: 
# num_bin_rad ..... number of bins in radial direction
# numb_bin_ang .... number of bins in [-1,1] interval of cos(polar angle)
# r_V_D_perp.dat .. file containing radii (in units of radial bin width deltax!!!,
#                .. radial potential of mean force, and the logarithm wperp=log(Dperp)
#                .. of the diffusion coefficient into the radial direction (organized
#                .. as 3 entries "r_i V_i wperp_i" per line
# lagtime1 ........ time interval at which transition counts were determined
# transition1.dat . file containing transition counts for lagtime1 (see below)
# lagtime2 transition2.dat etc.
#
# IMPORTANT: transition matrix files are stored one entry per line such that
#   fastest index: over initial radial bin ( 1 ... num_bin_rad )
#   intermediate index: over final radial bin ( 1 ... num_bin_rad )
#   slowest index: over final angular bin ( 1 ... num_bin_rad )
#  T(1,1<--1)
#  T(1,1<--2)
#  ...
#  T(1,1<--num_bin_rad)
#  T(1,2<--1)
#  ...
#  T(1,2<--num_bin_rad)
#  ...
#  T(num_bin_ang,num_bin_rad<--num_bin_rad)
#
# parameters to adjust:
#    dw = 0.5                      # MC move width for wpar=log(Dpar)
#    temp = .01                    # effective temperature (1 -> Bayesian; ~0 -> MaxLikelihood)
#    nmc = 10000                   # number of MC moves
#    num_MC_update = 100           # number of moves between MC step width adjustments
#                                  #     (0 if no adjustment)
#    lmax = 25                     # maximum order of legendre polynomials used in expansion of propagator
#    epsilon = 0.0                 # "smoothness" parameter for log(D) (if 0, smoothness prior is not used)
# usage:
#
# python ./diffusion_1D_rot.py num_bin_rad numb_bin_ang r_V_D_perp.dat lagtime1 transition1.dat lagtime2 transition2.dat lagtime3 transition3.dat ...
#
# copyright: Gerhard Hummer (NIH, August 2012)
#

import numpy as np
from numpy import matrix
import sys
import scipy
from scipy import linalg, special

def main():
    """main program to determine angular diffusion coefficient from transition counts"""
    
    print "# python program to determine angular diffusion coefficient from transition counts"
    print "# copyright: Gerhard Hummer (NIH, August 2012)\n"


# PARAMETERS

    lmax = 50                     # maximum order of legendre polynomials used in expansion of propagator
    
    # Monte Carlo parameters
    dw = 0.5                      # w=log(D) MC move width
    temp = .01                    # effective temperature used in MC sampling of posterior
    nmc = 1000                     # number of MC moves
    num_MC_update = 100           # number of moves between MC step width adjustments
                                  #     (0 if no adjustment)
    epsilon = 0.1                 # "smoothness" parameter for log(D) (if 0, smoothness prior is not used)

    print "# maximum order of legendre polynomials used in cos(theta) expansion:", lmax
    print "# initial MC move width of log(D):", dw
    print "# effective temperature used in Bayesian sampling:", temp
    print "# number of MC moves", nmc
    print "# smoothness parameter epsilon (if 0, no smoothness prior used):", epsilon
    if num_MC_update>0:
        print "# interval between updates of MC move width:", num_MC_update
    else:
        print "# no updates of MC move width"

# INPUT

    # number of radial histogram bins
    num_bin_rad = np.int32(sys.argv[1])
    print "# number of radial bins:", num_bin_rad

    # number of angular histogram bins
    num_bin_ang = np.int32(sys.argv[2])
    print "# number of angular bins:", num_bin_ang

    # number of lagtimes
    num_lag = ( len ( sys.argv ) - 4 ) / 2
    print "# number of lagtimes:", num_lag
    lagtimes = np.zeros((num_lag),dtype=np.float64)

    # transition matrix
    transition = np.zeros((num_bin_ang,num_bin_rad,num_bin_rad,num_lag),dtype=np.int)
    
# read radius, potential and diffusion coefficient
    rad = np.zeros(num_bin_rad,dtype=np.float64)        # radii at bin centers
    vperp = np.zeros(num_bin_rad,dtype=np.float64)      # potential of mean force in radial direction
    wperp = np.zeros(num_bin_rad,dtype=np.float64)      # log of diffusion coefficient in radial direction
    dvfile = open(sys.argv[3],"r")     # name of file containing rad, vperp and wperp
    i = 0
    for line in dvfile:
        (rad[i],vperp[i],wperp[i]) = line.split()
        i += 1
    dvfile.close()                     # close dvfile
    if i <> num_bin_rad:
        sys.stderr.write("error in V_D_perp.dat")
        quit()
        
# read transition matrices
    for i in range(4,len(sys.argv),2):
        j = i / 2 - 3
        print "# lag-time transition-matrix ", j, ":", sys.argv[i], sys.argv[i+1]
        lagtimes[j] = np.float64(sys.argv[i])
# read transition matrix line by line (one line per entry!!!)
#   fastest index: over starting radial bin ( 1 ... num_bin_rad )
#   intermediate index: over ending radial bin ( 1 ... num_bin_rad )
#   slowest index: over ending angular bin ( 1 ... num_bin_rad )
        tfile = open(sys.argv[i+1],"r")
        k = 0  # line counter 
        for line in tfile:
# extract bin indices from line counter
            ( k1, kk2 ) = divmod ( k , num_bin_rad ** 2 )
            ( k2, k3 )  = divmod ( kk2 , num_bin_rad )
            transition[k1,k2,k3,j] = np.float64(line)
            k += 1
        tfile.close()  # close transition matrix file
        if k <> num_bin_rad ** 2 * num_bin_ang:
            sys.stderr.write("wrong number of entries in %s\n"%(sys.argv[i+1]))
            quit()

# INITIALIZATION
    wpar = np.zeros(num_bin_rad,dtype=np.float64)  # log(Dpar)
    wt = np.zeros(num_bin_rad,dtype=np.float64)    # trial w
    naccw = np.float64(0.0)       # number accepted w moves
    naccw_update = np.float64(0.0)       # number accepted w moves between adjusts

# initialize rate matrix for radial diffusion (unchanged)
    rate = init_rate_matrix(num_bin_rad,vperp,wperp)

# initialize Legendre polynomials used in expansion of angular diffusion operator
    legendre_poly = setup_legendre_polynomials(lmax,num_bin_ang)

# initialize potential v[i] and wpar[i]=log(D(i)/deltax**2)
    wpar[:] = wperp[:]

# initial log-likelihood (before optimization)
    log_like = log_like_lag ( num_bin_rad, num_bin_ang, num_lag, rate, wpar, lagtimes, transition, rad, lmax, legendre_poly, epsilon )
    print "# initial log-likelihood: ", log_like
    sys.stderr.write("initial log-likelihood: %17.9e\n"%(log_like))

    sys.stderr.write("\n MC-move log-like acc(w)\n")

# loop over Monte Carlo moves
    for imc in range(nmc):

# diffusion move
        wt[:] = wpar[:]
        i = np.random.randint(0,num_bin_rad)
        wt[i] += dw * ( np.random.random() - 0.5 )
        log_like_try = log_like_lag ( num_bin_rad, num_bin_ang, num_lag, rate, wt, lagtimes, transition, rad, lmax, legendre_poly, epsilon )
        
        # Metropolis acceptance
        dlog = log_like_try - log_like
        r = np.random.random()
        if ( r < np.exp(dlog/temp) ):
            wpar[i] = wt[i]
            naccw += 1.
            naccw_update += 1.
            log_like = log_like_try
        
# print intermediate results
        if  (imc%100 == 0) | (imc == nmc-1):
             sys.stderr.write("%d %17.9e %13.5e\n"%(imc, log_like, naccw / ( imc + 1. )))
            
# update MC move widths
        if ( num_MC_update > 0 ):
            if ( ( imc + 1 ) % num_MC_update == 0 ):
                dw = dw * np.exp ( 0.1 * ( naccw_update / num_MC_update - 0.3 ) )
                naccw_update = 0.
                
                sys.stderr.write("new MC steps: %d %13.5e\n"%(imc, dw))

# print final results ( potential and diffusion coefficient)
    print "# final log-likelihood =", log_like
    print "# final MC move width of log(D):", dw
    print "#\n#     index  rad[i]/deltax  diffusion-coefficient"
    for i in range(num_bin_rad):
        sys.stdout.write("%8d %13.5e %13.5e\n"%(i,rad[i],np.exp(wpar[i])))
    return()

def init_rate_matrix(n,v,w):
    """initialize rate matrix from potential vector v and diffusion vector w = log(D(i)/deltax^2)"""
# initialize n x n matrix with zeros
    rate = np.zeros((n,n),dtype=np.float64)
# off-diagonal elements
    for i in range(n-1):
        rate[i+1,i] = np.exp(w[i]-0.5*(v[i+1]-v[i]))
        rate[i,i+1] = np.exp(w[i]-0.5*(v[i]-v[i+1]))
# diagonal elements: columns sum to zero        
    rate[0,0]     = - rate[1,0]
    rate[n-1,n-1] = - rate[n-2,n-1]
    for i in range(1,n-1):
        rate[i,i] = - rate[i-1,i] - rate[i+1,i]

    return(rate)

def log_like_lag ( num_bin_rad, num_bin_ang, num_lag, rate, w, lagtimes, transition, rad,lmax,legendre_poly, epsilon ):
    """calculate log-likelihood summed over differen lag times"""
    tiny = np.empty((num_bin_ang,num_bin_rad,num_bin_rad))
    tiny.fill(1.e-32) # lower bound of propagator (to avoid NaN's)
    log_like = np.float64(0.0)
# rate matrix for radial diffusion
    for ilag in range(num_lag):
# use elementwise maximum with tiny to avoid NaN errors
        propagator =  np.log(np.maximum(propagator_angular_diffusion(num_bin_rad,num_bin_ang,rate,w,rad,lagtimes[ilag],lmax,legendre_poly),tiny))
# sum up log likelihood
        log_like += np.sum( transition[:,:,:,ilag] * propagator[:,:,:] )
# smoothness prior for log(D)
    if (epsilon > 0.0):
        for i in range(num_bin_rad-1):
            log_like += - ( ( w[i+1] - w[i] ) / epsilon ) ** 2

    return ( log_like )

def setup_legendre_polynomials(lmax,num_bin_ang):
    """set up legendre polynomials P_l(x_k)"""
    legendre_poly = np.zeros((lmax+1,num_bin_ang),dtype=np.float64)
    x =  np.zeros((num_bin_ang),dtype=np.float64)
# set up cos(theta) grid
    for k in range(num_bin_ang):
        x[k] = -1. + 2. * ( k + 0.5 ) / np.float64(num_bin_ang) # evenly spaced bins in [-1,1]
# set up P_l(x_k) Legendre polynomial coefficients
    for l in range(lmax+1):
        legendre_poly[l,:] = 0.5 * ( 2. * l + 1. ) * scipy.special.legendre(l)(x[:])
    return(legendre_poly)

def propagator_angular_diffusion(num_bin_rad,num_bin_ang,rate,w,rad,lagtime,lmax,legendre_poly):
    """calculate propagator for angular diffusion as matrix exponential"""
# initialize arrays
    rate_l = np.zeros((num_bin_rad,num_bin_rad),dtype=np.float64)
    propagator = np.zeros((num_bin_ang,num_bin_rad,num_bin_rad),dtype=np.float64)
    sink = np.zeros((num_bin_rad),dtype=np.float64)
# set up sink
    sink[:] = np.exp(w[:]) / rad[:]**2
# loop over l (index of legendre polynomials)
    for l in range(lmax+1):             # sum over legendre polynomials 0,...,lmax
        rate_l[:,:] = rate[:,:]         # take rate matrix for radial diffusion
        for i in range(num_bin_rad):    # and add sink term D_par(i) * l * ( l + 1 ) / r(i)^2
            rate_l[i,i] = rate[i,i] - l * ( l + 1 ) * sink[i]
        mat_exp = scipy.linalg.expm2(lagtime*rate_l) # matrix exponential
# increment propagator by solution of sink eqution for each l
        for k in range(num_bin_ang):
            propagator[k,:,:] += legendre_poly[l,k] * mat_exp[:,:]

    return(propagator)

if __name__ == '__main__':
    main()
