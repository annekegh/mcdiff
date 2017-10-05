"""Script to compare different methods to extract profiles:
calc flux
AG, April 9, 2013
AG, April 25, 2013
AG, Jan 11, 2016
AG, Sept 14, 2016"""

import numpy as np
import mcdiff
from mcdiff.outreading import read_F_D_edges, read_Drad
import matplotlib.pyplot as plt

##### UNITS #####
# F -- in kBT
# D -- in angstrom**2/ps
# edges -- in angstrom

# prop -- n x n, with n=number of z-bins
# prop2d -- dim_rad x n x n, with dim_rad=number of r-bins


# TODO
# rate*lagtime
# lagtime in ps
# rate in 1/dt
# so do *dt or /dt somewhere 

#################################
##### PROPERTIES: MSD #####     #
#################################


# functions
#----------

def radialintegral_propagator_radial_diffusion(n,dim_rad,rate,wrad,lmax,bessel0_zeros,
           besselsintegral,power):

    rmax = np.float64(dim_rad)

    # initialize arrays
    rate_l = np.zeros((n,n),dtype=np.float64) # N x N
    matrix = np.zeros((n,n),dtype=np.float64) # N x N
    # set up sink
    sink = np.zeros((n),dtype=np.float64)   # N

    # loop over l (index of Bessel function)
    for l in range(lmax):
        sink = np.exp(wrad)*bessel0_zeros[l]**2/rmax**2 # sink term D_par(i) * (b_l)**2
        # in units np.exp(wrad) [dr**2/dt] / rmax**2 [dr**2], so in units [1/dt]

        rate_l[:,:] = rate[:,:]                 # take rate matrix for 1-D diffusion
        rate_l.ravel()[::n+1] -= sink           # and add sink term

        # with rate_l, construct now the inverse of this matrix
        # rate_l_power = rate**power  =>  unit is 1/dt**power
        # assume it has no zero eigenvalues
        if power == -1:
            rate_l_power = np.linalg.inv(rate_l)   # unit [dt]
        elif power == -2:
            rate_l1 = np.linalg.inv(rate_l)
            rate_l_power = np.dot(rate_l1,rate_l1)   # unit [dt**2]
        else:
            raise ValueError(power)

        # bessels no unit, is per r-bin
        # rate_l_power is in 1/dt**power
        # so matrix is in 1/dt**power
        matrix[:,:] += besselsintegral[l] * rate_l_power[:,:]   # in unit 1/time**power

    return matrix
    

def timeintegral_propagator_radial_diffusion(n,dim_rad,rate,wrad,
           lmax,bessel0_zeros,bessels,power=-1):
    """calculate matrices that are needed for the timeintegral of two-D propagator

    n  --  dim_trans, dimension transition matrix, usually number of bins in z-direction
    dim_rad  --  dimension transition matrix, always equal to len(redges)
    rate  --  rate matrix for 1-D diffusion in z-direction, in [1/dt]
    wrad  --  ln Drad, radial diffusion coefficient, dimension n
              Drad = exp(wrad), in [dr**2/dt]
    lmax  --  number of Bessels

    bessels0_zeros  --  first lmax zeros, no unit
    bessels  --  dimension lmax x dim_rad, no unit, is per r-bin
    power  --  rate matrix is taken to this power (as a consequence of integrating the propagator)

    rate_l  --  rate matrix including sink equation, in [1/dt]
    matrix  --  in 1/dt**power
            --  so power=-1 gives unit dt, power=-2 gives unit dt**2"""
    assert power in [-1,-2]  # for now, only these are implemented

    rmax = np.float64(dim_rad)
    # initialize arrays
    rate_l = np.zeros((n,n),dtype=np.float64)         # N x N
    matrix = np.zeros((dim_rad,n,n),dtype=np.float64) # dim_rad x N x N

    # set up sink
    sink = np.zeros((n),dtype=np.float64)   # N

    # loop over l (index of Bessel function)
    for l in range(lmax):
        sink = np.exp(wrad)*bessel0_zeros[l]**2/rmax**2 # sink term D_par(i) * (b_l)**2
        # in units np.exp(wrad) [dr**2/dt] / rmax**2 [dr**2], so in units [1/dt]

        rate_l[:,:] = rate[:,:]                 # take rate matrix for 1-D diffusion
        rate_l.ravel()[::n+1] -= sink           # and add sink term

        # with rate_l, construct now the inverse of this matrix
        # rate_l_power = rate**power  =>  unit is 1/dt**power
        # assume it has no zero eigenvalues
        if power == -1:
            rate_l_power = np.linalg.inv(rate_l)   # unit [dt]
        elif power == -2:
            rate_l1 = np.linalg.inv(rate_l)
            rate_l_power = np.dot(rate_l1,rate_l1)   # unit [dt**2]
        else:
            raise ValueError(power)

        # bessels no unit, is per r-bin
        # rate_l_power is in 1/dt**power
        # so matrix is in 1/dt**power
        for k in range(dim_rad):
            matrix[k,:,:] += bessels[l,k] * rate_l_power[:,:]   # in unit 1/time**power

    return matrix


def setup_bessel_functions_integral(lmax,dim_rad,analytical=True):
    """set up Bessel functions first type zero-th order J_0(b_l x)
    Input
      lmax  --  number of Bessel functions taken into account, each labeled
                by the zero by which the argument is rescaled
      dim_rad  --  len(redges), number of radial bins
    Output
      bessel0_zeros  --  no unit
      bessels  --  no unit: in units per-r-bin"""

    import scipy
    from scipy import linalg, special

    # first lmax zeros of 0th order Bessel first type
    bessel0_zeros = scipy.special.jn_zeros(0,lmax)   # no unit
    # 1st order Bessel first type in those zeros
    bessel1_inzeros = scipy.special.j1(bessel0_zeros)   # no unit
    # 2nd order Bessel first type in those zeros
    bessel2_inzeros = scipy.special.jv(2,bessel0_zeros)   # no unit


    # rmax = dim_rad in high precision ~ radius where all functions are zero
    rmax = np.float64(dim_rad)
    #print "rmax",rmax   # this is float(dim_rad)

    if analytical:
        besselsintegral = rmax**2 *(2./bessel0_zeros/bessel1_inzeros 
                             - 4*bessel2_inzeros/bessel0_zeros**2/bessel1_inzeros**2)
    # or maybe:
    #-------------
    #print scipy.special.besselpoly(....)
    else:
        # Verification
        #-------------
        # set up Bessel functions
        r = np.arange(dim_rad,dtype=np.float64)+0.5
        #print "r",r   # this is [0,1,2,...dim_rad-1]
        #actually, this is [0.5,1.5,2.5,...,dim_rad-1+0.5] so always the middle of the dim_rad bins
        bessels = np.zeros((lmax,dim_rad),dtype=np.float64)
        for l in range(lmax):
            bessels[l,:] = 2*r*scipy.special.j0(r/rmax*bessel0_zeros[l]) / bessel1_inzeros[l]**2 /rmax**2
            # in units r [dr] / rmax**2 [dr**2], so in units [1/dr]
            # so in units 'per r-bin' !
        besselsintegral = np.sum(bessels*r**2,-1)   # in units [dr**2]
    #print "besselsintegral     ",besselsintegral
    #print "besselsintegral sum (not relevant)   ", np.sum(besselsintegral)

    return besselsintegral


###########################################################################

def calc_theor_MSD_1D(prop,dz,donorm=True):
    """theoretical MSD 1D, from 1D propagator"""
    # assume that propagator has no periodic boundary conditions
    n = prop.shape[0]

    dis2 = np.zeros((n,n))
    for i in range(n):
        dis2[i,:] = (i-np.arange(n))**2 * dz**2   # in unit A^2

    a = prop*dis2

    if donorm:
        # - is not necessary if reflective boundaries
        # - is useful if somewhere absorbing boundaries
        norm = np.sum(prop,0)  # = sum of probs of all ending states
                               # = array for different initial positions
        a /= norm   # this will divide column a[:,i] by value norm[i]
                    # this is equivalent to normalizing the columns of prop
    # weighted sum
    msd = np.sum(a,axis=0)  # sum over all possible end states => gives array for different initial positions
    return msd   # in unit A^2

def calc_theor_MSD_2D(prop2d,dz,donorm=True):
    """theoretical MSD 1D, from 2D propagator normal/radial diffusion"""
    assert len(prop2d.shape)==3
    dim_rad = prop2d.shape[0]
    n = prop2d.shape[1]

    dis2 = np.zeros((n,n))
    for i in range(n):
        dis2[i,:] = (i-np.arange(n))**2 * dz**2  # in unit A^2

    a = prop2d * dis2
    #This means: for k in range(dim_rad): a[k,:,:] = prop2d[k,:,:]*dis2
    # dim a = dim_rad x n x n = radial-bin x end-z-state x initial-z-state

    if donorm:
        # - would not be really necessary if reflective boundaries and if prop2d were normalized better
        #   but prop2d is not perfectly normalized because of the boundary R(r=s)=0
        # - is useful if somewhere an absorbing boundary
        norm = np.sum(np.sum(prop2d,0),0)  # = sum of probs of all radial bins and ending-z-states
                                           # = array for different initial z-positions
        #a /= norm  # this will divide a[:,:,i] by value norm[i]
        # if norm <= 0: do not do anything   # TODO is this correct?????????????
        # if norm >0: divide by norm
        for i in range(n):
            if norm[i]>0:
                a[:,:,i] /= norm[i]
    # weighted sum
    msd = np.sum(np.sum(a,0),0)  # sum over all possible end states (radial-bins and ending-z-states)
                                 # this gives array for different initial positions
    return msd   # in unit A^2

def calc_theor_MSDrad(prop2d,dr,donorm=True):
    """theoretical MSD 1D, from 2D propagator normal/radial diffusion"""
    assert len(prop2d.shape)==3
    dim_rad = prop2d.shape[0]
    n = prop2d.shape[1]

    dr_midst = (np.arange(dim_rad)+0.5)*dr
    dis2 = dr_midst**2  # vector length dim_rad, # in unit A^2

    if True:   # the average square of the traveled radial distance <r**2>
        a = np.zeros(prop2d.shape)
        for i in range(n):
            for j in range(n):
                a[:,i,j] = prop2d[:,i,j]*dis2   # python numpy do something like a = dis2* prop2d ?  TODO 

        if donorm:
            # normalize?
            norm = np.sum(np.sum(prop2d,0),0)
            for i in range(n):
                 if norm[i]>0.:
                     a[:,:,i] /= norm[i]
        # weighted sum
        msd = np.sum(np.sum(a,0),0)

    # this gives the same  (start in bin 0)
    #---------------------
    #msd = 0.
    #for k in range(n):
    #    msd += np.sum(dr_midst**2 * prop[:,k,0]*dr)/np.sum(prop[:,k,0]*dr)
    #msd /= end-start

    return msd  # in unit A^2

def do_theor_MSD(rate,prop,prop0,dz,n,edges,redges,measured,lagtime,fig):
    # TODO clean up
    #prop0  -- 1D
    #prop  -- 2D
    nK = len(rate)
    mult = nK/n
    ##### PLOT #####
    # STUPID
    # BUT not stupid if I ALSO fold my trajectories...
    # stupid, because not unfolded...

    # MSD: theoretical
    a = calc_theor_MSD_1D(prop0,dz,n)
    # i did here an extension, with nK, neighbour, ...

    import matplotlib.pyplot as plt
    import scipy

    # MSD: theoretical
    plt.figure(1)
    plt.figure(2)
    #for l in [1e0,1e1,1e2,1e3,1e4,1e5,1e6]:#1e7,1e8]: #[256,512,1024,2048,4096,8192]: #[0.1,1,2,4,8,16,32,]:
    for l in [1e3,1e4,1e5,1e6,]: #[256,512,1024,2048,4096,8192]: #[0.1,1,2,4,8,16,32,]:
        prop0 = scipy.linalg.expm2(rate*l)
        msd = np.zeros(n)
        lim = n*(mult-1)/2 
        #lim = n*neighbour
        nrange = np.arange(n)
        #mm = 1e12
        #mM = 0
        #mr = 1
        for n0 in range(n):
            mid = (mult-1)*n/2+n0
            st = mid-lim
            end = mid+lim
            nrange = np.arange(st,end+1)
            #print n0, st, mid, end,lim, len(nrange)
            tosum = prop0[st:end+1,mid]*(nrange-mid)**2
            msd[n0] = np.sum(tosum)
            #mm = min(mm,min(tosum))
            #mM = max(mM,max(tosum))
            #mr = max(mm/mM,mr)
           # print n0,min(tosum),max(tosum) #tosum[0],tosum[-1]
        #print "msd",msd
        #print mm,mM,mr
        plt.figure(1)
        plt.plot(edges[:n],msd/l)
        plt.savefig("msd.png")
        plt.figure(2)
        mid = nK/2
        st = mid+n/2; st = min(nK-1,st)   # hack
        print mid,st,prop0.shape
        plt.plot(prop0[:,st],color='blue')
        plt.plot(prop0[:,mid],color='red')
        plt.ylim([0,0.15])
        plt.savefig("prop0.png")
    
    measured_D = measured[0]
    measured_edges = measured[2]

    plt.figure()
    plt.plot(edges[:-1],a)
    plt.plot(measured_edges[:-1],measured_D*lagtime)  #######
    plt.ylabel("MSD [A^2]")
    plt.savefig("%s.MSD-1dprop.%i.png"%(fig,lagtime))

    ##### PLOT #####
    # stupid, because not unfolded...
    # MSD
    a = calc_theor_MSD_2D_rad(prop,dz,n)

    plt.figure()
    plt.plot(edges[:-1],a)

    # why do I compute this???
    b = np.zeros(prop.shape)
    for m in range(b.shape[0]):
        b[m,:,:] = prop[m,:,:]*redges[m]**2
    b = b[:,:n,:n]
    #print b
    #plt.plot(edges[:-1],np.sum(np.sum(b.clip(min=0),0),1)/np.sum(np.sum(prop,0),1)) 

    plt.plot(measured_edges[1:],measured_D*lagtime)  #######
    plt.ylabel("MSD [A^2]")
    plt.savefig("%s.MSDb.%i.png"%(fig,lagtime))

