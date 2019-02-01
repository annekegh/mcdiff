"""Script to compare different methods to extract profiles:
calc flux
AG, April 9, 2013
AG, April 25, 2013
AG, Jan 11, 2016
AG, Sept 2016"""


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
##### PROPERTIES: MFPT #####
#################################

def calc_mfpt_szabo(F,D,dx,dt,st,end,getmax=False):
    """Compute the mean first passage time (MFPT) - formula Szabo 1980"""
    # TODO probably something still wrong with bins summation (too many, too few, D middle...)
    # dx -- in angstrom
    # dt -- in ps
    # D -- in angstrom**2/ps
    # mfpt -- in ps
    # numbering of bins starts with 0

    assert end<len(F) and st<len(F) and end>=0 and st>=0
    # select bins [st, st+1, ..., end-1, end], so this includes end-st+1 bins
    nbins = end-st+1
    v = F - min(F)
    part = np.exp(-v[st:end+1])
    part /= np.sum(part)
    partcum = np.zeros(nbins)
    mfpt = np.zeros(nbins)
    if st>0:
        d = 0.5*(D[st:end+1]+D[st-1:end])  # middle
    else:
        raise NotImplemented
    #d = D[st:end+1]

    # formula Szabo
    #---------------
    if st < end:
        #side = 
        # partcum
        for i in range(nbins):
            partcum[i] = np.sum(part[:i]) *dx   # integral
        integrand = partcum/(d*part)
        # mfpt
        for i in range(nbins):
            mfpt[i] = np.sum(integrand[i:]) *dx  # integral

    elif end < st:
        raise ValueError("not implemented")
        # partcum
        for i in range(nbins):
            partcum[i] = np.sum(part[i:]) *dx  # integral xmiddle to xfinal(=end=reflecting)
        integrand = partcum/(d*part)
        # mfpt
        for i in range(len(mfpt)):
            mfpt[i] = np.sum(integrand[:i]) *dx # integral xinit(=st=absorbing) to xmiddle

    maxmfpt = max(mfpt)   # mfpt[0] or mfpt[-1]
    tau = np.sum(mfpt*part)/np.sum(part)   # in ps, weighted mfpt

    if True:
        #print "partcum",partcum
        #print "tot partcum",np.sum(part)
        #print "mfpt",mfpt
        print("st,end",st,end, end=' ')
        print("tau",tau,"maxmfpt",maxmfpt)

    if getmax:
        return maxmfpt
    else:
        return tau

def calc_mfpt_othermethod(F,D,dx,dt,st,end,getmax=False):
    # not correct yet... this method is just wrong, probably
 
    # construct rate matrix and its pseudo-inverse
    #---------------------------------------------
    rate = construct_rate_matrix_from_F_D(v,D,dx,dt,pbc=False)  # NOPBC    # in units 1/dt
    # cut
    rate_small = rate[st:end+1,st:end+1]
    rate_small1 = np.linalg.pinv(rate_small)  # pseudo-inverse    # in units dt
    # compute mfpt
    mfpt_small = np.zeros(rate_small.shape,float)
    part = np.exp(-v)   # no unit, not normalized
    part_small = part[st:end+1]/np.sum(part[st:end+1])   # normalize
    for i in range(nbins):
        # fill up each row: MFPT(i,j) = starting in j, ending in i
        mfpt_small[i,:] = (-rate_small1[i,i]+rate_small1[i,:])/part_small[i]
    mfpt_small *= dt  # mfpt in ps

    #mfpt = np.zeros(rate.shape,float)
    #mfpt[a:b+1,a:b+1] = mfpt_small
    #tau = mfpt_small[0,-1]  # the longest time
    #tau = mfpt[st,end]
    if False:
        plt.figure()
        plt.plot(mfpt_small)
        plt.savefig("mfpt.matrix.png")

    if False:
        x3 = np.arange(len(mfpt))
        plt.figure()
        I = x3 #range(n)
        J = x3 #range(n)
        CS = plt.contourf(I,J,mfpt.T,30)
        plt.xlabel("ending bin i")
        plt.ylabel("starting bin j")
        # starting in bin 12, then ...
        #plt.plot([x3[12],x3[12]],color='k',linewidth=2)
        plt.plot([x3[0],x3[-1]],[x3[12],x3[12]],color='k',linewidth=2)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel('mean first passage time (ps)')

        #plt.figure()
        #plt.plot(mfpt_small)
        plt.savefig("mfpt.matrix.%i.%i.png" %(st,end))


# THIS ONE WORKS? BUT USE HUMMER INSTEAD
def calc_mfpt(F,D,dx,dt,st,end,getmax=False,side=None):
    """Compute the mean first passage time (MFPT)"""
    # dx -- in angstrom
    # dt -- in ps
    # D -- in angstrom**2/ps
    # mfpt -- in ps
    # numbering of bins starts with 0

    assert st>=0 and end>=0 and st<len(F) and end<len(F)
    assert side in ["both","left","right"]   # absorbing sides, otherwize: reflective

    assert st < end  # TODO #####################################

    # mfpt vector
    #-------------
    # construct rate matrix and inverse
    v = F-min(F)

    rate = construct_rate_matrix_from_F_D(v,D,dx,dt,pbc=False,st=st,end=end,side=side)  # in 1/dt
    print("st,end",st,end,rate.shape)
    rate1 = np.linalg.inv(rate)  # in units dt

    # construct mfpt vector
    mfpt = np.zeros(len(rate),float)
    for i in range(len(mfpt)):
        mfpt[i] = -np.sum(rate1[:,i]) # mfpt[i] is mfpt for starting position i
 
    mfpt *= dt  # in unit "ps"
    maxmfpt = max(mfpt)
    
    # global mfpt
    #------------
    F_small = v[st+1:end]  # submatrix
    part = np.exp(-F_small)   # no unit
    tau = np.sum(part*mfpt)/np.sum(part)    # normalized in interval

    # checking
    #----------
    if False:
        vals,vecs = np.linalg.eig(rate)
        print("vals rate (1/dt)",vals)
        print("mfpt (ps)",mfpt)
    if True:
        print("st,end,side",st,end,side, end=' ') 
        print("tau (ps)",tau, end=' ')
        print("maxmfpt (ps)", max(mfpt))

    if False: #True:
        plt.figure()
        if st<end:
            binrange = np.arange(st,end+1)
        else:
            binrange = np.arange(end,st+1)
        plt.plot(binrange,mfpt)
        plt.plot(binrange[0],[mfpt[0]],"o")
        plt.xlabel("bin")
        plt.ylabel("MFPT (ps)")
        plt.title("tau=%.1f ps, maxmfpt=%.1f ps" %(tau,max(mfpt)))
        plt.savefig("mfpt.vector.%i.%i.%s.png" %(st,end,side))
        plt.close()

    if getmax:
        return maxmfpt
    else:
        return tau

#-------------
#EXTRA August 22, 2016
# this is conditional MFPT
#-------------

# EXTRA functions
# mfpt vector
#-------------

def get_mfpt_from_rate_k12(rate,dt,k1,k2,part,):
    rate1 = np.linalg.inv(rate)    # R^(-1)
    rate12 = np.dot(rate1,rate1)   # R^(-2)
    rate13 = np.dot(rate12,rate1)  # R^(-3)

    totalflux = - k1*rate1[0,:] - k2*rate1[-1,:]  # flux on either side, to bin 0 or to bin N
    #print "totalflux"      # -k*R^(-1)
    #print totalflux  # this is 1 if k1 and k2 are both taken into account

    # Compute mean first passage time
    a = k1*rate12[0,:]+k2*rate12[-1,:]   # k*R^(-2)
    mfpt = a/totalflux                   # k*R^(-2) / (-k*R^(-1) )

    # Compute variance of first passage time
    b = - k1*rate13[0,:] - k2*rate13[-1,:]    ####   TODO correct?  # -k*R^(-3)
    mfpt2 = 2*b/totalflux           # mean      #  -k*R^(-3) / (-k*R^(-1) )

    var_mfpt = mfpt2-mfpt**2        # variance
    std_mfpt = np.sqrt(var_mfpt)    # standard deviation

    # weighted with Boltzmann distribution
    tau = np.sum(mfpt*part)/np.sum(part)
    #print "max",max(mfpt),"tau",tau

    return mfpt,tau, std_mfpt

def get_mfpt_from_rate_k12_before(rate,prop,dt,k1,k2,part,t):
    rate1 = np.linalg.inv(rate)    # R^(-1)
    rate12 = np.dot(rate1,rate1)   # R^(-2)
    rate13 = np.dot(rate12,rate1)  # R^(-3)
    rate1_prop = np.dot(rate1,prop)     # R^(-1)*e^(RT)
    rate12_prop = np.dot(rate12,prop)   # R^(-2)*e^(RT)

    totalflux = - k1*rate1[0,:] - k2*rate1[-1,:]  # flux on either side, to bin 0 or to bin N
    #print "totalflux"      # -k*R^(-1)
    #print totalflux  # this is 1 if k1 and k2 are both taken into account

    #flux before = probability to have left before certain time
    # flux until time T = int(k*exp(Rt), t=0..T)
    # = flux before T = k*R^(-1)*e^(RT) - k*R^(-1)
    flux_after = -k1*rate1_prop[0,:] - k2*rate1_prop[-1,:]  # -k*R^(-1)*e^(RT)
    flux_before = totalflux - flux_after
    #print "*** T",t
    #print "*** flux_before",flux_before

    # Compute mean first passage time
    a1 =   k1*rate12[0,:]       + k2*rate12[-1,:]        # k*R^(-2)
    a2 =   k1*t*rate1_prop[0,:] + k2*t*rate1_prop[-1,:]  # k*T*R^(-1)*e^(RT)
    a3 = - k1*rate12_prop[0,:]  - k2*rate12_prop[-1,:]   # k*R^(-2)*e^(RT)
    a = a1+a2+a3
    mfpt = a/flux_before

    # Compute variance of first passage time   # TODO not fixed yet
    b = - k1*rate13[0,:] - k2*rate13[-1,:]    ####   TODO correct?  # -k*R^(-3)
    mfpt2 = 2*b/totalflux           # mean      #  -k*R^(-3) / (-k*R^(-1) )

    var_mfpt = mfpt2-mfpt**2        # variance
    std_mfpt = np.sqrt(var_mfpt)    # standard deviation

    # weighted with Boltzmann distribution
    tau = np.sum(mfpt*part)/np.sum(part)
    #print "max",max(mfpt),"tau",tau

    return mfpt,tau, std_mfpt

def get_median_from_rate_k12(rate,prop,dt,k1,k2,part,):
    # TODO
    rate1 = np.linalg.inv(rate)    # R^(-1)
    rate12 = np.dot(rate1,rate1)   # R^(-2)
    rate13 = np.dot(rate12,rate1)  # R^(-3)
    rate1_prop = np.dot(rate1,prop)     # R^(-1)*e^(RT)
    rate12_prop = np.dot(rate12,prop)   # R^(-2)*e^(RT)

    totalflux = - k1*rate1[0,:] - k2*rate1[-1,:]  # flux on either side, to bin 0 or to bin N
    #print "totalflux"      # -k*R^(-1)
    #print totalflux  # this is 1 if k1 and k2 are both taken into account

    #flux before = probability to have left before certain time
    # flux until time T = int(k*exp(Rt), t=0..T)
    # = flux before T = k*R^(-1)*e^(RT) - k*R^(-1)
    flux_after = -k1*rate1_prop[0,:] - k2*rate1_prop[-1,:]  # -k*R^(-1)*e^(RT)
    flux_before = totalflux-flux_after


    # Compute mean first passage time
    a1 =   k1*rate12[0,:]       + k2*rate12[-1,:]        # k*R^(-2)
    a2 =   k1*t*rate1_prop[0,:] + k2*t*rate1_prop[-1,:]  # k*T*R^(-1)*e^(RT)
    a3 = - k1*rate12_prop[0,:]  - k2*rate12_prop[-1,:]   # k*R^(-2)*e^(RT)
    a = a1+a2+a3
    mfpt = -np.log(2.)*(rate1[0,:]+rate1[-1,:])
    print("median",mfpt)

    # Compute variance of first passage time   # TODO not fixed yet
    b = - k1*rate13[0,:] - k2*rate13[-1,:]    ####   TODO correct?  # -k*R^(-3)
    mfpt2 = 2*b/totalflux           # mean      #  -k*R^(-3) / (-k*R^(-1) )

    var_mfpt = mfpt2-mfpt**2        # variance
    std_mfpt = np.sqrt(var_mfpt)    # standard deviation

    # weighted with Boltzmann distribution
    tau = np.sum(mfpt*part)/np.sum(part)
    #print "max",max(mfpt),"tau",tau

    return mfpt,tau, std_mfpt


#-------------
def calc_mfpt_hummer(F,D,dx,dt,b1,b2,init=None,doprint=True,dofig=False,t=None,side="right",weigh=False):
    """Compute the mean first passage time (MFPT)"""
    # dx -- in angstrom
    # dt -- in ps
    # D -- in angstrom**2/ps
    # mfpt -- in ps
    # side -- right (default), left, both
    # numbering of bins starts with 0
    # t -- default is None to compute mean fpt,  
    assert side in ["right","left","both"]

    # init -- initial z-bin
    # b1 and b2 - absorbing bins
    assert b1>=-1 and b2<=len(F) and b1+1<b2
    if init is None:
        init = b1+1
    assert init>b1 and init<b2

    # rate matrix: k1,k2    
    v = F-min(F)
    rate = construct_rate_matrix_from_F_D(v,D,dx,dt,pbc=True)  # in 1/dt

    k1 = rate[b1,b1+1]  # at boundary b1   # if b1=-1, then b1+1=0, then k1=rate[-1,0], okay with PBC
    if b2 < len(F):
        k2 = rate[b2,b2-1]  # at boundary b2   # if b2=len(F), then b2-1=len(F)-1, then k2=rate[len(F),len(F)-1]
    elif b2 == len(F):
        k2 = rate[0,-1]  # if b2=len(F), then b2-1=len(F)-1,
                          # then k2=rate[len(F),len(F)-1] should be rate[0,-1]

    # construct rate matrix and inverse
    subrate = construct_rate_matrix_from_F_D(v,D,dx,dt,pbc=False,st=b1,end=b2,side="both")  # in 1/dt

    # for global mfpt
    F_small = v[b1+1:b2]  # submatrix
    part = np.exp(-(F_small-min(F_small)))   # no unit

    if t is None:
        # mfpt
        #---------
        mfpt1,tau1,std1 = get_mfpt_from_rate_k12(subrate,dt,0,k2,part)   # in ps   # this is idential to "exit to right"
        mfpt2,tau2,std2 = get_mfpt_from_rate_k12(subrate,dt,k1,0,part)   # in ps  # this is exit to the left
        mfpt3,tau3,std3 = get_mfpt_from_rate_k12(subrate,dt,k1,k2,part)   # in ps  # this is exit on both sides

    elif t == "median": # t is not None, t is not a float
        # mfpt median
        #------------
        import scipy
        prop = scipy.linalg.expm(subrate*t)
        mfpt1,tau1,std1 = get_median_from_rate_k12(subrate,prop,dt,0,k2,part)   # in ps   # this is idential to "exit to right"
        mfpt2,tau2,std2 = get_median_from_rate_k12(subrate,prop,dt,k1,0,part)   # in ps  # this is exit to the left
        mfpt3,tau3,std3 = get_median_from_rate_k12(subrate,prop,dt,k1,k2,part)   # in ps  # this is exit on both sides
    elif isinstance(t,float):   # t is a float
        # mfpt before time t
        #-------------------
        assert t>0.   # t should be a positive waiting time
        import scipy
        prop = scipy.linalg.expm(subrate*t)
        mfpt1,tau1,std1 = get_mfpt_from_rate_k12_before(subrate,prop,dt,0,k2,part,t)   # in ps   # this is idential to "exit to right"
        mfpt2,tau2,std2 = get_mfpt_from_rate_k12_before(subrate,prop,dt,k1,0,part,t)   # in ps  # this is exit to the left
        mfpt3,tau3,std3 = get_mfpt_from_rate_k12_before(subrate,prop,dt,k1,k2,part,t)   # in ps  # this is exit on both sides

    #print "check:",np.sum((mfpt1-mfpt3[::-1])**2)  # check symmetry between profiles

    # selected mfpt
    #--------------
    sel_mfpt1 = mfpt1[init-b1-1]  # right
    sel_mfpt2 = mfpt2[init-b1-1]  # left
    sel_mfpt3 = mfpt3[init-b1-1]  # both
    sel_std1 = std1[init-b1-1]
    sel_std2 = std2[init-b1-1]
    sel_std3 = std3[init-b1-1]

    # checking
    #----------
    if False:
        vals,vecs = np.linalg.eig(rate)
        print("vals rate (1/dt)",vals)
        print("mfpt (ps)",mfpt)
      
    if doprint:  #rlb = right left both    # This is nice printing!
        print("b1,b2,init",b1,b2,init, end=' ')
        print("tau-rlb",tau1,tau2,tau3, end=' ')
        print("maxmfpt-r", max(mfpt1), end=' ')
        print("sel-rlb", sel_mfpt1, sel_mfpt2, sel_mfpt3, end=' ')
        print("std", sel_std1,sel_std2,sel_std3)

    if dofig:
        plt.figure()
        binrange = np.arange(b1+1,b2)
        plt.plot(binrange,mfpt1,label="exit right")
        plt.plot(binrange,mfpt2,label="exit left")
        plt.plot(binrange,mfpt3,label="exit both sides")

        #plt.plot(binrange[0],[mfpt[0]],"o",)
        plt.xlabel("initial bin")
        plt.ylabel("MFPT (ps)")
        plt.legend(loc='best')
        plt.title("exit to the RIGHT: <tau>=%.1f ps, maxmfpt=%.1f ps" %(tau1,max(mfpt1)))
        plt.savefig("mfpt.hummer.vector.%i.%i.png" %(b1,b2,))
        plt.close()

    # choose side, for "return" value
    if side == "right":
        #mfpt = mfpt1
        sel_mfpt = sel_mfpt1
        if weigh: sel_mfpt = tau1
    elif side == "left":
        #mfpt = mfpt2
        sel_mfpt = sel_mfpt2
        if weigh: sel_mfpt = tau2
    elif side == "both":
        #mfpt = mfpt3
        sel_mfpt = sel_mfpt3
        if weigh: sel_mfpt = tau3
    else: raise ValueError


    return sel_mfpt


#############
# STATISTICS
#############
#############


#-------------

# distribution of Lparallel   # TODO
# L^2(t=t_exit) = mean (r^2, t=t_exit)
# distrib = p(r^2,t=t_exit)   #yes?????????????????

# distrib L^2(t|exit) = P()

# distribution of MFPT
def distribution_mfpt(F,D,dx,dt,b1,b2,trange):
    """Compute the mean first passage time (MFPT)"""
    # dx -- in angstrom
    # dt -- in ps
    # D -- in angstrom**2/ps
    # mfpt -- in ps
    # numbering of bins starts with 0

    # init -- initial z-bin
    # b1 and b2 - absorbing bins
    assert b1>=-1 and b2<=len(F) and b1+1<b2

    # rate matrix
    v = F-min(F)
    # construct rate matrix and inverse
    subrate = construct_rate_matrix_from_F_D(v,D,dx,dt,pbc=False,st=b1,end=b2,side="both")  # in 1/dt

    import scipy
    middle = len(subrate)/2   # this is initial bin "k", choose the middle   # TODO not necessary ??
    rho = []
    cdf = []
    for t in trange:
        mat = scipy.linalg.expm(subrate*t)   # needed in rho and cdf
        # (1) distribution rho(t)
        #------------------------
        # compute rho(t|k), probability to exit at time t
        rho_t = -np.sum(np.dot(subrate,mat),0)
        rho.append(rho_t[middle])
        #print "t",t,"rho",rho_t[middle]

        # (2) cumulative distribution int(rho(t'),t'=0..t)
        #-------------------------------------------------
        # compute CDF of rho(t|k), probability to exit < time t
        # int( R e^Rt, 0..t) = [e^Rt]_0..t = e^Rt-1
        cdf_t = np.diag(np.ones(len(subrate)))-mat
        cdf_t = np.sum(cdf_t,0)
        cdf.append(cdf_t[middle])
        #print "t",t,"cdf",cdf_t[middle]

    rho = np.array(rho)
    # rho contains rho(t|k) for all t values in trange
    print("trange",len(trange),"len",len(rho),"rho_t",rho_t.shape, "rho",rho.shape)

    cdf = np.array(cdf)
    # cdf contains cdf(t|k) for all t values in trange

    print("trange",len(trange),"len",len(cdf),"cdf_t",cdf_t.shape, "cdf",cdf.shape)

    # Writing data to file
    #----------------------
    filename = "t_rho.dat"
    data = np.zeros((len(trange),3))
    data[:,0] = trange
    data[:,1] = rho
    data[:,2] = cdf
    np.savetxt(filename,data)
    print("file written...", filename)

    if True:
        # Plot rho(t|k)   HERE k = MIDDLE
        #--------------------------------
        plt.figure()
        plt.subplot(2,1,1)
        plt.title("P(exit at time t)")
        plt.plot(trange,rho,"o")
        plt.xlabel("t [ps]")

        plt.subplot(2,1,2)
        plt.semilogx(trange,rho,"o")
        plt.xlabel("t [ps]")
        plt.savefig("t_rho.png")

    if True:
        # Plot CDF(t|k)   HERE k = MIDDLE
        #--------------------------------
        fig = plt.figure()
        plt.subplot(2,1,1)
        plt.title("P(exit < time t)")
        plt.plot(trange,cdf,"o")
        plt.xlabel("t [ps]")

        plt.subplot(2,1,2)
        plt.semilogx(trange,cdf,"o")
        plt.xlabel("t [ps]")
        plt.savefig("t_rho_cdf.png")

    calc_statistics_rho(trange,rho)

    return rho,cdf

def calc_statistics_rho(trange,rho):
    print("--")
    print("statistics")
    print("--")
    # mean
    dtrange = trange[1:]-trange[:-1]
    norm = np.sum(rho[1:]*dtrange)                  # not perfect, binning
    m = np.sum(((trange*rho)[1:])*dtrange)/norm
    # std
    m2 = np.sum(((trange**2*rho)[1:])*dtrange)/norm
    std = np.sqrt(m2-m**2)
    std2_a = np.sum((((trange-m)**2*rho)[1:])*dtrange)/norm
    std_a = np.sqrt(std2_a)
    print("norm", norm)
    print("mean", m)
    print("std", std,std_a)
    return m,std



def get_inverse_cdf(trange,cdf):
    """function that computes the inverse CDF based on the CDF
    Input
        trange -- array of t values
        cdf -- array of CDF(t) values,
               i.e. the cumulative distribution function (CDF)
    Output
        inverse_cdf -- array with Nrl times
                    inverse_cdf[0] minimum considered time, cdf is 0 = 0/Nrl
                    inverse_cdf[-1] largest considered time, cdf is 1 = Nrl/Nrl
                    cdf[inverse_cdf[i]] = i/Nrl, in interval [0,1]
        delta_inverse_cdf -- array with differences between subsequent elements
                             of inverse_cdf, so inverse_cdf[1:]-inverse_cdf[:-1]
    """
    Nrl = len(cdf)
    inverse_cdf = np.zeros(Nrl)
    inverse_cdf[0] = trange[0]
    k = 0
    # uniform grid
    #y = np.arange(0,1,0.00001)   # equal amount?
    Ny = Nrl    # equal amount???   # TODO
    y = np.arange(0,Ny)/float(Ny)
    for n in range(1,Ny):
        # check if value from uniform distr is in the next CDF interval
        # find first k where cdf is bigger
        while cdf[k] < y[n] and k < Nrl:
            k += 1
        # store corresponding time 
        inverse_cdf[n] = trange[k]
        # or interpolate
        #inverse_cdf[n] = trange[k-1] + (trange[k] - trange[k-1]) * (y[n] - cdf[k-1])/(cdf[k] - cdf[k-1]) 
        #if k >= Nrl:
        #    break
    delta_inverse_cdf = np.concatenate((np.diff(inverse_cdf), [0]))

    #print "cdf",cdf
    #print "inverse_cdf",inverse_cdf
    #print "delta_inverse_cdf",delta_inverse_cdf
    #print "0, 1/4, 1/2, 3/4, -1"
    #print inverse_cdf[0],inverse_cdf[Nrl/4],inverse_cdf[Nrl/2],inverse_cdf[3*Nrl/4],inverse_cdf[-1]
    #print trange[0],trange[Nrl/4],trange[Nrl/2],trange[3*Nrl/4],trange[-1]   # this is very different
    return inverse_cdf,delta_inverse_cdf

def draw_from_cdf(N,inverse_cdf,delta_inverse_cdf):
    """Draw from a ???? distribution
    Input
        N -- how many draws
        inverse_cdf -- array with Nrl times,
                        with cdf[inverse_cdf[i]] = i/Nrl, in interval [0,1]
        delta_inverse_cdf -- array with differences between subsequent elements
                             of inverse_cdf, so inverse_cdf[1:]-inverse_cdf[:-1]
    Output
        y -- array of N draws
    """
    Nrl = len(inverse_cdf)
    indices_float = np.random.uniform(size=N, high=Nrl-1)   # draw from uniform distrib
    indices = np.array(indices_float,'i')   # floor / ceil ?
    # interpolate around inverse_cdf[indices]
    y = inverse_cdf[indices] + (indices_float-indices)*delta_inverse_cdf[indices]
    # note: the simple formula would be:   y = inverse_cdf[indices]
    return y

def statistics_from_distribution(trange,rho,cdf):
    # distribution t by drawing myself
    #---------------------------------
    inverse_cdf,delta_inverse_cdf = get_inverse_cdf(trange,cdf)

    def get_quantile(inverse_cdf,percent):
        # 0% -> not in array
        # few% -> [0]
        # almost100% -> [-1]
        # 100% -> not in array
        # approximate as follows:
        # 0% -> [0]
        # 100% -> [-1]
        N = len(inverse_cdf)
        indexlow = np.floor((N-1)*percent)
        indexhigh = np.ceil((N-1)*percent)
        frac = N*percent - np.floor(N*percent)
        select = inverse_cdf[indexlow]*(1-frac) + frac*inverse_cdf[indexhigh]
        return select

    print("--")
    print("statistics")
    print("--")
    Nrl = len(inverse_cdf)
    for percent in np.arange(0.,1.001, 0.1):
        quant = get_quantile(inverse_cdf,percent)
        print(percent, quant)
    print("--")

    print("Getting the 90% quantile figured out")
    t90 = get_quantile(inverse_cdf,0.9)
    print("t90",t90)

    # where is t90 located?
    index = np.digitize([t90],trange)
    index = index[0]
    print("bin-t90, trange[bin-t90]:  ",index,trange[index])

    r = rho[index:]
    tr = trange[index:]
    # normalize rho[index:]
    #dtrange = tr[1:]-tr[:-1]
    #norm = np.sum(r[1:]*dtrange)                  # not perfect, binning
    #r /= norm
    # use this for statistics
    m90,std90 = calc_statistics_rho(tr,r)     # YES!!! THIS IS IT!!!



def draw_from_distribution(trange,cdf):

    print(done)
    inverse_cdf,delta_inverse_cdf = get_inverse_cdf(trange,cdf)

    ally = []
    for i in range(10000):
        N = 40000
        N = 1
        y = draw_from_cdf(N,inverse_cdf,delta_inverse_cdf)
        #print y
        #print "y",np.mean(y),np.std(y)
        ally.append(np.mean(y))

    print("statistics")
    print("mean, mean90")
    print(np.mean(ally), np.mean(sorted(ally)[int(len(ally)*0.9):]))


    # statistics of average of N=40
    hist,levels = np.histogram(ally,bins=30,)
    histcdf = np.zeros(hist.shape)
    for i in range(len(hist)):
        histcdf[i] = np.sum(hist[:i])
    histcdf /= np.sum(hist)
    mid = (levels[:-1]+levels[1:])/2.
    mid /= 1000   # from ps to ns
    print(len(hist),len(levels))

    #if True:
    if False:
        # Plot
        #--------------------------------
        plt.figure()
        plt.subplot(211)
        plt.plot(mid,hist)
        plt.xlabel("time [ns]")
        plt.ylabel("prob")
        plt.subplot(212)
        plt.plot(mid,histcdf)
        plt.plot(mid,np.ones(mid.shape)*0.25,'k')
        plt.plot(mid,np.ones(mid.shape)*0.75,'k')
        plt.plot(mid,np.ones(mid.shape)*0.975,'r')
        plt.plot(mid,np.ones(mid.shape)*0.025,'r')
        plt.plot(mid,np.ones(mid.shape)*0.5,'g')
        #plt.plot(np.ones(2)*mean, [0.,1.],'g')
        #plt.plot(np.ones(2)*(mean-2*ste), [0.,1.],'r')
        #plt.plot(np.ones(2)*(mean+2*ste), [0.,1.],'r')
        plt.ylim(ymin=0.)
        plt.xlabel("time [ns]")
        plt.ylabel("CDF")

        plt.savefig("t_rho_histogram_40000samples.png")




##############
def whatever_construct_propagator_from_F_D(F,D,Drad,dz,dr,dt,lagtime,lmax,dim_rad,pbc=True):
    # F -- in kBT
    # D -- in unit=angstrom**2/ps, D/unit has no dimension
    #      convert to w = ln(D/unit)
    # dz, dr -- angstrom
    # dt -- in ps
    # propagator  --  no unit, is per r-bin per z-bin
#    # normal
#    rate = construct_rate_matrix_from_F_D(F,D,dz,dt,pbc=pbc)   # in 1/dt


    # radial
    from mcdiff.twod import setup_bessel_functions, propagator_radial_diffusion
    wradunit = np.log(dr**2/dt)
    wrad = np.log(Drad)-wradunit

    bessel0_zeros,bessels = setup_bessel_functions(lmax,dim_rad)
    n = len(rate)
    propagator = propagator_radial_diffusion(n,dim_rad,rate,wrad,lagtime,
           lmax,bessel0_zeros,bessels,)   # no unit, is per r-bin per z-bin
    #return propagator





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
        print("file written...",figname)
    return P

