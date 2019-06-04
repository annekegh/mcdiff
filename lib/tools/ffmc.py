#!/usr/bin/env python
#
# usage:
#
# python ./diffusion_1D_new.py number_of_bins lagtime1 transition1.dat lagtime2 transition2.dat lagtime3 transition3.dat ...
#
# copyright: Gerhard Hummer (NIH, July 2012)
# An Ghysels (August 2012)
#

import numpy as np
import sys
import scipy
import scipy.linalg
import copy


#------------------------
# READING FUNCTIONS
#------------------------

def guess_dim_transition_linebyline(filename):
    """get number of histogram bins from transition matrix file (one line per entry!!!)"""
    #print "Reading...", filename
    f = file(filename)
    count = 0
    for line in f:   # skip lines starting with #
        if not line.startswith("#"):
            assert len(line.split())==1
            count += 1
    f.close()
    dim_trans = int(np.sqrt(count))
    return dim_trans

def guess_dim_transition_square(filename):
    """get number of histogram bins from transition matrix file"""
    #print "Reading...", filename
    f = file(filename)
    for line in f:   # skip lines starting with #
        if not line.startswith("#"):
            dim_trans = len(line.split())
    f.close()
    return dim_trans

def read_transition_linebyline(filename,dim_trans):
    """read transition matrix line by line (one line per entry!!!)
    comment lines starting with # are skipped
    comment line starting with #bins contains the bin edges"""

    transition = np.zeros((dim_trans,dim_trans))
    bins = []
    f = file(filename)
    for line in f:
        if line.startswith("#bins"):
            words = line.split()
            bins = [float(word) for word in words[1:]]  # skip first
            bins = np.array(bins)   # bin edges
            assert len(bins) == dim_trans+1
        if not line.startswith("#"):
            transition[0,0] = np.float64(line)   # put first line in [0,0] element
            k = 1
            break
    for line in f:  # read other lines
        if not line.startswith("#"):
            (k1, k2 ) = divmod ( k , dim_trans )
            transition[k1,k2] = np.float64(line)
            k += 1
    if k != dim_trans**2:
        print("wrong number of entries in ", filename)
        quit()
    f.close()

    if len(bins) == 0:  # use linear
        bins = np.arange(dim_trans+1)
    return transition,bins


def read_transition_square(filename,dim_trans):
    """read transition matrix line by line (whole row per line!!!)
    comment lines starting with # are skipped
    comment line starting with #bins contains the bin edges"""

    transition = np.zeros((dim_trans,dim_trans))
    bins = []
    Dt = 0
    f = file(filename)
    row = 0
    for line in f:
        if line.startswith("#Dt"):
            words = line.split()
            assert len(words)==2
            Dt = float(words[1])   # this is the lagtime!!
        if line.startswith("#bins"):
            words = line.split()
            bins = [float(word) for word in words[1:]]  # skip first
            bins = np.array(bins)   # bin edges

         # TODO   assert len(bins) == dim_trans+1
        if not line.startswith("#"):
            words = line.split()
            transition[row,:] = [int(word) for word in words]
            row += 1
    if row != dim_trans:
        print("wrong number of entries in ", filename)
        quit()
    f.close()

    if len(bins) == 0:  # use linear
        bins = np.arange(dim_trans+1)
    return transition,bins,Dt



def read_all_transitions(args,form="square"):
    """read all input: lagtimes and transitions"""
    print("args:", args)

    # number of lag times
    num_lag = len(args)/2
    lagtimes = np.zeros((num_lag),dtype=np.float64)
    print("number of lagtimes:", num_lag)

    # dimension of transition matrix
    if form == "square":
        dim_trans = guess_dim_transition_square(args[1])
    else:
        dim_trans = guess_dim_transition_linebyline(args[1])

    transition = np.zeros((num_lag,dim_trans,dim_trans),np.float64)
    print("dimension trans:", dim_trans)

    # read lagtimes and transition files
    for i in range(num_lag):
        lagtimes[i] = np.float64(args[2*i])

        filename = args[2*i+1]
        if form == "square":
            trans,bins,Dt = read_transition_square(filename,dim_trans)
            if Dt != 0:
                lagtimes[i] = Dt
        else:
            trans,bins = read_transition_linebyline(filename,dim_trans)
        transition[i,:,:] = trans
        print("window ", i, ": lagtime", lagtimes[i], "or",args[2*i], "from file", args[2*i+1])

    #### TODO let's hope I have the same bins every time!!!!
    ###  TODO add quality checks: bins, lagtime, dim_trans

    return lagtimes,transition,bins



#------------------------
# MONTE CARLO
#------------------------

class MCState(object):
    def __init__(self,pbc=True):
        self.pbc = pbc  # whether to use periodic boundary conditions

    def set_MC_params(self,dv,dw,temp,nmc,num_MC_update,dt,temp_end=None):
        self.dv = dv
        self.dw = dw
        self.temp = temp
        self.temp_start = temp
        self.temp_end = temp_end
        self.nmc = nmc
        self.num_MC_update = num_MC_update
        self.dt = dt

        if self.temp_end == None:
            self.temp_end = temp
            self.dtemp = 0.
        else:
            if self.nmc%self.num_MC_update == 0:
                nupdate = self.nmc/self.num_MC_update-1
            else:
                nupdate = self.nmc/self.num_MC_update
            self.dtemp = (self.temp_end-self.temp_start)/float(nupdate)  #if changing by adding
            #self.fdtemp = (self.temp_end/self.temp_start)**(1./nupdate)  #if changing by multiplying

    def print_MC_params(self):
        print("----- Settings MC -----")
        print("dv(MC-potential)=", self.dv)
        print("dw(MC-logD)=", self.dw)
        print("temp=", self.temp)
        print("n(MC)=", self.nmc)
        print("n(update)=", self.num_MC_update)
        print("dt=",self.dt)
        #print "k=", self.k
        print("-"*20)

    def init_MC(self):
        self.num_F = self.num_bin
        if self.pbc:
            self.num_D = self.num_bin
        else:
            self.num_D = self.num_bin-1
        self.naccv = 0              # number accepted v moves
        self.naccw = 0              # number accepted w moves
        self.naccv_update = 0       # number accepted v moves between adjusts
        self.naccw_update = 0       # number accepted w moves between adjusts
        # initialize potential v[i] and w[i]=log(D(i))
        self.v = np.float64(np.zeros(self.num_F))
        self.w = np.float64(np.zeros(self.num_D))
        if True:   # TODO improve - make better guess for v and w
            dx = self.bins[1]-self.bins[0]  # in angstrom
            unit = dx**2/self.dt            # in angstrom**2/ps
            self.w -= np.log(unit)          # exp(w) is in units dx**2/dt
            #print dx, 1./dx**2, unit, "D:", np.exp(self.w)*unit,"in angstrom**2/ps"
         
        # Smoothing
        #------------
        self.k = 10.*2./np.log(1.05)**2 * 0.01  # TODO make this a setting
        print("k, ik", self.k)
        # as in paper GH:
        #self.k = 1./0.05 # in units s/angstrom**2
        #self.k /= unit  # in units 'unit'
        #print "k,GH",self.k, "in unit"  # in angstrom**2/ps
        #print "k,GH",self.k * unit,"in angstrom**2/ps"  # in angstrom**2/ps

        # calc initial likelihood
        self.init_log_like()

    def init_log_like(self):
        # initialize log_like
        log_like = log_like_lag ( self.num_bin, self.num_lag, self.v, self.w, self.lagtimes, self.transition,pbc=self.pbc )
        if np.isnan(log_like):
            raise Error("Initial likelihood diverges")
        # add smoothing
        E_w = string_energy(self.w,self.k,self.pbc)
        self.log_like = log_like - E_w  # because
        print("initial log-likelihood:", self.log_like)
        self.all_log_like = np.zeros(self.nmc,float)

    def add_transitions(self,lagtimes,transition,bins):
        self.lagtimes = lagtimes
        self.transition = transition
        self.bins = bins   # bin edges
        # some handy dimensions
        self.num_lag   = len(lagtimes)
        self.num_bin   = len(bins)-1
        self.num_edges = len(bins)
        self.dim_trans = len(bins)-1

    def mcmove_potential(self):
        i = np.random.randint(0,self.num_F)
        vt = copy.deepcopy(self.v)          # temporary v
        vt[i] += self.dv * ( np.random.random() - 0.5 )
        log_like_try = log_like_lag ( self.num_F, self.num_lag, vt, self.w, self.lagtimes, self.transition, pbc=self.pbc )
        return vt,log_like_try

    def mcmetropolis_potential(self,vt,log_like_try):
        dlog = log_like_try - self.log_like
        r = np.random.random()
        if r < np.exp(dlog/self.temp):
            self.v = vt
            self.naccv += 1
            self.naccv_update += 1
            self.log_like = log_like_try

    def mcmove_diffusion(self):
        i = np.random.randint(0,self.num_D)
        wt = copy.deepcopy(self.w)
        wt[i] += self.dw * (np.random.random()-0.5)
        log_like_try = log_like_lag(self.num_bin, self.num_lag, self.v, wt, self.lagtimes, self.transition, pbc=self.pbc )
        # add restraints to smoothen
        E_wt = string_energy(wt,self.k,self.pbc)
        log_like_try -= E_wt  # because
        print("L,E", log_like_try, E_wt, log_like_try+E_wt)
        return wt,log_like_try

    def mcmetropolis_diffusion(self,wt,log_like_try):
        dlog = log_like_try - self.log_like
        r = np.random.random()  #in [0,1[
        if r < np.exp(dlog/self.temp): # accept if dlog increases, accept maybe if decreases
            self.w = wt
            self.naccw += 1
            self.naccw_update += 1
            self.log_like = log_like_try

    def print_intermediate(self,imc,printfreq):
        step = imc+1
        if  (imc%printfreq == 0) | (step == self.nmc):
            print(imc, self.log_like, float(self.naccv)/step, float(self.naccw)/step)

    def update_temp(self,imc):
        if self.num_MC_update > 0:
            if (imc+1)%self.num_MC_update == 0:
                self.temp += self.dtemp
                #self.temp *= self.fdtemp
                print("new MC temp:", imc, self.temp)

    def update_movewidth(self,imc):
        "adapt dv and dw such that acceptance ratio stays around 30 procent, or so"  # TODO
        if self.num_MC_update > 0:
            if ( ( imc + 1 ) % self.num_MC_update == 0 ):
                self.dv *= np.exp ( 0.1 * ( float(self.naccv_update) / self.num_MC_update - 0.3 ) )
                self.dw *= np.exp ( 0.1 * ( float(self.naccw_update) / self.num_MC_update - 0.3 ) )
                self.naccv_update = 0
                self.naccw_update = 0
                print("new MC steps:", imc, self.dv, self.dw)

    def print_log_like(self):
        print("===== log_like =====")
        for i in range(self.nmc/20):
             print(" ".join([str(val) for val in self.all_log_like[20*i:20*(i+1)]]))
        print("="*10)

    def print_statistics(self):
        print("===== Statistics =====")
        print("nmc       ", self.nmc)
        print("naccv     ", self.naccv)
        print("naccw     ", self.naccw)
        print("accv ratio", "%5.1f" %(float(self.naccv)/self.nmc*100),"%")
        print("accw ratio", "%5.1f" %(float(self.naccw)/self.nmc*100),"%")
        print("="*10)

    def print_final(self): 
        """print final results (potential and diffusion coefficient)"""
        #print "\n     index  bin  potential  diffusion-coefficient(shifted-by-half-bin)"
        print("\n%8s %8s %8s  %13s %s" % ("index","bin-str","bin-end","potential","diffusion-coefficient(shifted-by-half-bin)"))
        if self.pbc:
            for i in range(self.num_bin):
                sys.stdout.write("%8d %8.3f %8.3f  %13.5e %13.5e\n"%( i,self.bins[i],self.bins[i+1],self.v[i],np.exp(self.w[i]) ) )
        else:
            for i in range(self.num_bin-1):
                sys.stdout.write("%8d %8.3f %8.3f  %13.5e %13.5e\n"%( i,self.bins[i],self.bins[i+1],self.v[i],np.exp(self.w[i]) ) )
            sys.stdout.write("%8d %8.3f %8.3f  %13.5e\n"%( self.num_bin-1,self.bins[-2],self.bins[-1],self.v[-1] ) )

        sys.stdout.write("#Done")

#------------------------
# EXTRA FUNCTIONS
#------------------------

def init_rate_matrix(n,v,w,pbc=True):
    """initialize rate matrix from potential vector v and diffusion
    vector w = log(D(i)/delta^2)"""
    assert len(v) == n  # number of bins
    if pbc:
        assert len(w) == n
    else:
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

    # corners
    if pbc:  # periodic boundary conditions
        rate[0,-1]  = np.exp(w[-1]-0.5*(v[0]-v[-1]))
        rate[-1,0]  = np.exp(w[-1]-0.5*(v[-1]-v[0]))
        rate[0,0]   = - rate[1,0] - rate[-1,0]
        rate[-1,-1] = - rate[-2,-1] - rate[0,-1]
    else:  # reflecting boundaries (is equal to a hard wall)
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


def log_likelihood(n,ilag,transition,lagtime,rate):
    """calculate log-likelihood from rate matrix and transition matrix
    assuming time step lagtime"""
    # calc propagator as matrix exponential
    propagator = scipy.linalg.expm(lagtime*rate)
    # sum over all transitions
    # in case of numerical issues try: np.log(np.abs(propagator[i,j]))
    log_like = np.float64(0.0)
    a = np.where(transition[ilag,:,:]>0, transition[ilag,:,:] * np.log(abs(propagator)), 0)
    #print transition[ilag,:,:]
    #print propagator
    #print np.log(propagator)
    #print done
    log_like += a.sum()  #use high precision
    return log_like

def log_like_lag(num_bin, num_lag, v, w, lagtimes, transition,pbc=True):
    """calculate log-likelihood summed over all umbrella windows"""
    log_like = np.float64(0.0)
    rate = init_rate_matrix(num_bin,v,w,pbc=pbc)
    for ilag in range(num_lag):
        # add several distributions
        log_like += log_likelihood(num_bin,ilag,transition,lagtimes[ilag],rate)
    return log_like


