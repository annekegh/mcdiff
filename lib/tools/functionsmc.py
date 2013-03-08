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


def guess_num_bin(filename):
    """get number of histogram bins from transition matrix file"""
    print "Reading...", filename
    f = file(filename)
    count = 0
    for line in f:   # skip lines starting with #
        if not line.startswith("#"):
            assert len(line.split())==1
            count += 1
    f.close()
    num_bin = int(np.sqrt(count))
    #print "number of bins:", num_bin
    return num_bin

def read_transition(filename,num_bin):
    """read transition matrix line by line (one line per entry!!!)
    lines starting with # are skipped
    line starting with # bins is followed by a line with the number of bins"""

    transition = np.zeros((num_bin,num_bin))
    bins = []
    f = file(filename)
    for line in f:
        if line.startswith("#bins"):
            words = line.split()
            bins = [float(word) for word in words[1:]]  # skip first
            bins = np.array(bins)
        if not line.startswith("#"):
            transition[0,0] = np.float64(line)   # put first line in [0,0] element
            k = 1
            break
    for line in f:  # read other lines
        if not line.startswith("#"):
            (k1, k2 ) = divmod ( k , num_bin )
            transition[k1,k2] = np.float64(line)
            #print k1, k2, transition[k1,k2,j]
            k += 1
    if k != num_bin**2:
        print "wrong number of entries in ", filename
        quit()
    f.close()
    if len(bins) == 0:
        bins = np.arange(num_bin+1)
    return transition,bins

def read_all_transitions(argv):
    print argv
    """read all input from sys.argv: lagtimes and transitions"""
    num_lag = ( len ( sys.argv ) - 2 ) / 2
    print "number of lagtimes:", num_lag
    lagtimes = np.zeros((num_lag),dtype=np.float64)

    num_bin = guess_num_bin(argv[3])
    transition = np.zeros((num_lag,num_bin,num_bin),np.float64)
    # read lagtimes and transition files
    for i in range(2,len(sys.argv),2):
        j = i / 2 - 1
        print "window ", j, ":", sys.argv[i], sys.argv[i+1]
        lagtimes[j] = np.float64(sys.argv[i])

        filename = sys.argv[i+1]
        trans,bins = read_transition(filename,num_bin)
        transition[j,:,:] = trans

    #### TODO let's hope I have the same bins every time!!!!


    #print "transition len",len(transition),len(transition[0])
    #for i in range(len(transition)-1):          # assure that all have same number of bins
    #    assert len(transition[i])==len(transition[i+1])
    return lagtimes,transition,bins


class MCState(object):
    def __init__(self):
        pass
    def set_MC_params(self,dv,dw,temp,nmc,num_MC_update):
        self.dv = dv
        self.dw = dw
        self.temp = temp
        self.nmc = nmc
        self.num_MC_update = num_MC_update

    def print_MC_params(self):
        print "dv(MC-potential)=", self.dv
        print "dw(MC-logD)=", self.dw
        print "temp=", self.temp
        print "n(MC)=", self.nmc
        print "n(update)=", self.num_MC_update

    def init_MC(self,num_bin):
        self.num_bin = num_bin
        #self.vt = np.zeros(num_bin)        # trial v
        #self.wt = np.zeros(num_bin) #-1)      # trial w
        self.naccv = np.float64(0.0)       # number accepted v moves
        self.naccw = np.float64(0.0)       # number accepted w moves
        self.naccv_update = np.float64(0.0)       # number accepted v moves between adjusts
        self.naccw_update = np.float64(0.0)       # number accepted w moves between adjusts
        # initialize potential v[i] and w[i]=log(D(i))
        self.v = np.float64(np.zeros(num_bin))
        self.w = np.float64(np.zeros(num_bin)) #-1))    # TODO make this smarter !!!!!!!!!!!!!!!!!

        #dx = self.bins[1]-self.bins[0]
        #self.w -= np.log(dx**2)

        self.init_log_like()

    def init_log_like(self):
        # initialize log_like
        self.log_like = log_like_lag ( self.num_bin, self.num_lag, self.v, self.w, self.lagtimes, self.transition )
        if np.isnan(self.log_like):
            raise Error("Initial likelihood diverges")
        print "initial log-likelihood:", self.log_like

    def add_transitions(self,lagtimes,transition,bins):
        self.lagtimes = lagtimes
        self.num_lag = len(lagtimes)
        self.transition = transition
        self.bins = bins

    def mcmove_potential(self):
        i = np.random.randint(0,self.num_bin)
        vt = copy.deepcopy(self.v)          # temporary v
        vt[i] += self.dv * ( np.random.random() - 0.5 )
        log_like_try = log_like_lag ( self.num_bin, self.num_lag, vt, self.w, self.lagtimes, self.transition )
        return vt,log_like_try

    def mcmetropolis_potential(self,vt,log_like_try):
        dlog = log_like_try - self.log_like
        r = np.random.random()
        if ( r < np.exp(dlog/self.temp) ):
            self.v = vt
            self.naccv += 1.
            self.naccv_update += 1.
            self.log_like = log_like_try

    def mcmove_diffusion(self):
        i = np.random.randint(0,self.num_bin)  #-1)
        wt = copy.deepcopy(self.w)
        wt[i] += self.dw * ( np.random.random() - 0.5 )
        log_like_try = log_like_lag ( self.num_bin, self.num_lag, self.v, wt, self.lagtimes, self.transition )
        return wt,log_like_try

    def mcmetropolis_diffusion(self,wt,log_like_try):
        dlog = log_like_try - self.log_like
        r = np.random.random()
        if ( r < np.exp(dlog/self.temp) ):
            self.w = wt
            self.naccw += 1.
            self.naccw_update += 1.
            self.log_like = log_like_try

    def print_intermediate(self,imc,printfreq):
        if  (imc%printfreq == 0) | (imc == self.nmc-1):
            print imc, self.log_like, self.naccv / ( imc + 1. ), self.naccw / ( imc + 1. )

    def update_movewidth(self,imc):
        if ( self.num_MC_update > 0 ):
            if ( ( imc + 1 ) % self.num_MC_update == 0 ):
                self.dv *= np.exp ( 0.1 * ( self.naccv_update / self.num_MC_update - 0.3 ) )
                self.dw *= np.exp ( 0.1 * ( self.naccw_update / self.num_MC_update - 0.3 ) )
                self.naccv_update = 0.
                self.naccw_update = 0.
                print "new MC steps:", imc, self.dv, self.dw


    def print_final(self): 
        # print final results ( potential and diffusion coefficient)
        print "\n     index  potential  diffusion-coefficient(shifted-by-half-bin)  bin  binmiddle"
        for i in range(self.num_bin):
            #print i, v[i], np.exp(w[i])
            sys.stdout.write("%8d %13.5e %13.5e  %13.5f\n"%( i,self.v[i],np.exp(self.w[i]),self.bins[i], ) )
        #sys.stdout.write("%8d %13.5e  %13.5f\n"%(self.num_bin-1,self.v[-1],np.exp(self.w[i]),bins[-1],))
        sys.stdout.write("#Done")


def init_rate_matrix(n,v,w,pbc=True):
    """initialize rate matrix from potential vector v and diffusion 
    vector w = log(D(i)/delta^2)"""
    rate = np.float64(np.zeros((n,n)))  # high precision

    # off-diagonal elements
    diffv = v[1:]-v[:-1]  # diffv[i] = v[i+1]-v[i]
    exp1 = w[:-1]-0.5*diffv
    exp2 = w[:-1]+0.5*diffv
    rate.ravel()[n::n+1] = np.exp(exp1)
    rate.ravel()[1::n+1] = np.exp(exp2)
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

def log_likelihood(n,ilag,transition,lagtime,rate):
    """calculate log-likelihood from rate matrix and transition matrix
    assuming time step lagtime"""
    # calc propagator as matrix exponential
    propagator = scipy.linalg.expm2(lagtime*rate)
    # sum over all transitions
    # in case of numerical issues try: np.log(np.abs(propagator[i,j]))
    log_like = np.float64(0.0)
    a = np.where(transition[ilag,:,:]>0, transition[ilag,:,:] * np.log(propagator), 0)
    log_like += a.sum()  #use high precision
    return log_like

def log_like_lag(num_bin, num_lag, v, w, lagtimes, transition):
    """calculate log-likelihood summed over all umbrella windows"""
    log_like = np.float64(0.0)
    rate = init_rate_matrix(num_bin,v,w)
    for ilag in range(num_lag):
        # add several distributions
        log_like += log_likelihood(num_bin,ilag,transition,lagtimes[ilag],rate)
    return log_like

if __name__ == '__main__':
    main()
