#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# An Ghysels (August 2012)
#

import numpy as np
import copy

from .utils import init_rate_matrix, string_energy, string_vecs, log_likelihood, log_like_lag
from .twod import rad_log_like_lag, setup_bessel_functions

from .model import Model, RadModel
from .model import SinusCosinusModel,CosinusModel, RadCosinusModel
from .model import StepModel, OneStepModel
from .outreading import read_Fcoeffs, read_Dcoeffs, read_Dradcoeffs, read_dv_dw, read_F_D_edges


#------------------------
# MONTE CARLO
#------------------------

class MCState(object):
    def __init__(self,pbc,lmax=-1):
        self.pbc = pbc  # whether to use periodic boundary conditions
        self.lmax = lmax
        self.do_radial = (self.lmax>0)  # True or False

    def set_MC_params(self,dv,dw,dwrad,D0,dtimezero,temp,nmc,num_MC_update,move_timezero,k,temp_end=None,):
        self.dv = dv
        self.dw = dw
        self.dwrad = dwrad
        self.D0 = D0
        self.dtimezero = dtimezero
        self.temp = temp
        self.temp_start = temp
        self.temp_end = temp_end
        self.nmc = nmc
        self.num_MC_update = num_MC_update
        self.move_timezero = move_timezero

        if self.temp_end == None:
            self.temp_end = temp
            self.dtemp = 0.
        else:
            nupdate = self.nmc/self.num_MC_update-1
            if nupdate >0:
                self.dtemp = (self.temp_end-self.temp_start)/float(nupdate)  #if changing by adding
                #self.fdtemp = (self.temp_end/self.temp_start)**(1./nupdate)  #if changing by multiplying
            else:
                self.dtemp = 0.
        self.naccv = 0                 # number accepted v moves
        self.naccw = 0                 # number accepted w moves
        self.naccwrad = 0              # number accepted wrad moves
        self.nacctimezero = 0          # number accepted timezero moves
        self.naccv_update = 0          # number accepted v moves between adjusts
        self.naccw_update = 0          # number accepted w moves between adjusts
        self.naccwrad_update = 0       # number accepted wrad moves between adjusts
        self.nacctimezero_update = 0   # number accepted timezero moves between adjusts

        self.k = k  # spring constant in function spring

    def set_model(self,model,data,ncosF,ncosD,ncosDrad,pull):
        self.data = data   # transitions etc
        # convert pull (kBT/angstrom) to dF (kBT) between bins
        if pull is not None:
            self.pull = -pull*self.data.dz   # dF/dz = -pull
        else:
            self.pull = pull   #0.

        ncosP = 0

        # derive model
        if self.do_radial > 0:
            if model == "RadCosinusModel":
                self.model = RadCosinusModel(self.data,self.D0,ncosF,ncosD,ncosP,ncosDrad)
            elif model == "RadModel":
                self.model = RadModel(self.data,self.D0,ncosF,ncosD,ncosP)
            else:
                raise ValueError( "model %s not found" % model)
            bessel0_zeros,bessels = setup_bessel_functions(self.lmax,self.model.dim_rad,)
            self.model.bessels = bessels
            self.model.bessel0_zeros = bessel0_zeros
            self.model.rate = init_rate_matrix(self.model.dim_v,self.model.v,self.model.w,self.pbc,self.pull)
        else:
            if model == "CosinusModel":
                self.model = CosinusModel(self.data,self.D0,ncosF,ncosD,ncosP)
                # this will default to Model(self,data) if ncosF and ncosD are both 0
            elif model == "SinusCosinusModel":
                self.model = SinusCosinusModel(self.data,self.D0,ncosF,ncosD,ncosP)
            elif model == "StepModel":
                self.model = StepModel(self.data,self.D0,ncosF,ncosD,ncosP)
            elif model == "OneStepModel":
                self.model = OneStepModel(self.data,self.D0,ncosF,ncosD,ncosP)
            elif model == "Model":
                self.model = Model(self.data,self.D0)
            else:
                raise ValueError("model %s not found" % model)
        assert self.pbc == self.model.pbc  # make sure PBC for model and transition matrix are identical

    def init_log_like(self):
        # initialize log_like
        if self.do_radial:
            self.model.rate = init_rate_matrix(self.model.dim_v,self.model.v,self.model.w,self.pbc,self.pull)
            log_like = rad_log_like_lag(self.model.dim_v, self.model.dim_rad,
                  self.data.dim_lt, self.model.rate, self.model.wrad,
                  self.data.list_lt, self.data.list_trans, self.model.redges,
                  self.lmax,self.model.bessel0_zeros,self.model.bessels, 0.)
        else:
            log_like = log_like_lag(self.model.dim_v, self.data.dim_lt,
                  self.model.v, self.model.w, self.model.list_lt,
                  self.data.list_trans, self.pbc,self.pull)
        
        if log_like is None:
            raise ValueError("Initial propagator has non-positive elements")
        elif np.isnan(log_like):
            raise ValueError("Initial likelihood diverges")
        self.log_like = log_like

        # add smoothing to diffusion profile
        if self.k > 0.:
            E_w = string_energy(self.model.w,self.k,self.pbc)
            self.string_vecs = string_vecs(len(self.model.w),self.pbc)
            self.log_like = log_like - E_w  # minus sign because surface=log_like

        print("initial log-likelihood:", self.log_like)
        self.all_log_like = np.zeros(self.nmc,float)

        # TODO make nicer
        if self.model.ncosF > 0:
            self.naccv_coeff = np.zeros(self.model.ncosF,int)
        if self.model.ncosD > 0:
            self.naccw_coeff = np.zeros(self.model.ncosD,int)
        if self.do_radial:
            if self.model.ncosDrad > 0:
                 self.naccwrad_coeff = np.zeros(self.model.ncosDrad,int)
            else: self.model.ncosDrad = -1

    def use_initfile(self,initfile,final=True):
        if self.model.ncosF > 0:
            v_coeff = read_Fcoeffs(initfile,final=True) # unit: v_coeff[0] in kBT
            nc = len(v_coeff)
            if nc > 0:
                print("USING initfile for v_coeff",initfile,nc,"coeffs")
                n = min(nc,self.model.ncosF)
                self.model.v_coeff[:n] = v_coeff[:n]
                self.model.update_v()
        else:
            F,D,edges = read_F_D_edges(initfile) # unit: F in kBT
            nc = len(F)
            assert nc == len(self.model.v)
            print("USING initfile for v",initfile,nc,"values")
            self.model.v = F   # unit: always in kBT

        if self.model.ncosD > 0:
            w_coeff = read_Dcoeffs(initfile,final=True) # unit: w_coeff[0] in angstrom**2/ps
            nc = len(w_coeff)
            if nc > 0:
                print("USING initfile for w_coeff",initfile,nc,"coeffs")
                n = min(nc,self.model.ncosD)
                self.model.w_coeff[:n] = w_coeff[:n]
                self.model.w_coeff[0] -= self.model.wunit
                self.model.update_w()
        else:
            F,D,edges = read_F_D_edges(initfile) # unit: D in angstrom**2/ps
            nc = len(D)
            assert nc == len(self.model.w)
            print("USING initfile for w",initfile,nc,"values")
            self.model.w = np.log(D)-self.model.wunit

        if self.do_radial:
          if self.model.ncosDrad > 0:
            coeff = read_Dradcoeffs(initfile,final=True) # unit: wrad_coeff[0] in angstrom**2/ps
            nc = len(coeff)
            if nc > 0:
                    print("USING initfile for wrad_coeff",initfile,nc,"coeffs")
                    n = min(nc,self.model.ncosDrad)
                    self.model.wrad_coeff[:n] = coeff[:n]
                    self.model.wrad_coeff[0] -= self.model.wradunit
                    self.model.update_wrad()
            else:
                #print self.model.wrad_coeff
                coeff = read_Dcoeffs(initfile,final=True) # unit: w_coeff[0] in angstrom**2/ps
                nc = len(coeff)
                if nc > 0:
                    print("USING initfile for wrad_coeff",initfile,nc,"coeffs, using w_coeff!")
                    n = min(nc,self.model.ncosDrad)
                    self.model.wrad_coeff[:n] = coeff[:n]
                    self.model.wrad_coeff[0] -= self.model.wradunit
                    self.model.update_wrad()
                    #print self.model.wrad_coeff
          else:
            Drad,redges = read_Drad(initfile) # unit: Drad in angstrom**2/ps
            nc = len(Drad)
            assert nc == len(self.model.wrad)
            print("USING initfile for wrad",initfile,nc,"values")
            self.model.wrad = np.log(Drad)-self.model.wradunit

        dv,dw = read_dv_dw(initfile,final=True)
        #TODO   self.dv = dv
        #TODO   self.dw = dw

    #======== MONTE CARLO MOVES ========

    def mcmove_timezero(self):
        timezero_try = self.model.timezero + self.dtimezero * (np.random.random()-0.5)
        if timezero_try > -0.5*self.data.min_lt:     # ensure that shortest lagtime shrinks to no less than 1/2
            lagtimes_try = self.data.list_lt + timezero_try
            log_like_try = log_like_lag(self.model.dim_v, self.data.dim_lt,
                self.model.v, self.model.w, lagtimes_try, self.data.list_trans, self.pbc, self.pull)

            # Metropolis acceptance
            if log_like_try is not None and not np.isnan(log_like_try):  # propagator is well behaved
                dlog = log_like_try - self.log_like
                r = np.random.random()
                if ( r < np.exp(dlog/self.temp) ): # accept if dlog increases, accept maybe if decreases
                    self.model.timezero = timezero_try
                    self.model.list_lt = lagtimes_try
                    self.nacctimezero += 1.
                    self.nacctimezero_update += 1.
                    self.log_like = log_like_try

    def mcmove_potential(self):
        # propose temporary v vector: vt
        if self.model.ncosF == 1:
            # if by accident I try to update a flat basis function
            # but this should never happen
            index = 0
            vt = copy.deepcopy(self.model.v)
            coefft = copy.deepcopy(self.model.v_coeff)
            log_like_try = self.log_like

        elif self.model.ncosF <= 0:
            # FIRST
            index = np.random.randint(0,self.model.dim_v)
            vt = copy.deepcopy(self.model.v)          # temporary v
            vt[index] += self.dv * (np.random.random()-0.5)
            # SECOND   #TODO
            #index = np.random.randint(0,self.model.dim_v)
            #vt = self.model.v + self.dv * (np.random.random()-0.5) *self.string_vecs[:,index]
        else:
            index = np.random.randint(1,self.model.ncosF)  # important: I skip the first flat basis function
            coefft = copy.deepcopy(self.model.v_coeff)
            coefft[index] += self.dv * (np.random.random()-0.5)
            vt = self.model.calc_profile(coefft, self.model.v_basis)

        log_like_try = log_like_lag(self.model.dim_v, self.data.dim_lt,
                vt, self.model.w, self.model.list_lt, self.data.list_trans, self.pbc, self.pull)

        # Metropolis acceptance
        if log_like_try is not None and not np.isnan(log_like_try):  # propagator is well behaved
            dlog = log_like_try - self.log_like
            r = np.random.random()
            if r < np.exp(dlog/self.temp): # accept if dlog increases, accept maybe if decreases
                self.model.v[:] = vt[:]
                if self.model.ncosF > 0:
                    self.model.v_coeff[:] = coefft[:]
                    self.naccv_coeff[index] += 1
                self.naccv += 1
                self.naccv_update += 1
                self.log_like = log_like_try
        if False:
            self.check_propagator(self.model.list_lt[0])
            print("loglike",self.log_like)

    def mcmove_diffusion(self):
        # propose temporary w vector: wt
        if self.model.ncosD <= 0:
            if self.k > 0:
                index = np.random.randint(0,self.model.dim_w)   # TODO what if string_vecs has different dimension???
                wt = self.model.w + self.dw * (np.random.random()-0.5) *self.string_vecs[:,index]
            else:
                index = np.random.randint(0,self.model.dim_w)
                wt = copy.deepcopy(self.model.w)          # temporary w
                wt[index] += self.dw * (np.random.random()-0.5)
        else:
            index = np.random.randint(0,self.model.ncosD)
            coefft = copy.deepcopy(self.model.w_coeff)
            coefft[index] += self.dw * (np.random.random()-0.5) 
            wt = self.model.calc_profile(coefft, self.model.w_basis)

        log_like_try = log_like_lag(self.model.dim_v, self.data.dim_lt,
                self.model.v, wt, self.model.list_lt, self.data.list_trans, self.pbc, self.pull)

        if log_like_try is not None and not np.isnan(log_like_try):   # propagator is well behaved
            # add restraints to smoothen
            if self.k > 0.:
                E_wt = string_energy(wt,self.k,self.pbc)
                log_like_try -= E_wt  # minus sign because surface=log_like
    
            # Metropolis acceptance
            dlog = log_like_try - self.log_like
            r = np.random.random()  #in [0,1[
            if r < np.exp(dlog/self.temp): # accept if dlog increases, accept maybe if decreases
                self.model.w[:] = wt[:]
                if self.model.ncosD > 0:
                    self.model.w_coeff[:] = coefft[:]
                    self.naccw_coeff[index] += 1
                self.naccw += 1
                self.naccw_update += 1
                self.log_like = log_like_try
        if False:
            self.check_propagator(self.model.list_lt[0])
            print("loglike",self.log_like)

    def mcmove_diffusion_radial(self):
        # propose temporary wrad
        if self.model.ncosDrad <= 0:
            index = np.random.randint(0,self.model.dim_wrad)
            wradt = copy.deepcopy(self.model.wrad)          # temporary wrad
            wradt[index] += self.dwrad * (np.random.random()-0.5)

        else:
            index = np.random.randint(0,self.model.ncosDrad)
            coefft = copy.deepcopy(self.model.wrad_coeff)
            coefft[index] += self.dwrad * (np.random.random()-0.5)
            wradt = self.model.calc_profile(coefft, self.model.wrad_basis)

        log_like_try = rad_log_like_lag(self.model.dim_v, self.model.dim_rad, self.data.dim_lt, self.model.rate, 
                 wradt, self.data.list_lt, self.data.list_trans, self.model.redges,self.lmax,self.model.bessel0_zeros,self.model.bessels, 0. )

        #print "dlog",log_like_try - self.log_like
        # Metropolis acceptance
        if log_like_try is not None and not np.isnan(log_like_try):   # propagator is well behaved  TODO implement
            dlog = log_like_try - self.log_like
            #print "dlog",dlog,log_like_try
            r = np.random.random()  #in [0,1[
            #if dlog > 0: print "aha",
            #print "dlog",dlog,self.log_like,log_like_try 
            if r < np.exp(dlog/self.temp): # accept if dlog increases, accept maybe if decreases
                #print "accpet"
                self.model.wrad = wradt
                if self.model.ncosDrad > 0:
                    #print self.model.wrad_coeff
                    self.model.wrad_coeff = coefft
                    self.naccwrad_coeff[index] += 1
                self.naccwrad += 1
                self.naccwrad_update += 1
                self.log_like = log_like_try
        else: print("WARNING: log_like_try behaves badly")


    def check_propagator(self,lagtime):
        import scipy
        rate = init_rate_matrix(self.model.dim_v,self.model.v,self.model.w,self.model.pbc,self.pull)
        vals,vecs = np.linalg.eig(rate)
        line = ""
        for v in vals:
            if v.imag < 1e-10: VAL=v.real
            else: VAL = v
            line += str(VAL)+" "

        propagator = scipy.linalg.expm(lagtime*rate)
        vals,vecs = np.linalg.eig(propagator)
        line2 = ""
        for v in vals:
            if v.imag < 1e-10: VAL=v.real
            else: VAL = v
            line2 += str(VAL)+" "
        tiny = 1e-10
        count = np.sum(propagator<tiny)
        #log_like = np.float64(0.0)  # use high precision
        #b = transition[ilag,:,:]*np.log(propagator.clip(tiny))
        #log_like += np.sum(b)
        print("count",count)
        print("ratematrix",line)
        print("propagatormatrix",line2)


    #======== UPDATE MC PARAMS ========
    def update_temp(self,imc):
        if self.dtemp != 0.:
          if self.num_MC_update > 0:
            if (imc+1)%self.num_MC_update == 0:
                self.temp += self.dtemp
                #self.temp *= self.fdtemp
                print("new MC temp:", imc, self.temp)

    def update_movewidth(self,imc):
        """adapt dv and dw such that acceptance ratio stays around 30 procent, or so"""  # TODO
        if self.num_MC_update > 0:
            if ( (imc+1) % self.num_MC_update == 0 ):
                if self.do_radial:
                    self.dwrad *= np.exp ( 0.1 * ( float(self.naccwrad_update) / self.num_MC_update - 0.3 ) )
                    #print "R",float(self.naccwrad_update) / self.num_MC_update
                    self.naccwrad_update = 0
                else:
                    if self.model.ncosF != 1:  # if I am not sampling one flat basisfunction
                        self.dv *= np.exp ( 0.1 * ( float(self.naccv_update) / self.num_MC_update - 0.3 ) )
                    self.dw *= np.exp ( 0.1 * ( float(self.naccw_update) / self.num_MC_update - 0.3 ) )
                    self.naccv_update = 0
                    self.naccw_update = 0
                print("new MC steps:", imc, self.dv, self.dw, self.dwrad)


    #======== PRINTING ========
    def print_MC_params(self,f=None,final=False):
        if f is None:
            import sys
            f = sys.stdout
        if final: print("----- final Settings MC -----", file=f)
        else:     print("----- Settings MC -----", file=f)
        print("dv(MC-potential)=", self.dv, file=f)
        print("dw(MC-logD)=", self.dw, file=f)
        print("dwrad(MC-logDrad)=", self.dwrad, file=f)
        print("temp=", self.temp, file=f)
        print("n(MC)=", self.nmc, file=f)
        print("n(update)=", self.num_MC_update, file=f)
        print("k=", self.k, file=f)
        print("-"*20, file=f)

    def print_intermediate(self,imc,printfreq):
        step = imc+1
        if  (imc%printfreq == 0) | (step == self.nmc):
            print(imc, self.log_like, float(self.naccv)/step, float(self.naccw)/step, float(self.naccwrad)/step)

    def print_log_like(self):
        # will only work if I actually filled it in
        print("===== log_like =====")
        for i in range(self.nmc/20):
             print(" ".join([str(val) for val in self.all_log_like[20*i:20*(i+1)]]))
        print("="*10)

    def print_statistics(self,f=None,):
        if f is None:
           import sys
           f = sys.stdout
        print("===== Statistics =====", file=f)
        print("nmc         ", self.nmc, file=f)
        print("naccv       ", self.naccv, file=f)
        print("naccw       ", self.naccw, file=f)
        print("naccwrad    ", self.naccwrad, file=f)
        print("nacctimezero", self.nacctimezero, file=f)
        print("accv ratio       ", "%5.1f" %(float(self.naccv)/self.nmc*100),"%", file=f)
        print("accw ratio       ", "%5.1f" %(float(self.naccw)/self.nmc*100),"%", file=f)
        print("accwrad ratio    ", "%5.1f" %(float(self.naccwrad)/self.nmc*100),"%", file=f)
        print("acctimezero ratio", "%5.1f" %(float(self.nacctimezero)/self.nmc*100),"%", file=f)
        print("="*10, file=f)
        if self.model.ncosF > 0:
            tot = max(1,np.sum(self.naccv_coeff))  # if all val are zero and sum is zero, then take 1
            print("naccv_coeff", file=f)
            for i,val in enumerate(self.naccv_coeff):
                print("%8d %8d %5.1f %s %5.1f %s" %(i,val,float(val)/tot*100,"%",float(val)/self.nmc*100,"%"), file=f)
        if self.model.ncosD > 0:
            tot = max(1,np.sum(self.naccw_coeff))
            print("naccw_coeff", file=f)
            for i,val in enumerate(self.naccw_coeff):
                print("%8d %8d %5.1f %s %5.1f %s" %(i,val,float(val)/tot*100,"%",float(val)/self.nmc*100,"%"), file=f)
        if self.do_radial:
          if self.model.ncosDrad > 0:
            tot = max(1,np.sum(self.naccwrad_coeff))
            print("naccwrad_coeff", file=f)
            for i,val in enumerate(self.naccwrad_coeff):
                print("%8d %8d %5.1f %s %5.1f %s" %(i,val,float(val)/tot*100,"%",float(val)/self.nmc*100,"%"), file=f)

    def print_coeffs_laststate(self,f,final=False):
        """print basis functions and other model parameters"""

        from mcdiff.log import print_coeffs
        v_coeff = None
        w_coeff = None
        wrad_coeff = None
        timezero = None

        if self.model.ncosF>0: v_coeff = self.model.v_coeff
        if self.model.ncosD>0: w_coeff = self.model.w_coeff
        if self.do_radial:
            if self.model.ncosDrad>0: wrad_coeff = self.model.wrad_coeff        
        if self.move_timezero>0: timezero = self.model.timezero

        print_coeffs(f,self.model,v_coeff,w_coeff,wrad_coeff,timezero,final=final)


    def print_laststate(self,f,final=False): 
        """print final results (potential and diffusion coefficient)
        f is a writable object"""

        self.print_MC_params(f,final=final)
        self.print_coeffs_laststate(f,final=final)

        from mcdiff.log import print_profiles

        v = self.model.v
        w = self.model.w
        if self.do_radial: wrad = self.model.wrad
        else: wrad = None
        if self.move_timezero: timezero = self.model.timezero
        else: timezero = None

        print_profiles(f,self.model,v,w,wrad=wrad,final=final)

