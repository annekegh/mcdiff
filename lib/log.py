#!/bin/python

import numpy as np

class Logger(object):
    def __init__(self,MC):
        self.nmc = MC.nmc
        self.freq = 100
        nf = MC.nmc/100+1  # if n=1000, then I want 11 stores / if n=999, then I want 10 stores
        self.nf = nf

        # arrays
        self.log_like  = np.zeros((nf),float)
        self.timezero  = np.zeros((nf),float)
        self.dv        = np.zeros((nf),float)
        self.dw        = np.zeros((nf),float)
        self.dtimezero = np.zeros((nf),float)
        #self.Ew        = np.zeros((nf),float)
        if MC.model.ncosF <= 0:
            self.v       = np.zeros((nf,MC.model.dim_v),float)
        else:
            self.v_coeff = np.zeros((nf,MC.model.ncosF),float)
        if MC.model.ncosD <= 0:
            self.w       = np.zeros((nf,MC.model.dim_w),float)
        else:
            self.w_coeff = np.zeros((nf,MC.model.ncosD),float)
        if MC.do_radial:
            self.dwrad = np.zeros((nf),float)
            if MC.model.ncosDrad <= 0:
                self.wrad = np.zeros((nf,MC.model.dim_wrad),float)
            else:
                self.wrad_coeff = np.zeros((nf,MC.model.ncosDrad),float)

    def log(self,j,MC):  # j is counter, counting starts with 1, ends with nmc
        i = j/self.freq   # if n=1000, freq=100, then store 0,100,...,900,1000 
        if i*self.freq == j:
            #print i,j,self.freq,self.v.shape,MC.model.v.shape
            self.log_like[i]  = MC.log_like
            self.timezero[i]  = MC.model.timezero
            self.dv[i]        = MC.dv
            self.dw[i]        = MC.dw
            self.dtimezero[i] = MC.dtimezero
            #self.Ew[i]        = MC.Ew
            if MC.model.ncosF > 0:
                self.v_coeff[i,:] = MC.model.v_coeff
            else:
                self.v[i,:]       = MC.model.v
            if MC.model.ncosD > 0:
                self.w_coeff[i,:] = MC.model.w_coeff
            else:
                self.w[i,:]       = MC.model.w
            if MC.do_radial:
                self.dwrad[i] = MC.dwrad
                if MC.model.ncosDrad > 0:
                    self.wrad_coeff[i,:] = MC.model.wrad_coeff
                else:
                    self.wrad[i,:] = MC.model.wrad

    def prettyprint(self,f):
        #f = file(filename+"2","w+")
        from pprint import pprint
        pprint (vars(self))
        for attr in dir(self):
            print >> f, "obj.%s = %s" % (attr, getattr(self, attr))
        #f.close()

    def dump(self,filename):
        f = file(filename,"w+")
        import pickle
        pickle.dump(self,f)
        f.close()

    def statistics(self,MC,st=0) : #f):
        s = st/self.freq
        # st  --  start (cutting out the first MC steps)
        if s >= self.nf:
            print "WARNING: supposed to skip %i MC steps, i.e. %i frames, but skipped none" %(self.nmc,s)

        def print_vector(vec,s):
            print "VEC",vec.shape
            assert len(vec.shape) == 2
            for i in range(vec.shape[1]):
                print i, np.mean(vec[s:,i]),np.std(vec[s:,i])

        print "===== stat v ====="
        if MC.model.ncosF <= 0:
            print_vector(self.v,s)
        else:
            vec = np.zeros((self.nf,MC.model.dim_v),float)
            for i in range(len(vec)):
                vec[i,:] = MC.model.calc_profile(self.v_coeff[i,:],MC.model.v_basis)
            print_vector(vec,s)
            print "===== stat v_coeff ====="
            print_vector(self.v_coeff,s)

        print "===== stat w ====="
        if MC.model.ncosD <= 0:
            print_vector(self.w,s)
        else:
            vec = np.zeros((self.nf,MC.model.dim_w),float)
            for i in range(len(vec)):
                vec[i,:] = MC.model.calc_profile(self.w_coeff[i,:],MC.model.w_basis)
            print_vector(vec,s)
            print "===== stat w_coeff ====="
            print_vector(self.w_coeff,s)

        if MC.do_radial:
            print "===== stat wrad ====="
            if MC.model.ncosDrad <= 0:
                print_vector(self.wrad,s)
            else:
                vec = np.zeros((self.nf,MC.model.dim_wrad),float)
                for i in range(len(vec)):
                    vec[i,:] = MC.model.calc_profile(self.wrad_coeff[i,:],MC.model.wrad_basis)
                print_vector(vec,s)
                print "===== stat wrad_coeff ====="
                print_vector(self.wrad_coeff)


def load_logger(filename):
    """Load an object that I dumped before"""
    f = open(filename)
    import pickle
    pic = pickle.load(f)
    f.close()
    return pic

