#!/bin/python

import numpy as np

class Logger(object):
    def __init__(self,MC):
        n = MC.nmc
        self.nmc = n
        self.freq = 100
        nf = int(n/100)
        self.nf = nf

        self.edges = MC.model.edges
        self.wunit = MC.model.wunit
        # arrays
        self.v         = np.zeros((nf,MC.model.dim_v),float)
        self.w         = np.zeros((nf,MC.model.dim_w),float)
        self.log_like  = np.zeros((nf),float)
        self.timezero  = np.zeros((nf),float)
        self.dv        = np.zeros((nf),float)
        self.dw        = np.zeros((nf),float)
        self.dtimezero = np.zeros((nf),float)
        #self.Ew        = np.zeros((nf),float)
        if MC.model.ncosF > 0:
            self.v_coeff   = np.zeros((nf,MC.model.ncosF),float)
        if MC.model.ncosD > 0:
            self.w_coeff   = np.zeros((nf,MC.model.ncosD),float)
        if MC.do_radial:
          if MC.model.ncosDrad > 0:
            self.wrad_coeff   = np.zeros((nf,MC.model.ncosDrad),float)


    def log(self,j,MC):
        i = j/self.freq
        self.v[i,:]       = MC.model.v
        self.w[i,:]       = MC.model.w
        self.log_like[i]  = MC.log_like
        self.timezero[i]  = MC.model.timezero
        self.dv[i]        = MC.dv
        self.dw[i]        = MC.dw
        self.dtimezero[i] = MC.dtimezero
        #self.Ew[i]        = MC.Ew
        if MC.model.ncosF > 0:
            self.v_coeff[i,:] = MC.model.v_coeff
        if MC.model.ncosD > 0:
            self.w_coeff[i,:] = MC.model.w_coeff
        if MC.do_radial:
          if MC.model.ncosDrad > 0:
            self.wrad_coeff[i,:] = MC.model.wrad_coeff

    def dump(self,filename):
        #f = file(filename+"2","w+")
        #from pprint import pprint
        #pprint (vars(self))
        #for attr in dir(self):
        #    print >> f, "obj.%s = %s" % (attr, getattr(self, attr))
        #f.close()
        f = file(filename,"w+")
        import pickle
        pickle.dump(self,f)
        f.close()

    def statistics(self,MC,st=0) : #f):
        s = st/self.nf
        # st  --  start (cutting out the first MC steps)
        print "===== stat v ====="
        for i in range(self.v.shape[1]):
            print i, np.mean(self.v[s:,i]),np.std(self.v[s:,i])
        print "===== stat w ====="
        for i in range(self.w.shape[1]):
            print i, np.mean(self.w[s:,i]),np.std(self.w[s:,i])
        if MC.model.ncosF > 0:
            print "===== stat v_coeff ====="
            for i in range(self.v_coeff.shape[1]):
                print i, np.mean(self.v_coeff[s:,i]),np.std(self.v_coeff[s:,i])
        if MC.model.ncosD > 0:
            print "===== stat v_coeff ====="
            for i in range(self.w_coeff.shape[1]):
                print i, np.mean(self.w_coeff[s:,i]),np.std(self.w_coeff[s:,i])
        if MC.do_radial:
          if MC.model.ncosDrad > 0:
            print "===== stat w_coeff ====="
            for i in range(self.wrad_coeff.shape[1]):
                print i, np.mean(self.wrad_coeff[s:,i]),np.std(self.wrad_coeff[s:,i])



def load_logger(filename):
    """Load an object that I dumped before"""
    f = open(filename)
    import pickle
    pic = pickle.load(f)
    f.close()
    return pic


