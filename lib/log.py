#!/bin/python

import numpy as np

class Logger(object):
    def __init__(self,MC):
        self.nmc = MC.nmc
        self.freq = 100
        nf = MC.nmc/self.freq+1  # if n=1000, then I want 11 stores / if n=999, then I want 10 stores
        self.nf = nf

        # arrays
        self.log_like  = np.zeros((nf),float)
        self.timezero  = np.zeros((nf),float)
        self.dv        = np.zeros((nf),float)
        self.dw        = np.zeros((nf),float)
        self.dwrad     = np.zeros((nf),float)
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
            self.dwrad[i]     = MC.dwrad
            self.dtimezero[i] = MC.dtimezero
            #self.Ew[i]        = MC.Ew
            if MC.model.ncosF > 0:
                self.v_coeff[i,:] = MC.model.v_coeff[:]
            else:
                self.v[i,:]       = MC.model.v[:]
            if MC.model.ncosD > 0:
                self.w_coeff[i,:] = MC.model.w_coeff[:]
            else:
                self.w[i,:]       = MC.model.w[:]
            if MC.do_radial:
                if MC.model.ncosDrad > 0:
                    self.wrad_coeff[i,:] = MC.model.wrad_coeff[:]
                else:
                    self.wrad[i,:] = MC.model.wrad[:]

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

    def print_MC_params(self,f=None,final=False):
        if f is None:
            import sys
            f = sys.stdout
        if final: print >>f, "----- final Settings MC -----"
        else:     print >>f, "----- Settings MC -----"
        print >>f, "dv(MC-potential)=", self.dv[-1]
        print >>f, "dw(MC-logD)=", self.dw[-1]
        if hasattr(self,"dwrad"):   # for older versions
            print >>f, "dwrad(MC-logDrad)=", self.dwrad[-1]
        #print >>f, "temp=", self.temp
        print >>f, "n(MC)=", self.nmc
        #print >>f, "n(update)=", self.num_MC_update
        #print >>f, "k=", self.k
        print >>f, "-"*20

    def get_profiles_average(self,model,st=0):
        # st  --  start (cutting out the first MC steps)
        # v,w,wrad  --  profiles
        # vst,wst,wradst  --  errors on profiles
        # v = F in kBT
        # w => D = exp(w*unitw)
        # wrad => Drad = exp(wrad*unitwrad)
        s = st/self.freq
        if s >= self.nf:
            print "WARNING: supposed to skip %i MC steps, i.e. %i frames, but skipped none" %(self.nmc,s)

        # Free energy
        if model.ncosF <= 0:    # I CHANGED THIS??????????????
            v_coeff = None
            F = np.mean(self.v[s:,:],0)
            Fst = np.std(self.v[s:,:],0)
        else:
            v_coeff = np.mean(self.v_coeff[s:,:],0)
            # compute profile
            vec = np.zeros((self.nf,model.dim_v),float)
            for i in range(len(vec)):
                vec[i,:] = model.calc_profile(self.v_coeff[i,:],model.v_basis)
            F = np.mean(vec[s:,:],0)
            Fst = np.std(vec[s:,:],0)

        # Diffusion profile
        if model.ncosD <= 0:
            w_coeff = None
            D0 = np.exp(self.w+model.wunit)
        else:
            w_coeff = np.mean(self.w_coeff[s:,:],0)
            # compute profile
            vec = np.zeros((self.nf,model.dim_w),float)
            for i in range(len(vec)):
                vec[i,:] = model.calc_profile(self.w_coeff[i,:],model.w_basis)            
            D0 = np.exp(vec+model.wunit)
        D = np.mean(D0[s:,:],0)
        Dst = np.std(D0[s:,:],0)

        # Radial diffusion profile
        if hasattr(model,"ncosDrad"):
            if model.ncosDrad <= 0:
                wrad_coeff = None
                D0 = np.exp(self.wrad+model.wradunit)
            else:
                wrad_coeff = np.mean(self.wrad_coeff[s:,:],0)
                # compute profile
                vec = np.zeros((self.nf,model.dim_wrad),float)
                for i in range(len(vec)):
                    vec[i,:] = model.calc_profile(self.wrad_coeff[i,:],model.wrad_basis)
                D0 = np.exp(vec+model.wradunit)
            Drad = np.mean(D0[s:,:],0)
            Dradst = np.std(D0[s:,:],0)

        else:
            wrad_coeff = None
            Drad = None
            Dradst = None

        error = [Fst,Dst,Dradst]

        if hasattr(model,"timezero"): timezero = model.timezero
        else: timezero = None
        return F,D,Drad, error, v_coeff,w_coeff,wrad_coeff, timezero

    def print_average(self,model,st=0):

        F,D,Drad, error, v_coeff,w_coeff,wrad_coeff, timezero = self.get_profiles_average(model,st=st)

        #print self.__dict__
        #v = np.mean(self.v,0)

        import sys
        self.print_MC_params(f=sys.stdout,final=True)
        print_coeffs(sys.stdout,model,v_coeff,w_coeff,wrad_coeff,timezero,final=True,)
        print_profiles(sys.stdout,model,F,D,Drad,final=True,error=error,unit="notinternal")


    def statistics(self,MC,st=0):
        # st  --  start (cutting out the first MC steps)
        s = st/self.freq
        if s >= self.nf:
            print "WARNING: supposed to skip %i MC steps, i.e. %i frames, but skipped none" %(self.nmc,s)

        def print_vector(vec,s):
            # vec has dimension  nsamples-in-MC x len(profile)
            print "VEC",vec.shape
            assert len(vec.shape) == 2
            for i in range(vec.shape[1]):
                print i, np.mean(vec[s:,i]),np.std(vec[s:,i])

        # Free energy
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

        # Diffusion profile
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

        # Radial diffusion profile
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
                print_vector(self.wrad_coeff,s)


#================================================

def load_logger(filename):
    """Load an object that I dumped before"""
    f = open(filename)
    import pickle
    pic = pickle.load(f)
    f.close()
    return pic

#================================================

def print_coeffs(f,model,v_coeff=None,w_coeff=None,wrad_coeff=None,timezero=None,final=False):
    """print basis functions and other model parameters"""
    if model.ncosF>0:
        if final: print >>f,"===== final v_coeff ====="
        else:     print >>f,"===== v_coeff ====="
        for i,val in enumerate(v_coeff):
            print >>f, "%8d %13.5e" %(i,val)
    if model.ncosD>0:
        if final: print >>f,"===== final w_coeff ====="
        else:     print >>f,"===== w_coeff ====="
        print >>f, "%8d %13.5e" %(0,w_coeff[0]+model.wunit)  # only the first needs to be shifted
        for i,val in enumerate(w_coeff[1:]):
            print >>f, "%8d %13.5e" %(i+1,val)
    if timezero is not None:
        if final: print >>f,"===== final timezero ====="
        else:     print >>f,"===== timezero ====="
        print >>f, "%13.5e" %(timezero)
    if wrad_coeff is not None:
      if model.ncosDrad > 0:
        if final: print >>f,"===== final wrad_coeff ====="
        else:     print >>f,"===== wrad_coeff ====="
        print >>f, "%8d %13.5e" %(0,wrad_coeff[0]+model.wradunit)  # only the first needs to be shifted
        for i,val in enumerate(wrad_coeff[1:]):
            print >>f, "%8d %13.5e" %(i+1,val)
    print >>f, "="*10



def print_profiles(f,model,v,w,wrad=None,final=False,error=None,unit="internal"): 
    """print profiles (potential and diffusion coefficients)
    f is a writable object"""
    # error is a list of arrays [Fst,Dst,Dradst]


    if final: print >>f,"===== final F D ====="
    else:     print >>f,"===== F D ====="

    # units:
    F = v  # in kBT
    edges = model.edges  # in angstrom
    if unit is "internal":
        D = np.exp(w+model.wunit)  # in angstrom**2 per [unit-lag-times]
        if wrad is not None:
            Drad = np.exp(wrad+model.wradunit)  # in angstrom**2 per [unit-lag-times]
             # TODO what with 2 pi r
    else:
        D = w
        Drad = wrad

    if error is not None:     # units were already okay!
        Fst = error[0] # in kBT
        Dst = error[1]
        Dradst = error[2]
        if unit is "internal":
            Dst = np.exp(error[1]+model.wunit)
            if wrad is not None:
                Dradst = np.exp(error[2]+model.wradunit)

    f.write("%8s %8s %8s  %13s %s" % ("index","bin-str","bin-end","potential","diffusion-coefficient(shifted-by-half-bin)\n"))
    if model.pbc: maxi = model.dim_v
    else: maxi = model.dim_v-1
    for i in range(maxi):
        if error is None:
            f.write("%8d %8.3f %8.3f  %13.5e %13.5e%s\n" %(i,edges[i],edges[i+1],F[i],D[i]) )
        else:
            f.write("%8d %8.3f %8.3f  %13.5e %13.5e  %13.5e %13.5e\n" %(
                    i,edges[i],edges[i+1],F[i],D[i],Fst[i],Dst[i] ) )
    if not model.pbc:  # some other type of line: only F, no D
        if error is None:
            f.write("%8d %8.3f %8.3f  %13.5e\n" %(model.dim_v-1,edges[-2],edges[-1],F[-1] ) )
        else:
            f.write("%8d %8.3f %8.3f  %13.5e  %13.5e\n" %(
                    i,edges[i],edges[i+1],F[-1],Fst[i],) )
    f.write("#Done\n")

    if wrad is not None:
        """print final results (radial diffusion coefficient)
        f is a writable object"""
        if final: print >>f,"===== final Drad ====="
        else:     print >>f,"===== Drad ====="
        edges = model.edges   # same as bins in 1-D z-direction
        f.write("%8s %8s %8s  %s" % ("index","bin-str","bin-end","diffusion-coefficient-at[i]\n"))
        for i in range(model.dim_wrad):
            #print "i,dim_wrad,len(edges)",i, self.model.dim_wrad, len(edges)
            if error is None:
                f.write("%8d %8.3f %8.3f  %13.5e\n" %(i,edges[i],edges[i+1],Drad[i] ) )
            else:
                f.write("%8d %8.3f %8.3f  %13.5e  %13.5e\n" %(i,edges[i],edges[i+1],Drad[i],Dradst[i] ) )
        f.write("#Done\n")

#================================================

def write_average_from_pic(picfilename,datfilename):
    from outreading import read_F_D_edges_logger

    logger = load_logger(picfilename)
    F,D,edges,Fst,Dst = read_F_D_edges_logger(logger)

    import sys
    sys.stdout = file(datfilename,"w+")
    logger.print_average(logger.model)
    sys.stdout = sys.__stdout__

