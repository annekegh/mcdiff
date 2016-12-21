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


    def average_profile_v_from_coeff(self,):
        a = np.zeros((self.nf,self.model.dim_v))
        for i in xrange(len(a)):
            a[i,:] = self.model.calc_profile(self.v_coeff[i,:],self.model.v_basis)
        v = np.mean(a,0)
        vst = np.std(a,0)
        return v,vst

    def average_profile_w_from_coeff(self,):
        a = np.zeros((self.nf,self.model.dim_w))
        for i in xrange(len(a)):
            a[i,:] = self.model.calc_profile(self.w_coeff[i,:],self.model.w_basis)
        w = np.mean(a,0)
        wst = np.std(a,0)
        d = np.mean(np.exp(a),0)        # unit is missing
        dst = np.std(np.exp(a),0)       # unit is missing
        return w,wst,d,dst

    def average_profile_wrad_from_coeff(self,):
        a = np.zeros((self.nf,self.model.dim_wrad))
        for i in xrange(len(a)):
            a[i,:] = self.model.calc_profile(self.wrad_coeff[i,:],self.model.wrad_basis)
        wrad = np.mean(a,0)
        wradst = np.std(a,0)
        drad = np.mean(np.exp(a),0)     # unit is missing
        dradst = np.std(np.exp(a),0)    # unit is missing
        return wrad,wradst,drad,dradst

    def print_average(self,model,st=0):
        from outreading import read_F_D_edges_logger, read_Drad_logger, read_coeff_logger
        F,D,edges,Fst,Dst = read_F_D_edges_logger(self)
        Drad,redges,Dradst = read_Drad_logger(self)
        v_coeff,w_coeff,wrad_coeff,v_coeff_st,w_coeff_st,wrad_coeff_st, timezero,timezero_st = read_coeff_logger(self)
        error = [Fst,Dst,Dradst]

        import sys
        self.print_MC_params(f=sys.stdout,final=True)
        print_coeffs(sys.stdout,model,v_coeff,w_coeff,wrad_coeff,timezero,final=True,)
        print_profiles(sys.stdout,model,F,D,Drad,final=True,error=error,unit="notinternal")

    def statistics(self,MC,st=0):

        # TODO I should change the function statistics
        # statistics_print versus statistics no print

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


    def plot_likelihood(self,figname):
        import matplotlib.pyplot as plt
        plt.figure()
        r = np.arange(len(self.log_like))*self.freq
        plt.plot(r,self.log_like)
        plt.savefig(figname)


class CombinedLogger(Logger):
    def __init__(self, log1, log2, do_radial=False):
        self.do_radial = do_radial
        self.nmc = log1.nmc + log2.nmc
        assert log1.freq==log2.freq, 'CombinedLogger requires Logger instances with equal frequency'
        self.freq = log1.freq
        nf = self.nmc/self.freq + 1
        self.nf = nf
        # arrays
        self.log_like  = np.concatenate([log1.log_like,log2.log_like])
        self.timezero  = np.concatenate([log1.timezero,log2.timezero])
        self.dv        = np.concatenate([log1.dv,log2.dv])
        self.dw        = np.concatenate([log1.dw,log2.dw])
        self.dwrad     = np.concatenate([log1.dwrad,log2.dwrad])
        self.dtimezero = np.concatenate([log1.dtimezero,log2.dtimezero])


        if hasattr(log1,"v"):
            self.v = np.concatenate([log1.v,log2.v])
        else:
            self.v_coeff = np.concatenate([log1.v_coeff,log2.v_coeff])
        if hasattr(log1,"w"):
            self.w = np.concatenate([log1.w,log2.w])
        else:
            self.w_coeff = np.concatenate([log1.w_coeff,log2.w_coeff])
        if do_radial:
            if hasattr(log1,"wrad"):
                self.wrad = np.concatenate([log1.wrad,log2.wrad])
            else:
                self.wrad_coeff = np.concatenate([log1.wrad_coeff,log2.wrad_coeff])


    def statistics(self,model):
        #... MC.model.blabla ==> model.blabla
        # TODO
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
        if model.ncosF <= 0:
            print_vector(self.v,s)
        else:
            vec = np.zeros((self.nf,model.dim_v),float)
            for i in range(len(vec)):
                vec[i,:] = model.calc_profile(self.v_coeff[i,:],model.v_basis)
            print_vector(vec,s)
            print "===== stat v_coeff ====="
            print_vector(self.v_coeff,s)

        # Diffusion profile
        print "===== stat w ====="
        if model.ncosD <= 0:
            print_vector(self.w,s)
        else:
            vec = np.zeros((self.nf,model.dim_w),float)
            for i in range(len(vec)):
                vec[i,:] = model.calc_profile(self.w_coeff[i,:],model.w_basis)
            print_vector(vec,s)
            print "===== stat w_coeff ====="
            print_vector(self.w_coeff,s)

        # Radial diffusion profile
        if self.do_radial:
            print "===== stat wrad ====="
            if model.ncosDrad <= 0:
                print_vector(self.wrad,s)
            else:
                vec = np.zeros((self.nf,model.dim_wrad),float)
                for i in range(len(vec)):
                    vec[i,:] = model.calc_profile(self.wrad_coeff[i,:],model.wrad_basis)
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
    f -- is a writable object
    error -- is a list of arrays [Fst,Dst,Dradst]
    final -- whether flag 'final' should be written
    """

    if final: print >>f,"===== final F D ====="
    else:     print >>f,"===== F D ====="

    # units:
    F = v  # in kBT
    edges = model.edges  # in angstrom
    if unit is "internal":
        D = np.exp(w+model.wunit)  # in angstrom**2 per [unit-lag-times]
        if wrad is not None:
            Drad = np.exp(wrad+model.wradunit)  # in angstrom**2 per [unit-lag-times]
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
                Dradst = np.exp(error[2]+model.wradunit)   #### let's hope this never happens, because this is not really okay ####

    f.write("%8s %8s %8s  %13s %s" % ("index","bin-str","bin-end","potential","diffusion-coefficient(shifted-by-half-bin)\n"))
    if model.pbc: maxi = model.dim_v
    else: maxi = model.dim_v-1
    for i in range(maxi):
        if error is None:
            f.write("%8d %8.3f %8.3f  %13.5e %13.5e\n" %(i,edges[i],edges[i+1],F[i],D[i]) )
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
    # TODO this can become function of logger?
    logger = load_logger(picfilename)

    import sys   # redirect output to chosen file
    sys.stdout = file(datfilename,"w+")
    logger.print_average(logger.model)
    sys.stdout = sys.__stdout__

def combine_picfiles(picfilename1,picfilename2,filename=None,do_radial=False):
    # take two picfiles, make loggers, and combine the loggers
    # if newpicfile is set to a filename, then a new picfile will be written
    # the function returns the combined logger
    logger1 = load_logger(picfilename1)
    logger2 = load_logger(picfilename2)
    logger = CombinedLogger(logger1,logger2,do_radial=do_radial)
    if filename is not None:
        logger.dump(filename)
    return logger

