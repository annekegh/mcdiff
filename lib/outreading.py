#!/bin/python

import numpy as np

def read_F_D_edges(filename):
    """File contains F and D in the following format
    bin-number  start-of-bin  end-of-bin  F(bin)  D(bin-to-next-bin)
    Lines starting with # are skipped.
    Assume F in units [kBT] and D in units [angstrom**2/ps]."""

    with open(filename,"r") as f:
        F = []
        D = []
        bins_str = []   # beginning of bin
        bins_end = []   # end of bin
    
        startline = "   index  bin-str  bin-end"
        for line in f:
            if line.startswith(startline): 
                break
        for line in f:
            if line.startswith("="):
                break
            if not line.startswith("#"):
                words = line.split()
                bins_str.append(float(words[1]))
                bins_end.append(float(words[2]))
                F.append(float(words[3]))
                if len(words) >= 5:
                    D.append(float(words[4]))
                elif len(words)!=4 and len(words)!=6:
                    raise ValueError("error in the format line"+line)
    edges = bins_str+[bins_end[-1]]  # last bin edge is added
    # this returns three vectors: F values, D values, edges values
    return np.array(F),np.array(D),np.array(edges)

def read_Drad(filename):
    """File contains Drad in the following format
    bin-number  start-of-bin  end-of-bin  Drad(bin)
    Lines starting with # are skipped.
    Assume Drad in units [angstrom**2/ps]."""

    D = []
    bins_str = []   # beginning of bin
    bins_end = []   # end of bin

    with open(filename,"r") as f:
        startline = "   index  bin-str  bin-end  diffusion-coefficient-at[i]"
        for line in f:
            if line.startswith(startline):
                break
        for line in f:
            if line.startswith("="):
                break
            if not line.startswith("#"):
                words = line.split()
                bins_str.append(float(words[1]))
                bins_end.append(float(words[2]))
                if len(words) == 4 or len(words)==5:
                    D.append(float(words[3]))
              #  else:  #if len(words)!=4:
              #      raise ValueError("error in the format line"+line)
    edges = bins_str+[bins_end[-1]]  # last bin edge is added
    return np.array(D),np.array(edges)

def read_many_profiles(list_filename,pic=False):
    """Read from a list of filenames: F,D"""
    if pic:
        F = []
        D = []
        E = []
        Fst = []
        Dst = []
        from mcdiff.log import load_logger
        for filename in list_filename:
            logger = load_logger(filename)
            f,d,edges,fst,dst = read_F_D_edges_logger(logger)
            F.append(f-min(f))
            D.append(d)
            E.append(edges)
            Fst.append(fst)
            Dst.append(dst)
        return F,D,E,Fst,Dst
 
    else:
        F = []
        D = []
        E = []  # edges
        for filename in list_filename:
            f,d,edges = read_F_D_edges(filename)  # three 1D np arrays
            F.append(f-min(f))
            D.append(d)
            E.append(edges)
        return F,D,E,None,None

def read_many_profiles_Drad(list_filename,pic=False):
    """Read from a list of filenames: Drad"""
    if pic:
        Drad   = []
        RE     = []
        Dradst = []
        from mcdiff.log import load_logger
        for filename in list_filename:
            logger = load_logger(filename)
            drad,redges,dradst = read_Drad_logger(logger)
            Drad.append(drad)
            RE.append(redges)
            Dradst.append(dradst)
        return Drad,RE,Dradst
 
    else:
        Drad   = []
        RE     = []
        for filename in list_filename:
            drad,redges = read_Drad(filename)
            Drad.append(drad)
            RE.append(redges)
        return Drad,RE,None


#==============================

def average_profile(prof):
    """average over list prof = [array1, array2,...] or over array prof = a-profile-on-every-line
    each array is a one-dimensional numpy array
    """
    if type(prof) is list:
        if len(prof) == 1:   # no need to take average
            mean = prof[0]
            if mean is not None:
                std  = np.zeros(len(mean),float)   # no std
            else:
                std = None
        else:
            for i in range(1,len(prof)): assert len(prof[i])==len(prof[0])   # all the same length
            profarr  = np.array(prof)         # convert
            mean = np.mean(profarr,0)
            std  = np.std(profarr,0)
    elif len(prof.shape) == 2:  # a numpy array, each row is a profile
        mean = np.mean(prof,0)
        std  = np.std(prof,0)
    else:   # just one profile, no need to take average
        mean = np.zeros(prof.shape)
        mean[:] = prof[:]   # make a copy
        std = np.zeros(len(mean),float)   # no std
    return mean,std

def average_profiles(F,D,Drad,E):
    """average over different profiles
    average profiles are taken (not the coeff)

    F, D, E are lists of profiles, or arrays where every line is a profile
    """
    assert len(F) == len(D)
    #assert len(F) == len(Drad)
    Fmean,Fst = average_profile(F)
    Dmean,Dst = average_profile(D)
    Dradmean,Dradst = average_profile(Drad)
    edges = E[0]   # just the first
    return Fmean,Dmean,Dradmean,edges,Fst,Dst,Dradst

    #for i in range(len(E)-1):
    #        assert len(E[i])==len(E[i+1])  # each file has same number of bins
    #        assert (E[i]==E[i+1]).all()  # check if edges are really identical
    #edges = np.array(E[0])   # and now just choose the first one
    #if False:#True:    # TODO do always?
    #    if len(F.shape) == 2:
    #        for i in range(F.shape[1]):
    #            F[:,i] += (-min(F[:,i]) ) #+ i)

    # keep edges a 1D vector ?????XXXX TODO 



def read_Fcoeffs(filename,final=False):
    startline = "===== v_coeff ====="
    if final: startline = "===== final v_coeff ====="

    coeff = []
    with open(filename,"r") as f:
        for line in f:
            if line.startswith(startline):
                break
        for line in f:
            if line.startswith("="):
                break
            if not line.startswith("#"):
                words = line.split()
                coeff.append(float(words[1]))
    return np.array(coeff)

def read_Dcoeffs(filename,final=False):
    startline = "===== w_coeff ====="
    if final: startline = "===== final w_coeff ====="

    coeff = []
    with open(filename,"r") as f:
        for line in f:
            if line.startswith(startline):
                break
        for line in f:
            if line.startswith("="):
                break
            if not line.startswith("#"):
                words = line.split()
                coeff.append(float(words[1]))

    return np.array(coeff)

def read_Dradcoeffs(filename,final=False):
    startline = "===== wrad_coeff ====="
    if final: startline = "===== final wrad_coeff ====="

    coeff = []
    with open(filename,"r") as f:
        for line in f:
            if line.startswith(startline):
                break
        for line in f:
            if line.startswith("="):
                break
            if not line.startswith("#"):
                words = line.split()
                coeff.append(float(words[1]))

    return np.array(coeff)

def read_dv_dw(filename,final=False):
    startline = "----- Settings MC -----"
    if final: startline = "----- final Settings MC -----"

    with open(filename,"r") as f:
        for line in f:
            if line.startswith(startline):
                break
        for line in f:
            if line.startswith("dv"):
                words = line.split()
                dv = float(words[1])
                break
        for line in f:
            if line.startswith("dw"):
                words = line.split()
                dw = float(words[1])
                break
    return dv,dw


#==============================
# using loggers
#==============================

def read_coeff_logger(logger):
    """read coeffs from logger and determine the average
    """
    if logger.model.ncosF > 0:
        v_coeff,v_coeff_st = average_profile(logger.v_coeff)
    else:
        v_coeff = None
        v_coeff_st = None

    if logger.model.ncosD > 0:
        w_coeff,w_coeff_st = average_profile(logger.w_coeff)
    else:
        w_coeff = None
        w_coeff_st = None

    if hasattr(logger.model,"ncosDrad"):
        if logger.model.ncosDrad > 0:
            wrad_coeff,wrad_coeff_st = average_profile(logger.wrad_coeff)
    else:
        wrad_coeff = None
        wrad_coeff_st = None

    # average of timezero
    timezero = np.mean(logger.timezero)
    timezero_st = np.std(logger.timezero)

    return v_coeff,w_coeff,wrad_coeff,v_coeff_st,w_coeff_st,wrad_coeff_st,timezero,timezero_st

def read_F_D_edges_logger(logger):
    """read F,D,edges from logger and determine the average
    average is taken from profiles (not coeffs)
    """
    if logger.model.ncosF <= 0:
        v,vst = average_profile(logger.v)
    else:
        v,vst = logger.average_profile_v_from_coeff()

    if logger.model.ncosD <= 0:
        w,wst = average_profile(logger.w)
        d,dst = average_profile(np.exp(logger.w))   # units are missing
    else:
        w,wst,d,dst = logger.average_profile_w_from_coeff()   # units are missing

    F = v*logger.model.vunit
    D = d*np.exp(logger.model.wunit)
    edges = logger.model.edges
    Fst = vst*logger.model.vunit
    Dst = dst*np.exp(logger.model.wunit)
    return F,D,edges,Fst,Dst

def read_Drad_logger(logger):
    """read Drad from logger and determine the average
    average is taken from profiles (not coeffs)
    """
    if hasattr(logger.model,"ncosDrad"):
        if logger.model.ncosDrad <= 0:
            w,wst = average_profile(logger.wrad)
            d,dst = average_profile(np.exp(logger.wrad))   # units are missing
        else:
            wrad,wradst,drad,dradst = logger.average_profile_wrad_from_coeff()   # units are missing

        Drad = drad*np.exp(logger.model.wradunit)
        redges = logger.model.redges
        Dradst = dradst*np.exp(logger.model.wradunit)
    else:
        Drad = None; redges = None; Dradst = None
    return Drad,redges,Dradst

def read_F_D_edges_logger_individualprofiles(logger):
    if logger.model.ncosF <= 0:
        F = logger.v*logger.model.vunit
        #print "Fshape",F.shape
    else:
        a = np.zeros((logger.nf,logger.model.dim_v))
        for i in range(len(a)):
            a[i,:] = logger.model.calc_profile(logger.v_coeff[i,:],logger.model.v_basis)
        F = a*logger.model.vunit
        #print "Fshape",F.shape

    if logger.model.ncosD <= 0:
        W = logger.w+logger.model.wunit
        D = np.exp(W)   # in angstrom**2/ps
    else:
        a = np.zeros((logger.nf,logger.model.dim_w))
        for i in range(len(a)):
            a[i,:] = logger.model.calc_profile(logger.w_coeff[i,:],logger.model.w_basis)
        W = a+logger.model.wunit
        D = np.exp(W)   # in angstrom**2/ps

    edges = logger.model.edges
    return F,D,edges

def read_Drad_logger_individualprofiles(logger):
    if logger.model.ncosDrad <= 0:
        Wrad = logger.wrad+logger.model.wradunit
        Drad = np.exp(Wrad)   # in angstrom**2/ps
    else:
        a = np.zeros((logger.nf,logger.model.dim_wrad))
        for i in range(len(a)):
            a[i,:] = logger.model.calc_profile(logger.wrad_coeff[i,:],logger.model.wrad_basis)
        Wrad = a+logger.model.wradunit
        Drad = np.exp(Wrad)   # in angstrom**2/ps

    edges = logger.model.edges
    return Drad,edges

