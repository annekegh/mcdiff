#!/bin/python

import numpy as np


def read_F_D_edges(filename):
    """File contains F and D in the following format
    bin-number  start-of-bin  end-of-bin  F(bin)  D(bin-to-next-bin)
    Lines starting with # are skipped.
    Assume F in units [kBT] and D in units [angstrom**2/ps]."""

    f = file(filename)
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
            elif len(words)!=4:
                raise ValueError("error in the format line"+line)
    f.close()
    edges = bins_str+[bins_end[-1]]  # last bin edge is added
    # this returns three vectors: F values, D values, edges values
    return np.array(F),np.array(D),np.array(edges)

def read_Drad(filename):
    """File contains Drad in the following format
    bin-number  start-of-bin  end-of-bin  Drad(bin)
    Lines starting with # are skipped.
    Assume Drad in units [angstrom**2/ps]."""

    f = file(filename)
    D = []
    bins_str = []   # beginning of bin
    bins_end = []   # end of bin

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
            if len(words) == 4:
                D.append(float(words[3]))
            else:  #if len(words)!=4:
                raise ValueError("error in the format line"+line)
    f.close()
    edges = bins_str+[bins_end[-1]]  # last bin edge is added
    return np.array(D),np.array(edges)

def read_many_profiles_Drad(list_filename,ave=False,pic=False):
  if pic:
    from mcdiff.log import load_logger
    if len(list_filename) == 1:
        logger = load_logger(list_filename[0])
    else:
        raise ValueError("list_filename should contain one pickled object in current implementation")
    F,D,edges,Fst,Dst = read_F_D_edges_logger(logger)
    Drad,redges,Dradst = read_Drad_logger(logger)
    if ave:
        return F,D,Drad,edges,Fst,Dst,Dradst
    else:
        return F,D,Drad,edges

  else:
 
    F = []
    D = []
    Drad = []
    E = []  # edges
    for filename in list_filename:
        drad,e = read_Drad(filename)
        f,d,e = read_F_D_edges(filename)
        F.append(f)
        D.append(d)
        Drad.append(drad)
        E.append(e)
    F = np.array(F).transpose()
    D = np.array(D).transpose()
    Drad = np.array(Drad).transpose()
    E = np.array(E).transpose()  ###
    #for i in range(len(E)-1):
    #    assert len(E[i])==len(E[i+1])  # each file has same number of bins
    #    assert (E[i]==E[i+1]).all()  # check if edges are really identical
    #edges = np.array(E[0])   # and now just choose the first one
    edges = E
    if True:
        if len(F.shape) == 2:
            for i in range(F.shape[1]):
                F[:,i] += (-min(F[:,i]) ) #+ i)

    if ave:   # average over the different files
        Fst = np.std(F,-1)
        Dst = np.std(D,-1)
        Dradst = np.std(Drad,-1)
        F = np.mean(F,-1)
        D = np.mean(D,-1)
        Drad = np.mean(Drad,-1)
        return F,D,Drad,edges,Fst,Dst,Dradst
    else:
        return F,D,Drad,edges

def read_many_profiles(list_filename,pic=False):
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

def average_profiles(F,D,E):
    "average over different profiles"
    # F, D, E are lists of profiles
    assert type(F) is list
    assert type(D) is list
    assert len(F) == len(D)
    if len(F) == 1:
        Fst = np.zeros(len(F[0]),float)
        Fmean = F[0]
    else:
        farr  = np.array(F)
        Fmean = np.mean(farr,0)
        Fst   = np.std(farr,0)
    if len(D) == 1:
        Dmean = D[0]
        Dst = np.zeros(len(D[0]),float)
    else:
        darr = np.array(D)
        Dmean = np.mean(darr,0)
        Dst = np.std(darr,0)
    edges = E[0]
    #for i in range(len(E)-1):
    #        assert len(E[i])==len(E[i+1])  # each file has same number of bins
    #        assert (E[i]==E[i+1]).all()  # check if edges are really identical
    #edges = np.array(E[0])   # and now just choose the first one
    #if False:#True:    # TODO do always?
    #    if len(F.shape) == 2:
    #        for i in range(F.shape[1]):
    #            F[:,i] += (-min(F[:,i]) ) #+ i)

    # keep edges a 1D vector ?????XXXX TODO 
    return Fmean,Dmean,edges,Fst,Dst


def read_Fcoeffs(filename,final=False):
    startline = "===== v_coeff ====="
    if final: startline = "===== final v_coeff ====="

    coeff = []
    f = file(filename)
    for line in f:
        if line.startswith(startline):
            break
    for line in f:
        if line.startswith("="):
            break
        if not line.startswith("#"):
            words = line.split()
            coeff.append(float(words[1]))
    f.close()
    return np.array(coeff)

def read_Dcoeffs(filename,final=False):
    startline = "===== w_coeff ====="
    if final: startline = "===== final w_coeff ====="

    coeff = []
    f = file(filename)
    for line in f:
        if line.startswith(startline):
            break
    for line in f:
        if line.startswith("="):
            break
        if not line.startswith("#"):
            words = line.split()
            coeff.append(float(words[1]))
    f.close()

    return np.array(coeff)

def read_Dradcoeffs(filename,final=False):
    startline = "===== wrad_coeff ====="
    if final: startline = "===== final wrad_coeff ====="

    coeff = []
    f = file(filename)
    for line in f:
        if line.startswith(startline):
            break
    for line in f:
        if line.startswith("="):
            break
        if not line.startswith("#"):
            words = line.split()
            coeff.append(float(words[1]))
    f.close()

    return np.array(coeff)

def read_dv_dw(filename,final=False):
    startline = "----- Settings MC -----"
    if final: startline = "----- final Settings MC -----"

    f = file(filename)
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
    f.close()
    return dv,dw

def read_F_D_edges_logger(logger):
    if logger.model.ncosF <= 0:
        v = np.mean(logger.v,0)
        vst = np.std(logger.v,0)
    else:
        a = np.zeros((logger.nf,logger.model.dim_v))
        for i in xrange(len(a)):
            a[i,:] = logger.model.calc_profile(logger.v_coeff[i,:],logger.model.v_basis)
        v = np.mean(a,0)
        vst = np.std(a,0)

    if logger.model.ncosD <= 0:
        w = np.mean(logger.w,0)
        wst = np.std(logger.w,0)
        dst = np.std(np.exp(logger.w),0)
    else:
        a = np.zeros((logger.nf,logger.model.dim_w))
        for i in xrange(len(a)):
            a[i,:] = logger.model.calc_profile(logger.w_coeff[i,:],logger.model.w_basis)
        w = np.mean(a,0)
        wst = np.std(a,0)
        dst = np.std(np.exp(a),0)

    F = v*logger.model.vunit
    W = w+logger.model.wunit
    D = np.exp(W)   # in angstrom**2/ps
    edges = logger.model.edges
    Fst = vst*logger.model.vunit
    Wst = wst+logger.model.wunit
    #Dst = D*(np.exp(Wst)-1)
    Dst = dst*np.exp(logger.model.wunit)
    return F,D,edges,Fst,Dst

def read_Drad_logger(logger):
    if logger.model.ncosDrad <= 0:
        wrad = np.mean(logger.wrad,0)
        wradst = np.std(logger.wrad,0)
    else:
        a = np.zeros((logger.nf,logger.model.dim_wrad))
        for i in xrange(len(a)):
            a[i,:] = logger.model.calc_profile(logger.wrad_coeff[i,:],logger.model.wrad_basis)
        wrad = np.mean(a,0)
        wradst = np.std(a,0)
        dradst = np.std(np.exp(a),0)

    Wrad = wrad+logger.model.wradunit
    Drad = np.exp(Wrad)   # in angstrom**2/ps
    redges = logger.model.redges
    Wradst = wradst+logger.model.wradunit
    #Dst = D*(np.exp(Wst)-1)
    Dradst = dradst*np.exp(logger.model.wradunit)
    return Drad,redges,Dradst

