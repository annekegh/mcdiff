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
  print "pic",pic
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
    for i in range(len(E)-1):
        assert len(E[i])==len(E[i+1])  # each file has same number of bins
        assert (E[i]==E[i+1]).all()  # check if edges are really identical
    edges = np.array(E[0])   # and now just choose the first one
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

def read_many_profiles(list_filename,ave=False,pic=False):
  if pic:
    from mcdiff.log import load_logger
    if len(list_filename) == 1:
        logger = load_logger(list_filename[0])
    else:
        raise ValueError("list_filename should contain one pickled object in current implementation")
    F,D,edges,Fst,Dst = read_F_D_edges_logger(logger)
    F -= min(F)
    if ave:
        print "doing average"
        return F,D,edges,Fst,Dst
    else:
        return F,D,edges

  else:
    F = []
    D = []
    E = []  # edges
    for filename in list_filename:
        f,d,e = read_F_D_edges(filename)
        F.append(f)
        D.append(d)
        E.append(e)
    F = np.array(F).transpose()
    D = np.array(D).transpose()
    for i in range(len(E)-1):
        assert len(E[i])==len(E[i+1])  # each file has same number of bins
        assert (E[i]==E[i+1]).all()  # check if edges are really identical
    edges = np.array(E[0])   # and now just choose the first one
    if True:
        if len(F.shape) == 2:
            for i in range(F.shape[1]):
                F[:,i] += (-min(F[:,i]) ) #+ i)

    if ave:   # average over the different files
        Fst = np.std(F,-1)
        Dst = np.std(D,-1)
        F = np.mean(F,-1)
        D = np.mean(D,-1)
        return F,D,edges,Fst,Dst
    else:       
        return F,D,edges


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
    F = np.mean(logger.v,0)
    W = np.mean(logger.w,0)
    wunit = logger.model.wunit
    D = np.exp(W+wunit)   # in angstrom**2/ps
    edges = logger.model.edges
    Fst = np.std(logger.v,0)
    Wst = np.std(logger.w,0)
    Dst = D*(np.exp(Wst)-1)
    return F,D,edges,Fst,Dst

def read_Drad_logger(logger):
    wradunit = logger.model.wradunit
    Wrad = np.mean(logger.wrad,0)
    Drad = np.exp(Wrad+wradunit)
    redges = logger.model.redges
    Wradst = np.std(logger.wrad,0)
    Dradst = Drad*(np.exp(Wradst)-1)

#  TODO  fix
    wrad_coeff = np.mean(logger.wrad_coeff,0)
    Wrad = logger.model.calc_profile(wrad_coeff,logger.model.wrad_basis)
    Drad = np.exp(Wrad+wradunit)

    wrad_coeffst = np.std(logger.wrad_coeff,0)
    wradpluswradst = logger.model.calc_profile(wrad_coeff+wrad_coeffst,logger.model.wrad_basis)
    Wradst = wradpluswradst - Wrad
    Dradst = Drad*(np.exp(Wradst)-1)

    return Drad,redges,Dradst

