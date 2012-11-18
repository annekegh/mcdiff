"""read trajectory and compare
directly with non-discretized propagator"""

import numpy as np
import scipy
import scipy.linalg

def read_traj(filename):
    traj = []
    f = file(filename)
    for line in f:
        if not line.startswith("#"):
            traj.append(float(line))
    f.close()
    return np.array(traj)

def calc_propagator_basis(time,coeffsF,coeffsD,coeffsP,freqsF,freqsD,freqsP):
    # actually, I know the freqs because of symmetry
    zpbc = 50.
    freqs = 2*np.pi/zpbc
    # assuming cos basisfunctions
    # freqs are basisset frequencies
    coeffs_dF_dz = -freqsF*coeffsF
    coeffs_dP_dz = -freqsD*coeffsP
    return prop

# TODO
def calc_overlap_basis(p_basis):
    L = len(p_basis)
    overlap = np.zeros((L,L),float)
    for i in xrange(L):
        for j in np.arange(j,L):
            a = np.sum(p_basis[:,i]*p_basis[:,j])
            overlap[i,j] = a
            overlap[j,i] = a
    return overlap

def calc_rhs_diffusionequation_basisfunctions(p_basis,rate):
    L = len(p_basis)
    rhs = np.zeros((L,L),float)
    for i in xrange(L):
        for j in xrange(L):
            a = np.sum(p_basis[:,i]*p_basis[:,j]*rate[:,j])  # I can do this better  TODO might be wrong
            rhs[i,j] = a
    return rhs


def log_likelihood_basis(transition,p_basis,rate,lagtime):
    overlap = calc_overlap_basis(p_basis)
    rhs = calc_rhs_diffusionequation_basisfunctions(p_basis,rate)
    vals,vecs = np.linalg.eig(np.dot(np.inv(overlap),rhs))
    # lagtime:
    prop = np.sum(vecs*np.exp(lagtime*vals),axis=X)
    pass
    # TODO

#---------------------------------------------------------------

def add_Tmat_cos(ncos,T,coor,lt,):
    # coor in units L
    l = len(coor)
    assert lt < l
    coor1 = coor[:-lt]
    coor2 = coor[lt:]
    ran = np.arange(-nc,nc+1)
    for i in range(len(coor1)):
        vec1 = ( 2*np.pi*coor1[i]*ran ).reshape(1,-1)
        vec2 = ( 2*np.pi*coor2[i]*ran ).reshape(-1,1)
        #print "ran",len(ran),"vec1",vec1.shape,"vec2",vec2.shape
        #Tre += np.dot( np.cos(vec2), np.cos(vec1))
        #Tim += np.dot( np.sin(vec2), np.sin(vec1))
    #T += Tre+Tim*1j

    T += np.dot(np.exp(-vec2*1j),np.exp(-vec1*1j))
    return T

def init_rate_cos(nc,ncos,F,D):
    # assume F = [F0,F1,F2,...,Fnc]
    # assume D = [D0,D1,D2,...,Dnc]
    rate = np.zeros((ncos,ncos),complex)
    for i,k in enumerate(range(-nc,nc+1,1)):
        for j,l in enumerate(range(-nc,nc+1,1)):
            #print "i,k",i,k,"j,l",j,l
            d = abs(l-k)
            if d <= nc:
                rate[i,j] += -k*l*D[i-j]   #D[d]
            for r,m in enumerate(range(-nc,nc+1)):
                d = abs(l-m+k)
                if d <= nc:  #we already know that abs(m) is <= nc 
                    rate[i,j] += -k*m*F[r]*D[j-r+i]   #D[d]
    # rate in units (2 pi / L)**2
    return rate

def log_likelihood_cos(ncos,rate,lt,T):
    prop = scipy.linalg.expm2(rate*lt)
    print "prop",prop
    tiny = 1e-20
    log_like = np.sum( np.where(abs(prop)>tiny, T * np.log(prop),0 ))  # Tim
    return log_like

def maketraj(filename):
    dx = 0.1  # in units L
    nmc = 100000
    x0 = 0.
    f = file(filename,"w+")
    for i in range(nmc):
        print >> f, x0-np.floor(x0)  # pbc: in [0,1]   # save every ps
        x0 += (np.random.rand()-0.5)*dx
    f.close()

filename = "dyn.z"
#maketraj(filename)
coor = read_traj(filename)
nc = 1
ncos = 2*nc+1
lt = 20  # in ps
print "Settings"
print "nc",nc,"ncos",ncos
print "lt",lt
print "len coor",len(coor)

T = np.zeros((ncos,ncos),complex)
add_Tmat_cos(ncos,T,coor,lt)


#F = np.zeros(nc+1)
#D = np.zeros(nc+1)
F = np.zeros(ncos)  # in units kBT times L  # the Fourier transform
D = np.zeros(ncos)  # in units D   times L
logs = []
vals = [0,1e-10,0.00001,0.0001,0.001,0.002,0.005,0.01,0.02,0.03,0.04,0.05,0.1,1.,2.,3.]
for val in vals:
    D[nc+1] = val
    F[nc+1]=0.1
    
    
    rate = init_rate_cos(nc,ncos,F,D)
    print "rate",rate
    print "eig rate", np.linalg.eig(rate)[0]
    
    log_like = log_likelihood_cos(ncos,rate,lt,T)
    logs.append(log_like)
    print "log_like",log_like

from pprint import pprint
pprint(vals)
pprint( logs)
print abs(np.array(logs))
