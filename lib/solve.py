#!/bin/python

import numpy as np
import scipy

from .outreading import read_F_D_edges,read_many_profiles
from .utils import init_rate_matrix

def propagate_with_sink(time,F,D,edges,sink_index):
    n = len(F)
    rate = init_rate_matrix(n,F,D,pbc=True)
    # undo PBC
    rate[0,-1] = 0.
    rate[-1,0] = 0.
    rate[0,0] = - rate[1,0]
    rate[-1,-1] = -rate[-2,-1]
    # sink
    rate[-1,-1] = -1000000
    rate[-2,-1] = 0.
    #rate[-1,-2] = 0.
    # fixed
    #rate[:,0] = 0.
    rate[0,:] = 0.
    #rate[0,0] = 1.
    # test
    #print "rate",rate
    #print "test",np.sum(rate,axis=0)

    prop = scipy.linalg.expm2(time*rate)
    #print prop
    init = np.ones(len(prop))
    init[0] = 1.

    profile = np.dot(prop,init)

    return profile

    #f = file("out.time.txt","w+")
    
    #f.close()

import sys
import matplotlib.pyplot as plt
filename = sys.argv[1]
F,D,edges = read_F_D_edges(filename)

def extend_vec(F,):
    n = len(F)
    newF = np.zeros((3*n),float)
    newF[:n] = F[0]
    newF[2*n:] = F[-1]
    newF[n:2*n] = F
    return newF

#F = extend_vec(F)
#D = extend_vec(D)

sink_index = len(F)-1
plt.figure()
for time in np.arange(0,1,0.2):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='blue')
for time in np.arange(1,10,1):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='green')
for time in np.arange(10,100,10.):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='red')
for time in np.arange(100,1000,100):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='black')
for time in np.arange(1000,10000,1000):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='blue')
for time in np.arange(10000,100000,10000):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='green')
for time in np.arange(100000,1000000,100000):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='red')
for time in np.arange(1000000,10000000,1000000):
    profile = propagate_with_sink(time,F,D,edges,sink_index)
    plt.plot(profile,color='black')


plt.savefig("fig_prop.png")




# get flux

def Deff(F,D):
    assert len(F) == len(D)
    a = np.sum(np.exp(-F))
    b = np.sum(np.exp(F)/D)
    return 1./(a*b)

def calc_unitD(edges):
    dx = edges[1]-edges[0]
    dt = 1.  # ps
    unitD = dx**2/dt  # in angstrom**2/ps
    return unitD

print("Diffusion constants")
unitD = calc_unitD(edges)
n = len(D)
print("unitD",unitD)
print("D[0]",D[0]*unitD)
print("D[n/2]",D[n/2]*unitD)
print("Deff",Deff(F,D)*unitD)
