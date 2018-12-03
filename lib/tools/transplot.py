#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')

def plot_transition(trans,shift,figname=None):
    if figname==None:
        figname = "contour.%i.png"%shift

    import matplotlib.pyplot as plt
    n = len(trans)
    x = list(range(n))
    X, Y = np.meshgrid(x, x)
    #print "plotting",X.shape, Y.shape, trans.shape
    plt.figure()
    plt.contourf(X,Y,trans)
    plt.axis('equal')
    plt.xlabel("start")
    plt.ylabel("end")
    plt.title("dn=%i" %shift)
    plt.savefig(figname)
    print("file written:",figname)

    #print "transition",transition.shape
    # in case it is needed to convert the format:
    #trans = transition[0,:,:]
    #plt.title("lt=%f, dt=%f, dn=%i" %(shift*dtc,dtc,shift))
