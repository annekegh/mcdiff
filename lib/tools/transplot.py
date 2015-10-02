#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# adapted by An Ghysels (August 2012)
#

import numpy as np
import matplotlib
matplotlib.use('Agg')

from ffmc import read_all_transitions



def plot_transition(transition,shift,dtc,filename=None):
    import matplotlib.pyplot as plt
    print "transition",transition.shape
    trans = transition[0,:,:]
    n = len(trans)
    x = range(n)
    X, Y = np.meshgrid(x, x)
    print "plotting",X.shape, Y.shape, trans.shape
    plt.figure()
    plt.contourf(X,Y,trans)
    plt.axis('equal')
    plt.title("lt=%f, dt=%f, dn=%i" %(shift*dtc,dtc,shift))
    if filename==None:
        filename = "contour.%i.png"%shift
    plt.savefig(filename)
    print "file written:",filename

#for shift in [1,10,100,1000,5000,]:
#    argv = ["%i" %shift,"../HexdWat/up1/pbctrans/transitions.nbins48.%i.pbc.dat" %shift]
#    dtc = 1.
#    print "Plotting...", argv[1]
#    print "lt=%f, dt=%f, dn=%i" %(shift*dtc,dtc,shift)
#    lagtimes,transition,bins = read_all_transitions(argv)
#    plot_transition(transition,shift,dtc)

