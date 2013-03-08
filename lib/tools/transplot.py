#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# adapted by An Ghysels (August 2012)
#

import numpy as np
import matplotlib
matplotlib.use('Agg')

from ffmc import read_all_transitions



def plot_transition(transition,shift,dtc):
    import matplotlib.pyplot as plt
    print transition.shape
    trans = transition[0,:,:]
    n = len(trans)
    x = range(n)
    X, Y = np.meshgrid(x, x)
    print X.shape, Y.shape, trans.shape
    plt.figure()
    plt.contour(X,Y,trans)
    plt.axis('equal')
    plt.title("lt=%f, dt=%f, dn=%i" %(shift*dtc,dtc,shift))
    plt.savefig("contour.%i.png"%shift)

for shift in [1,10,100,1000,5000,]:
    argv = ["%i" %shift,"../HexdWat/up1/pbctrans/transitions.nbins48.%i.pbc.dat" %shift]
    dtc = 1.
    print "Plotting...", argv[1]
    print "lt=%f, dt=%f, dn=%i" %(shift*dtc,dtc,shift)
    lagtimes,transition,bins = read_all_transitions(argv)
    plot_transition(transition,shift,dtc)

