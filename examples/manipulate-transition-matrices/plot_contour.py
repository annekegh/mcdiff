"""How to plot the transition matrix (square)

usage:
  python plot_contour.py
"""

from mcdiff.tools.ffmc import read_transition_square
from mcdiff.tools.transplot import plot_transition

shift = 20  # lag time, or "dn"
nbins = 100
filename = "../data/HexdWat/transitions.nbins100.20.pbc.A.dat"

print("Plotting...", filename)
trans,bins,Dt = read_transition_square(filename,nbins)
plot_transition(trans,shift)

if True:
    # EXTRA
    def plot_transition_histogram(trans,shift):
        import matplotlib.pyplot as plt
        import numpy as np
        plt.figure()
        x = np.arange(len(trans))
        y = np.sum(trans,axis=0)
        y2 = np.sum(trans,axis=1)
        plt.plot(x,y,"o")
        plt.plot(x,y2)
        print(trans[:6,:6])
        print(y)
        print(y2)
        plt.savefig("contour.%i.hist.png"%shift)

    plot_transition_histogram(trans,shift)

