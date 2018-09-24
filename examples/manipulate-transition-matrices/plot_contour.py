"""How to plot the transition matrix (square)"""

from mcdiff.tools.ffmc import read_transition_square
from mcdiff.tools.transplot import plot_transition

shift = 20  # lag time, or "dn"
nbins = 100
filename = "../data/HexdWat/transitions.nbins100.20.pbc.A.dat"

print "Plotting...", filename
trans,bins,Dt = read_transition_square(filename,nbins)
plot_transition(trans,shift)

