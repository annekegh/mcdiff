"""Test how to merge two transition matrices
AG, March 14, 2012"""

from mcdiff.transitions import Transitions, merge_transitions

filename1 = "../data/HexdWat/transitions.nbins100.20.pbc.A.dat"
filename2 = "../data/HexdWat/transitions.nbins100.20.pbc.B.dat"

trans = Transitions([filename1,filename2])  # the list may be longer if you want to merge more matrices

outfile = "transitions.merge.dat"
merge_transitions(trans,outfile)

