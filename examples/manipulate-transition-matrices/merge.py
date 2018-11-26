"""Test how to merge two transition matrices
AG, Nov 26, 2018

usage:
   python merge.py

(or: python merge.py > out.merge)
"""

from mcdiff.transitions import Transitions, merge_transitions

print("--------------")
print("merge only two")
print("--------------")
filename1 = "../data/HexdWat/transitions.nbins100.20.pbc.A.dat"
filename2 = "../data/HexdWat/transitions.nbins100.20.pbc.B.dat"
trans = Transitions([filename1,filename2])  # the list may be longer if you want to merge more matrices
outfile = "transitions.merge.AB.dat"
merge_transitions(trans,outfile)

print("----------------------------")
print("merge all replicates A,B,C,D")
print("----------------------------")
list_filenames = ["../data/HexdWat/transitions.nbins100.20.pbc.%s.dat"%letter for letter in "ABCD"]
print list_filenames

trans = Transitions(list_filenames)
outfile = "transitions.merge.dat"
merge_transitions(trans,outfile)

