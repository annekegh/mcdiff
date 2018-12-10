"""Extract transition matrix from trajectory, e.g. in ucomo2.dat format
AG, August 19, 2012
AG, December 10, 2018, cleaned up to be added to examples."""

from mcdiff.tools.extract import write_Tmat_square, transition_matrix_add1
from mcdiff.tools.histogram import read_coor
import numpy as np

# setup
#-------
indir = "./"
outdir = "pbctrans/"
dt = 1.  # ps  storage between coord
zpbc = 67.92547  # depends on the system
nbins = 100
dx = zpbc/nbins
edges = np.arange(-nbins*dx/2,(nbins+1)*dx/2.,dx)
assert len(edges) == nbins+1


# read trajectory
#-----------------
z = read_coor(indir+"ucomo2.dat",rv=True,axis=2) # com=True
# axis = 0,1,2 refers to x,y,z coordinate
# rv = whether the ucomo2.dat file format is used, see below.
# com = True (default) means that the COM-of-lipid columns are skipped

# In case you have a different format, you could use your own
# reading function here.
# This example trajectory only has 200 timesteps, very short.

# construct transition matrix for different lag times
#----------------------------------------------------
for shift in [1,2,3,4,5,10,15,20,30,40,50]:
    lt = shift*dt  # lag time in ps
    print("*"*20+"\nshift "+str(shift)+" lagtime "+str(lt)+" ps")

    A = np.zeros((nbins+2,nbins+2),int)
    print("nbins, len(edges), shape(A), z.shape")
    print(nbins, len(edges), A.shape, z.shape)
  
    # the N centers of mass of the N O2 molecules,
    # or the 2*N atoms of the N O2 molecules
    for i in range(z.shape[1]):
        traj = z[:,i]-zpbc*np.floor(z[:,i]/zpbc+0.5) # put z=0 in the center
        A = transition_matrix_add1(A,traj,edges,shift=shift)

    count = "pbc"  
    Tmatfile = "%s/transitions.nbins%i.%i.%s.dat"%(outdir,nbins,shift,count)
    write_Tmat_square(A[1:-1,1:-1],Tmatfile,lt,count,edges=edges,dt=dt,dn=shift)

#----------------------------------------------------

"""format of ucomo2.dat file

e.g., for a system with 2 O2 molecules (or 2 O atoms),
every line contains:

timestep  COM[x] COM[y] COM[z] O[x] O[y] O[z] O[x] O[y] O[z]

where
timestep = timestep (in ps or ..., this is not read currently)
COM = center of mass of the lipids (assumed centered close to z=0)
O[x] O[y] O[z] = coords of center of mass of an oxygen molecule,
                 with respect to COM (!), so assume that shifting
                 has already happened

or with  (this amounts to the same thing)
e.g., for a system with 2 oxygen atoms, the meaning could be,

O[x] O[y] O[z] = coords of an oxygen atom,
                 with respect to COM (!), so assume that shifting
                 has already happened
"""

