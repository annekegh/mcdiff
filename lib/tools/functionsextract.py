"""Take trajectory and create transition matrix
AG, August 19, 2012"""

import numpy as np

def count_2D(B,X,Y,Z,edges,redges,shift=1):
    print B.shape
    assert len(B.shape) == 3
    assert B.shape[0] == len(redges)
    assert B.shape[1] == len(edges)+1
    assert B.shape[2] == len(edges)+1
    
    #Nz = len(edges)-1
    #Nr = len(redges)-1
    digitized = np.digitize(Z,edges)

    zstart = digitized[:-shift]
    zend = digitized[shift:]
    dX = X[shift:] - X[:-shift]
    dY = Y[shift:] - Y[:-shift]
    dR = np.sqrt(dX**2+dY**2)
    dr = np.digitize(dR,redges)

    #print max(digitized)
    #print max(dr)

    assert len(dr) == len(zstart)
    assert len(dr) == len(zend)

#    i=0
    for start,end,r in zip(zstart,zend,dr):
        B[r-1,end,start] += 1 

#        if r > 10:
        #nbins = len(edges)-1
        #D = min( abs(end-start), abs(end-start-nbins), abs(end-start+nbins) )
        #if D >= 5:
#           print "step", i
#           print shift, start, end
           #print stop
#        i+=1

    return B


def read_traj0000(filename):
    f = file(filename)
    data = []
    for line in f:
      if not line.startswith("#"):
      ####  data.append(float(line))
        data.append(float(line.split()[0]))
    return np.array(data)


def write_Tmat_square(A,filename,lt,count,edges=None,dt=None,dn=None):
    """Write the transition matrix counts in fixed format
    lt  --  lag time in ps
    edges  --  bin edges
    """
    f = file(filename,"w+")
    L = len(A)  # number of bins + 1
    print >> f, "#lt   ", lt  # lag time
    print >> f, "#count", count # how transitions were counted

    if dt != None:
        print >> f, "#dt   ", dt
    if dn != None:
        print >> f, "#dn   ", dn
    if edges != None:
        print >> f, "#edges ", " ".join([str(b) for b in edges])

    for i in xrange(L):
        print >> f, " ".join([str(val) for val in A[i,:]])
    f.close()

def write_Tmat_cube(B,filename,lt,count,edges=None,redges=None,dt=None,dn=None):
    """Write the transition matrix counts in fixed format
    lt  --  lag time in ps
    edges  --  bin edges
    redges  --  bin edges of radial bins
    """
    f = file(filename,"w+")
    size_r = B.shape[0]  # number of bins + 1
    size_z = B.shape[1]
    assert B.shape[2] == B.shape[1]
    
    print >> f, "#lt   ", lt  # lag time
    print >> f, "#count", count # how transitions were counted

    if dt != None:
        print >> f, "#dt   ", dt
    if dn != None:
        print >> f, "#dn   ", dn
    if edges != None:
        print >> f, "#edges ", " ".join([str(b) for b in edges])
    if redges != None:
        print >> f, "#redges ", " ".join([str(b) for b in redges])

    for i in xrange(size_r):
      for j in xrange(size_z):
        print >> f, " ".join([str(val) for val in B[i,j,:]])
      print >> f, "-"
    f.close()


def transition_matrix_add1(A,x,edges,shift=1):
    assert len(x.shape) == 1
    assert shift < len(x)
    assert len(A) == len(edges)+1
    digitized = np.digitize(x,edges)
    #if False:  # pbc
    #  for i in digitized:
    #    assert i > 0
    #    assert i < len(edges)
    #print "min,max,A.shape"
    #print min(digitized),max(digitized),A.shape
    #print "min,max",
    #print min(x),max(x)
    # periodic boundary conditions: just checking
    print "check boundary", sum(A[0,:]), sum(A[-1,:]), sum(A[:,0]), sum(A[:,-1])
    #i = 0
    for start,end in zip(digitized[:-shift],digitized[shift:]):
        #nbins = len(edges)-1
        #D = min( abs(end-start), abs(end-start-nbins), abs(end-start+nbins) )
        #if D >= 5:
        #   print "step", i
        #   print shift, start, end
        #   #print stop
        A[end,start]+=1
        #i += 1
    return A

def transition_matrix_add2(A,x,edges,shift=1):
    assert len(x.shape) == 1

    assert (x>edges[0]).all()
    assert (x<edges[-1]).all()
    N = len(edges)-1
    L = edges[-1]-edges[0]
    print edges[0], edges[-1], "L",L, "L/2",L/2.
    digitized = [min(int(N*(L/2.+xi)/L),N-1) for xi in x]
    #print digitized
    print x
    for i in range(len(x)-shift):
        A[digitized[i+shift],digitized[i]] += 1
    return A

# other counting
   #A[-1,1:] += A[0,1:]
    #A[1:,-1] += A[1:,0]
    #A[-1,-1] += A[0,0]
    #if filename is not None:
    #    write_Tmat_square(A[1:-1,1:-1],filename+"."+str(shift)+".pbc.dat")

def cut_transition_square(filename,dim_trans,start,end,outfile):
    """Read transitions from file and cut out the piece start->end,
    so size NxN with N-end-start
    counting starts from 0"""
    from mcdiff.reading import read_transition_square
    transition = read_transition_square(filename,dim_trans)
    header = read_transition_header(filename)

    A = transition[start:end+1,start:end+1]

    write_Tmat_square(A,filename,lt,count,edges=None,dt=None,dn=None)


#=========== OLDER =========
def write_Tmat_linebyline(A,filename,edges=None):
    f = file(filename,"w+")
    L = len(A)  # number of bins + 1
    if edges != None:
        print >> f, "#edges ", " ".join([str(b) for b in edges])
    for i in xrange(L):
        for j in xrange(L):
            print >> f, A[i,j]
    f.close()


