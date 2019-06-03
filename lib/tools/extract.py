"""Take trajectory and create transition matrix
AG, August 19, 2012

Get survival probability
AG, August 21, 2013"""

import numpy as np

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
    L = len(A)  # number of bins + 1
    with open(filename,"w+") as f:
        f.write("#lt    {}\n".format(lt))  # lag time
        f.write("#count {}\n".format(count)) # how transitions were counted
    
        if dt != None:
            f.write("#dt    {}\n".format(dt))
        if dn != None:
            f.write("#dn    {}\n".format(dn))
        if edges is not None:
            f.write("#edges  "+" ".join([str(b) for b in edges])+"\n")
        for i in range(L):
            f.write(" ".join([str(val) for val in A[i,:]])+"\n")

def write_Tmat_cube(B,filename,lt,count,edges=None,redges=None,dt=None,dn=None):
    """Write the transition matrix counts in fixed format
    lt  --  lag time in ps
    edges  --  bin edges
    redges  --  bin edges of radial bins
    """
    size_r = B.shape[0]  # number of bins + 1
    size_z = B.shape[1]
    assert B.shape[2] == B.shape[1]

    with open(filename,"w+") as f:

        f.write("#lt    {}\n".format(lt))  # lag time
        f.write("#count {}\n".format(count)) # how transitions were counted

        if dt != None:
            f.write("#dt    {}\n".format(dt))
        if dn != None:
            f.write("#dn    {}\n".format(dn))
        if edges is not None:
            f.write("#edges  "+" ".join([str(b) for b in edges])+"\n")
        if redges is not None:
            f.write("#redges  "+" ".join([str(b) for b in redges])+"\n")

        for i in range(size_r):
            for j in range(size_z):
                f.write(" ".join([str(val) for val in B[i,j,:]])+"\n")
            f.write("-\n")

#===============================================
#  transition matrix: count the transitions in trajectory
#===============================================

def transition_matrix_add1(A,x,edges,shift=1):
    assert len(x.shape) == 1
    assert shift < len(x)
    assert len(A) == len(edges)+1
    digitized = np.digitize(x,edges)
    #if False:  # pbc
    #  for i in digitized:
    #    assert i > 0
    #    assert i < len(edges)
    #print("min,max,A.shape")
    #print(min(digitized),max(digitized),A.shape)
    #print("min,max",min(x),max(x))
    # periodic boundary conditions: just checking
    print(("check boundary", sum(A[0,:]), sum(A[-1,:]), sum(A[:,0]), sum(A[:,-1])))

    for start,end in zip(digitized[:-shift],digitized[shift:]):
        A[end,start]+=1
    return A

def transition_matrix_add1_npt(A,x,zpbc,nbins,shift=1):
    assert len(x.shape) == 1
    assert shift < len(x)
    assert len(x) == len(zpbc)
    assert len(A) == nbins+2
    arr = np.arange(nbins+1)-nbins/2.
    digitized = np.zeros(x.shape)
    for i in range(len(digitized)):
        digitized[i] = np.digitize([x[i]],(arr*zpbc[i]/nbins))
    # periodic boundary conditions: just checking
    print(("check boundary", sum(A[0,:]), sum(A[-1,:]), sum(A[:,0]), sum(A[:,-1])))

    for start,end in zip(digitized[:-shift],digitized[shift:]):
        A[end,start]+=1
    return A


def transition_matrix_add2(A,x,edges,shift=1):
    assert len(x.shape) == 1

    assert (x>edges[0]).all()
    assert (x<edges[-1]).all()
    N = len(edges)-1
    L = edges[-1]-edges[0]
    print((edges[0], edges[-1], "L",L, "L/2",L/2.))
    digitized = [min(int(N*(L/2.+xi)/L),N-1) for xi in x]
    #print(digitized)
    print(x)
    for i in range(len(x)-shift):
        A[digitized[i+shift],digitized[i]] += 1
    return A

#===============================================
#  transition matrix 2D: count the transitions in trajectory,
#  both normal and radial distance
#===============================================

def count_2D(B,X,Y,Z,edges,redges,shift=1):
    print(B.shape)
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

    for start,end,r in zip(zstart,zend,dr):
        B[r-1,end,start] += 1

    return B

#===============================================
# exit times: detect mean first passage times in the trajectory
#===============================================


def get_exittime(traj,cut,cutcenter,par="exit"):

    if par == "cross-lr":
        return get_crosstime(traj,cut,cutcenter,side="lr")   # kind of bad code here... TODO
    elif par == "cross-rl":
        return get_crosstime(traj,cut,cutcenter,side="rl")
    elif par != "exit":
        raise ValueError("this parameter is not found:",par)

    assert len(traj)>0
    assert cut>cutcenter

    # False is 0, True is 1
    is_out = abs(traj)>cut
    is_incenter = abs(traj)<cutcenter

    #detect inside
    try:
        countin = is_incenter.tolist().index(True)
    except ValueError:
        countin = None
	#print "warning-1, particle never within cutoff region"

    # detect outside
    try:
        if countin is not None:
            countout = is_out.tolist().index(True,countin+1)
        else:
            countout = is_out.tolist().index(True)
    except ValueError:
        countout = None
        #print "warning-2, particle never leaves the region"

    if countin != None and countout != None:
        # store time difference
        diff = countout-countin
        #print "exit time: shift,letter",shift,letter,"count",countin,countout,diff, diff*lt , "ps"
        #print "countin ",countin
        #print "countout",countout
    else: diff = None
    return countin,countout,diff

def get_exittime_old(traj,cut,cutcenter,):
    assert len(traj)>0
    assert cut>cutcenter
    count = 0

    # False is 0, True is 1
    is_out = abs(traj)>cut
    is_incenter = abs(traj)<cutcenter

    #detect inside
    for i in range(count,len(traj)):
        if is_incenter[i]: break
        else: count += 1
    countin = count

    for i in range(count,len(traj)):
        if is_out[i]: break
        else: count += 1
    countout = count
    #print "countin ",countin
    #print "countout",countout

    if countin >= len(traj)-1:
        #particle never within cutoff region
        diff = None
        #print "warning-1"
    elif countout >= len(traj)-1:
        #particle never leaves the region
        diff = None
        #print "warning-2"
    else:
        # store time difference
        diff = countout-countin
        #print diff
        #print "exit time: shift,letter",shift,letter,"count",countin,countout,diff, diff*lt , "ps"
    return countin,countout,diff


def fpt_add(exittimes1,exittimes2,x,edges,b1,b2,):
    """Collect first passage times"""
    assert len(x.shape) == 1
    nbins = len(edges)-1
    assert b1 >= -1
    assert b2 <= nbins+1   # TODO think about
    digitized = np.digitize(x,edges)

    isout_b1 = (digitized<=b1).tolist()
    isout_b2 = (digitized>=b2).tolist()
    inside = (digitized>b1) * (digitized<b2)

    for i,start in enumerate(digitized):
        # start is the bin at time i
        if inside[i]:  # if at time i inside [b1,b2]
            #print "yes",i,start,inside[i]
            try:
                exit_b1 = isout_b1.index(True,i)  # first exit time to b1 after i
                #print "1.",exit_b1
                exittimes1[start].append(exit_b1-i)
            except ValueError:
                pass
            try:
                exit_b2 = isout_b2.index(True,i)  # first exit to b2 after i
                #print "2.",exit_b2
                exittimes2[start].append(exit_b2-i)
            except ValueError:
                pass
        #else:
            #print "no",i,start,inside[i]

def fpt_all_trajs(trajs,edges,b1,b2):
    # b1 = left bin
    # b2 = right bin
    nbins = len(edges)-1
    exittimes1 = [[] for i in range(nbins)]  # for every initial position, exit to b1
    exittimes2 = [[] for i in range(nbins)]  # for every initial position, exit to b2
    for i,traj in enumerate(trajs):
        #print "-"*5,"traj",i
        fpt_add(exittimes1,exittimes2,traj,edges,b1,b2,)
    return exittimes1, exittimes2

def get_crosstime(traj,cut1,cut2,side):
    assert len(traj)>0
    assert cut1<=cut2
    if cut1 == cut2:
        return countin,countout,diff
    count = 0

    # False is 0, True is 1
    is_above_cut1 = traj>cut1
    is_below_cut2 = traj<cut2
    is_inside = np.logical_and(is_above_cut1,is_below_cut2)
    is_out = np.logical_not(is_inside)
    is_out_left = np.logical_not(is_above_cut1)
    is_out_right = np.logical_not(is_below_cut2)

#    #detect outside
#    try:
#        countout_init = is_out.tolist().index(True)
#        countout_init_left  = is_out_left.tolist().index(True)
#        countout_init_right = is_out_right.tolist().index(True)
#
#        if is_out_left[countout_init]: region = 1
#        elif is_out_right[countout_init]: region = 3
#        else: raise ValueError("what's wrong here??")
#    except ValueError:
#        countout_init = None
#        #print "warning-1, particle never within cutoff region"

    # detect outside left
    try:
        countout_init_left  = is_out_left.tolist().index(True)
    except ValueError:
        countout_init_left  = None
    #detect crossing from left to inside
    if countout_init_left is not None:
        try:
            countin_left = is_inside.tolist().index(True,countout_init_left)
        except ValueError:
            countin_left = None
    else:
        countin_left = None
    # detect crossing to outside on other side
    if countin_left is not None:
        try:
            countout_right = is_out_right.tolist().index(True,countin_left+1)
        except ValueError:
            countout_right = None
    else:
        countout_right = None


    # detect outside right
    try:
        countout_init_right = is_out_right.tolist().index(True)
    except ValueError:
        countout_init_right  = None
    #detect crossing from right to inside
    if countout_init_right is not None:
        try:
            countin_right = is_inside.tolist().index(True,countout_init_right)
        except ValueError:
            countin_right = None
    else:
        countin_right = None
    # detect crossing to outside on other side
    if countin_right is not None:
        try:
            countout_left = is_out_left.tolist().index(True,countin_right+1)
        except ValueError:
            countout_left = None
    else:
            countout_left = None

    # either recrossing to the same region, either crossing to other side

    if countin_left != None and countout_right != None:
        # store time difference
        diff_lr = countout_right-countin_left
        #print "countin ",countin_left
        #print "countout",countout_right
    else: diff_lr = None
    if countin_right != None and countout_left != None:
        # store time difference
        diff_rl = countout_left-countin_right
        #print "countin ",countin_right
        #print "countout",countout_left
    else: diff_rl = None

    if side == "lr":
        return countin_left,countout_right,diff_lr
    elif side == "rl":
        return countin_right,countout_left, diff_rl
    else:
        raise ValueError("side not known")


#===============================================
# survival: detect survival of individual bins in the trajectory
#===============================================

def calc_survival_probability(list_coor,edges,shift=1):
    """Calculate the survival probability
    without counting particles that left a bin during the lag time"""
    initial = np.zeros(len(edges)+1,float)
    survived = np.zeros(len(edges)+1,float)
    for x in list_coor:
        calc_survival_probability_add(initial,survived,x,edges,shift=shift) 
    return survived/initial, initial, survived

def calc_survival_probability_add(initial,survived,x,edges,shift=1):
    assert len(x.shape)==1
    assert shift < len(x)
    assert len(initial)==len(survived)
    assert len(initial)==len(edges)+1

    digitized = np.digitize(x,edges)
    for i in range(len(digitized)-shift):
        start = digitized[i]
        same = (digitized[i:i+shift+1]==start).all()
        #print(start,same,int(same),digitized[i:i+shift+1])
        initial[start] += 1.    # should be float
        survived[start] += float(same)   # should be float

    return survived/initial, initial, survived

def indices_survived(x,edges,shift=1):
    """check whether a coordinate is still in the same bin after the lag time"""
    assert len(x.shape) == 1
    assert shift < len(x)

    indices = []
    digitized = np.digitize(x,edges)
    for i in range(len(digitized)-shift):
        start = digitized[i]
        same = (digitized[i:i+shift+1]==start).all()
        if same:
            indices.append(i)
    return indices

#=========== OLDER =========
# older - instead of square - not much used anymore
def write_Tmat_linebyline(A,filename,edges=None):
    L = len(A)  # number of bins + 1
    with open(filename,"w+") as f:
        if edges is not None:
            f.write("#edges  "+" ".join([str(b) for b in edges])+"\n")
        for i in range(L):
            for j in range(L):
                f.write(str(A[i,j]))

