#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# An Ghysels (August 2012)
#

import numpy as np

#------------------------
# READING FUNCTIONS
#------------------------


def guess_dim_transition_square(filename):
    """get number of histogram bins from transition matrix file"""
    #print "Reading...", filename
    with open(filename,"r") as f:
        for line in f:   # skip lines starting with #
            if not line.startswith("#"):
                dim_trans = len(line.split())
                break
    return dim_trans

def guess_dim_transition_cube(filename):
    """get number of histogram bins from transition matrix file"""
    #print "Reading...", filename
    # Stored format: ...
    #   ...
    row = 0      # not really necessary
    dim_rad = 0  # dimension of transition matrix
    lenrow = 0   # dimension of transition matrix
    with open(filename,"r") as f:
        for line in f:
            if line.startswith("-"):
                assert lenrow == row
                dim_rad += 1  # I reached next radial point
                row = 0
            elif not line.startswith("#"):
                lenrow = len(line.split())
                row += 1
    return dim_rad,lenrow

def read_transition_header(filename):
    header = {}
    with open(filename,"r") as f:
        for line in f:
            if line.startswith("#"):
                words = line.split()
                if words[0] == "#edges":
                    edges = [float(word) for word in words[1:]]  # skip first
                    header['edges'] = np.array(edges)   # bin edges
                    header['dz'] = (edges[-1]-edges[0])/float(len(edges)-1)  # bin width
                elif words[0] == "#dt":
                    header['dt'] = float(words[1])
                elif words[0] == "#dn":
                    header['dn'] = int(words[1])
                elif words[0] == "#lt":
                    header['lt'] = float(words[1])
                elif words[0] == "#file":
                    filename = words[1]
                    header['file'] = filename
                elif words[0] == "#count":
                    count = words[1]
                    header['count'] = count
                elif words[0] == "#abc":
                    a = float(words[1])
                    b = float(words[2])
                    c = float(words[3])
                    header.update({'a':a,'b':b,'c':c})
                elif words[0] == "#redges":
                    redges = [float(word) for word in words[1:]]  # skip first
                    header['redges'] = np.array(redges)   # bin redges

    print("HEADER")
    print(header)
    check_content_header(header)
    return header

def check_content_header(header):
    if 'dt' in header and 'dn' in header:
        #print header['lt'], header['dn'], header['dt'], header['dn']*header['dt']
        #if 'lt' in header: assert header['lt'] == (header['dn']*header['dt'])
        if 'lt' in header: assert abs(header['lt'] - (header['dn']*header['dt']))<1e-5
        else: header['lt'] = (header['dn']*header['dt'])
    elif 'dt' in header and 'lt' in header:
        if 'dn' in header: assert header['lt'] == (header['dn']*header['dt'])
        else: header['dn'] = (header['lt']*header['dt'])
    elif 'dn' in header and 'lt' in header:
        if 'dt' in header: assert header['lt'] == (header['dn']*header['dt'])
        else: header['dt'] = (header['lt']*header['dn'])
    else:
        print("Two out of the following have to be specified: dt, dn, lt.")
        raise ValueError("can not construct lag time lt")

    print("lt", header['lt'])
    if 'count' in header:
        assert header['count'] in ["pbc","cut"]  #so far what I know
    else:
        print("counting method count is not specified correctly")
        raise ValueError("count should be specified in transition matrix file")

    # TODO if not 'edges' in header:
       # print "Bin edges should be specified in transition matrix file"
       # raise ValueError("could not find bin edges")
       # header['edges'] = None


def read_transition_square(filename,dim_trans):
    """read transition matrix line by line (whole row per line!!!)
    comment lines starting with # are skipped"""

    transition = np.zeros((dim_trans,dim_trans),int)
    with open(filename,"r") as f:
        row = 0
        for line in f:
            if not line.startswith("#"):
                words = line.split()
                assert len(words) == dim_trans
                transition[row,:] = [int(word) for word in words]
                row += 1
        if row != dim_trans:
            print("wrong number of entries in ", filename)
            quit()
    return transition

def read_transition_cube(filename,dim_rad,dim_trans):
    """Read transitions from N1 x N2 x N2 matrix"""
    # dim_rad is len(redges), is dimension of transition matrix
    # dim_trans is dimension of transition matrix
    transition = np.zeros((dim_rad,dim_trans,dim_trans))
    row = 0
    radbin = 0
    with open(filename,"r") as f:
        for line in f:
            if line.startswith("-"):
                if row != dim_trans:
                    print("wrong number of entries in ", filename)
                    quit()
                radbin += 1  # I reached next radial point
                row = 0
            elif not line.startswith("#"):
                words = line.split()
                assert len(words) == dim_trans
                transition[radbin,row,:] = [int(word) for word in words]
                row += 1
    return transition

#=========================== NOT MUCH USED/NOT UPDATED ============================

def guess_dim_transition_linebyline(filename):
    """get number of histogram bins from transition matrix file (one line per entry!!!)"""
    #print "Reading...", filename
    count = 0
    with open(filename,"r") as f:
        for line in f:   # skip lines starting with #
            if not line.startswith("#"):
                assert len(line.split())==1
                count += 1
    dim_trans = int(np.sqrt(count))
    return dim_trans


def read_transition_linebyline(filename,dim_trans):
    """read transition matrix line by line (one line per entry!!!)
    comment lines starting with # are skipped
    comment line starting with #edges contains the bin edges"""

    transition = np.zeros((dim_trans,dim_trans))
    edges = []
    with open(filename,"r") as f:
        for line in f:
            if line.startswith("#edges"):
                words = line.split()
                edges = [float(word) for word in words[1:]]  # skip first
                edges = np.array(edges)   # bin edges
                assert len(edges) == dim_trans+1
            if not line.startswith("#"):
                transition[0,0] = np.float64(line)   # put first line in [0,0] element
                k = 1
                break
        for line in f:  # read other lines
            if not line.startswith("#"):
                (k1, k2 ) = divmod ( k , dim_trans )
                transition[k1,k2] = np.float64(line)
                k += 1
        if k != dim_trans**2:
            print("wrong number of entries in ", filename)
            quit()

    if len(edges) == 0:  # use linear
        edges = np.arange(dim_trans+1)
    return transition,edges


#=========================== BACKUP ============================

def read_transition_square_bkp(filename,dim_trans):
    """read transition matrix line by line (whole row per line!!!)
    comment lines starting with # are skipped
    comment line starting with #edges contains the bin edges"""

    transition = np.zeros((dim_trans,dim_trans))
    edges = []
    Dt = 0
    row = 0
    with open(filename,"r") as f:
        for line in f:
            if line.startswith("#Dt"):
                words = line.split()
                assert len(words)==2
                Dt = float(words[1])   # this is the lagtime!!
            if line.startswith("#edges"):
                words = line.split()
                edges = [float(word) for word in words[1:]]  # skip first
                edges = np.array(edges)   # bin edges
    
             # TODO   assert len(edges) == dim_trans+1
            if not line.startswith("#"):
                words = line.split()
                transition[row,:] = [int(word) for word in words]
                row += 1
        if row != dim_trans:
            print("wrong number of entries in ", filename)
            quit()

    if len(edges) == 0:  # use linear
        edges = np.arange(dim_trans+1)
    return transition,edges,Dt

def read_all_transitions_bkp(args,form="square"):
    """read all input: lagtimes and transitions"""
    print("args:", args)

    # number of lag times
    num_lag = len(args)/2
    lagtimes = np.zeros((num_lag),dtype=np.float64)
    print("number of lagtimes:", num_lag)

    # dimension of transition matrix
    if form == "square":
        dim_trans = guess_dim_transition_square(args[1])
    else:
        dim_trans = guess_dim_transition_linebyline(args[1])

    transition = np.zeros((num_lag,dim_trans,dim_trans),np.float64)
    print("dimension trans:", dim_trans)

    # read lagtimes and transition files
    for i in range(num_lag):
        lagtimes[i] = np.float64(args[2*i])

        filename = args[2*i+1]
        if form == "square":
            trans,edges,Dt = read_transition_square(filename,dim_trans)
            if Dt != 0:
                lagtimes[i] = Dt
        else:
            trans,edges = read_transition_linebyline(filename,dim_trans)
        transition[i,:,:] = trans
        print("window ", i, ": lagtime", lagtimes[i], "or",args[2*i], "from file", args[2*i+1])

    #### let's hope I have the same edges every time!!!!
    ###  add quality checks: edges, lagtime, dim_trans

    return lagtimes,transition,edges

