"""Script to plot histogram of xyz coords
AG, August 13, 2012"""

import numpy as np
import matplotlib.pyplot as plt


def read_coor_x_or_y_or_z(filename):
    #print "Reading...", filename
    data = []
    with open(filename,"r") as f:
        for line in f:
            if not line.startswith("#"):
                data.append([float(word) for word in line.split()])
    return np.array(data)

def read_data_rv(filename):
    """Extract data from files of Rick Venable"""
    data = []
    with open(filename,"r") as f:
        for line in f:
            words = line.split()
            data.append([float(word) for word in words])

    data = np.array(data)
    #print "data",data.shape
    #print data[:10,:10]
    return data

def read_coor(filename,rv=False,axis=2,com=False):
    # axis: 0 (x), 1 (y), or 2 (z)
    if rv:
        data = read_data_rv(filename)
        #times = data[:,0]
        if com:
            coms = data[:,1:4]
            return coms[:,axis::3]
        else:
            coords = data[:,4:]
            return coords[:,axis::3]
    else:
        data = read_coor_x_or_y_or_z(filename)
        return data   # this is an array

def shift_wrt_layer(data1,data2):
    print("shifting")
    d2 = np.array(data2)[:,0]   # format: mu sd

    print("d2",d2.shape)
    print("data1",data1.shape,data2.shape)
    coor = np.zeros(data1.shape,float)
    for i in range(data1.shape[1]):
        coor[:,i] = data1[:,i] - d2  # shift
    return coor

def store_hist(bin_edges,hist,filename,):
    assert len(bin_edges) == len(hist)+1
    bin_midst = [(bin_edges[i]+bin_edges[i+1])/2. for i in range(len(bin_edges)-1)]
    with open(filename,"w+") as f:
        for i in range(len(hist)):
            print(bin_edges[i], bin_edges[i+1], bin_midst[i], hist[i], file=f)
    print("file written...",filename)

def plot_histogram_pbc(coor,zpbc,figname,):
    nbins = 100
    low = -zpbc/2
    high = zpbc/2
    bins = np.linspace(low,high,nbins)
    hist, bin_edges = np.histogram(coor,bins,normed=True)
    bin_midst = [(bin_edges[i]+bin_edges[i+1])/2. for i in range(len(bin_edges)-1)]

    plt.figure()
    plt.plot(bin_midst,hist)
    plt.bar(np.array(bin_midst)+zpbc,hist)
    plt.savefig(figname+".png")
    print("file written...",figname+".png")

    plt.figure()
    maxloghist = max(np.log(hist))
    plt.plot(bin_midst,-np.log(hist)+maxloghist)
    plt.ylim(0,3.5)   # hard coded ########### TODO
    plt.ylabel("F in kBT")
    plt.savefig(figname+".one.log.png")
    print("file written...",figname+".one.log.png")
    
    # another period
    plt.plot(np.array(bin_midst)+zpbc,-np.log(hist)+maxloghist,color='blue')
    plt.savefig(figname+".log.png")
    print("file written...",figname+".log.png")

    # print to a file for later use
    with open(figname+".log.txt","w+") as f:
        print("#freeenergy: binmidst F", file=f)
        print("#zpbc", zpbc, file=f)
        for i in range(len(bin_midst)):
            print(bin_midst[i], -np.log(hist[i])+maxloghist, file=f) 
    print("file written...",figname+".log.txt")


def plot_histogram(coor,figname,ymin=None,ymax=None):
    nbins = 100
    hist, bin_edges = np.histogram(coor, nbins, normed=True)
    bin_midst = [(bin_edges[i]+bin_edges[i+1])/2. for i in range(len(bin_edges)-1)]
 
    plt.figure()
    plt.plot(bin_midst,hist)
    plt.savefig(figname+".png")
    print("file written...",figname+".png")

    plt.figure()
    maxloghist = max(np.log(hist))
    plt.plot(bin_midst,-np.log(hist)+maxloghist)
    plt.ylim(ymin,ymax)
    plt.ylabel("F in kBT")
    plt.savefig(figname+".log.png")
    print("file written...",figname+".log.png")

def add_constant_tolist(list_vecs,list_c):
    list_tot = [vec for vec in list_vecs]  # a copy
    L = len(list_vecs[0])
    for c in list_c:
        cvec = c*np.ones(L,float)
        list_tot.append(cvec)
    return list_tot

def add_constant_toarray(array_vecs,list_c):
    L = len(array_vecs)
    array_tot = np.zeros(array_vecs.shape,float)
    array_tot[:] = array_vecs[:]
    for c in list_c:
        cvec = c*np.ones((L,1),float)
        array_tot = np.append(array_tot,cvec,axis=1)        
    return array_tot

def add_vec_toarray(array_vecs,list_vecs):
    L = len(array_vecs)
    array_tot = np.zeros(array_vecs.shape,float)
    array_tot[:] = array_vecs[:]
    for vec in list_vecs:
        array_tot = np.append(array_tot,vec.reshape(-1,1),axis=1)
    return array_tot

def plotvecs(list_vecs):
    # for instance: 3 plots, each 1000 time steps [1000data,1000data,1000data]
    # then this returns something of size 1000x3, which will give 3 lines when plotted
    for i in range(len(list_vecs)):
        assert len(list_vecs[i]) == len(list_vecs[0])
    return np.array(list_vecs).transpose()


def plot_crd_z(data1,data2,dtc,outdir,zpbc):
    print("data", data1.shape, data2.shape)
    from .distance import plot_vs_time

    toplot = add_vec_toarray(data1[:,::2],[data2[:,0],data2[:,0]+zpbc/2.,data2[:,0]-zpbc/2.])
    plot_vs_time(toplot,dtc,outdir+"/fig_crd-z.png")

    coor = shift_wrt_layer(data1,data2)
    toplot = add_constant_toarray(coor[:,::2], [zpbc/2.,-zpbc/2.])
    plot_vs_time(toplot,dtc,outdir+"/fig_crd-z.shift.png")

    coor -= zpbc*np.floor(coor/zpbc+0.5)
    toplot = add_constant_toarray( coor[:,::2], [zpbc/2.,-zpbc/2.])
    kwargs = {"linewidth":0,"marker":"o","markersize":1,"markeredgewidth":0}
    plot_vs_time(toplot,dtc,outdir+"/fig_pbccrd-z.shift.png",**kwargs)

def plot_crd_z_hexd(data1,dtc,outdir,zpbc):
    from .distance import plot_vs_time

    print(data1.shape)
    coor = data1[:,::2]
    toplot = coor
    plot_vs_time(toplot,dtc,outdir+"/fig_crd-z.png")

    coor -= zpbc*np.floor(data1[:,::2]/zpbc+0.5)
    toplot = add_constant_toarray( coor[:,::2], [zpbc/2.,-zpbc/2.])
    kwargs = {"linewidth":0,"marker":"o","markersize":1,"markeredgewidth":0}
    plot_vs_time(toplot,dtc,outdir+"/fig_pbccrd-z.png",**kwargs)

