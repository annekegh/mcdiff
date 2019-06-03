#!/bin/python

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plotsettings():
    plt.rc(('xtick','ytick','axes'), labelsize=24.0)
    plt.rcParams['lines.linewidth'] = 3
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.minor.size'] = 4
    #plt.rcParams['xtick.major.width'] = 2  # default 0.5
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.minor.size'] = 4
    #plt.rcParams['ytick.major.width'] = 2  # default 0.5
    plt.rcParams['figure.subplot.left']=0.15
    plt.rcParams['figure.subplot.bottom']=0.14
    plt.rcParams['legend.fontsize'] = 18
def plotsettingsax(ax):
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        #line.set_color('green')
        line.set_markersize(10)
        line.set_markeredgewidth(2)
plotsettings()
linestyles = ["-","--"]

def plot_F(F,filename,edges,title="free energy",pbc=True,legend=None,grey=False,
           error=None,transparent=False):
    # F in units kBT
    # F is a list of arrays, each array is a profile
    # error is a list of arrays or just an array
    plt.figure()

    if grey and len(F) > 1:
        for i,f in enumerate(F):  #each column
            dx = edges[i][1]-edges[i][0]
            x = edges[i][:len(f)]+dx/2.
            L = edges[i][-1]-edges[i][0]    ########## TODO
            plt.plot(x,f,color="grey")
            if pbc:
                plt.plot(x+L,f,color="grey")
        plt.plot(x,f,color="black")  # overprint the last one
        if pbc:
            plt.plot(x+L,f,color="black")

    else:
        for i,f in enumerate(F):  #each column
            dx = edges[i][1]-edges[i][0]
            x = edges[i][:len(f)]+dx/2.
            L = edges[i][-1]-edges[i][0]   ################################### TODO
            # error bars
            if error is None:
                plt.plot(x,f)
            else:
                if type(error) is list:
                    assert len(error) == len(F)
                    std=error[i]
                else:
                    std=error
                assert len(std) == len(f)
                plt.errorbar(x,f,yerr=std)
            if pbc:
                plt.plot(x+L,f)
    plt.xlabel(r"$z$ ($\AA$)")
    plt.ylabel(r"$F$ ($k_BT$)")
    plt.title(title)
    plt.ylim(0,4.5)
    #plt.ylim(xmin=0)
    if legend is not None: plt.legend(legend)
    plt.grid(True)
    plt.savefig(filename,transparent=transparent)

def plot_D(D,filename,edges,title="diffusion",pbc=True,legend=None,grey=False,
           error=None,transparent=False):
    # D in units angstrom**2/ps
    # D is a list of arrays, each array is a profile
    # error is a list of arrays or an array
    plt.figure()

    if grey:
        for i,d in enumerate(D):  #each column
            dx = edges[i][1]-edges[i][0]
            x = edges[i][:len(d)]+dx
            L = edges[i][-1]-edges[i][0]   ################################### TODO
            plt.plot(x,d,color="grey")
            if pbc:
                plt.plot(x+L,d,color="grey")
        plt.plot(x,d,color="black")  # overprint the last one
        if pbc:
            plt.plot(x+L,d,color="black")
    else:
        for i,d in enumerate(D):
            dx = edges[i][1]-edges[i][0]
            x = edges[i][:len(d)]+dx
            L = edges[i][-1]-edges[i][0]   ################################### TODO
            # error bars
            if error is None:
                plt.plot(x,d)
            else:
                if type(error) is list:
                    assert len(error) == len(D)
                    std=error[i]
                else:
                    std=error
                assert len(std) == len(d)
                plt.errorbar(x,d,yerr=std)
            if pbc:
                plt.plot(x+L,d)
    #plt.xlabel("z [A]")
    #plt.ylabel("D [A^2/ps]")
    #plt.ylabel(r"$D_\mathrm{n}$ ($\AA^2/$ps)")
    plt.xlabel(r"$z$ ($\AA$)")
    plt.ylabel(r"$D_\perp$ ($\AA^2/$ps)")
    #plt.ylim(0.,0.6)
    plt.ylim(ymin=0)
    plt.title(title)
    if legend is not None:  # and error is None:    # not sure why this worked like that??
        plt.legend(legend)
    #plt.tight_layout()
    plt.grid(True)
    plt.savefig(filename,transparent=transparent)

def plot_Drad(Drad,filename,edges,title="rad-diffusion",pbc=True,legend=None,grey=False,
              error=None,transparent=False):
    # Drad in units angstrom**2/ps
    # Drad is a list of arrays, each array is a profile
    # error is a list of arrays or an array
    plt.figure()

    if grey:
        for i,d in enumerate(Drad):  #each column
            dx = edges[i][1]-edges[i][0]
            x = edges[i][:len(d)]+dx
            L = edges[i][-1]-edges[i][0]   ################################### TODO
            plt.plot(x,d,color="grey")
            if pbc:
                plt.plot(x+L,d,color="grey")
        plt.plot(x,d,color="black")  # overprint the last one
        if pbc:
            plt.plot(x+L,d,color="black")
    else:
        for i,d in enumerate(Drad):
            dx = edges[i][1]-edges[i][0]
            x = edges[i][:len(d)]+dx
            L = edges[i][-1]-edges[i][0]   ################################### TODO
            # error bars
            if error is None:
                plt.plot(x,d)
            else:
                if type(error) is list:
                    assert len(error) == len(Drad)
                    std=error[i]
                else:
                    std=error
                assert len(std) == len(d)
                plt.errorbar(x,d,yerr=std)
            if pbc:
                plt.plot(x+L,d)
    plt.xlabel(r"$z$ ($\AA$)")
    plt.ylabel(r"$D_{||}$ ($\AA^2/$ps)")
    #plt.ylim(0.,0.7)
    plt.title(title)
    if legend is not None:  # and error is None:    # not sure why this worked like that??
        plt.legend(legend)
    #plt.tight_layout()
    plt.grid(True)
    plt.savefig(filename,transparent=transparent)


def plot_both(F,D,filename,edges,transparent=False):
    # assume F,D,edges are lists of arrays, each array is a profile
    plt.figure()
    for i in range(len(F)):
        f = F[i]
        d = D[i]
        dx = edges[i][1]-edges[i][0]
        x_F = edges[i][:len(f)]+dx/2.
        x_D = edges[i][:len(d)]+dx
    
        plt.subplot(2,1,1)
        plt.plot(x_F,f)
        plt.ylabel("F [kBT]")
    
        plt.subplot(2,1,2)
        plt.plot(x_D,d)
        plt.ylabel("D [A^2/ps]")
        plt.xlabel("z [A]")

    plt.savefig(filename,transparent=transparent)

def plot_three(F,D,Drad,filename,edges,transparent=False):
    # assume F,D,edges are lists of arrays, each array is a profile
    plt.figure()
    for i in range(len(F)):
        f = F[i]
        d = D[i]
        drad = Drad[i]
        dx = edges[i][1]-edges[i][0]
        x_F    = edges[i][:len(f)]+dx/2.   # these are the middle points of the bins
        x_D    = edges[i][:len(d)]+dx      # if periodic, then one more than if non-periodic
        x_Drad = edges[i][:len(f)]+dx/2.   # these are the middle points of the bins


        plt.subplot(3,1,1)
        plt.plot(x_F,f)
        plt.ylabel(r"$F$ ($k_BT$)")

        plt.subplot(3,1,2)
        plt.plot(x_D,d)
        #plt.ylabel(r"$D_\mathrm{n}$ ($\AA^2/$ps)")
        plt.ylabel(r"$D_\perp$ ($\AA^2/$ps)")

        plt.subplot(3,1,3)
        plt.plot(x_Drad,drad)
        #plt.ylabel(r"$D_\mathrm{r}$ ($\AA^2/$ps)")
        plt.ylabel(r"$D_{||}$ ($\AA^2/$ps)")

    plt.xlabel(r"$z$ ($\AA$)")
    plt.savefig(filename,transparent=transparent)

def plot_ratio(D,Drad,filename,edges,title="anisotropy",transparent=False):
    # assume F,D,edges are lists of arrays, each list is a profile
    plt.figure()
    plt.title(title)
    for i in range(len(D)):
        d = D[i]
        drad = Drad[i]
        dx = edges[i][1]-edges[i][0]
        x_D    = edges[i][:len(d)]+dx      # if periodic, then one more than if non-periodic
        x_Drad = edges[i][:len(drad)]+dx/2.   # these are the middle points of the bins

     #   plt.subplot(2,1,1)
     #   plt.plot(x_D,d)
     #   plt.plot(x_Drad,drad,'o-')
     #   plt.ylabel("D,Drad [A^2/ps]")

     #   plt.subplot(2,1,2)
        if len(d) != len(drad):  pass   # TODO XXXXXXXXXXXXXXX  ratio = Drad[-nc:]/D[-nc:]
        nc = min(len(d),len(drad))   # number of components
        x = edges[i][:nc]

        dx = x[1]-x[0]
        x_ratio = x+3*dx/4.   # halfway
        dave = (d[:-1]+d[1:])/2.  # interpolate between points
        dave = np.array([(d[0]+d[-1])/2.]+dave.tolist())
        #ratio = drad[:nc]/d[:nc]
        #plt.plot(x_ratio,ratio)
        ratio = drad[:nc]/dave[:nc]
        plt.plot(x,ratio)
        plt.plot([x[0],x[-1]],[1.,1.],color='k',lw=1)
        #plt.ylabel("Drad/D")

        plt.xlabel(r"$z$ ($\AA$)")
        #plt.ylabel(r"$D_\mathrm{r}/D_\mathrm{n}$")
        plt.ylabel(r"$D_{||}/D_\perp$")


        #plt.subplot(3,1,3)
        #plt.plot(x,drad[:nc]-dave[:nc])
        #plt.plot([x[0],x[-1]],[0.,0.],color='k',lw=1)
        #plt.ylabel("Drad-D")
        #plt.xlabel("z [A]")

    plt.savefig(filename,transparent=transparent)


def make_plots(F,D,Drad,edges,filename,pbc=True,legend=None,grey=False,
               title=None,error=None,ave=False,transparent=False):
    # assume F in units kBT, is a list
    # assume D, Drad in units angstrom**2/ps, is a list
    # assume error is a list of [list of arrays or a list] for F and D and Drad
    outF = filename+"_F.png"
    outD = filename+"_D.png"
    outDrad = filename+"_Drad.png"
    outDratio = filename+"_Dratio.png"
    outboth = filename+"_both.png"

    if error is not None:
        Fst = error[0]
        Dst = error[1]
        Dradst = error[2]
    else:
        Fst = None
        Dst = None
        Dradst = None

    if ave:
        from .outreading import average_profiles
        F_ave_mean, D_ave_mean, Drad_ave_mean, edges_mean, F_ave_st, D_ave_st, Drad_ave_st = average_profiles(F,D,Drad,edges) 
        f,d,drad,ed  = ([F_ave_mean],[D_ave_mean],[Drad_ave_mean],[edges_mean],)  #make lists of these arrays
        fst,dst,dradst = (F_ave_st,D_ave_st,Drad_ave_st)
    else:
        f,d,drad,ed = (F,D,Drad,edges,)  # these are arrays
        fst,dst,dradst = (Fst,Dst,Dradst)      # these are lists of arrays or lists

    plot_F(f,outF,ed,pbc=pbc,grey=grey,title=title,error=fst,transparent=transparent)
    plot_D(d,outD,ed,pbc=pbc,grey=grey,title=title,error=dst,transparent=transparent,legend=legend)
    #if not None in Drad:
    if all(isinstance(d,np.ndarray) for d in drad):
        plot_Drad(drad,outDrad,edges,pbc=pbc,grey=grey,title=title,error=dradst,transparent=transparent,legend=legend)
        plot_three(f,d,drad,outboth,edges,transparent=transparent)
        plot_ratio(d,drad,outDratio,edges,transparent=transparent)
    else:
        plot_both(f,d,outboth,edges,)
