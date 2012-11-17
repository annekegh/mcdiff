#!/bin/python

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plotsettings():
    plt.rc(('xtick','ytick','axes'), labelsize=24.0)
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['figure.subplot.left']=0.15
    plt.rcParams['figure.subplot.bottom']=0.14
    plt.rcParams['legend.fontsize'] = 18
def plotsettingsax(ax):
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        #line.set_color('green')
        line.set_markersize(10)
        line.set_markeredgewidth(2)
plotsettings()


def plot_F(F,filename,edges,title="free energy",pbc=True,grey=False,
           error=None):
    # assume F in units kBT
    dx = edges[1]-edges[0]
    x = edges[:len(F)]+dx/2.
    L = edges[-1]-edges[0]
    plt.figure()

    if grey and len(F.shape) > 1:
        for f in F.transpose():  #each column
            plt.plot(x,f,color="grey")
            if pbc:
                plt.plot(x+L,f,color="grey")
        plt.plot(x,f,color="black")  # overprint the last one
        if pbc:
            plt.plot(x+L,f,color="black")

    else:
        if error is None:
            plt.plot(x,F)  # plot in middle of bin
        else:
            plt.errorbar(x,F,yerr=error)
        if pbc:
           plt.plot(x+L, F)
    plt.xlabel("z [A]")
    plt.ylabel("F [kBT]")
    plt.title(title)
    plt.ylim(0,5)
    plt.savefig(filename)

def plot_D(D,filename,edges,title="diffusion",pbc=True,legend=None,grey=False,error=None):
    # assume D in units angstrom**2/ps
    dx = edges[1]-edges[0]
    x = edges[:len(D)]+dx    # if periodic, then one more than if non-periodic
    L = edges[-1]-edges[0]
    plt.figure()

    if grey and len(D.shape) > 1:
        for d in D.transpose():  #each column
            plt.plot(x,d,color="grey")
            if pbc:
                plt.plot(x+L,d,color="grey")
        plt.plot(x,d,color="black")  # overprint the last one
        if pbc:
            plt.plot(x+L,d,color="black")
    else:
        if error is None:
            plt.plot(x,D)
        else:
            plt.errorbar(x,D,yerr=error)
        if pbc:
            plt.plot(x+L, D)
    plt.xlabel("z [A]")
    plt.ylabel("D [A^2/ps]")
    plt.title(title)
    if legend is not None and error is None:
        plt.legend(legend)
    plt.savefig(filename)

def plot_Drad(D,filename,edges,title="rad-diffusion",pbc=True,legend=None,grey=False,error=None):
    # assume D in units angstrom**2/ps
    dx = edges[1]-edges[0]
    x = edges[:len(D)]+dx/2.    # if periodic, then one more than if non-periodic
    L = edges[-1]-edges[0]
    plt.figure()

    if grey and len(D.shape) > 1:
        for d in D.transpose():  #each column
            plt.plot(x,d,color="grey")
            if pbc:
                plt.plot(x+L,d,color="grey")
        plt.plot(x,d,color="black")  # overprint the last one
        if pbc:
            plt.plot(x+L,d,color="black")
    else:
        if error is None:
            plt.plot(x,D)
        else:
            plt.errorbar(x,D,yerr=error)
        if pbc:
            plt.plot(x+L, D)
    plt.xlabel("z [A]")
    plt.ylabel("Drad [A^2/ps]")
    plt.title(title)
    if legend is not None and error is None:
        plt.legend(legend)
    plt.savefig(filename)


def plot_both(F,D,filename,edges):
    dx = edges[1]-edges[0]
    x_F = edges[:len(F)] + dx/2.  # these are the middle points of the bins
    x_D = edges[:len(D)] + dx  # if periodic, then one more than if non-periodic

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(x_F,F)
    plt.ylabel("F [kBT]")

    plt.subplot(2,1,2)
    plt.plot(x_D,D)
    plt.ylabel("D [A^2/ps]")
    plt.xlabel("z [A]")

    plt.savefig(filename)

def plot_three(F,D,Drad,filename,edges):
    dx = edges[1]-edges[0]
    x_F = edges[:len(F)] + dx/2.  # these are the middle points of the bins
    x_D = edges[:len(D)] + dx  # if periodic, then one more than if non-periodic
    x_Drad = edges[:len(Drad)] + dx/2.  # these are the middle points of the bins

    plt.figure()
    plt.subplot(3,1,1)
    plt.plot(x_F,F)
    plt.ylabel("F [kBT]")

    plt.subplot(3,1,2)
    plt.plot(x_D,D)
    plt.ylabel("D [A^2/ps]")

    plt.subplot(3,1,3)
    plt.plot(x_Drad,Drad)
    plt.ylabel("Drad [A^2/ps]")
    plt.xlabel("z [A]")

    plt.savefig(filename)

def plot_ratio(D,Drad,filename,edges):
    dx = edges[1]-edges[0]
    #x_D = edges[:len(D)] + dx  # if periodic, then one more than if non-periodic
    #x_Drad = edges[:len(Drad)] + dx/2.  # these are the middle points of the bins
    nc = min(len(D),len(Drad))  # number of components

    x = edges[:nc]  #halfway

    plt.figure()
    ratio = Drad[:nc]/D[:nc]
    plt.plot(x+3*dx/4.,ratio)  # halfway
    if len(D) != len(Drad):  # if not periodic...
        ratio = Drad[-nc:]/D[-nc:]
        # TODO
    plt.savefig(filename)

def make_plots(F,D,Drad,edges,filename,pbc=True,legend=None,grey=False,title=None,error=None,):
    outF = filename+"_F.png"
    outD = filename+"_D.png"
    outDrad = filename+"_Drad.png"
    outDratio = filename+"_Dratio.png"
    outboth = filename+"_both.png"

    if error is not None:
        Fst = error[0]
        Dst = error[1]
    else:
        Fst = None
        Dst = None

    plot_F(F,outF,edges,pbc=pbc,grey=grey,title=title,error=Fst)
    plot_D(D,outD,edges,pbc=pbc,legend=legend,grey=grey,title=title,error=Dst)
    if Drad is not None:
        plot_Drad(Drad,outDrad,edges,pbc=pbc,legend=legend,grey=grey,title=title,)
        plot_three(F,D,Drad,outboth,edges,)
        plot_ratio(D,Drad,outDratio,edges)
    else:
        plot_both(F,D,outboth,edges,)

