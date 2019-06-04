"""Script to compare different methods to extract profiles:
calc flux
AG, April 9, 2013
AG, April 25, 2013
AG, Jan 11, 2016
AG, Sept 2016"""

#lt = 10

import numpy as np
import matplotlib.pyplot as plt
import mcdiff
from mcdiff.outreading import read_F_D_edges, read_Drad
from mcdiff.utils import construct_rate_matrix_from_F_D
from mcdiff.permeability.deff import calc_Dave_midF, calc_Dave_notmidF

##### UNITS #####
# F -- in kBT
# D -- in angstrom**2/ps
# edges -- in angstrom

# st, end -> using bins [st,st+1,...,end-1]
# this means that number of used bins = end-1-st+1 = end-st

# TODO
# rate*lagtime
# lagtime in ps
# rate in 1/dt
# so do *dt or /dt somewhere 

#################################
##### PERMEABILITIES        #####
#################################

### 1) NORMAL PERMEABILITY

def calc_permeability(F,D,dx,st,end,edges=None,ref=0,doprint=False):
    # F in units kBT, D in units A^2/ps, dx in units A
    # choose reference:
    # ref = 0 is default #### ASSUME BULK value at bin 0
    # ref = st would be another plausible assumption

    """ from bin st until bin end, using bins [st,st+1,...,end-1]
    number of used bins = end-1 - st + 1 = end-st (typical Python way of slicing)
    edges[st] is between bin [st-1] and bin [st]
    edges[end] is between bin [end-1] and bin [end]
    h = edges[end]-edges[st] = length of covered bins = height
    e.g. all bins: st=0, end=len(F)
    """
    if edges is None:
        edges = np.arange(len(F)+1)*dx   # is accurate if dx is accurate...
    h = edges[end]-edges[st]  # angstrom, is accurate
    #print "st,end",st,end

    aveD = [(D[0]+D[-1])/2.] + ((D[1:]+D[:-1])/2.).tolist()  # assume PBC
    aveD = np.array(aveD)    # angstrom**2/ps

    Fref = F[ref]
    part = np.exp(-(F-Fref))      # no units, partition function, not normalized
    dRdx = 1./(part*aveD)         # permeation resistance per unit length, unit ps/A^2
    R = np.sum(dRdx[st:end]) *dx  # resistance to permeability, unit ps/A
                                  # integrate from x=-25.5=edges[st] to x=25.5=edges[end]
    P = 1./R                      # permeability, unit A/ps
    Deff = h*P                    # effective D, in A**2/ps

    # effective length
    #heff = D[ref]/P  # in units A
    heff = aveD[ref]/P  # in units A   TODO
    diff_h = heff-h  # in units A     # this looks weird on graph??
    ratio_h = heff/h
    if doprint:
        ##print "st,end %3i %3i"%(st,end),
        print("st,end,h %7.2f %7.2f %7.2f"%(edges[st],edges[end],h), end=' ')
        print("P",P, "Deff",Deff, "heff",heff,"R",R, "dRdx",dRdx[st])
    return P,Deff,heff,diff_h

### 2) RADIAL PERMEABILITY

def calc_permeability_radial_h(F,Drad,dx,radius,st,end,edges=None,doprint=True):
    # F in units kBT, D in units A^2/ps, dx in units A
    """ from bin st until bin end
    see function calc_permeability
    """
    if edges is None:
        edges = np.arange(len(F)+1)*dx
    h = edges[end]-edges[st]  # angstrom

    # naieve:   (I might not have the correct end bin)
    #P,Deff,R = calc_permeability_radial(F[st:end],Drad[st:end],dx,radius)

    #Deff = calc_Dave_midF(F,Drad,st=st,end=end)   # average Drad, unit A^2/ps  # TODO
    Deff = calc_Dave_notmidF(F,Drad,st=st,end=end)   # average Drad, unit A^2/ps
    P = Deff/radius                               # permeability, unit A/ps
    R = radius/Deff                               # permeability resistance, unit ps/A

    if doprint:
        print("st,end,h %7.2f %7.2f %7.2f"%(edges[st],edges[end],h), end=' ')
        print("P",P, "Deff",Deff, "R",R)

    return P,Deff,R

def calc_permeability_radial(F,Drad,dx,radius):
    #Deff = calc_Dave_midF(F,Drad)   # average Drad, unit A^2/ps    # TODO
    Deff = calc_Dave_notmidF(F,Drad)   # average Drad, unit A^2/ps
    P = Deff/radius                 # permeability, unit A/ps
    R = radius/Deff                 # permeability resistance, unit ps/A
    return P,Deff,R

########
### 3) OXYGEN TRANSPORT PARAMETER

def calc_oxygentransportparameter(F,D,edges):
    """compute W ~ D*exp(-F)"""
    # F in units kBT, D in units A^2/ps
    aveD = [(D[0]+D[-1])/2.] + ((D[1:]+D[:-1])/2.).tolist()
    aveD = np.array(aveD)
    part = np.exp(-(F-min(F)))
    edges_mid = (edges[:-1]+edges[1:])/2.
    product = aveD*part     # in unit A^2/ps
    return product,edges_mid


########
### 4) EQUILIBRIUM DISTRIBUTION PROFILE WHEN MEASURING PERMEABILITY

def calc_flux(rate,prob,dt):
    """compute flux through bin borders

    fluxp -- fluxp[i] = flux from bin i to bin i+1
    fluxn -- fluxn[i] = flux from bin i+1 to bin i (absolute value)
    flux  -- flux[i] = nett flux = fluxp - fluxn

    len(flux) vector is len(prob)-1
    rate matrix can be larger, but first bin of prob and first
    bin of rate matrix should match"""
    fluxp = np.zeros(len(prob)-1,float)  # positive, left to right
    fluxn = np.zeros(len(prob)-1,float)  # negative, right to left
    for i in range(len(fluxp)):
        # flux through border between bin i and bin i+1
        # two parts: (1) from bin i to i+1   minus  (2) from bin i+1 to i
        fluxp[i] = rate[i+1,i]*prob[i]
        fluxn[i] = rate[i,i+1]*prob[i+1]
    fluxp /= dt   # dt in ps, so flux in 1/ps
    fluxn /= dt   # dt in ps, so flux in 1/ps
    flux = fluxp-fluxn
    return flux,fluxp,fluxn  # in 1/ps

def calc_permeability_distribution(F,D,dx,dt,st,end,A,B,figname=None,ref=0,doprint=True):
    """Compute the permeability from the rate matrix and the steady-state solution

    ref -- reference bin, standard bin 0 (ASSUME BULK water)
    doprint -- whether to print partitioning

    prob[st]  = A.  (no units)
    prob[end] = B.  (no units)

    so difference is Delta = A-B (no units) 
    """
    # F in units kBT, D in units A^2/ps, dx in units A

    # check
    assert st>=0
    assert end<len(F)
    assert ref>-1 and ref<len(F)
    #assert A!=0. or B!=0.

    rate = construct_rate_matrix_from_F_D(F,D,dx,dt)   # in 1/dt, PBC

    # compute steady-state solution prob
    #-----------------------------------
    # prob (no units) is not normalized
    # solve subrate*p = vec
    subrate = rate[st+1:end,st+1:end]    # covers bins [st+1,...,end-1]
    r1 = np.linalg.inv(subrate)

    # rhs of equation 
    vec1 = np.zeros(len(subrate),float)
    vec2 = np.zeros(len(subrate),float)
    vec1[0]  = -rate[st+1,st]    #vec[0]  = -rate[1,0]
    vec2[-1] = -rate[end-1,end]  #vec[-1] = -rate[N,N+1]

    # p = subrate^-1 * vec
    p1 = np.dot(r1,vec1)
    p2 = np.dot(r1,vec2)

    # plug in vector of original size, with boundaries
    # end points of vector prob are kept constant at A and B
    prob1 = np.zeros(len(subrate)+2,float)
    prob2 = np.zeros(len(subrate)+2,float)
    prob1[1:-1] = p1
    prob2[1:-1] = p2
    prob1[0]  = 1.
    prob2[-1] = 1.
    prob = A*prob1 + B*prob2

    # compute flux (here constant)
    #-----------------------------
    flux1,fluxp1,fluxn1 = calc_flux(rate[st:,st:],prob1,dt)  # in 1/ps
    flux2,fluxp2,fluxn2 = calc_flux(rate[st:,st:],prob2,dt)  # in 1/ps
    flux, fluxp, fluxn  = calc_flux(rate[st:,st:],prob, dt)  # in 1/ps   # with A,B
    # check flux is constant in steady-state with fixed boundaries
    #print np.all(flux==flux[0])  # not True, because of very small differences
    assert np.all(np.isclose(flux1, np.ones(len(flux1))*flux1[0], rtol=1e-04,) )
    assert np.all(np.isclose(flux2, np.ones(len(flux2))*flux2[0], rtol=1e-04,) )
    assert np.all(np.isclose(flux,  np.ones(len(flux)) *flux[0],  rtol=1e-04,) )
    #print "flux1 (in 1/ps)",flux1[0]
    #print "flux2 (in 1/ps)",flux2[0]
    #print "flux  (in 1/ps)",flux[0]

    # compute permeability
    #---------------------
    # particles/time = flux J = permeability * Delta-concentration-per-length
    # solve for permeability:   P = flux*dx
    P =  flux1[0]*dx     # in A/ps, assume flux is constant
    # reference is naturally in bin st, put in bin ref
    P *= np.exp(F[ref]-F[st])

 #   # this is the same as the following:
 #   P2 = -flux2[0]*dx     # in A/ps, assume flux is constant
 #   # reference is naturally in bin end, put in bin end
 #   P2 *= np.exp(F[ref]-F[end])
 #   # verify
 #   assert np.isclose(P,P2,rtol=1e-06,)

    h = dx*(len(prob)-1)  # membrane thickness, boundaries bin 0 and bin -1, in A
    R = 1./P              # resistance to permeability, in ps/A
    Deff = h*P            # effective D, in A**2/ps

    if doprint:
        #print "permeability (in A/ps)",P
        print("h %7.2f  P %10.5f  Deff %10.5f  R %10.5f"%(h,P,Deff,R))

    # Now consider steady-state, but actually imitate equilibrium
    #------------------------------------------------------------
    # by choosing A and B wisely

    # recalculate prob_equil for plotting
    # prob1: ref is st, prob2: ref is end, prob: not really ref
    # prob_equil: reference in bin ref, i.e. prob_equil[ref]=1.  # this is a choice
    prob_equil = prob1*np.exp(-F[st]+F[ref]) + prob2*np.exp(-F[end]+F[ref])

    if doprint:
        analyze_partitioning_steadystate(F,st,end,prob1,prob2,prob,ref,)

    if figname is not None:
        plt.figure()
        #plt.subplot(2,1,1)
        plt.plot(prob1,label='left')
        plt.plot(prob2,label='right')
        rescale = max(max(prob1),max(prob2))/max(prob_equil)
        #rescale=1.
        plt.plot(prob_equil*rescale,'-',color='grey',label='equil')
        plt.plot(prob,'--',color='k',lw='2',label="A=%.2f B=%.2f"%(A,B))
        plt.xlabel("bin number")
        plt.ylabel("prob (no unit)")
        plt.legend(loc='best')
        plt.title("steady-state A*left+B*right")

        plt.savefig(figname)
        print("file written...",figname)

        plt.figure()
        #plt.subplot(2,1,2)
        plt.plot(fluxp,label='pos')
        plt.plot(fluxn,'--',label='neg')
        plt.xlabel("bin number")
        plt.ylabel("flux (1/ps)")
        plt.title("flux[0]=%f.4"%(flux[0]))
        plt.legend(loc='best')
        plt.savefig(figname+"flux.png")
        print("file written...",figname+"flux.png")

    return P,Deff,R,h,prob


########
### 5) PARTITIONING
###    in steady-state, equilibrium, transient regime

def calc_mean_prob(prob):
    # membrane = everything except first and last bin
    prob_mem = np.sum(prob[1:-1])/(len(prob)-2)    # = mean condentration in membrane
    prob_wat = (prob[0]+prob[-1])/2.               # = mean concentration in water (borders)
    return prob_mem,prob_wat

def analyze_partitioning_steadystate(F,st,end,prob1,prob2,prob,ref,):
    """Compare the membrane/water partitioning beteen steady-state and equilibrium"""

    # partitioning in steady-state
    #-----------------------------
    # steady-state, but actually imitates equilibrium
    # prob1: ref is st, prob2: ref is end, prob: not really ref
    # probE: put probE[ref]=1.  # this is a choice
    probE = prob1*np.exp(-F[st]+F[ref]) + prob2*np.exp(-F[end]+F[ref])

    prob1_mem,prob1_wat = calc_mean_prob(prob1)
    prob2_mem,prob2_wat = calc_mean_prob(prob2)
    prob_mem, prob_wat  = calc_mean_prob(prob)   # with A,B
    probE_mem,probE_wat = calc_mean_prob(probE)  # steady-state imitates equilibrium

    # partitioning in equilibrium
    #-----------------------------
    part = np.exp(-(F-F[ref]))   # put reference in bin ref

    part_mem,part_wat = calc_mean_prob(part[st:end+1])

    print("--- analyze partitioning steady state ---")
    print("st  end    prob[st]   prob[end]    prob_mem   prob_wat prob_mem/prob_wat")
    print("%i  %i       %7.2f     %7.2f     %7.4f    %7.4f    %7.4f        %s"%(st,end,prob[0],
          prob[-1],prob_mem,prob_wat,prob_mem/prob_wat, "prob-A-B"))
    print("%i  %i       %7.2f     %7.2f     %7.4f    %7.4f    %7.4f        %s"%(st,end,prob1[0],
          prob1[-1],prob1_mem,prob1_wat,prob1_mem/prob1_wat, "prob-right"))
    print("%i  %i       %7.2f     %7.2f     %7.4f    %7.4f    %7.4f        %s"%(st,end,prob2[0],
          prob2[-1],prob2_mem,prob1_wat,prob2_mem/prob1_wat, "prob-left"))
    print("%i  %i       %7.2f     %7.2f     %7.4f    %7.4f    %7.4f        %s"%(st,end,probE[0],
          probE[-1],probE_mem,probE_wat,probE_mem/probE_wat, "prob-steady-state-equil"))
    print("%i  %i       %7.2f     %7.2f     %7.4f    %7.4f    %7.4f        %s"%(st,end,part[st],
          part[end],part_mem,part_wat,part_mem/part_wat, "prob-equilibr"))
    print("-"*3)


# TODO type comments for this function
def calc_partition_coefficient(F,st,end,doprint=True):
    assert st>0 or end < len(F)-1
    F0 = F[0]
    F1 = np.mean(F[:st])     # left
    F2 = np.mean(F[st:end])  # center
    F3 = np.mean(F[end:])    # right
    DeltaF_10 = F1-F0   #
    DeltaF_20 = F2-F0   #
    K = np.exp(-DeltaF_20)

    v = F-min(F)
    c0 = np.exp(-v[0])
    c1 = np.mean(np.exp(-v[:st]))
    c2 = np.mean(np.exp(-v[st:end]))
    c3 = np.mean(np.exp(-v[end:]))
    K_10 = c1/c0
    K_20 = c2/c0

    if doprint:
        print("st,end", end=' ')
        print("F0-left-cent", end=' ')
        print("DF_left-0", end=' ')
        print("DF_cent-0", end=' ')
        print("K=exp(-DF)", end=' ')
        print("K=left/0", end=' ')
        print("K=cent/0")

    if True:
        print(st,end, end=' ')
        print(F0,F1,F2, end=' ')
        print(DeltaF_10, end=' ')
        print(DeltaF_20, end=' ')
        print(K, end=' ')
        print(K_10, end=' ')
        print(K_20)

    if False:
        print("st,end",st,end, end=' ')
        print("F0-left-cent",F0,F1,F2, end=' ')
        print("DF_left-0",DeltaF_10, end=' ')
        print("DF_cent-0",DeltaF_20, end=' ')
        print("K=exp(-DF)",K, end=' ')
        print("K=left/0",K_10, end=' ')
        print("K=cent/0",K_20)

    return K_20

def analyze_partitioning_transient(F,D,pbc,dt,st,end,times,p_0,figname):
    from mcdiff.utils import init_rate_matrix
    import scipy
    import scipy.linalg
    import numpy.linalg

    # store profiles
    plt.figure()
    ax = plt.subplot(211)
    plt.plot(F,"o")
    plt.plot([st,st],[min(F)-0.1,max(F)+0.1],color='r',lw='2')
    plt.plot([end,end],[min(F)-0.1,max(F)+0.1],color='r')
    plt.ylabel("F (kBT)")
    ax = plt.subplot(212)
    plt.plot(np.arange(len(D))+0.5,D,"o")
    plt.plot([st,st],[min(D)-0.1,max(D)+0.1],color='r')
    plt.plot([end,end],[min(D)-0.1,max(D)+0.1],color='r')
    plt.ylabel("D (A/ps)")
    plt.xlabel("bin number")
    plt.savefig(figname+".FD.png")
    plt.close()

    data = np.zeros((len(times),7))
    plt.figure(1)
    plt.figure(2)
    for i,lagtime in enumerate(times):
        rate = init_rate_matrix(len(F),F,np.log(D),pbc,)
        prop = scipy.linalg.expm(lagtime*rate)
        prob = np.dot(prop,p_0)
        plt.figure(1)
        plt.plot(prob)
        p_mem,p_wat = calc_mean_prob(prob[st:end+1])   # p_wat is only border
        p_wat2 = (np.sum(prob[:st+1]) + np.sum(prob[end:]))/(len(F)-end+st+1)  # all water
        flux,fluxp,fluxn = calc_flux(rate,prob,dt)
        plt.figure(2)
        plt.plot(flux)
        data[i,0] = p_mem
        data[i,1] = p_wat
        data[i,2] = prob[st]
        data[i,3] = prob[end]
        data[i,4] = flux[st]    # st to st+1
        data[i,5] = flux[end-1] # end-1 to end
        data[i,6] = p_wat2

    plt.figure(1)
    plt.xlabel("bin number")
    plt.ylabel("concentration")
    plt.savefig(figname+".profiles.png")
    plt.figure(2)
    plt.xlabel("bin number")
    plt.ylabel("flux (1/ps)")
    plt.savefig(figname+".fluxes.png")

    # equilibrium
    equil = np.exp(-(F-min(F)))
    # rescale to match
    equil *= np.sum(p_0)/np.sum(equil)
    equil_mem,equil_wat = calc_mean_prob(equil[st:end+1])
    equil_wat2 = (np.sum(equil[:st+1]) + np.sum(equil[end:]))/(len(F)-end+st+1)
    equil_st  = equil[st]
    equil_end = equil[end]
    equil_memwat  = equil_mem / equil_wat
    equil_memwat2 = equil_mem / equil_wat2

    # for plots
    p_mem = data[:,0]
    p_wat = data[:,1]
    p_st  = data[:,2]
    p_end = data[:,3]
    flux_st  = data[:,4]
    flux_end = data[:,5]
    p_wat2= data[:,6]

    plt.figure(3)
    ax = plt.subplot(111)
    plt.plot(times,p_mem,'s-',label='m')
    plt.plot(times,p_wat2,'s-',label='w')
    plt.plot(times,p_st, 'o',label='b1')
    plt.plot(times,p_end,'o',label='b2')
    plt.plot(times,p_wat,'o-',label='b12')
    plt.plot(times,np.ones(len(times))*equil_wat2,'-',color='grey')
    plt.plot(times,np.ones(len(times))*equil_wat,'-',color='grey')
    plt.plot(times,np.ones(len(times))*equil_mem,'-',color='grey')
    plt.legend(loc='best')
    plt.xlabel("lag times (ps)")
    plt.ylabel("concentration")
    plt.title("membrane, water, borders")
    ax.set_xscale("log", nonposx='clip')
    #ax.set_yscale("log", nonposy='clip')
    plt.savefig(figname+".conc.png")

    plt.figure(4)
    ax = plt.subplot(111)
    plt.plot(times,p_mem/p_wat2,'s-',label='m/w')
    plt.plot(times,p_mem/p_wat,'o-',label='m/b12')
    plt.plot(times,np.ones(len(times))*equil_memwat2,'-',color='grey') #label='equil2')
    plt.plot(times,np.ones(len(times))*equil_memwat,'-',color='grey') #label='equil')
    plt.legend(loc='best')
    plt.xlabel("lag times (ps)")
    plt.ylabel("ratio conc")
    ax.set_xscale("log", nonposx='clip')
    #ax.set_yscale("log", nonposy='clip')
    plt.title("ratio membrane, water, borders")
    plt.savefig(figname+".concratio.png")

    plt.figure(5)
    ax = plt.subplot(111)
    plt.plot(times,flux_st, label='fl-st (1/ps)')
    plt.plot(times,flux_end,label='fl-end (1/ps)')
    #plt.plot(times,p_end-p_st, label='Deltap')
    plt.legend(loc='best')
    plt.xlabel("lag times (ps)")
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    plt.savefig(figname+".fluxDeltap.png")
    plt.close('all')

################################
# LINK WITH COUNTING CROSSINGS
################################

def convert_P_cref_to_rate(P,cref):
    return P*cref*2

def convert_P_cref_to_crossings(P,cref,T):
    return P*cref*2*T

def convert_F_D_to_rate(F,D,dx,st,end,npart,ref=None):
    if ref is None:
        ref = st
    # TODO st and end, height, in cref
    # npart is number of particles in the system
    P, Deff, heff, diff_h = calc_permeability(F,D,dx,st,end,ref=ref,doprint=False)   # P in angstrom/ps
    v = F-min(F)
    cref = np.exp(-v[ref])/np.sum(np.exp(-v))   # fraction in bin ref (no unit)
    zpbc = len(F)*dx         # periodic box size
    cref = cref/zpbc*npart   # number of particles per z-length (in # per angstrom)

    print("P (in A/ps)",P)
    print("P (in cm/s)",P*10000)  # TODO wild guess double check

    return P*cref*2



# Andreas: dcma crossings
def expected_number_of_crossings_per_area(P,number_of_permeant_molecules,
                                     cross_section_area, simulation_time,
                                     skip=None, threshold=0.05):
        """
        Given a certain number of permeant molecules in the system, the
        area of the membrane's cross-section, and the length of the simulation,
        how many molecules will cross the membrane on average?
        P = r / (2c_w), where r is the crossing rate per area and unit time.

        Args:
            number_of_permeant_molecules
            cross_section_area: average area of the bilayer mid-plane in A^2
            simulation_time: length of the simulation in ns
            skip: same as in calculate_permeability
            threshold: same as in calculate_permeability

        Returns:
            A tuple (int, int, int)
                - expected number of crossing
                - minimum of a 90% confidence interval
                - maximum of a 90% confidence interval
        """
        p, resist, min_i, max_i = self.calculate_permeability(skip, threshold)
        distribution = np.exp(-self.free_energy)
        distribution *= (number_of_permeant_molecules/sum(distribution))
        # distribution in (average # molecules / min_i bins)
        # min_i bins ^= sum(bin_widths)
        concentration_in_water = (
                sum(distribution[:min_i]) /
                (sum(self.bin_widths[:min_i]) * cross_section_area)
        )  # in #molecules/A^3
        simulation_time *= 1e3 # in ps
        p *= 1e-4  # now in Angstrom/ps = 10^-8 cm / 10^-12 s = 10^4 cm/s
        r = p * 2.0 * concentration_in_water  # in 1 / (A^2 ps)
        n_crossings = r * simulation_time * cross_section_area
        min_n, max_n = poisson.interval(0.9, n_crossings)
        return n_crossings, min_n, max_n
