"""analyze traj
AG, July 15, 2012"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import struct


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

class RunData(object):
    akma_time = 48.8882143060371e-3  # in ps
    def __init__(self):
        self.started = False
        self.startedout = False
        self.distdone = False
        
    def read_charmm_traj(self,filename,dtc=None):
        """Read rom CHARMM trajectory"""
        import pychm.future.io.charmm as iocharmm
        import pychm.scripts.getprop as getprop
        taco = iocharmm.open_dcd(filename)

        data = taco.get_massive_dump()

        #from pprint import pprint
        #pprint (vars(taco))
        #print "frame_dt",taco.frame_dt
 
        x = data['x']
        y = data['y']
        z = data['z']
        xtl = data['xtl']

        if not self.started:
            self.x = x
            self.y = y
            self.z = z
            self.xtl = xtl
            self.started = True
            self.size = x.shape[1]
            self.acdone = [False]*self.size

            self.nstep = taco.nsteps
            self.nsavv = taco.nsavv  # freq to save velocities
            self.del_t = taco.del_t  # TODO  WHATT??? Timestep in AKMA-units. Bit-copy from the 32-bit real number
            self.del_t = struct.unpack('f', struct.pack('i', taco.del_t))[0] * RunData.akma_time  # should be identical to self.dt
            # nprint = freq to print energy      500
            # iprfrq = freq to give statistics  5000
            # ntrfrq = freq to remove trans/rot  100
            self.dt = 0.001   # ps  # make sure to be float!!!                # TODO hard coded
            if dtc == None:
                self.dtc = self.dt * taco.nsavc
                self.nsavc = taco.nsavc  # freq to save coords
            else:
                self.dtc = float(dtc)
                self.nsavc = float(dtc)/self.dt

        else:

            self.x = np.append(self.x,x,0)
            self.y = np.append(self.y,y,0)
            self.z = np.append(self.z,z,0)
            self.xtl = np.append(self.xtl,xtl,0)
            if dtc == None:
                dtc = self.dt * taco.nsavc
                nsavc = taco.nsavc
            else:
                nsavc = float(dtc) / self.dt

            assert self.nsavc == nsavc
            assert self.nsavv == taco.nsavv
            #assert self.dt == taco.dt
            assert self.dtc == dtc

        self.xpbc = self.xtl[:,0]
        self.ypbc = self.xtl[:,2]
        self.zpbc = self.xtl[:,-1]

        self.nstep = self.x.shape[0]
        print "Added: %15i Total: %15i CHARMM: nsteps/nsavc = %i/%i=%i" %(x.shape[0],
                self.nstep,taco.nsteps,taco.nsavc,taco.nsteps/taco.nsavc)

    def read_charmm_out(self,filename):

        outDict = getprop.getProp(open(filename),'dynatime','dynaener',
                'dynapresse','dynapressi','dynatemp')  #'averener','dynavdw')
        if not self.startedout:
            self.ener = np.array(outDict['dynaener'])
            self.time = np.array(outDict['dynatime'])
            self.presse = np.array(outDict['dynapresse'])
            self.pressi = np.array(outDict['dynapressi'])
            self.temp   = np.array(outDict['dynatemp'])
            self.startedout = True
        else:
            self.ener = np.append(self.ener,outDict['dynaener'])
            self.time = np.append(self.time,outDict['dynatime'])
            self.presse = np.append(self.presse,outDict['dynapresse'])
            self.pressi = np.append(self.pressi,outDict['dynapressi'])
            self.temp   = np.append(self.temp,outDict['dynatemp'])

    def print_settings(self):
        print "Settings"
        for attr in ["started","size","distdone",
                    "nstep","nsavc","nsavv",
                    "dt","dtc","del_t"]:
            print "%12s " %attr, self.__dict__[attr]
            

    def update(self,minstep=0):
        self.x_ave = np.mean(self.x[minstep:],0)
        self.y_ave = np.mean(self.y[minstep:],0)
        self.z_ave = np.mean(self.z[minstep:],0)

        self.xcom = np.mean(self.x,1)
        self.ycom = np.mean(self.y,1)
        self.zcom = np.mean(self.z,1)
        self.xcom_ave = np.mean(self.xcom)
        self.ycom_ave = np.mean(self.ycom)
        self.zcom_ave = np.mean(self.zcom)

        #self.xpbc = self.xtl[:,0]
        #self.ypbc = self.xtl[:,2]
        #self.zpbc = self.xtl[:,-1]

    def remove_com_motion(self):
        #print self.x.shape
        #print self.xcom.shape
        self.x[:,:] = (self.x[:,:].transpose()-self.xcom[:]).transpose()
        self.y[:,:] = (self.y[:,:].transpose()-self.ycom[:]).transpose()
        self.z[:,:] = (self.z[:,:].transpose()-self.zcom[:]).transpose()
        
        self.xcom = np.mean(self.x,1)
        self.ycom = np.mean(self.y,1)
        self.zcom = np.mean(self.z,1)
        self.xcom_ave = np.mean(self.xcom)
        self.ycom_ave = np.mean(self.ycom)
        self.zcom_ave = np.mean(self.zcom)

    def find_middle_layer(self,outdir,moltypes=None,outfile=None):
        if moltypes == None:
            moltypes = rd.select
        for moltype in moltypes:
            print "=====\n%s" % moltype
            cor = self.z[:,self.select[moltype]]
            mu = np.mean(cor,1)
            sd = np.std(cor,1)
            # check if lower when shifted
            mask = cor < 0.
            maskshift = (mask.transpose()*self.xpbc).transpose()
            zshift = cor + maskshift
            mu1 = np.mean(zshift,1)
            sd1 = np.std(zshift,1)
            # shift other direction: this gives the same standard deviation

            print "average of last 200 timesteps"
            print "not shifted   ",
            print "mean:", np.mean(mu[-100:]),  "sd", np.mean(sd[-100:])
            print "when shifted  ",
            print "mean:", np.mean(mu1[-100:]), "sd", np.mean(sd1[-100:])

            if outfile is not None:
                f = file(outfile+".%s.z" %moltype,"w+")
                print >> f, "# average z-coor %s" %moltype
                for i in xrange(len(mu)):
                    print >> f, mu[i], sd[i]    
                f.close()

            #plot_vs_time(np.array([mu,mu-sd,mu+sd,mu1,mu1-sd1,mu1+sd1]).transpose(),self.dtc,outdir+"/fig_musd.%s.png"%moltype)


    def smooth_traj_manually(self):
        for i,cor in enumerate([self.x, self.y, self.z]):
            pbc = [self.xpbc, self.ypbc, self.zpbc][i]
            for at in range(self.size): #[0,]:  #xrange(self.size):    # HARDCODED TODO
                diff = cor[1:,at]-cor[:-1,at]
                plusmask = diff > (pbc[:-1]*0.75)
                minmask = diff < -(pbc[:-1]*0.75)
                plusmask = plusmask.tolist()
                minmask = minmask.tolist()

                start = -1
                while True:
                    try:
                        start = start + plusmask[start+1:].index(True) +1
                        cor[start+1:,at] -= pbc[start]  # shift
                        #print start, plusmask[start]
                    except ValueError:
                        break
                start = -1
                while True:
                    try:
                        start = start + minmask[start+1:].index(True) +1
                        cor[start+1:,at] += pbc[start]  # shift
                        #print start, minmask[start]
                    except ValueError:
                        break


    def preliminary_plots(self,outdir):
        histogram(self.x[:,0],outdir+"/fig_hist.x0.png")
        histogram(self.z[:,0],outdir+"/fig_hist.z0.png")
        plot_vs_time(self.x[:,range(10)],self.dtc,outdir+"/fig_x0.png")
        plot_vs_time(self.ener,self.dt,outdir+"/fig_ener.png")
        plot_vs_time(self.temp,self.dt,outdir+"/fig_temp.png")
        plot_vs_time(self.pressi,self.dt,outdir+"/fig_pressi.png")
        plot_vs_time(self.presse,self.dt,outdir+"/fig_presse.png")
        plot_vs_time(np.array([self.xcom,self.ycom,self.zcom]).transpose(),self.dtc,outdir+"/fig_com.png")

    def write_pbc(self,outfile):
        f = file(outfile,"w+")
        print >> f, "# box parameters a b c"
        for t in range(self.nstep):
            print >> f, self.xpbc[t], self.ypbc[t], self.zpbc[t]
        f.close()

    def write_coor(self,index,outfile,moltype,pbc=False):
        if index == 0:
            coor = self.x[:,self.select[moltype]]
        elif index == 1:
            coor = self.y[:,self.select[moltype]]
        elif index == 2:
            coor = self.z[:,self.select[moltype]]
        else:
            raise Error("index should be 0, 1, or 2")
        if pbc:
            if index == 0:
                pbccoor = self.xpbc
            elif index == 1:
                pbccoor = self.ypbc
            elif index == 2:
                pbccoor = self.zpbc

        f = file(outfile,"w+")
        print >> f, "# coor[%i]" %index
        for step in range(coor.shape[0]):
            vec = coor[step,:]
            if pbc:
                vec -= pbccoor[step]*np.floor(vec/pbccoor[step]+0.5)
            line = " ".join([str(vec[at]) for at in range(len(vec))])
            print >> f, line
        f.close()

    def calc_dist_normal(self):
        if self.distdone == True: return
        xy2 = (self.x[:,:] - self.x[0,:])**2 + (self.y[:,:] - self.y[0,:])**2
        z = abs(self.z[:,:] - self.z[0,:])
        r2 = xy2 + z**2
        self.dist_xy = np.sqrt(xy2)
        self.dist_z  = z
        self.dist_r  = np.sqrt(r2)
        self.distdone = True

    def calc_dist(self):   # kind of oversampling
        if self.distdone == True: return
        alldist_xy = np.zeros((self.nstep,self.size),float)
        alldist_z  = np.zeros((self.nstep,self.size),float)

        for lt in range(1,self.nstep):  # lt is lagtime
            # fill in alldist[i]
            diff_xy = (self.x[:-lt,:] - self.x[lt:,:])**2 + (self.y[:-lt,:] - self.y[lt:,:])**2
            diff_z  = abs(self.z[:-lt,:] - self.z[lt:,:])
            #xy2 = 0.
            #z = 0.
            #for j in range(self.nstep-lt):
            #    xy2 += (self.x[j+lt,:] - self.x[j,:])**2 + (self.y[j+lt,:] - self.y[j,:])**2
            #    z += self.z[j+lt,:] - self.z[j,:]
            alldist_xy[lt,:] = np.mean(diff_xy,0)  #xy2/(self.nstep-lt)
            alldist_z[lt,:]  = np.mean(diff_z,0)   #z/(self.nstep-lt)

        self.dist_xy = np.sqrt(alldist_xy)
        self.dist_z  = alldist_z
        self.dist_r  = np.sqrt(alldist_xy+alldist_z**2)
        self.distdone = True

    def autocorr_ave(self,select):
        if np.sum(self.acdone) == 0:
            self.ac = np.zeros((self.nstep,self.size,3),float)
        for i,at in enumerate(xrange(self.size)):
            if not self.acdone[at]:
                self.ac[:,i,0] = autocorr(self.x[:,at])
                self.ac[:,i,1] = autocorr(self.y[:,at])
                self.ac[:,i,2] = autocorr(self.z[:,at])
                self.acdone[at] = True
        AC = np.mean(self.ac[:,select,:],1)
        print "AC", AC.shape
        return AC


    def fit_r_vs_time(self,outdir):
        if not self.distdone: self.calc_dist()
        for moltype in self.select:
            print "=====\n%s" % moltype
            a = self.dist_r[:,self.select[moltype]]
            average = np.mean(a,1)
            #plot_vs_time(average,self.dtc,outdir+"/fig_distmean.%s.png"%moltype)
            fit_sqrt_vs_time(average,self.dtc,outdir+"/fig_distmean.%s.png"%moltype,title=moltype)

    def fit_z_vs_time(self,outdir):
        if not self.distdone: self.calc_dist()
        for moltype in self.select:
            print "=====\n%s" % moltype
            a = self.dist_z[:,self.select[moltype]]
            average = np.mean(a,1)
            #plot_vs_time(average,self.dtc,outdir+"/fig_distz.%s.png"%moltype)
            fit_sqrt_vs_time(average,self.dtc,outdir+"/fig_distz.%s.png"%moltype,title=moltype)

    def fit_xy_vs_time(self,outdir):
        if not self.distdone: self.calc_dist()
        for moltype in self.select:
            print "=====\n%s" % moltype
            a = self.dist_xy[:,self.select[moltype]]
            average = np.mean(a,1)
            #plot_vs_time(average,self.dtc,outdir+"/fig_distxy.%s.png"%moltype)
            fit_sqrt_vs_time(average,self.dtc,outdir+"/fig_distxy.%s.png"%moltype,title=moltype)

    def analyze_dist(self,outdir,outfile=None):
        self.calc_dist()
        if outfile is not None:
            f = open(outdir+"/"+outfile+".dist.obj",'wb')
            pickle.dump(self,f)
            f.close()
        self.fit_r_vs_time(outdir)
        self.fit_z_vs_time(outdir)
        self.fit_xy_vs_time(outdir)
        plot_vs_time(self.dist_r[:,:10],self.dtc,outdir+"/fig_dist0.png")


    def create_histograms(self):
        """Assume used on trajectories that have NOT been expanded"""
        for moltype in self.select:
            print "=====\n%s" % moltype
            kwargs = {"range":(-40,40), "normed":True}
            histogram(self.x[:,self.select[moltype]].ravel(),"fig_hist.x.%s.png"%moltype,**kwargs)
            histogram(self.y[:,self.select[moltype]].ravel(),"fig_hist.y.%s.png"%moltype,**kwargs)
            histogram(self.z[:,self.select[moltype]].ravel(),"fig_hist.z.%s.png"%moltype,**kwargs)
            # log scale
            kwargs = {"range":(-40,40), "normed":True, "log":True}
            histogram(self.z[:,self.select[moltype]].ravel(),"fig_hist.zlog.%s.png"%moltype,**kwargs)

    def transition_matrix(self,bins,filename=None,shift=1):
        #self.Tmat = []
        for moltype in self.select:
          if moltype == "O2":    #XXXXXXXXXXXXX TODO FIX !!!!!!!!!!!!!!!!!!!!!!
            for i,allcor in enumerate([self.x, self.y, self.z]):
                label = ["x","y","z"][i] + "." + moltype
                cor = allcor[:,self.select[moltype]].transpose()

                nb = len(bins)                  # number of bins
                A = np.zeros((nb+1,nb+1),int)   # also those that are too high/too low
                for traj in cor:                # add from several atoms
                    fill_transition_matrix(A,traj,bins,shift=shift)
                #fill_transition_matrix2(A,cor,bins,shift=shift)

                # assume external ones are identical transitions, therefore add ?
                # or put to zero?    # TODO
                #A[0,-1] = 0
                #A[-1,0] = 0

                #self.Tmat.append(A)                
                if filename is not None:
                    write_Tmat(A,filename+"."+str(shift)+"."+label+".dat")
                    write_Tmat(A[1:-1,1:-1],filename+"."+str(shift)+"."+label+".cut.dat")

                # periodic boundary conditions
                A[-1,1:] += A[0,1:]
                A[1:,-1] += A[1:,0]
                A[-1,-1] += A[0,0]
                if filename is not None:
                    write_Tmat(A[1:,1:],filename+"."+str(shift)+"."+label+".pbc.dat")


def autocorr(x,shift=True):
    if shift is True:
        vec = x - np.mean(x)
    else:
        vec = x
    result = np.correlate(vec, vec, mode='full')
    return result[result.size/2:]

def fill_transition_matrix(A,x,bins,shift=1):
    digitized = np.digitize(x,bins)
    #ind = np.argsort(x)
    #digitized = np.searchsorted(x[ind], bins, "right")
    #digitized = 
    #print x[:20]
    #print digitized[:20]
    #print bins

    start = digitized[0]
    for i,val in enumerate(digitized[shift:]):
        start = digitized[i]
        end = val
        #print i,start,end
        A[start,end] += 1
        #start = end

def fill_transition_matrix2(A,cor,bins,shift=1):
    """This method is equivalent but much too slow,
    because of the the np.unique and the classIds.index(i) commands."""
    # project trajectory on bins
    class_ids = np.zeros(cor.shape,int)
    for i,traj in enumerate(cor):
        class_ids[i,:] = np.digitize(traj,bins)
    n, t = class_ids.shape
    js = range(t-1)

    # all possible bin values
    k = len(bins)+1
    classIds = range(k)

    # fill matrix
    transitions = np.zeros((k,k))
    for state_0 in js:
            state_1 = state_0+shift
            if state_1 >= t: break
            state_0 = class_ids[:,state_0]
            state_1 = class_ids[:,state_1]
            initial = np.unique(state_0)
            for i in initial:
                ending = state_1[state_0 == i]
                uending = np.unique(ending)
                row = classIds.index(i)
                for j in uending:
                    col = classIds.index(j)
                    transitions[row, col] += sum(ending == j)
        #row_sum = transitions.sum(axis=1)
        #p = np.dot(np.diag(1/(row_sum+(row_sum == 0))), transitions)
        #self.p = np.matrix(p)
    A[:,:] = transitions[:,:]

def write_Tmat(A,filename):
    f = file(filename,"w+")
    L = len(A)  # number of bins + 1
    for i in xrange(L):
        for j in xrange(L):
            print >> f, A[i,j]
    f.close()



#=====================
# PLOTS
#=====================

def histogram(x,filename,**kwargs):

    plt.figure()
    plt.hist(x,bins=100,**kwargs)
    plt.savefig(filename)

    #hist(x, bins=10, range=None, normed=False, weights=None,
    #   cumulative=False, bottom=None, histtype='bar', align='mid',
    #   orientation='vertical', rwidth=None, log=False,
    #   color=None, label=None, **kwargs)

def plot_vs_time(x,dt,filename,**kwargs):
    plt.figure()
    t = np.arange(len(x))*dt
    plt.plot(t,x,**kwargs)
    plt.xlabel("t [ps]")
    plt.savefig(filename)

def fit_sqrt_vs_time(rdata,dt,figname,title=None,verbose=False,t0=0.,std=None,sqrt=True):
    # use input data correctly
    if sqrt:
        r = rdata
        r2 = rdata**2
    else:
        r = np.sqrt(abs(rdata))
        r2 = rdata
    if std == None:
        std2 = None
    else:
        if sqrt: std2 = std**2
        else: std2 = std

    t = np.arange(len(r))*dt+t0

    # fit distance
    p = np.polyfit(t[:50],r[:50],1)
    if verbose:
        print "linear fit (50)"
        print p

    p = np.polyfit(t,r2,1)   # this fit is plotted
    if verbose:
        print "linear fit square dist"
        print p
    p = np.poly1d(p)
    fitted_1 = p(t)
    a_1 = p.c[0]  # this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                # = 1e-8 meter**2/second = 1e-4 cm**2/s
    b_1 = p.c[1]

    A = np.power.outer(t,np.arange(0.5,1.,1.))
    outputs = np.linalg.lstsq(A,r)[0]
    if verbose:
        print "fit root time"
        print outputs
    fitted_2 = outputs[0]**2*t
    a_2 = outputs[0]

    plt.figure()
    if std == None:
        plt.plot(t,r2)
    else:
        plt.errorbar(t,r2,yerr=std2)
    plt.plot(t,fitted_1)
    plt.plot(t,fitted_2)
    plt.xlabel("t [ps]")
    plt.ylabel("MSD [A^2]")
    plt.legend(["data","a*t+b","a*t"])
    if title is None: title = ""
    title += ", a=%f [e-4cm**2/s]  b=%f" %(a_1,b_1)
    plt.title(title)
    plt.savefig(figname)
    return a_1,b_1,a_2


def analyze_dist_multipleruns(list_rd,outdir):
    assert len(list_rd) > 0
    # take some settings from first run
    dtc = list_rd[0].dtc
    select = list_rd[0].select
    nstep = list_rd[0].nstep
    nrd = len(list_rd)
    t = np.arange(0,nstep*dtc,dtc)

    allD = np.zeros((nrd,3,len(select)))
    print "="*5
    print "Results"

    for it, moltype in enumerate(select):
        print moltype
        tot_xy = np.zeros((nstep,nrd),float)
        tot_z  = np.zeros((nstep,nrd),float)
        tot_r  = np.zeros((nstep,nrd),float)
        for i,rd in enumerate(list_rd):
            print i
            a = rd.dist_xy[:,rd.select[moltype]]
            average = np.mean(a,1)
            tot_xy[:,i] = average
            p = np.polyfit(t,average**2,1)
            allD[i,0,it] = p[0]  # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                        # = 1e-8 meter**2/second = 1e-4 cm**2/s

            a = rd.dist_z[:,rd.select[moltype]]
            average = np.mean(a,1)
            tot_z[:,i] = average
            p = np.polyfit(t,average**2,1)
            allD[i,1,it] = p[0]  # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                        # = 1e-8 meter**2/second = 1e-4 cm**2/s

            a = rd.dist_r[:,rd.select[moltype]]
            average = np.mean(a,1)
            tot_r[:,i] = average
            p = np.polyfit(t,average**2,1)
            allD[i,2,it] = p[0]  # a_1 this is in angstrom**2/ps = 1e-20/1e-12 meter**2/second
                        # = 1e-8 meter**2/second = 1e-4 cm**2/s

        m_xy = np.mean(tot_xy,1)
        #s_xy = np.std(tot_xy,1)
        m_z  = np.mean(tot_z,1)
        #s_z  = np.std(tot_z,1)
        m_r  = np.mean(tot_r,1)
        #s_r  = np.std(tot_r,1)
        fit_sqrt_vs_time(m_xy,dtc,outdir+"/fig_distxy.%s.average.png"%moltype,title=moltype)
        fit_sqrt_vs_time(m_z, dtc,outdir+"/fig_distz.%s.average.png"%moltype,title=moltype)
        fit_sqrt_vs_time(m_r, dtc,outdir+"/fig_distmean.%s.average.png"%moltype,title=moltype)

    
    print "="*20
    print "Diffusion estimates"
    print allD.shape
    print allD

    print "="*5
    for it,moltype in enumerate(select):
        print moltype
        print "x"
        print allD[:,0,it]
        print "z"
        print allD[:,1,it]
        print "r"
        print allD[:,2,it]

    print "="*5
    print "Diffusion constant"
    for it,moltype in enumerate(select):
        print moltype
        print "   %20s   %20s" %("Dmean[e-4cm^2/s]","Dstd")
        print "x  %20.10f   %20.10f" %(np.mean(allD[:,0,it]), np.std(allD[:,0,it]))
        print "z  %20.10f   %20.10f" %(np.mean(allD[:,1,it]), np.std(allD[:,1,it]))
        print "r  %20.10f   %20.10f" %(np.mean(allD[:,2,it]), np.std(allD[:,2,it]))
    print "="*5

   

#=====================
# RUN
#=====================

def collect_data(indir,chunkstart,chunkend,outdir,outfile,smooth=True,ext='trj',dtc=None):
    rd = RunData()

    for i in range(chunkstart,chunkend):
        if smooth:
            filename = indir+"/dyn"+str(i)+".nojump."+ext
        else:
            filename = indir+"/dyn"+str(i)+"."+ext
        rd.read_charmm_traj(filename,dtc=dtc)
        print "read:", filename

# TODO make optional
#        filename = indir+"/dyn"+str(i)+".out"
#        rd.read_charmm_out(filename)
#        print "read:", filename

        #print rd.x.shape, len(rd.time), len(rd.ener)

    #plot_vs_time(rd.x[:,0],rd.dtc,outdir+"/fig_x0.noshift.png")
    #plot_vs_time(rd.x[:,90],rd.dtc,outdir+"/fig_x90.noshift.png")

  #TODO  #rd.update()
    #if smooth:
    #    rd.smooth_traj_manually()

    # store
    if outfile is not None:
        if smooth: filename = outdir+"/"+outfile+".smooth.obj"
        else: filename = outdir+"/"+outfile+".nosmooth.obj"
        f = open(filename,'wb')
        pickle.dump(rd,f)
        f.close()
    return rd


