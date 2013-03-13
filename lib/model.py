#!/usr/bin/env python
#
# copyright: Gerhard Hummer (NIH, July 2012)
# An Ghysels (August 2012)
#

import numpy as np

"""
count  --  how transitions were counted [pbc, cut, ...]
edges  --  edges of the bins
nbins  --  number of bins
dim_v  --  dimension of v vector (discretized F profile: free energy)
dim_w  --  dimension of w vector (discretized ln D profile: diffusion)
dim_trans  --  dimension of transition matrix

list_lt  --  array with lagtimes
timezero  --  time offset t0

v[i]  --  potential in bin i [in unit kT]
w[i]  --  log of diffusion coefficient between bin i and bin i+1 [in unit dx**2/dt]

Radial
redges  --  edges of radial bins
dim_rad  --  len(redges), number of transitions in radial direction
wrad[i]  --  log of diffusion coefficient in radial direction in bin i (put in middle of bin)
"""

class Model(object):
    def __init__(self,trans,D0):
        # trans is a Transitions object
        self.edges = trans.edges
        self.count = trans.count

        if trans.count == "pbc":
            self.nbins = len(self.edges)-1
            self.dim_w = self.nbins

        elif trans.count == "cut":
            self.nbins = len(self.edges)-1
            self.dim_w = self.nbins-1

        elif trans.count == "pbccut":
            self.nbins = len(self.edges)
            self.dim_w = self.nbins

        assert trans.dim_trans == self.nbins
        self.dim_trans = self.nbins
        self.dim_v = self.nbins

        if "pbc" in trans.count: self.pbc = True  # w has equal length as v
        else: self.pbc = False   # w is 1 shorter than v

        import copy
        self.list_lt = copy.deepcopy(trans.list_lt)  # in the model we may have "effective" lagtimes

        self.init_model(D0)

    def init_model(self,D0):
        print "INIT MODEL"
        # initialize potential v[i] and w[i]=log(D(i))
        # initialize F and D (v and w) to a constant value
        self.v = np.float64(np.zeros(self.dim_v))  # in kBT
        self.w = np.float64(np.zeros(self.dim_w))  # in dx**2/dt
        # units
        dt = 1. # ps
        dx = self.edges[1]-self.edges[0]  # in angstrom
        unit = dx**2/dt            # in angstrom**2/ps
        self.vunit = 1.  # in kBT
        self.wunit = np.log(unit)  # in angstrom**2/ps

        self.w += (np.log(D0)-self.wunit)  # initialize to D0 in A**2/ps

        # timezero
        self.timezero = 0.


class CosinusModel(Model):
    """Derived class, making use of basis functions"""
    def __init__(self,trans,D0,ncosF,ncosD,ncosP):
        assert ncosF >= 0 # number of basis functions
        assert ncosD >= 0 # number of basis functions
        self.ncosF = ncosF
        self.ncosD = ncosD
        self.ncosP = ncosP
        Model.__init__(self,trans,D0)
        if self.ncosF > 0:
            self.init_model_cosF()
        if self.ncosD > 0:
            self.init_model_cosD(D0)
        if self.ncosP > 0:
            self.init_model_cosP()

    def create_basis_center(self,dim,ncos,type="cos"):
        x = np.arange(dim)
        basis = [np.cos(2*k*np.pi*(x+0.5)/dim)/(k+1) for k in range(ncos)]
        return np.array(basis).transpose()
    def create_basis_border(self,dim,dim2,ncos):
        x = np.arange(dim)
        basis = [np.cos(2*k*np.pi*(x+1.)/dim2)/(k+1) for k in range(ncos)]
        #print "basis",basis
        return np.array(basis).transpose()

    def init_model_cosF(self):
        print "INIT COS MODEL F", self.ncosF
        self.v_coeff = np.zeros((self.ncosF),float)
        self.v_basis = self.create_basis_center(self.dim_v,self.ncosF)
        self.update_v()

    def init_model_cosD(self,D0):
        print "INIT COS MODEL D", self.ncosD
        self.w_coeff = np.zeros((self.ncosD),float)
        ###
        self.w_coeff[0] += (np.log(D0)-self.wunit)   # initialize to D = D0 A**2/ps

        # important: periodicity self.v_dim
        self.w_basis = self.create_basis_border(self.dim_w,self.dim_v,self.ncosD)
        self.update_w()

    def init_model_cosP(self):
        print "INIT COS MODEL P", self.ncosP
        self.p_coeff = np.zeros((self.ncosP),float)
        self.p_coeff[0] = 1.  # normalization TODO
        self.p_basis = self.create_basis_center(self.dim_v,self.ncosP)
        self.update_p()


    def calc_profile(self,coeff,basis):
        a = coeff * basis
        return np.sum(a,axis=1)

    def update_v(self,coeff=None):
        if coeff is not None:
            assert len(coeff)==len(self.v_coeff)
            self.v_coeff = coeff
        self.v = self.calc_profile(self.v_coeff,self.v_basis)

    def update_w(self,coeff=None):
        if coeff is not None:
            assert len(coeff)==len(self.w_coeff)
            self.w_coeff = coeff
        self.w = self.calc_profile(self.w_coeff,self.w_basis)

    def update_p(self,coeff=None):
        if coeff is not None:
            assert len(coeff)==len(self.p_coeff)
            self.p_coeff = coeff
        self.p = self.calc_profile(self.p_coeff,self.p_basis)


#class HarrModel(Model):
class StepModel(Model):
    """Derived class, making use of basis functions"""
    def __init__(self,trans,D0,ncosF,ncosD,ncosP):
        assert ncosF >= 0 # number of basis functions
        assert ncosD >= 0 # number of basis functions
        self.ncosF = ncosF
        self.ncosD = ncosD
        self.ncosP = ncosP
        Model.__init__(self,trans,D0)
        if self.ncosF > 0:
            self.init_model_F()
        if self.ncosD > 0:
            self.init_model_D(D0)
        if self.ncosP > 0:
            self.init_model_P()

    def init_model_F(self,):
        print "INIT STEP MODEL F", self.ncosF
        self.v_coeff = np.zeros((self.ncosF),float)  # heights
        dx = self.dim_v /2. /self.ncosF
        self.v_x0 = np.arange(0,self.ncosF*dx,dx)  # location
        self.v_basis = self.create_basis_center(self.dim_v,self.v_x0)
        self.update_v()

    def init_model_D(self,D0):
        print "INIT STEP MODEL D", self.ncosD
        self.w_coeff = np.zeros((self.ncosD),float)  # heights
        dx = self.dim_v /2. /self.ncosD   # use length of v vector 
        self.w_x0 = np.arange(0,self.ncosD*dx,dx)
        self.w_basis = self.create_basis_border(self.dim_w,self.dim_v,self.w_x0)
        self.w_coeff += (np.log(D0)-self.wunit)   # initialize to D0 A**2/ps
        self.update_w()

    def create_basis_center(self,dim,allx0):
        x = np.arange(dim)+0.5
        basis = [np.where((x>=x0)& (x<=dim-x0),1.,0.) for x0 in allx0]
        return np.array(basis).transpose()
    def create_basis_border(self,dim,dim2,allx0):
        x = np.arange(dim)+1.
       # dx = allx0[1]-allx0[0]
        # important: periodicity of dim2
       # basis = []
       # for i,x0 in enumerate(allx0):
       #     b = np.where((x>=x0)& (x<=dim2-x0),1.,0.)
       #     bstart = int(np.floor(x0))
       #     bend = min(int(np.ceil(x0)+dx), dim/2)
       #  
       #     print bstart,bend,x0,x0+dx,dx,
       #     slope = np.arange(bstart-x0,bend-x0)/dx
       #     b[bstart:bend] = slope
       #     print slope
       #     if not self.pbc:
       #         b[-bend:dim-bstart] = slope[::-1]
       #     print b
       #     basis.append(b)
        basis = [np.where((x>=x0)& (x<=dim2-x0),1.,0.) for x0 in allx0]   # ik vermoed hier eenf fout TODO
        #basis = [np.where((x>=x0)& (x<=dim2-x0),np.minimum((x-x0)/dx,1),0.) for x0 in allx0]
        return np.array(basis).transpose()

    def calc_profile(self,coeff,basis):
        a = coeff * basis
        return np.sum(a,axis=1)

    def update_v(self,coeff=None):
        if coeff is not None:
            assert len(coeff)==len(self.v_coeff)
            self.v_coeff = coeff
        self.v = self.calc_profile(self.v_coeff,self.v_basis)

    def update_w(self,coeff=None):
        if coeff is not None:
            assert len(coeff)==len(self.w_coeff)
            self.w_coeff = coeff
        self.w = self.calc_profile(self.w_coeff,self.w_basis)


class OneStepModel(StepModel):

    def __init__(self,trans,D0,ncosF,ncosD,ncosP):
        StepModel.__init__(self,trans,D0,ncosF,ncosD,ncosP)

    def init_model_F(self,):
        print "INIT ONESTEP MODEL F", self.ncosF
        self.v_coeff = np.zeros((self.ncosF),float)  # heights
        dx = float(self.dim_v) /self.ncosF
        self.v_x0 = np.arange(0,self.ncosF*dx,dx)  # location
        self.v_basis = self.create_basis_center(self.dim_v,self.v_x0)
        self.update_v()

    def init_model_D(self,D0):
        print "INIT ONESTEP MODEL D", self.ncosD
        self.w_coeff = np.zeros((self.ncosD),float)  # heights
        dx = float(self.dim_v) /self.ncosD   # use length of v vector 
        self.w_x0 = np.arange(0,self.ncosD*dx,dx)
        self.w_basis = self.create_basis_border(self.dim_w,self.dim_v,self.w_x0)
        self.w_coeff[0] += (np.log(D0)-self.wunit)   # initialize to D0 A**2/ps
        self.update_w()

    def create_basis_center(self,dim,allx0):
        x = np.arange(dim)+0.5
        basis = [np.where((x>=x0),1.,0.) for x0 in allx0]
        return np.array(basis).transpose()
    def create_basis_border(self,dim,dim2,allx0):
        x = np.arange(dim)+1.
        basis = [np.where((x>=x0),1.,0.) for x0 in allx0]
        return np.array(basis).transpose()


class RadModel(CosinusModel):
    def __init__(self,trans,D0,ncosF,ncosD,ncosP):
        # trans is a RadTransitions object

        # setting the normal Model things
        CosinusModel.__init__(self,trans,D0,ncosF,ncosD,ncosP)

        # adding the radial model things
        self.redges = trans.redges
        self.dim_rad = trans.dim_rad  # len(redges)
        self.init_model_rad()

    def init_model_rad(self,):
        print "INIT RADMODEL"
        self.wrad = np.float64(np.zeros(self.dim_v))  # one for each bin, in dx**2/dt
        self.dim_wrad = len(self.wrad)
        # initialize wrad to a constant value of D = 1 angstrom**2/ps
        dt = 1. # ps
        dr = self.redges[1]-self.redges[0]  # in angstrom
        unit = dr**2/dt                     # in angstrom**2/ps
        self.wrad -= np.log(unit)           # exp(wrad) is in units dx**2/dt

        # units
        self.wradunit = np.log(unit)    # in angstrom**2/ps

class RadCosinusModel(RadModel):
    """Derived class, making use of basis functions"""
    def __init__(self,trans,D0,ncosF,ncosD,ncosP,ncosDrad):
        # usual things
        RadModel.__init__(self,trans,D0,ncosF,ncosD,ncosP)
        # using ncosDrad basis functions
        assert ncosDrad >= 0
        self.ncosDrad = ncosDrad
        if self.ncosDrad > 0:
            self.init_model_cosDrad()

    def init_model_cosDrad(self):
        print "INIT COS MODEL Drad", self.ncosDrad
        self.wrad_coeff = np.zeros((self.ncosDrad),float)
        self.wrad_coeff[0] = -self.wradunit  # the constant term

        assert self.dim_v == self.dim_wrad
        # important: periodicity self.dim_v = self.dim_rad!!!
        self.wrad_basis = self.create_basis_center(self.dim_v,self.ncosDrad)
        self.update_wrad()

    def update_wrad(self,coeff=None):
        if coeff is not None:
            assert len(coeff)==len(self.wrad_coeff)
            self.wrad_coeff = coeff
        self.wrad = self.calc_profile(self.wrad_coeff,self.wrad_basis)

