from __future__ import absolute_import

import numpy as np
from ._tmodtom import transitmodel

class LCModel(object):

    def __init__(self):
        self.nplanets = 0
        self.T0 = []
        self.period = []
        self.impact = []
        self.rprs = []
        self.ecosw = []
        self.esinw = []
        self.rvamp = []
        self.occ = []
        self.ell = []
        self.alb = []
        self.rvtime = []

    def add_star(self, rho=1.5, ld1=0.2,
        ld2=0.4, ld3=0.0, ld4=0.0, dil=0.0,
        veloffset=0.0, zpt=0.0):
        """
        add details of the star about which the planet(s) orbit
        ldX are the limb darkening parameters
        if only ld1 and and ld2 are non-zero a quadrativ limb darkening
        model is used
        if ld3 and ld4 are also non-zoro then a 4-parameter limb
        darkening law is used

        rho = mean stellar density in g/cc**3
        dil is the proportion of the total light not coming from the
            target star 0.5 means thatyou ahve two stars of equal
            brightness
        """
        self.rho = rho
        self.ldp = [ld1,ld2,ld3,ld4]
        self.dil = dil
        self.veloffset = veloffset
        self.zpt = zpt

    def update_star(self,**kwargs):
        [setattr(self,k,v) for k,v in kwargs.items()]

    def update_planet(self,pnum,**kwargs):
        for k,v in kwargs.items():
            valarr = getattr(self,k)
            valarr[pnum] = v
            setattr(self,k,valarr)


    def add_planet(self,replace=None,
        T0=1.0, period=1.0, impact=0.1,
        rprs=0.1, ecosw=0.0, esinw=0.0,
        rvamp=0.0, occ=0.0, ell=0.0,
        alb=0.0):
        if replace == None:
            self.nplanets = self.nplanets + 1
            pnum = self.nplanets - 1
            self.add_dimention_to_planet_params()
        else:
            pnum = replace

        self.T0[pnum] = T0
        self.period[pnum] = period
        self.impact[pnum] = impact
        self.rprs[pnum] = rprs
        self.ecosw[pnum] = ecosw
        self.esinw[pnum] = esinw
        self.rvamp[pnum] = rvamp
        self.occ[pnum] = occ
        self.ell[pnum] = ell
        self.alb[pnum] = alb

    def add_dimention_to_planet_params(self):
        self.T0 = np.r_[self.T0,0.0]
        self.period = np.r_[self.period,0.0]
        self.impact = np.r_[self.impact,0.0]
        self.rprs = np.r_[self.rprs,0.0]
        self.ecosw = np.r_[self.ecosw,0.0]
        self.esinw = np.r_[self.esinw,0.0]
        self.rvamp = np.r_[self.rvamp,0.0]
        self.occ = np.r_[self.occ,0.0]
        self.ell = np.r_[self.ell,0.0]
        self.alb = np.r_[self.alb,0.0]

    def add_data(self, time=None, itime=None,
        ntt=None, tobs=None, omc=None,
        datatype=None):
        """
        Add data after all the planets are added!!
        """
        if time is None:
            self.time = np.arange(0,10,0.0188)
        else:
            self.time = time
        npt = len(self.time)
        nmax = 1500000

        if itime is None:
            default_cadence = 1625.3 / 86400.
            self.itime = np.zeros(npt) + default_cadence
        else:
            self.itime = itime

        if ntt is None:
            self.ntt = np.zeros(self.nplanets)
        else:
            self.ntt = ntt

        if tobs is None:
            self.tobs = np.empty([self.nplanets,nmax])
        else:
            self.tobs = tobs

        if omc is None:
            self.omc = np.empty([self.nplanets,nmax])
        else:
            self.omc = omc

        if datatype is None:
            self.datatype = np.zeros(npt)
        else:
            self.datatype = datatype

    def add_rv(self,rvtime=None, rvitime=None):
        if rvtime is None:
            self.rvtime = np.arange(
                self.T0[0],4.*self.period[0],1.0)
        else:
            self.rvtime = rvtime

        if rvitime is None:
            default_cadence = 30. / 1440. #30 mins
            self.rvitime = np.zeros_like(self.rvtime) + default_cadence
        else:
            self.rvitime = rvitime


    @property
    def transitmodel(self):
        """
        return a transit model
        calling of model is

        transitmodel(nplanets,sol,time,itime,
        ntt,tobs,omc,datatype)

        sol is [rho,ld1,ld2,ld3,ld4,
        dil,veloffset,zpt,T0,per,b,rprs,ecosw,esinw,
        rvamp,occ,ell,alb]
        """
        ld1,ld2,ld3,ld4 = self.ldp
        sol = np.zeros(8 + self.nplanets*10)

        sol[0:8] = [self.rho,ld1,ld2,ld3,ld4,self.dil,
            self.veloffset,self.zpt]

        sol[8:] = np.array([self.T0,self.period,
            self.impact,self.rprs,
            self.ecosw,self.esinw,
            self.rvamp,self.occ,
            self.ell,self.alb]).T.flatten()

        if np.shape(self.rvtime)[0] != 0:
            time = np.r_[self.time,self.rvtime]
            itime = np.r_[self.itime,self.rvitime]
            datatype = np.r_[
                self.datatype,np.ones_like(self.rvtime)]

            fmod = transitmodel(self.nplanets,
                sol,time,itime,self.ntt,
                self.tobs,self.omc,datatype)
            nrv = np.shape(self.rvtime)[0]
            self._transitmodel = fmod[:-nrv]
            self._rvmodel = fmod[-nrv:]

        else:
            time = self.time
            itime = self.itime
            datatype = self.datatype

            self._transitmodel = transitmodel(self.nplanets,
                sol,time,itime,self.ntt,
                self.tobs,self.omc,datatype)
            self._rvmodel = None



        return self._transitmodel - 1.0

    @property
    def rvmodel(self):
        try:
            return self._rvmodel
        except AttributeError:
            tmod = self.transitmodel
            return self._rvmodel



    def get_ancil_vals(self):
        npt = len(self.time)
        itime = np.zeros(self._npt) + (self.cadence)
        ntt = np.zeros(self._nplanets)
        tobs = np.empty([self._nplanets,npt])
        omc = np.empty([self._nplanets,npt])
        datatype = np.zeros(npt)
        return [itime,ntt,tobs,omc,datatype]


def give_me_earth():
    M = LCModel()
    M.add_star()
    M.add_planet(rprs=0.009155,period=365.25)
    M.add_data(time=np.arange(0,1000,0.0188))
    return (M.time, M.transitmodel)















