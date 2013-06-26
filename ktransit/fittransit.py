import ktransit 
from scipy import optimize
import numpy as np 

class FitTransit():

    def __init__(self):
        self.mod = ktransit.LCModel()
    #    self.nplanets = self.mod.nplanets
        self.free_parameters(fitparstar=[],
        fitparplanet=[])

    @property 
    def nplanets(self):
        self._nplanets = self.mod.nplanets
        return self._nplanets

    def add_guess_star(self,rho=1.5, ld1=0.2, 
        ld2=0.4, ld3=0.0, ld4=0.0, dil=0.0,
        veloffset=0.0, zpt=0.0):
        self.starguess_d={
            'rho': rho,
            'ld1': ld1,
            'ld2': ld2,
            'ld3': ld3,
            'ld4': ld4,
            'dil': dil,
            'veloffset': veloffset,
            'zpt': zpt}
        self.mod.add_star(rho, ld1, 
        ld2, ld3, ld4, dil,
        veloffset, zpt)

    def update_star(self,**kwargs):
        goodlist = ['rho', 'ld1',
        'ld2', 'ld3', 'ld4', 'dil',
        'veloffset', 'zpt']
        kwnew = dict([(k,v) for k,v in kwargs.iteritems() if k in goodlist])
        self.mod.add_star(**kwnew)

    def update_planet(self,pnum,**kwargs):
        goodlist = ['T0',
            'period',
            'impact',
            'rprs',
            'ecosw', 
            'esinw',
            'rvamp',
            'occ',
            'ell', 
            'alb']
        kwnew = dict([(k,v) for k,v in kwargs.iteritems() if k in goodlist])
        self.mod.add_planet(replace=pnum,**kwnew)

    def add_guess_planet(self, 
        T0=1.0, period=1.0, impact=0.1,
        rprs=0.1, ecosw=0.0, esinw=0.0,
        rvamp=0.0, occ=0.0, ell=0.0, 
        alb=0.0):
        pname = 'pnum' + str(self.nplanets)
        planet_d = {
            'T0': T0,
            'period': period,
            'impact': impact,
            'rprs': rprs,
            'ecosw': ecosw, 
            'esinw': esinw,
            'rvamp': rvamp,
            'occ': occ,
            'ell': ell, 
            'alb': alb}
        self.planetguess_d = {
            pname: planet_d}
        self.mod.add_planet(
            T0=T0, period=period, impact=impact,
            rprs=rprs, ecosw=ecosw, esinw=ecosw,
            rvamp=rvamp, occ=occ, ell=ell, 
            alb=alb)

    def add_data(self,time=None, flux=None,
        ferr=None, itime=None, 
        ntt=None, tobs=None, omc=None,
        datatype=None):

        if time is None:
            self.time = np.arange(0,10,0.0188)
        else:
            self.time = time
        npt = len(self.time)
        nmax = 1500000

        if flux is None:
            self.flux = np.zeros(npt)
        else:
            self.flux = flux

        if ferr is None:
            self.ferr = np.zeros(npt) + 0.1
        else:
            self.ferr = ferr

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

        self.mod.add_data(time=time, itime=itime, 
            ntt=ntt, tobs=tobs, omc=omc,
            datatype=datatype)

    def free_parameters(self,
        fitparstar=['rho','zpt'],
        fitparplanet=['T0','period','impact','rprs']):
        self.fitparstar = fitparstar
        self.fitparplanet = fitparplanet

    def calc_model(self):
        self._tmod = self.mod.transitmodel
        return self._tmod

    def residuals(self):
        self._res =  (self.flux - self.calc_model()) / self.ferr
        return self._res


    def ret_lstsq(self,fitpars):
        """

        """
        fitstar_d = dict(zip(self.fitparstar,fitpars[:len(self.fitparstar)]))
        self.update_star(**fitstar_d)

        nps = len(self.fitparstar)
        npp = len(self.fitparplanet) / self.nplanets
        for i in xrange(self.nplanets):
            fitplanet_d = dict(
                zip(self.fitparplanet,
                fitpars[nps+i*npp:nps+(i+1)*npp
                ]))

            self.update_planet(pnum=i,**fitplanet_d)

        return self.residuals()

    def do_fit(self):
        fitvalstar = [self.starguess_d[k] for k in self.fitparstar]
        fitvalplanet = []
        for i in xrange(self.nplanets):
            pdic = self.planetguess_d['pnum' + str(i)]
            fitvalplanet = np.r_[fitvalplanet,
            [pdic[k] for k in self.fitparplanet]]
        fitpars = np.r_[fitvalstar,fitvalplanet]
        self.fitout = optimize.leastsq(self.ret_lstsq,fitpars,full_output=True)
        self.make_fitoutdicts()

        self.transitmodel = self.mod.transitmodel

    def make_fitoutdicts(self):
        self.fitresult = self.fitout[0]
        nps = len(self.fitparstar)
        npp = len(self.fitparplanet)
        stellarout = self.fitresult[0:nps]
        self.fitresultstellar = dict(zip(self.fitparstar,stellarout))
        self.fitresultplanets = {}
        for i in xrange(self.nplanets):
            planetout = self.fitresult[nps + i*npp:nps + (i+1)*npp]
            self.fitresultplanets['pnum' + str(i)] = dict(
                zip(self.fitparplanet,planetout))

    class Parameter:
        def __init__(self, value):
                self.value = value

        def set(self, value):
                self.value = value

        def __call__(self):
                return self.value

    def fit(function, parameters, time, flux):
        def f(params):
            i = 0
            for p in parameters:
                p.set(params[i])
                i += 1
            return y - function(x)

        p = [param() for param in parameters]
        optimize.leastsq(f, p)

def make_plot(time,obsf,model):
    #import matplotlib
    #matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    fig, (ax1,ax2) = plt.subplots(2, 1, sharex=True, sharey=False)
    ax1.scatter(time,obsf,s=2,color='k',label='obs data')
    ax1.plot(time,model,color='r',label='model fit')
    ax1.legend()
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Norm Flux')
    ax1.set_xlim([np.min(time)-1,np.max(time)+1])
    yrng = np.max(model) - np.min(model)
    ax1.set_ylim([np.min(obsf) - 0.2*yrng,
        np.max(obsf)+0.2*yrng])

    ax2.scatter(time,obsf - model,s=2,color='k',
        label='residuals')
    ax2.legend()
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Residuals')
    yrng2 = np.max(obsf - model) - np.min(obsf - model)
    ax2.set_ylim([np.min(obsf - model) - 0.2*yrng2,
        np.max(obsf - model)+0.2*yrng2])


    fig.set_figwidth(10)
    fig.set_figheight(8)

    return fig

if __name__ == '__main__':
    ##make fake data
    M = ktransit.LCModel()
    M.add_star()
    M.add_planet(period=1.001,impact=0.1,T0=1.51,rprs=0.05)
    M.add_data()

    tmod = M.transitmodel
    #addnoise
    from numpy.random import normal
    tmod_noisy = tmod + normal(0.0,0.00005,size=len(tmod))

    freeparstar = ['rho','zpt']
    freeparplanet = ['period','rprs','T0','impact']

    fitT = FitTransit()
    fitT.add_guess_star(rho=5.0)
    fitT.add_guess_planet(period=1.0,impact=0.7,T0=1.5,rprs=0.1)
    fitT.add_data(time=M.time,flux=tmod_noisy)
    fitT.free_parameters(freeparstar,freeparplanet)
    fitT.do_fit()

    #update model
    
    M.add_star(**fitT.fitresultstellar)
    for i in xrange(fitT.nplanets):
        use_d = fitT.fitresultplanets['pnum' + str(i)]
        M.add_planet(replace=i,**use_d)

    fig = make_plot(M.time,tmod_noisy+1.0,fitT.transitmodel+1.0)
    fig.savefig('fitresult.pdf')

