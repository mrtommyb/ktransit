from ktransit import FitTransit, LCModel, plot_results
import kplr
import numpy as np

def med_filt(x, y, dt=4.):
    """
    De-trend a light curve using a windowed median.

    """
    x, y = np.atleast_1d(x), np.atleast_1d(y)
    assert len(x) == len(y)
    r = np.empty(len(y))
    for i, t in enumerate(x):
        inds = (x >= t - 0.5 * dt) * (x <= t + 0.5 * dt)
        r[i] = np.median(y[inds])
    return r

def norm_by_quarter(flux,ferr,quarter):
    for i in np.unique(quarter):
        ferr[quarter == i] /= np.median(
            flux[quarter == i])
        flux[quarter == i] /= np.median(
            flux[quarter == i])
    return flux, ferr


def byebyebaddata(flux,ferr,quality):
    finite = np.isfinite(flux)
    qflag = quality == 0
    mask = np.logical_and(finite,qflag)
    flux[~mask] = np.nan
    ferr[~mask] = np.nan
    return flux, ferr

def i_hate_nans(time,flux,ferr, quarter, quality):
    finite = np.isfinite(flux)
    return (time[finite],flux[finite],
            ferr[finite], quarter[finite], quality[finite])

def main(kic):
    client = kplr.API()
    star = client.star(kic)
    lcs = star.get_light_curves(short_cadence=False)
    time, flux, ferr, quality,quarter = [], [], [], [], []
    for lc in lcs:
            with lc.open() as f:
                # The lightcurve data are in the first FITS HDU.
                hdu_data = f[1].data
                time = np.r_[time,hdu_data["time"]]
                #fluxval = hdu_data["sap_flux"]
                #fluxmed = np.median(fluxval)
                #flux = np.r_[flux,fluxval / fluxmed]
                #ferr = np.r_[ferr,hdu_data["sap_flux_err"] / fluxmed]
                flux = np.r_[flux,hdu_data["sap_flux"]]
                ferr = np.r_[ferr,hdu_data["sap_flux_err"]]
                quality = np.r_[quality,hdu_data["sap_quality"]]
                quarter = np.r_[quarter,f[0].header["QUARTER"] +
                    np.zeros(len(hdu_data["time"]))]


    flux, ferr = byebyebaddata(flux,ferr,quality)
    time, flux, ferr, quarter, quality = i_hate_nans(time,
                                                     flux, ferr,
                                                     quarter, quality)
    flux, ferr = norm_by_quarter(flux, ferr,quarter)
    cflux = (flux / med_filt(time,flux,dt=1.0)) - 1.0



    fitT = FitTransit()
    fitT.add_guess_star(rho=2.45)
    k0 = star.kois[0]
    fitT.add_guess_planet(
        period=k0.koi_period,
        impact=k0.koi_impact,
        T0=k0.koi_time0bk,
        rprs=k0.koi_ror)

    k1 = star.kois[1]
    fitT.add_guess_planet(
        period=k1.koi_period,
        impact=k1.koi_impact,
        T0=k1.koi_time0bk,
        rprs=k1.koi_ror)

    k2 = star.kois[2]
    fitT.add_guess_planet(
        period=k2.koi_period,
        impact=k2.koi_impact,
        T0=k2.koi_time0bk,
        rprs=k2.koi_ror)

    fitT.add_data(time=time,flux=cflux,ferr=ferr)
    freeparstar = ['rho','zpt']
    freeparplanet = [
    'period','T0','impact','rprs']
    fitT.free_parameters(freeparstar,freeparplanet)

    fitT.do_fit()

    return (time,cflux,ferr),fitT




if __name__ == '__main__':
    kic = 8478994
    (time,cflux,ferr),fitT = main(kic)
    fitT.print_results()
    fig = ktransit.plot_results(time,cflux,fitT.transitmodel)
    fig.savefig('ktransitfit.png')
