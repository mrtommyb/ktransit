import numpy as np 
import kplr

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
        flux[quarter == i] /= np.median(
            flux[quarter == i])
        ferr[quarter == i] /= np.median(
            flux[quarter == i])
    return flux, ferr

def byebyebaddata(flux,ferr,quality):
    finite = np.isfinite(flux)
    qflag = quality == 0
    mask = np.logical_and(finite,qflag)
    flux[~mask] = np.nan
    ferr[~mask] = np.nan
    return flux, ferr

def i_hate_nans(time,flux,ferr):
    finite = np.isfinite(flux)
    return time[finite],flux[finite],ferr[finite]

def get_clean_data(kic):
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
    flux, ferr = norm_by_quarter(flux, ferr,quarter)
    cflux = (flux / med_filt(time,flux,dt=1.0)) - 1.0

    time,cflux,ferr = i_hate_nans(time,cflux,ferr)

    return (time,cflux,ferr), star

