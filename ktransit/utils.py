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

