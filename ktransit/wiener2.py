from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import numpy as np
from scipy import optimize, fftpack, signal


"""
This function is a modified version of the astroML wiener_filter function
https://github.com/astroML/astroML

Copyright (c) 2012-2013, Jacob Vanderplas All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer. Redistributions in binary form must
reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the
distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""




def wiener2(t, h, signal='Cauchy', noise='flat', return_PSDs=False,
                  signal_params=None, noise_params=None,lowcut=None):
    """Compute a Wiener-filtered time-series
    This is a hacked version of the astroML wiener_filter to fit a
    Cauchy distribution rather than a Gaussian

    Also does gaussian and a two parameter model with cauchy + gaussian

    Parameters
    ----------
    t : array_like
        evenly-sampled time series, length N
    h : array_like
        observations at each t
    signal : str (optional)
        currently only 'Cauchy' is supported
    noise : str (optional)
        currently only 'flat' is supported
    return_PSDs : bool (optional)
        if True, then return (PSD, P_S, P_N)
    signal_guess : tuple (optional)
        A starting guess at the parameters for the signal.  If not specified,
        a suitable guess will be estimated from the data itself. (see Notes
        below)
    noise_guess : tuple (optional)
        A starting guess at the parameters for the noise.  If not specified,
        a suitable guess will be estimated from the data itself. (see Notes
        below)

    Returns
    -------
    h_smooth : ndarray
        a smoothed version of h, length N

    Notes
    -----
    The Wiener filter operates by fitting a functional form to the PSD::

       PSD = P_S + P_N

    The resulting frequency-space filter is given by::

       Phi = P_S / (P_S + P_N)

    This entire operation is equivalent to a kernel smoothing by a
    kernel whose Fourier transform is Phi.

    Choosing Signal/Noise Parameters
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    the arguments ``signal_guess`` and ``noise_guess`` specify the initial
    guess for the characteristics of signal and noise used in the minimization.
    They are generally expected to be tuples, and the meaning varies depending
    on the form of signal and noise used.  For ``gaussian``, the params are
    (amplitude, width).  For ``flat``, the params are (amplitude,).

    See Also
    --------
    scipy.signal.wiener : a static (non-adaptive) wiener filter
    """
    # Validate signal
    if signal != 'Cauchy':
        raise ValueError("only signal='Cauchy' is supported")
    if signal_params is not None and len(signal_params) != 2:
        raise ValueError("signal_params should be length 2")

    # Validate noise
    if noise != 'flat':
        raise ValueError("only noise='flat' is supported")
    if noise_params is not None and len(noise_params) != 1:
        raise ValueError("noise_params should be length 1")

    # Validate t and hd
    t = np.asarray(t)
    h = np.asarray(h)

    if (t.ndim != 1) or (t.shape != h.shape):
        raise ValueError('t and h must be equal-length 1-dimensional arrays')

    # compute the PSD of the input
    N = len(t)
    Df = 1. / N / (t[1] - t[0])
    f = fftpack.ifftshift(Df * (np.arange(N) - N / 2))

    H = fftpack.fft(h)
    PSD = abs(H) ** 2

    # fit signal/noise params if necessary
    if signal_params is None:
        amp_guess = np.max(PSD[1:])
        width_guess = np.min(np.abs(f[PSD[1:] < np.mean(PSD[1:])]))
        signal_params = (amp_guess, width_guess)
    if noise_params is None:
        noise_params = (np.mean(PSD[1:]),)

    # Set up the Wiener filter:
    #  fit a model to the PSD: sum of signal form and noise form

    def signalL2(x, A, width):
        """ Lorentzian """
        width = abs(width) + 1E-99  # prevent divide-by-zero errors
        L1 = A * (width / ((x)**2 +width**2))
        return L1

    def noise(x, n):
        return n * np.ones(x.shape)

    # use [1:] here to remove the zero-frequency term: we don't want to
    # fit to this for data with an offset.
    min_func = lambda v: np.sum((PSD[1:] - signalL2(
        f[1:], v[0], v[1]) - noise(f[1:], v[2])) ** 2)
    v0 = tuple(signal_params) + tuple(noise_params)
    v,d1,d2 = optimize.fmin_l_bfgs_b(min_func, v0, approx_grad=True)
    P_S = signalL2(f, v[0], v[1])
    P_N = noise(f, v[2])

    #shall cutoff our filter at low frequency?
    #obviously not, that would be silly
    if lowcut != None:
            cutoff_freq = lowcut
            mask = f[1:] < cutoff_freq
            P_S[mask] = 0.0

    Phi = P_S / (P_S + P_N)
    Phi[0] = 1  # correct for DC offset

    # Use Phi to filter and smooth the values
    h_smooth = fftpack.ifft(Phi * H)

    if not np.iscomplexobj(h):
        h_smooth = h_smooth.real

    if return_PSDs:
        return h_smooth, PSD, P_S, P_N, Phi
    else:
        return h_smooth

def wienerLG(t, h, signal='LG', noise='flat', return_PSDs=False,
                  signal_params=None, noise_params=None,lowcut=None,
                  Gauss_bounds=None):
    """Compute a Wiener-filtered time-series
    This is a hacked version of the astroML wiener_filter to fit a
    Cauchy + Gaussian distribution rather than a Gaussian

    Parameters
    ----------
    t : array_like
        evenly-sampled time series, length N
    h : array_like
        observations at each t
    signal : str (optional)
        only 'LG' is supported
    noise : str (optional)
        currently only 'flat' is supported
    return_PSDs : bool (optional)
        if True, then return (PSD, P_S, P_N)
    signal_guess : tuple (optional)
        A starting guess at the parameters for the signal.  If not specified,
        a suitable guess will be estimated from the data itself. (see Notes
        below)
    noise_guess : tuple (optional)
        A starting guess at the parameters for the noise.  If not specified,
        a suitable guess will be estimated from the data itself. (see Notes
        below)

    Returns
    -------
    h_smooth : ndarray
        a smoothed version of h, length N

    Notes
    -----
    The Wiener filter operates by fitting a functional form to the PSD::

       PSD = P_S + P_N

    The resulting frequency-space filter is given by::

       Phi = P_S / (P_S + P_N)

    This entire operation is equivalent to a kernel smoothing by a
    kernel whose Fourier transform is Phi.

    Choosing Signal/Noise Parameters
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    the arguments ``signal_guess`` and ``noise_guess`` specify the initial
    guess for the characteristics of signal and noise used in the minimization.
    They are generally expected to be tuples, and the meaning varies depending
    on the form of signal and noise used.  For ``gaussian``, the params are
    (amplitude, width).  For ``flat``, the params are (amplitude,).

    See Also
    --------
    scipy.signal.wiener : a static (non-adaptive) wiener filter
    """
    # Validate signal
    if signal != 'LG':
        raise ValueError("only signal='LG' is supported")
    if signal_params is not None and len(signal_params) != 5:
        raise ValueError("signal_params should be length 2")

    # Validate noise
    if noise != 'flat':
        raise ValueError("only noise='flat' is supported")
    if noise_params is not None and len(noise_params) != 1:
        raise ValueError("noise_params should be length 1")

    # Validate t and hd
    t = np.asarray(t)
    h = np.asarray(h)

    if (t.ndim != 1) or (t.shape != h.shape):
        raise ValueError('t and h must be equal-length 1-dimensional arrays')

    # compute the PSD of the input
    N = len(t)
    Df = 1. / N / (t[1] - t[0])
    f = fftpack.ifftshift(Df * (np.arange(N) - N / 2))

    H = fftpack.fft(h)
    PSD = abs(H) ** 2

    # fit signal/noise params if necessary
    if signal_params is None:
        amp_guess = np.max(PSD[1:])
        width_guess = np.min(np.abs(f[PSD[1:] < np.mean(PSD[1:])]))
        amp_guess2 = 0.1
        width_guess2 = 0.5
        mu_guess2 = 10.
        signal_params = (amp_guess, width_guess,
            amp_guess2,width_guess2,mu_guess2)
    if noise_params is None:
        noise_params = (np.mean(PSD[1:]),)

    # Set up the Wiener filter:
    #  fit a model to the PSD: sum of signal form and noise form

    def signalLG(x, A, width, A2, width2, mu2):
        """ Lorentzian + Gaussian """
        width = abs(width) + 1E-99  # prevent divide-by-zero errors
        L1 = A * (width / ((x)**2 +width**2))
        L2 = A2 * np.exp(-0.5 * ((x-mu2) / width2) ** 2)
        L3 = A2 * np.exp(-0.5 * ((x+mu2) / width2) ** 2)
        return L1 + L2 + L3
    def noise(x, n):
        return n * np.ones(x.shape)

    # use [1:] here to remove the zero-frequency term: we don't want to
    # fit to this for data with an offset.
    min_func = lambda v: np.sum((PSD[1:] - signalLG(
        f[1:], v[0], v[1], v[2], v[3], v[4]) - noise(f[1:], v[5])) ** 2)
    v0 = tuple(signal_params) + tuple(noise_params)
    if Gauss_bounds is not None:
        bounds = bounds = [(None,None), (None,None), (None,None),
                (None,None), Gauss_bounds, (None,None)]
    else:
        bounds = None
    v,d1,d2 = optimize.fmin_l_bfgs_b(min_func, v0, approx_grad=True,
        bounds=bounds)
    P_S = signalLG(f, v[0], v[1], v[2], v[3], v[4])
    P_N = noise(f, v[5])

    #shall cutoff our filter at low frequency?
    #obviously not, that would be silly
    if lowcut != None:
            cutoff_freq = lowcut
            mask = f[1:] < cutoff_freq
            P_S[mask] = 0.0

    Phi = P_S / (P_S + P_N)
    Phi[0] = 1  # correct for DC offset

    # Use Phi to filter and smooth the values
    h_smooth = fftpack.ifft(Phi * H)

    if not np.iscomplexobj(h):
        h_smooth = h_smooth.real

    if return_PSDs:
        return h_smooth, PSD, P_S, P_N, Phi
    else:
        return h_smooth
