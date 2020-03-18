import numpy as np

from scipy.stats import binned_statistic

def model_rise(a, b, c, d, t):
    """Daveport 2014 flare rise phase model."""
    return 1. + a * t + b * t**2 + c * t**3 + d * t**4

def model_decay(a, b, c, d, t):
    """Daveport 2014 flare decay phase model."""
    return a * np.exp(b * t) + c * np.exp(d * t)


def aflare(t, tpeak, dur, ampl, upsample=True, uptime=10):
    '''
    The Analytic Flare Model evaluated for a single-peak (classical).
    Reference Davenport et al. (2014) http://arxiv.org/abs/1411.3723
    Use this function for fitting classical flares with most curve_fit
    tools.
    Note: this model assumes the flux before the flare is zero centered
    Parameters
    ----------
    t : 1-d array
        The time array to evaluate the flare over
    tpeak : float
        The time of the flare peak
    dur : float
        The duration of the flare
    ampl : float
        The amplitude of the flare
    upsample : bool
        If True up-sample the model flare to ensure more precise energies.
    uptime : float
        How many times to up-sample the data (Default is 10)
    Returns
    -------
    flare : 1-d array
        The flux of the flare model evaluated at each time
    '''
    _fr = [1.00000, 1.94053, -0.175084, -2.24588, -1.12498]
    _fd = [0.689008, -1.60053, 0.302963, -0.278318]

    fwhm = dur # crude approximation for a triangle shape would be dur/2.
    
    if upsample:
        dt = np.nanmedian(np.diff(t))
        timeup = np.linspace(min(t)-dt, max(t)+dt, t.size * uptime)

        flareup = np.piecewise(timeup, [(timeup<= tpeak) * (timeup-tpeak)/fwhm > -1.,
                                        (timeup > tpeak)],
                                    [lambda x: (_fr[0]+                       # 0th order
                                                _fr[1]*((x-tpeak)/fwhm)+      # 1st order
                                                _fr[2]*((x-tpeak)/fwhm)**2.+  # 2nd order
                                                _fr[3]*((x-tpeak)/fwhm)**3.+  # 3rd order
                                                _fr[4]*((x-tpeak)/fwhm)**4. ),# 4th order
                                     lambda x: (_fd[0]*np.exp( ((x-tpeak)/fwhm)*_fd[1] ) +
                                                _fd[2]*np.exp( ((x-tpeak)/fwhm)*_fd[3] ))]
                                    ) * np.abs(ampl) # amplitude

        # and now downsample back to the original time...
        ## this way might be better, but makes assumption of uniform time bins
        # flare = np.nanmean(flareup.reshape(-1, uptime), axis=1)

        ## This way does linear interp. back to any input time grid
        # flare = np.interp(t, timeup, flareup)

        ## this was uses "binned statistic"
        downbins = np.concatenate((t-dt/2.,[max(t)+dt/2.]))
        flare,_,_ = binned_statistic(timeup, flareup, statistic='mean',
                                     bins=downbins)

    else:
        flare = np.piecewise(t, [(t<= tpeak) * (t-tpeak)/fwhm > -1.,
                                 (t > tpeak)],
                                [lambda x: (_fr[0]+                       # 0th order
                                            _fr[1]*((x-tpeak)/fwhm)+      # 1st order
                                            _fr[2]*((x-tpeak)/fwhm)**2.+  # 2nd order
                                            _fr[3]*((x-tpeak)/fwhm)**3.+  # 3rd order
                                            _fr[4]*((x-tpeak)/fwhm)**4. ),# 4th order
                                 lambda x: (_fd[0]*np.exp( ((x-tpeak)/fwhm)*_fd[1] ) +
                                            _fd[2]*np.exp( ((x-tpeak)/fwhm)*_fd[3] ))]
                                ) * np.abs(ampl) # amplitude
      

    return flare

# ------------------- priors -------------------------------------

def logit(function):
    '''Make a probability distribution
    a log probability distribution.'''
    def wrapper(*args, **kwargs):
        result = function(*args, **kwargs)
        # ignore division by zero because you want to have the -np.inf results
        np.seterr(divide='ignore') 
        result = np.log(result)
        return result
    return wrapper


@logit
def gaussian_prior(x, mu, sigma):
    '''Evaluate a normalized Gaussian function
    with mu and sigma at latitude x.
    
    Parameters:
    ------------
    x : float
        latitude between -pi/2 and pi/2
    '''
    if np.isnan(x):
        #print("O")
        return 0
    else:
       # print(1 / (sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2)))

        return  1 / (sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))


# ------------- Fit FWHM with uniform prior

def lnprior_fwhm(theta):
    dur = theta
    if 0.0 < dur:
        return 0.0
    else:
        return -np.inf

def loglikelihood_fwhm(theta, time, y, yerr):
    dur = theta
    m = aflare(time, 0, dur, 1.)
    l1 = (y - m)**2 / yerr**2 + np.log(2 * np.pi * yerr**2)
    lsum = -.5 * np.nansum(l1) 
    return lsum

# ------------- Fit a, FWHM, and tpeak with uniform priors

def lnprior_full(theta):
    tpeak, dur, ampl = theta
    if -5. < tpeak < 100 and 0.0 < dur and 0.0 < ampl:
        return 0.0
    return -np.inf

def loglikelihood_full(theta, time, y, yerr):
    tpeak, dur, ampl = theta
    m = aflare(time, tpeak, dur, ampl)
    l1 = (y - m)**2 / yerr**2 + np.log(2 * np.pi * yerr**2)
    lsum = -.5 * np.nansum(l1) 
    return lsum




# ------------- Fit the rise phase only -----------------

def loglikelihood_rise(theta, time, y, yerr):
    a, b, c, d = theta
    m = model_rise(a, b, c, d, time)
    l1 = (y - m)**2 / yerr**2 + np.log(2 * np.pi * yerr**2)
    lsum = -.5 * np.nansum(l1) 
    return lsum

    
def lnprior_rise(theta):
    """
    
    Parameters:
    ------------
    theta : tuple
        start time, duration, amplitude
    x : array
        time array to constrain start time
    """
    a, b, c, d =  theta
    prior = (gaussian_prior(a, 1.941, 0.008) +
             gaussian_prior(b, -0.175, 0.032) +
             gaussian_prior(c, -2.246, 0.039) +
             gaussian_prior(d, -1.125, 0.016) 
            )
    return calculate_posterior_value_that_can_be_passed_to_mcmc(prior)


def lnprior_rise_uniform(theta):
    
    if np.isfinite(theta).all(0):
        return 0.0
    else:
        return -np.inf


# ---------------- Fit decay phase only


def lnprior_decay(theta):
    """
    Parameters:
    ------------
    theta : tuple
        start time, duration, amplitude
    x : array
        time array to constrain start time
    """
    truths = [0.6890, -1.600, 0.3030, -0.2783] #a,b,c,d
    e_truths = [0.0008, 0.003, 0.0009, 0.0007]
    a, b, c, d =  theta
    prior = (gaussian_prior(a, truths[0], e_truths[0]) +
             gaussian_prior(b, truths[1], e_truths[1]) +
             gaussian_prior(c, truths[2], e_truths[2]) +
             gaussian_prior(d, truths[3], e_truths[3]) 
            )
    return calculate_posterior_value_that_can_be_passed_to_mcmc(prior)

def loglikelihood_decay(theta, time, y, yerr):
    a, b, c, d = theta
    m = model_decay(a, b, c, d, time)
    l1 = (y - m)**2 / yerr**2 + np.log(2 * np.pi * yerr**2)
    lsum = -.5 * np.nansum(l1) 
    return lsum


# -----------------------------------------------------------


def calculate_posterior_value_that_can_be_passed_to_mcmc(lp):
    '''Do some checks to make sure MCMC will work.'''
    if not np.isfinite(lp):
        return -np.inf
    if np.isnan(lp):
        return -np.inf
    else:
        return lp

    
# -----------------------------------------------------------

def choose_lnprob(lnprior, loglikelihood):
    
    def lnprob(theta, x, y, yerr):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + loglikelihood(theta, x, y, yerr)
        
    return lnprob



