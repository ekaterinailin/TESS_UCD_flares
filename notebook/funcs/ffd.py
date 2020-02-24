import numpy as np
from scipy.optimize import fmin

def ML_powerlaw_estimator(EDs, a):
    '''
    Power law ML estimator from
    Maschberger and Kroupa (2009),
    formula (9).

    Parameters:
    -----------
    data : Series or np.array
        data that is suspected to follow
        a power law relation
    '''
    if np.array(a <= 1.).any():
        raise ValueError('Power law exponent must be >1.')
    n = len(EDs)
    if n == 0:
        raise ValueError('No data.')
    Y = EDs.min()
    if Y < 0:
        raise ValueError('Negative value encountered in data.')
    a = de_bias_alpha(n, a)
    Yexp = (np.power(Y,1-a))
    T = np.log(EDs).sum()
    Z = de_biased_upper_limit(EDs, a)
    Zexp = (np.power(Z,1-a))

    return n / (a - 1) + n * ((Zexp * np.log(Z) - Yexp * np.log(Y)) / (Zexp - Yexp)) - T


def de_biased_upper_limit(data, a):
    '''
    De-biases the upper limits for a
    ML power law exponent estimator.
    Uses formular (13) (and (14)) from
    Maschberger and Kroupa (2009).

    Parameters:
    -----------
    data : Series or array
        data that is suspected to follow
        a power law relation
    a : float or array of floats
        quasi de-biased ML estimator for alpha
        (de_bias_alpha before inserting here!)

    Returns:
    ---------
    Quasi de-biased upper limit.
    '''
    if len(data) == 0:
        raise ValueError('No data.')
    if (data < 0).any():
        raise ValueError('Negative values '
                         'encountered in data.')
    Xn = data.max()
    X1 = data.min()
    if Xn == X1:
        raise ValueError('Data range is zero.')
    n = len(data)
    G = (1. - a) * np.log(Xn / X1)#(14)
    base = 1.+ (np.exp(G) - 1.) / n
    exponent = 1. / (1. - a)
    return Xn * np.power(base, exponent)

def de_bias_alpha(n, alpha):
    '''
    De-biases the power law value
    according to Maschberger and Kroupa (2009),
    formula (12).

    Paramaters:
    ------------
    n : int
        Size of the data
    alpha : float or array of floats
        Power law exponent value from ML estimator

    Returns:
    -----------
    quasi de-biased ML estimator for alpha
    '''
    if np.array(np.isnan(n) | np.isnan(np.array(alpha)).any()):
        raise ValueError('de_bias_alpha: one or '
                         'both arg(s) is/are NaN')
    return (alpha - 1.) * n / (n - 2) + 1.
    
def plot_percentile_percentile(self, ax, sig_level=0.05, **kwargs):
    '''
    Plot the percentile-percentile, or
    probability-probability distribution, as
    suggested by Maschberger and Kroupa 2009.

    Parameters:
    --------------
    ax : Axes object
        panel to plot to
    sig_level : 0 < float < 1
        significance level for acceptance region
    '''
    if self.cutoff_ED_lower is not None:
        data = self.ED
    if self.cutoff_energy_lower is not None:
        data = self.energy
    alpha = self.alpha
    if alpha is None:
        raise ValueError('Compute power law exponent first.')
    sorted_data = np.sort(data)
    pp = calculate_cumulative_powerlaw_distribution(sorted_data, alpha)
    y = (np.arange(1, len(pp) + 1) - .5) / len(pp)
    limit = calculate_KS_acceptance_limit(len(data), sig_level=sig_level)
    ax.plot(pp, y, **kwargs)
    ax.plot(pp, pp + limit, c='k', label='$p = {}$'.format(sig_level))
    ax.plot(pp, pp - limit, c='k')
    ax.set_xlabel(r'$P_i$')
    ax.set_ylabel(r'S')
    ax.set_title('Percentile-percentile plot with acceptance region')
    return
        
def fit_beta_to_powerlaw(a, freq, alpha, alpha_err, tot_obs_time, mode="ED"):
    '''
    Fit beta via non-linear least squares to a power
    law with given alpha using the cumulative
    FFD. Generate uncertainty using jackknife algorithm.
    '''
    def LSQ(x0,a,freq,alpha):
        zw = ((x0 / (np.power(a,alpha-1.) * (alpha-1.))-freq)**2).sum()
        return np.sqrt(zw)

    N = len(a)
    if N==0:
        raise ValueError('No data.')
        
    #jackknife uncertainty
    x0starts = {'ED' : 10, 'energy' : 1e25}
    _beta = np.array([fmin(LSQ,x0=x0starts[mode],
                          args=(np.delete(a,i, None),np.delete(freq,i,None),alpha),
                          disp=0)[0] for i in range(N)])
                          
    #cumulative beta = beta_cum
    print(_beta.mean())
    beta = _beta.mean() / tot_obs_time
    beta_err = np.sqrt( (N-1) / N * ( (_beta / tot_obs_time - beta)**2 ).sum() )
    
    #power law beta = beta_cum * |alpha-1|
    beta = beta * np.abs(alpha - 1.)
    
    #propagate errors on alpha to beta
    beta_err = np.sqrt(beta_err**2 * (alpha - 1.)**2 + beta**2 * alpha_err**2)
    
    return beta, beta_err

def plot_powerlaw(ax, a, alpha, beta, mode, custom_xlim=None, **kwargs):
    '''
    Plot the power law fit to the FFD.

    Parameters:
    -----------
    ax : matplotlibe Axes object
        plot to insert the power law in to
    kwargs : dict
        Keyword arguments to pass to plt.plot()

    Return:
    --------
    3 power law points to construct a line
    in log-log representation.
    '''
    if custom_xlim is None:
        x = np.linspace(np.nanmin(a),np.nanmax(a),3)
    else:
        mi, ma = custom_xlim
        x = np.linspace(mi, ma, 3)
    y = beta / np.abs(alpha - 1.) * np.power(x,-alpha+1.)
    ax.plot(x, y,  **kwargs)
    return