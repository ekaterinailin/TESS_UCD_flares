import copy
import pandas as pd
import numpy as np

from scipy.optimize import fmin

class FFD(object):
    """Flare frequency distribution.
    
    
    """
    def __init__(self, dataframe, alpha=None, alpha_err=None,
                 beta=None, beta_err=None, total_obs_time=1.,
                 ID=None):
        
        self.f = dataframe
        self.alpha = alpha
        self.alpha_err = alpha_err
        self.beta = beta
        self.beta_err = beta_err
        self.total_obs_time = total_obs_time
        self.ID = ID
    
    def ed_and_freq(self, energy_correction=False,
                    recovery_probability_correction=False,
                    multiple_stars=False):
        
        if ((energy_correction==False) & (recovery_probability_correction==False)):
            key = "no_corr"
        
        elif ((energy_correction==True) & (recovery_probability_correction==False)):
            key = "ed_corr"
        
        elif ((energy_correction==True) & (recovery_probability_correction==True)):
            key = "edrecprob_corr"
            
        else:
            raise KeyError("This set of parameters for energy correction, recovery" \
                           " probability correction is not implemented.")
            
        ed, f = self._ed_and_counts(key, multiple_stars)
        
        return ed, f / self.total_obs_time
    
    def _ed_and_counts(self, key, multiple_stars):
        
        def cum_dist(df, col, ID):
            
            return np.cumsum(np.ones_like(df[col].values))
        
        def getfreq_cum_dist(df, col, ID):
            
            freq = _get_freq(df, ID, col)
            return np.cumsum(1 / freq)
        
        def cum_dist_rec_prob(df, col, ID):
            
            return np.cumsum(1. / df.recovery_probability.values)
        
        def getfreq_cum_dist_rec_prob(df, col, ID):
            
            freq = _get_freq(df, ID, col)
            return np.cumsum(1. / df.recovery_probability.values / freq)
            
        
        vals = {"no_corr":{False: ["ed_rec", cum_dist],
                           True: ["ed_rec", getfreq_cum_dist]},
                "ed_corr":{False: ["ed_corr", cum_dist],
                           True: ["ed_corr", getfreq_cum_dist]},
                "edrecprob_corr":{False: ["ed_corr", cum_dist_rec_prob],
                           True: ["ed_corr", getfreq_cum_dist_rec_prob]}}
        
        df = self.f.copy(deep=True)
        col, func = vals[key][multiple_stars]
        df = df.sort_values(by=col, ascending=False)
        
        ed = df[col].values
        freq = func(df, col, self.ID)
        
        return ed, freq
    
        
    def fit_beta_to_powerlaw(self, ed, freq, mode="ED"):
        '''
        Fit beta via non-linear least squares to a power
        law with given alpha using the cumulative
        FFD. Generate uncertainty using jackknife algorithm.
        '''
        def LSQ(x0, ed, freq, alpha):
            zw = ((x0 / (np.power(ed,alpha-1.) * (alpha-1.))-freq)**2).sum()
            return np.sqrt(zw)

        N = len(ed)
        if N==0:
            raise ValueError('No data.')
            
        #jackknife uncertainty
        x0starts = {'ED' : 10, 'energy' : 1e25}
        _beta = np.array([fmin(LSQ,x0=x0starts[mode],
                              args=(np.delete(ed,i),np.delete(freq,i),self.alpha),
                              disp=0)[0] for i in range(N)])

        #cumulative beta = beta_cum
        beta = _beta.mean() / self.tot_obs_time
        beta_err = np.sqrt( (N-1) / N * ( (_beta / self.tot_obs_time - beta)**2 ).sum() )

        #power law beta = beta_cum * |alpha-1|
        beta = beta * np.abs(self.alpha - 1.)

        #propagate errors on alpha to beta
        beta_err = np.sqrt(beta_err**2 * (self.alpha - 1.)**2 + beta**2 * self.alpha_err**2)

        #set attributes
        self.beta = beta
        self.beta_err = beta_err
        
        return _beta * np.abs(self.alpha - 1.), beta, beta_err
    
    def plot_powerlaw(self, ax, ed, custom_xlim=None, **kwargs):
        '''
        Plot the power law fit to the FFD. [No tests]

        Parameters:
        -----------
        ax : matplotlibe Axes object
            plot to insert the power law in to
        ed : array
            array of ED or energy values in the FFD
        custom_xlim : 2-tuple
            minimum, maximum ED/energy value for power law
        kwargs : dict
            Keyword arguments to pass to plt.plot()

        Return:
        --------
        3 power law points to construct a line
        in log-log representation.
        '''
        if custom_xlim is None:
            x = np.linspace(np.nanmin(ed),np.nanmax(ed),3)
        else:
            mi, ma = custom_xlim
            x = np.linspace(mi, ma, 3)
        y = self.beta / np.abs(self.alpha - 1.) * np.power(x,-self.alpha + 1.)
        ax.plot(x, y,  **kwargs)
        return

    def fit_powerlaw(self, ed, alims=[1.01,3.]):
        '''
        Calculate the un-biased ML power law estimator
        from Maschberger and Kroupa (2009), sections
        3.1.4. and 3.1.5.

        Parameters:
        ------------
        ed: array
            EDs or energies that supposedly folow a power law
        alims:
            parameter range for power law exponent
            
        Return:
        -------
        alpha, alpha_err - float, float
            power law exponent and its jackknife uncertainty
        '''
        
        # solve eq. 9 using scipy.fmin, define jacknife uncertainty
        N = len(ed)
        _alpha = np.array([fmin(ML_powerlaw_estimator, x0=2.,
                               args=(np.delete(ed,i),), disp=0)[0]
                          for i in range(N)])
        
        # alpha is the mean value
        alpha = _alpha.mean()
        
        # uncertainty is the standard deviation
        sig_alpha = np.sqrt( (N-1) / N * ( (_alpha - alpha)**2 ).sum() )
        
        return alpha, sig_alpha
    
def _get_freq(dataframe, ID, sort):
    
    freq = []
    df = dataframe.copy(deep=True)
    try:
        df = df.sort_values(by=sort, ascending=True)
    except:
        raise KeyError(f"The flare table needs a {sort} column.")
    
    for i in range(df.shape[0]):
        try:
            f = df.iloc[:i+1]
            freq.append(len(set(f[ID].values)))
        except KeyError:
            raise KeyError("Pass the column name of target IDs to the FFD constructor: ID = ???.")

    return np.array(freq[::-1]) / freq[-1]

def ML_powerlaw_estimator(x0, EDs):
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
    if np.array(x0 <= 1.).any():
        raise ValueError('Power law exponent must be >1.')
    n = len(EDs)
    if n == 0:
        raise ValueError('No data.')
    Y = EDs.min()
    if Y < 0:
        raise ValueError('Negative value encountered in data.')
    x0 = de_bias_alpha(n, x0)
    Yexp = (np.power(Y,1-x0))
    T = np.log(EDs).sum()
    Z = de_biased_upper_limit(EDs, x0)
    Zexp = (np.power(Z,1-x0))

    return np.abs(n / (x0 - 1) + n * ((Zexp * np.log(Z) - Yexp * np.log(Y)) / (Zexp - Yexp)) - T)
    

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
