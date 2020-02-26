import pytest
import copy

import numpy as np
import pandas as pd

from ..ffd import (FFD,
                   _get_freq,
                   ML_powerlaw_estimator,
                   de_biased_upper_limit,
                   de_bias_alpha)


def generate_random_power_law_distribution(a, b, g, size=1, seed=None):
    """Power-law generator for pdf(x)\propto x^{g-1}
    for a<=x<=b
    """
    if seed is not None:
        np.random.seed(seed)
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)



def test_init_FFD():
    
    # Generate a flare table
    a, b, g, size = 10, 1e3, -1, 200
    pwl = generate_random_power_law_distribution(a, b, g, size=size, seed=80)
    df = pd.DataFrame({"ed_rec":pwl,
                   "ed_corr":pwl*1.2,
                   "recovery_probability":.8,
                   "TIC":list(np.arange(1,21))*10
                    })
    
    # init an FFD object
    ffd = FFD(df)
    
    # check correct initialisation
    assert ffd.f.shape[0] == 200
    assert (ffd.f.columns.values == ['ed_rec', 'ed_corr',
                                     'recovery_probability', 
                                     'TIC']).all()

def test_ed_and_freq():
    """Tests _ed_and_counts, too."""
    
    # Generate a flare table
    a, b, g, size = 10, 1e3, -1, 200
    pwl = generate_random_power_law_distribution(a, b, g, size=size, seed=80)
    df = pd.DataFrame({"ed_rec":pwl,
                   "ed_corr":pwl*1.2,
                   "recovery_probability":.8,
                   "TIC":list(np.arange(1,21))*10
                    })
    
    # Init an FFD object
    ffd = FFD(df)

    # Check if the results are correct for all combinations of parameters
    # single star
    # --------------------------------------------------------------------

    # no correction
    ed, freq = ffd.ed_and_freq(energy_correction=False,
                               recovery_probability_correction=False,
                                multiple_stars=False)
    assert (freq == np.arange(1,201)).all()
    assert (ed == df.ed_rec.sort_values(ascending=False)).all()

    # only energy correction
    ed, freq = ffd.ed_and_freq(energy_correction=True,
                               recovery_probability_correction=False,
                                multiple_stars=False)
    assert (freq == np.arange(1,201)).all()
    assert (ed == df.ed_corr.sort_values(ascending=False)).all()

    # energy and frequency correction
    ed, freq = ffd.ed_and_freq(energy_correction=True,
                               recovery_probability_correction=True,
                                multiple_stars=False)

    assert (ed == df.ed_corr.sort_values(ascending=False)).all()
    assert (freq == np.arange(1,201)/.8).all()

    # multiple stars
    # -------------------------------------------------------------------

    ffd = FFD(df)

    # You must pass a Key to ID
    with pytest.raises(KeyError):
        ed, freq = ffd.ed_and_freq(energy_correction=False,
                               recovery_probability_correction=False,
                                multiple_stars=True)

    ffd = FFD(df,ID="TIC")

    # no correction
    ed, freq = ffd.ed_and_freq(energy_correction=False,
                               recovery_probability_correction=False,
                                multiple_stars=True)

    assert (np.diff(freq) > 0.).all()
    assert (ed == df.ed_rec.sort_values(ascending=False)).all()

    # only energy correction
    ed, freq = ffd.ed_and_freq(energy_correction=True,
                               recovery_probability_correction=False,
                                multiple_stars=True)

    assert (np.diff(freq) > 0.).all()
    assert (ed == df.ed_corr.sort_values(ascending=False)).all()
    assert freq[0] == 1
    _f = copy.copy(freq)

    # energy and frequency correction
    ed, freq = ffd.ed_and_freq(energy_correction=True,
                               recovery_probability_correction=True,
                                multiple_stars=True)

    assert (ed == df.ed_corr.sort_values(ascending=False)).all()
    assert (np.diff(freq) > 0.).all()
    assert _f/.8 == pytest.approx(freq)
    
    # Check failing case:
    # -------------------------------------------------------------
    
    with pytest.raises(KeyError):
        ed, freq = ffd.ed_and_freq(energy_correction=False,
                               recovery_probability_correction=True,
                                multiple_stars=True)
        
    with pytest.raises(KeyError):
        ed, freq = ffd.ed_and_freq(energy_correction=False,
                               recovery_probability_correction=True,
                                multiple_stars=False)

def test__get_freq():
    
    # Generate a flare table
    N = 20
    testdf = pd.DataFrame({"ID":np.arange(N)%4,
                           "sortcol":np.arange(200, 200-N, -1)})
    
    # call _get_freq
    f = _get_freq(testdf, "ID", "sortcol")
    
    # Check if the result is as expected
    assert (f == np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                           1., 1., 1., 1., 1., 1., 0.75, 0.5, 0.25])).all()

    # If the required keys don't exist throw errors
    with pytest.raises(KeyError):
        testdf = pd.DataFrame({"ID":np.arange(N)%4})
        f = _get_freq(testdf, "ID", "sortcol")

    with pytest.raises(KeyError):
        testdf = pd.DataFrame({"sortcol":np.arange(200, 200-N, -1)})
        f = _get_freq(testdf, "ID", "sortcol")
        
def test_fit_powerlaw():
    # Generate a flare table
    a, b, g, size = 10, 1e3, -1, 200
    pwl = generate_random_power_law_distribution(a, b, g, size=size, seed=80)
    df = pd.DataFrame({"ed_rec":pwl,
                       "ed_corr":pwl*1.2,
                       "recovery_probability":.8,
                       "TIC":list(np.arange(1,21))*(size//20)
                       })

    # init an FFD object
    ffd = FFD(df)
    
    # check the result
    assert (1.963553855895996, 0.08012203082491737) == ffd.fit_powerlaw(df.ed_rec.values)
    
    
def test_ML_powerlaw_estimator():
    dataformat = [np.array, pd.Series]
    for df in dataformat:
        ed = df([1,1,1,1,2,2,4])
        x0 = 1.
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(x0, ed)
        x0 = 1.5
        assert ML_powerlaw_estimator(x0, ed) == pytest.approx(1.6190804181576444)
        ed = df([])
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(x0, ed)
        ed = df([-3,2,2,2])
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(x0, ed)
        ed = df([1,1,1,1])
        with pytest.raises(ValueError):
            ML_powerlaw_estimator(x0, ed)

            
def test_de_biased_upper_limit():
    dataformat = [np.array, pd.Series]
    for df in dataformat:
        data = df([])
        with pytest.raises(ValueError):
            de_biased_upper_limit(data, 2.)
        data = df([1,10, 100])
        assert de_biased_upper_limit(data, 1000000.) == pytest.approx(100., rel=1e-4)
        data = df([1,1,1])
        with pytest.raises(ValueError):
            de_biased_upper_limit(data, 3.)
        data = df([-1,1,1])
        with pytest.raises(ValueError):
            de_biased_upper_limit(data, 3.)

            
def test_de_bias_alpha():
    assert de_bias_alpha(200,1) == 1.
    with pytest.raises(ValueError):
        de_bias_alpha(np.nan,2)
    with pytest.raises(ValueError):
        de_bias_alpha(30,np.nan)
    with pytest.raises(ValueError):
        de_bias_alpha(np.nan,np.nan)
    with pytest.raises(ZeroDivisionError):
        de_bias_alpha(2,2)
