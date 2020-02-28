import pytest
import copy

import numpy as np
import pandas as pd

from ..ffd import (FFD,
                   _get_multistar_factors,
                   ML_powerlaw_estimator,
                   de_biased_upper_limit,
                   de_bias_alpha,
                   calculate_average_number_of_exceeding_values,
                   _calculate_number_of_exceeding_values,
                   stabilised_KS_statistic,
                   calculate_KS_acceptance_limit,
                   _apply_stabilising_transformation
                    )


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

def test__get_multistar_factors():
    
    # Generate a flare table
    N = 20
    testdf = pd.DataFrame({"ID":np.arange(N)%4,
                           "sortcol":np.arange(200, 200-N, -1)})
    
    # call _get_multistar_factors
    f = _get_multistar_factors(testdf, "ID", "sortcol")
    
    # Check if the result is as expected
    assert (f == np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                           1., 1., 1., 1., 1., 1., 0.75, 0.5, 0.25])).all()

    # If the required keys don't exist throw errors
    with pytest.raises(KeyError):
        testdf = pd.DataFrame({"ID":np.arange(N)%4})
        f = _get_multistar_factors(testdf, "ID", "sortcol")

    with pytest.raises(KeyError):
        testdf = pd.DataFrame({"sortcol":np.arange(200, 200-N, -1)})
        f = _get_multistar_factors(testdf, "ID", "sortcol")
        
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
    
def test_is_powerlaw_truncated():
    # Generate a flare table
    a, b, g, size = 10, 1e3, -1, 200
    pwl = generate_random_power_law_distribution(a, b, g, size=size, seed=80)
    ffd = FFD()
    ffd.alpha = 2.

    # sort in ascending order
    sortpwl = np.sort(pwl)

    # remove highest energies to see where it starts to become truncated
    assert ffd.is_powerlaw_truncated(pwl) == False
    assert ffd.is_powerlaw_truncated(sortpwl[:-8]) == False
    assert ffd.is_powerlaw_truncated(sortpwl[:-9]) == True
    assert ffd.is_powerlaw_truncated(sortpwl[:-20]) == True
    assert ffd.is_powerlaw_truncated(sortpwl[:-100]) == True
    
    
def test_is_powerlaw():

    a, b, g, size = 10, 1e3, -1., 200
    pwl = generate_random_power_law_distribution(a, b, g, size=size, seed=80)

    # init an FFD object
    ffd = FFD()

    ffd.alpha = 2.
    # pwl is a power law with exponent 2
    assert ffd.is_powerlaw(pwl)

    ffd.alpha = 2.3
    # pwl is not a power law with exponent 2
    assert not ffd.is_powerlaw(pwl)


    a, b, g, size = 10, 1e3, -1., 20
    pwl = generate_random_power_law_distribution(a, b, g, size=size, seed=80)

    # init an FFD object
    ffd = FFD()

    ffd.alpha = 2.
    # pwl is a power law with exponent 2
    assert ffd.is_powerlaw(pwl)

    ffd.alpha = 2.3
    # pwl is not a power law with exponent 2 but 20 is too small of a sample
    assert ffd.is_powerlaw(pwl)

    a, b, g, size = 10, 1e3, -1., 50
    pwl = generate_random_power_law_distribution(a, b, g, size=size, seed=80)

    with pytest.raises(TypeError):
        ffd.is_powerlaw()
    
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
        

def test_calculate_average_number_of_exceeding_values():
    data = np.linspace(10,1e4, 300)
    alpha = 2.
    mean, std = calculate_average_number_of_exceeding_values(data, alpha, 1000, seed=2311)
    assert mean == 3.
    assert std == 0.
    mean, std = calculate_average_number_of_exceeding_values(data, alpha, 1000, seed=10)
    assert mean == 0.
    assert std == 0.


def test__calculate_number_of_exceeding_values():
    data = np.linspace(10,1e4, 300)
    assert _calculate_number_of_exceeding_values(data, 2., seed=10) == 0
    assert _calculate_number_of_exceeding_values(data, 2., seed=2311) == 3
    with pytest.raises(ValueError):
        _calculate_number_of_exceeding_values(np.arange(3), 2., seed=2311)
        

def test_stabilised_KS_statistic():
    sizes = [1e2,1e3,1e4]
    minval, maxval = 10, 1e4
    datas = [generate_random_power_law_distribution(minval, maxval, -1., size=int(size), seed=10) for size in sizes]
    
    KSlist = [stabilised_KS_statistic(data, 2., False) for data in datas]
    assert KSlist[0] > KSlist[1]
    assert KSlist[1] > KSlist[2]
    

def test_calculate_KS_acceptance_limit():
    with pytest.raises(ValueError):
        calculate_KS_acceptance_limit(0, sig_level=0.05)
    with pytest.raises(ValueError):
        calculate_KS_acceptance_limit(0, sig_level=-0.05)
    with pytest.raises(ValueError):
        calculate_KS_acceptance_limit(0, sig_level=1.05)
    assert (calculate_KS_acceptance_limit(100, sig_level=0.05)
            == pytest.approx(0.13581, rel=1e-4))
    assert (calculate_KS_acceptance_limit(100, sig_level=0.01)
            == pytest.approx(0.16276, rel=1e-4))
    assert (calculate_KS_acceptance_limit(100, sig_level=0.01)
            > calculate_KS_acceptance_limit(1000, sig_level=0.01))
    

def test__apply_stabilising_transformation():
    u = [.1,.2,.3]
    assert _apply_stabilising_transformation(u).shape[0] == 3
    assert (np.diff(u) > 0).all()

    # giving an empty arrays returns an empty array
    assert (_apply_stabilising_transformation([]) == np.array([])).all()

    # Passing negative values throws an error
    with pytest.raises(ValueError):
        _apply_stabilising_transformation([-1.,.2,.4])