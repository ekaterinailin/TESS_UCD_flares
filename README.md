# TESS_UCD_flares - Find and analyse flares in TESS UCD light curves

Schmidt et al. (2020) (in prep.) studied flare in ultracool dwarfs in TESS Cycle 1 light curves. Here we present the de-trending prescription,
and the procedure of injection-recovery of synthetic flares. 

## Contents

`custom_apertures/`

`flare_tables/`

`notebook/`

    - `sample_injrec.py`: Use this script to sample injection and recovery of flares in a light curves
    - `NB4_Vet_flares.ipynb`: A notebook that allows you to look at the light curve and its flares at different stages of the analysis process. **NEEDS CLEANING** 
    - `NB3_....ipynb`: A notebook in which you can use your injection-recovery sample to characterize flares that you found. **NEEDS CLEANING**
    - `NB1....ipynb`: A notebook to simply show how to de-trend and find flares using a custom detrending procedure.
    - `funcs/`
        - `helper.py`
        - `custom_detrending.py`
        - `tests/` tests to run with `pytest funcs/`
    

## Installation requirements

### AltaiPony

Follow the installation instructions for AltaiPony [here](https://altaipony.readthedocs.io/en/latest/install.html#installation).

### Other requirements

- numpy
- pandas
- scipy
- astropy

## Custom de-trending

TESS ultracool dwarf light curves exhibit various trends, both astrophysical and systematic that we need to remove before searching the residual time series for flares. The main trends are

- movement of the spacecraft resulting in both linear and non-linear trends,
- stellar rotational modulation caused by spots that is sometimes strictly periodic and sinusoidal.

The latter can be an extremely fast modulation of the order of a few hours. 

To account for these effects we adopted the following procedure:

1. If the difference in flux from the beginning to the end of the light curve is >20% we fitted a 3rd order spline to the binned light curve to remove strong global trends that take effect on time scales of multiple days. 
2. When no strong global trends are present, strong rotational modulation can be removed iteratively using sine fits with periods obtained from Lomb-Scargle periodograms until not strong periods are left in the light curve.
3. Next, we remove remaining trends with using a rolling median on time scales of 10 h.
4. To remove trends on shorter time scales, we mask and pad outliers, and apply a Savitzky-Golay filter with window sizes that are three times shorter that the dominant residual frequency. The minimum window size is 2.5 h and is determined using a FFT of the residual light curve.
5. In a final step we filter the light curve with the minimum window size  and estimate the scatter in the light curve using a rolling standard deviation (window size= 30 min) with masked and padded positive outliers, i.e. flare candidates. The masked and padded values are assigned the mean uncertainty from the rest of the light curve.


