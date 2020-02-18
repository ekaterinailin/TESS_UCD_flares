# TESS_UCD_flares - Find and analyse flares in TESS UCD light curves

Schmidt et al. (2020) (in prep.) studied flare in ultracool dwarfs in TESS Cycle 1 light curves. Here we present the de-trending prescription,
and the procedure of injection-recovery of synthetic flares. 

## Contents

- `custom_apertures/`: custom aperture light curve fits files
- `flare_tables/`: final results 
- `notebook/`
 - `sample_injrec.py`: Use this script to sample injection and recovery of flares in a set of light curves
 - `analyse_injrec.py`: Use this script to combine stored injection recovery tables with light curves to characterize them and store final  results in `flare_tables/`
 - `NB1_Find_flares_in_custom_aperture_light_curves.ipynb`: A notebook to simply show how to de-trend and find flares using a custom detrending procedure, adn then characterize the candidates with injection-recovery of synthetic flares.
 - `NB4_Vet_flares.ipynb`: A notebook that allows you to look at the light curve and its flares at different stages of the analysis process. **NEEDS CLEANING** 
 - `NB3_....ipynb`: A notebook in which you can use your injection-recovery sample to characterize flares that you found. **NEEDS CLEANING**
 - `funcs/`
   - `helper.py` reading files, masking light curves, etc.
   - `custom_detrending.py` custom_detrending function lives here with all the functions that it calls
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

## Injection-recovery

When an iterative and complex but deterministic procedure is used to de-trend a light curve there is no analytic way to determine

- the recovery probability or
- the true energy 

of a flare. Every light curve is different, and recovery probability and measured energy largely depend on its amplitude and duration convolved with time sampling, and de-trending effects. Injecting synthetic flares spanning the parameter space of amplitude and duration into the not yet de-trended light curve resolves both problems at once. The light curve with synthetic flares injected can undergo the same de-trending procedure as the one without. Recovering the injected flares results in a measure for recovery probability as a function of amplitude and duration, and the difference between the injected and recovered flare energies can be used as a correction factor for the flares found in the original light curve. 

Injection-recovery involves the following seven steps:

1. Randomly generate flare amplitude, duration and peak time for a synthetic flare. Mask the flare candidates from the original light curve to avoid overlap.
2. Add the synthetic flare to the original light curve.
3. De-trend that light curve using the same procedure that was used to detect the original flare candidates.
4. Determine the properties for the recovered flare and save them together with the injected ones.
5. Repeat 1.-4. such that the parameter space covers the properties of the original flares.
6. Apply the correction factor to amplitude and duration of the original flare to find its intrinsic properties.
7. Find the recovery probability of flares with these intrinsic properties.

Injection recovery must be performed on every light curve individually, and computational cost scale with the complexity of the de-trending prescription.




This procedure assumes that the flares follow a certain flare shape, here as it is described in Davenport 2014.


