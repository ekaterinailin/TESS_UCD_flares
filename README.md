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



