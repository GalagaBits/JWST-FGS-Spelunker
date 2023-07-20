# Spelunker — NIRISS FGS quicklook pipeline 

![](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/spelunker.png)

-------------------------------------------------------------------------------------

`spelunker` is a package that assists and finds technical anomalies and stellar properties from guidestars. 

**Authors:** Derod Deal (dealderod@ufl.edu), Néstor Espinoza (nespinoza@stsci.edu)

## Statement of need

Every time JWST observes an object, it simultaneously observes a nearby star --- a so-called "guide star" --- with the NIRISS Fine Guidance Sensor (FGS) that is used to keep the telescope locked on the target of interest. While researchers typically focus on their science targets, the guide star data can be extremely interesting on its own right. On the one hand, telescope-level anomalies could be detected (and, in principle, corrected) using this guide star data. On the other, this data also provides a "free" sky survey in the infrared (0.6 to 5 microns), on which short (~hours to days) time series of stars are recorded --- which one could "mine" if a pipeline existed for it to search for, e.g., stellar variability or even exoplanet transits: a true treasure chest. Here we present a first version of an automated, public quick-look time-series data processing pipeline for NIRISS FGS data. The pipeline is able to generate time-series for several metrics of the FGS data in an automated fashion, including fluxes and PSF variations, along with derived products from those such as periodograms that can aid on their analysis given only a JWST program ID number.

## Installation

To install `spelunker`, use `pip install`.

```bash
pip install spelunker
```

## Using the library

Get started with `spelunker` with only two lines of code.

```python
import spelunker

spk = spelunker.load('/Users/ddeal/JWST-Treasure-Chest/', pid=1534)
```
With our object `spk`, we can start interpeting data from guidestars, such as making a plot of the tracked guidestars within a Program ID.

```python
spk.guidestar_plot()
```
<p align='center'>
    <img src="https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/plots/guidestar_positions.png"  width=80% height=80%>
</p>

We can also plot the timeseries of fitted gaussian parameters (use `spk.gauss2d_fit` to apply gaussian fitting to all guidestar frames) or the flux timeseries. Mnemonics from JWST technical events can be overplotted on any timeseries, such as high-gain antenna (HGA) movement or the FGS tracks a new guidestar.

```python
import matplotlib.pyplot as plt

spk.mast_api_token = 'insert a token from auth.MAST here'

fig, ax = plt.subplots(figsize=(12,4),dpi=200)

ax = spk.mnemonics_local('GUIDESTAR')
ax = spk.mnemonics('SA_ZHGAUPST', 60067.84, 60067.9) 
ax.plot(spk.fg_time, spk.fg_flux)
plt.legend(loc=3)
plt.xlim(60067.84, 60067.9)
plt.show()
```
<img src="https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/plots/mnemonics.png">

For more information on the tools under `spelunker` and how to get started, visit the [quickstart guide](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/notebooks/fgs-spelunker_quickstart.ipynb). Get acquainted with `spelunker` with the following example notebooks:

- [Guidestar Targets](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/notebooks/examples/guidestar_targets.ipynb)
- [Pixel centroid changes and mnemonics](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/notebooks/examples/pixel_centroid_mnemonics.ipynb)


## Licence and attribution

This project is under the MIT License, which can be viewed [here](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/LICENSE).