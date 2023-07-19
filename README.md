# Spelunker — NIRISS FGS quicklook pipeline 

![](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/spelunker.png)

JWST FGS Spelunker is a package that assists and finds technical anomalies from the JWST and stellar properties from guidestars. 

### Abstract

Every time JWST observes an object, it simultaneously observes a nearby star --- a so-called "guide star" --- with the NIRISS Fine Guidance Sensor (FGS) that is used to keep the telescope locked on the target of interest. While researchers typically focus on their science targets, the guide star data can be extremely interesting on its own right. On the one hand, telescope-level anomalies could be detected (and, in principle, corrected) using this guide star data. On the other, this data also provides a "free" sky survey in the infrared (0.6 to 5 microns), on which short (~hours to days) time series of stars are recorded --- which one could "mine" if a pipeline existed for it to search for, e.g., stellar variability or even exoplanet transits: a true treasure chest. Here we present a first version of an automated, public quick-look time-series data processing pipeline for NIRISS FGS data. The pipeline is able to generate time-series for several metrics of the FGS data in an automated fashion, including fluxes and PSF variations, along with derived products from those such as periodograms that can aid on their analysis given only a JWST program ID number. We present preliminary analyses on a handful of the longest FGS time series to date, highlighting some of the properties of the data, how it can help JWST users and the prospects of using this data for search for astrophysical signals in the archive.

## FGS Spelunker Quickstart

To install and use FGS Spelunker, see the [quickstart guide](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/notebooks/fgs-spelunker_quickstart.ipynb).


![](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/plots/1541movie.gif)

![](https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/plots/guidestar_positions.png)