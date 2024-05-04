---
title: 'Spelunker: A quick-look Python pipeline for JWST NIRISS FGS Guide Star Data'
tags:
  - Python
  - JWST
  - time domian astronomy
  - NIRISS/FGS
authors:
  - name: Derod Deal
    orcid: 0009-0006-6758-4751
    affiliation: 1
  - name: Néstor Espinoza
    orcid: 0000-0001-9513-1449
    affiliation: 2
affiliations:
 - name: Department of Astronomy, University of Florida P.O. Box 112055, Gainesville, FL, USA
   index: 1
 - name: Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
   index: 2
date: 14 November 2023
bibliography: paper.bib
---


# Summary

``spelunker`` is a Python library that provides several tools for analyzing and visualizing JWST NIRISS FGS guide star data products. The library can download guide star data from the Mikulski Archive for Space Telescopes [@marston_overview_2018] for any given Program ID in a single line of code. Through an efficient parallelization process, the pipeline can use this data to extract photometry in seconds and point-spread function information from all frames within a selected Program ID in minutes. Data loaded within the package allows for many possible analyses such as Gaussian fitting, periodograms of each fitted Gaussian parameter, flux time series, and photometry optimization.

In addition, ``spelunker`` provides visualization and analysis tools to study these products in detail, including the incorporation of JWST engineering telemetry, which can be used to put the JWST FGS data products in context with other observatory variables that might help explain data patterns both in the primary science data products and in the JWST FGS data itself. The library also cross-references, tracks and stores guide star metadata for the user. This metadata includes information such as the GAIA ID of the guide star, coordinates, and magnitudes, among others. This pipeline empowers users to study guide star data with applications like technical anomaly detection and time-domain astronomy at a high 64 ms cadence. Existing pipeline such as ``jwst`` provides the structure for calibrataing NIRISS/FGS, NIRSpec, NIRCam, and MIRI data, such as with JWST ``Datamodels``. ``spelunker``  tools are engineered to explore FGS guide stars in two main categories: as a method to detect technical anomalies in JWST data, and as a pipeline to study time-domain astronomy with guide stars.

# Statement of need

The James Webb Space Telescope [@gardner_james_2023] produces some of the highest sensitivity imaging of the cosmos across all instruments. One of them, the NIRISS Fine Guidance Sensor [@doyon_jwst_2012], provides guide star imaging with a passband of 0.6 to 5 microns through two separate channels, each with a $2.3’ \times 2.3’$ field of view (FOV) and a sampling rate of 64 ms—data that is taken in parallel and is thus available for every JWST observing program. While the onboard system uses guide stars to provide information to the attitude control system (ACS) which stabilizes the observatory, the astronomical community can also use the data products associated with these 64 ms cadence images as science products. Usages range from studying guide star photometry in search of transient phenomena to using these data to identify and investigate technical anomalies that might occur during scientific observations with the rest of the JWST instruments. Despite this wide range of possible usages, these data products are not straightforward to manipulate and analyze, and there is no publicly available package to download, investigate, and research guide star data. ``spelunker`` is a Python library that was developed to enable access to these guide star data products and their analysis.

![There are seven parameters `gauss2d_fit` measures: amplitude (counts of the guide star), x pixel coordinate, y pixel coordinate, the x and y standard deviations, theta (orientation of the Gaussian model), and the offset (the background counts). This diagram visualizes what each parameter represents on the Gaussian model. \label{fig:Gaussian_diagram}](Gaussian_diagram.png)


# Overview of Spelunker

``spelunker`` allows anyone to download and utilize fine guidence guides tar data from the JWST. Users can use the following lines of code to download GS-FG data into an object:

```python
spk = spelunker.load(pid=1534)
```

``spelunker`` uses ``astroquery`` and MAST [@marston_overview_2018] to find and download GS-FG FITS files. There are several functions that can manipulate and analyze guide star data:

- **``gauss2d_fit``** A spatial Gaussian fit is applied to each of the frames loaded into the `spk` object. The Gaussian will fit the amplitude, pixel coordinates, pixel standard deviations, the model orientation, and the background offset. A diagram of Gaussian fitting parameters is shown in \autoref{fig:Gaussian_diagram}. The Gaussian measurements are then stored in ``spk.gaussfit_results`` as an astropy table. Fitting spatial Gaussians to your guide star data can allow you to reveal technical anomalies that might be confused with science data from NIRISS or other JWST instruments and within the guide star flux time series (see \autoref{fig:guidestar_1803}).


![A  snippet from the guide star time series from Cycle 1 GO Program ID 1803, observation 1 and visit 1. **Top** — The guide star time series of PID 1803 after loading it into ``spelunker`` using ``timeseries_binned_plot``. The time series uses the sum of counts in each guide star fine guidance (GS-FG) frame. The data has no significant features. **Middle** — The same time series after applying ``optimize_photometry`` to the guide star light curve. There are now prominent drops in the flux, which were previously unseen with the raw time series data. **Bottom** — Gaussian fitted x pixel coordinate and y pixel coordinate for each frame in this section of time series data. The guide star shifts around in this time series, highlighting the core function of the ACS. \label{fig:guidestar_1803}](timeseries_plot.png)

- ** ``mnemonics``** Users can access JWST engineering telemetry and mnemonics using the ``mnemonics`` function. With a MAST API token, any mnemonic is accessible. High-gain antenna (HGA) movement and NIRISS filter wheel current are two examples of events that overplot science data to identify technical events on the telescope. Anomalies detection in guide star data or data from NIRISS, NIRCAM, NIRSpec, and MIRI is one of the primary capabilities of this function. Using ``mnemonics`` requires the additional installation of the ``jwstuser`` library.

- **``periodogram``** This function uses the Lomb-Scarle periodogram [@lomb_least-squares_1976;@scargle_studies_1982] to detect periodicities in guide star Gaussian fits. Periods in Gaussian fitted parameters like x and y pixel coordinates highlight systematics for an entire PID.

- **``optimize_photometry``** ``optimize_photometry`` extracts the highest SNR pixels to optimize raw guide star photometry loaded from ``spelunker``. \autoref{fig:guidestar_1803} demonstrates that ``optimize_photometry`` reveals more information from guide star time series than that produced by the sum of counts in each frame.

With the mentioned tools, ``spelunker`` utilizes object oriented programming (OOP) to store handy variables and its outputs, for instance, 1D and 2D time series, guide star time arrays, and JWST data models. Running ``gauss2d_fit``, ``periodogram``, and ``mnemonics`` will store their outputs in accessible attributes. Useful properties of the guide star are stored in these attributes (for instance, guide star galactic coordinates, GAIA ID, and stellar magnitudes). 

There are various plotting and visualization tools integrated into ``spelunker``'s workflow. One useful function is ``timeseries_binned_plot``, which automatically plots a binned time series. The functions ``gauss2d_fit``, ``periodogram``, and ``mnemonics`` have ``matplotlib`` axes returned for straightforward plotting. Animations of spatial time series are another visualization tool covered under ``spelunker``. Of course, users have the option to use the guide star data and results from attributes to generate plots. 


# Acknowledgements

We would like to thank the Space Telescope Science Institute and the National Astronomy Consortium for the opportunity to develop this project. In particular, we acknowledge funding and support from the 2023 version of the Space Astronomy Summer Program (SASP) at STScI that made it possible for the authors to work together on this project. 

# References
