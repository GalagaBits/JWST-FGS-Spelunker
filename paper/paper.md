---
title: 'Spelunker: A quick-look Python pipeline for JWST NIRISS FGS Guidestar Data'
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

``spelunker`` is a Python library that provides several tools for analyzing and visualizing JWST NIRISS FGS guide star data products. The library can download guide star data from the Mikulski Archive for Space Telescopes [@marston_overview_2018] for any given Program ID in a single line of code. Through efficient parallelization processes, the pipeline can use this data to extract photometry in seconds and point-spread function information from all frames within a selected Program ID in minutes. Data loaded within the package allows for many possible analyses such as Gaussian fitting, periodograms of each fitted Gaussian parameter and flux timeseries, and photometry optimization.

In addition, ``spelunker`` provides visualization and analysis tools to study these products in detail, including the incorporation of JWST engineering telemetry, which can be used to put the JWST FGS data products in context with other observatory variables that might help explain data patterns both in the primary science data products and on the JWST FGS data itself. The library also cross-references, tracks and stores guide star metadata for the user. This metadata includes information such as the GAIA ID of the guide star, coordinates, and magnitudes, among others. This pipeline empowers users to study guidestar data with applications like technical anomaly detection and time-domain astronomy at a high 64 ms cadence.

# Statement of need

The James Webb Space Telescope [@gardner_james_2023] produces some of the highest sensitivity imaging of the cosmos across all instruments. One of them, the NIRISS Fine Guidance Sensor [@doyon_jwst_2012], provides guide star imaging with a passband of 0.6 to 5 microns through two separate channels, each with a $2.3’ \times 2.3’$ field of view (FOV) and a sampling rate of 64 ms—data that is taken in parallel and is thus available for every JWST observing program. While the onboard system uses guidestars to guide the attitude control system (ACS) which stabilizes the observatory, the astronomical community can also use the data products associated with these 64 ms cadence images as science products. Usages range from studying guide star photometry in search of transient phenomena to using this data to identify and investigate technical anomalies that might occur during scientific observations with the rest of the JWST instruments. Despite this wide range of possible usages, these data products are not straightforward to manipulate and analyze, and there is no publicly available package to download, investigate, and research guidestar data. ``spelunker`` is a Python library that was developed to enable access to these guide star data products and their analysis.

![A snippet from the guidestar timeseries from Cycle 1 GO Program ID 1803, observation 1, and visit 1. **Top** — The guidestar timeseries of PID 1803 after loading it into ``spelunker`` using ``timeseries_binned_plot``. The timeseries uses the sum of counts in each guidestar fine guidence (GS-FG) frame. The data has no significant features. **Middle** — The same timeseries after applying ``optimize_photometry`` to the guidestar lightcurve. There are now prominent drops in the flux, which were previously unseen with the raw timeseries data. **Bottom** — Gaussian fitted x pixel coordinate and y pixel coordinate for each frame in this section of timeseries data. The guidestar shifts around in this timeseries, likely highlighting the functions of the ACS. \label{fig:guidestar_1803}](timeseries_plot.png)


# Overview of Spelunker

``spelunker`` allows anyone to download and ultilze fine guidence guidestar data from the JWST. Users can use the following lines of code to download GS-FG data into an object:

```python
spk = spelunker.load(pid=1534)
```

``spelunker`` uses ``astroquery`` and MAST [@marston_overview_2018] to find and download GS-FG FITS files. There are several functions that can manipulate and analyze guidestar data:

- **``gauss2d_fit``** A spatial Gaussian fit is applied to each of the frames loaded into the `spk` object. The Gaussian will fit the amplitude, pixel coordinates, pixel standard deviations, the model's orientation, and the background offset. A diagram of Gaussian fitting parameters is shown in \autoref{fig:Gaussian_diagram}. The Gaussian measurements are then stored in ``spk.gaussfit_results`` as an astropy table. Fitting spatial Gaussians to your guidestar data can allow you to reveal technical anomalies that might be confused with science data from NIRISS or other JWST instruments and within the Guidestar flux timeseries (see \autoref{fig:guidestar_1803}).


![There are seven parameters `gauss2d_fit` measures: amplitude (counts of the guidestar), x pixel coordinate, y pixel coordinate, the x and y standard deviations, theta (orientation of the Gaussian model), and the offset (the background counts). This diagram visualizes what each parameter represents on the Gaussian model. \label{fig:Gaussian_diagram}](Gaussian_diagram.png)

- **``mnemonics``** Users can access JWST engineering telemetry and mnemonics using the ``mnemonics`` function. With a MAST API token, any mnemonic is accessible. High-gain antenna (HGA) movement and NIRISS filter wheel current are two examples of events that overplot science data to identify technical events on the telescope. Anomalies detection in guidestar data or data from NIRISS, NIRCAM, NIRSpec, and MIRI is one of the primary abilities of this function.

- **``periodogram``** This function uses the Lomb-Scarle periodogram [@lomb_least-squares_1976;@scargle_studies_1982] to detect periodicities in guidestar Gaussian fits. Periods in Gaussian fitted parameters like x and y pixel coordinates highlight systematics for an entire PID.

- **``optimize_photometry``** ``optimize_photometry`` extracts the highest SNR pixels to optimize raw guidestar photometry loaded from ``spelunker``. \autoref{fig:guidestar_1803} demonstrates that ``optimize_photometry`` reveals more information from guidestar timeseries than a guidestar timeseries produced by the sum of counts in each frame.

With the mentioned tools, ``spelunker`` utilizes object oriented programming (OOP) to store handy variables and ``spelunker`` outputs, for instance, 1D and 2D timeseries, guidestar time arrays, and JWST datamodels. Running ``gauss2d_fit``, ``periodogram``, and ``mnemonics`` will store their outputs in accessible attributes. Useful properties of the guidestar are stored in these attributes (for instance, guidestar galactic coordinates, GAIA ID, and stellar magnitudes).

There are various plotting and visualization tools integrated into ``spelunker``'s workflow. One useful function is ``timeseries_binned_plot``, which automatically plots a binned timeseries. The functions ``gauss2d_fit``, ``periodogram``, and ``mnemonics`` have ``matplotlib`` axes returned for straightforward plotting. Animations of spatial timeseries are another visualization tool covered under ``spelunker``. Of course, users have the option to use the guidestar data and results from attributes to generate plots. 


# Acknowledgements

We would like to thank the Space Telescope Science Institute and the National Astronomy Consortium for the opportunity to develop this project. In particular, we acknowledge funding and support from the 2023 version of the Space Astronomy Summer Program (SASP) at STScI that made it possible for the authors to work together on this project. 

# References
