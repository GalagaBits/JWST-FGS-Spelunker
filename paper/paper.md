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

# Statement of need

The James Webb Space Telescope [@gardner_james_2023] produces some of the highest sensitivity imaging of the cosmos across all instruments. One of them, the NIRISS Fine Guidance Sensor [@doyon_jwst_2012], provides guide star imaging with a passband of 0.6 to 5 microns through two separate channels, each with a $2.3’ \times 2.3’$ field of view (FOV) and a sampling rate of 64 ms—data that is taken in parallel and is thus available for every JWST observing program. While the onboard system uses guidestars to guide the attitude control system (ACS) which stabilizes the observatory, the astronomical community can also use the data products associated with these 64 ms cadence images as science products. Usages range from studying guide star photometry in search of transient phenomena to using this data to identify and investigate technical anomalies that might occur during scientific observations with the rest of the JWST instruments. Despite this wide range of possible usages, these data products are not straightforward to manipulate and analyze, and there is no publicly available package to download, investigate, and research guidestar data. ``spelunker`` is a Python library that was developed to enable access to these guide star data products and their analysis.

# Summary

``spelunker`` is a Python library that provides several tools for analyzing and visualizing JWST NIRISS FGS guide star data products. The library can download guide star data from the Mikulski Archive for Space Telescopes [@marston_overview_2018] for any given Program ID in a single line of code:

```python
spk = spelunker.load(pid=1534)
```

Through efficient parallelization processes, the pipeline can use this data to extract photometry in seconds and point-spread function information from all frames within a selected Program ID in minutes. Data loaded within the package allows for many possible analyses relevant to existing science programs and to perform science with the JWST FGS guide star data. In addition, ``spelunker`` provides visualization and analysis tools to study these products in detail, including the incorporation of JWST engineering telemetry, which can be used to put the JWST FGS data products in context with other observatory variables that might help explain data patterns both in the primary science data products and on the JWST FGS data itself. The library also cross-references, tracks and stores guide star metadata for the user. This metadata includes information such as the GAIA ID of the guide star, coordinates, and magnitudes, among others.

# Acknowledgements

We would like to thank the Space Telescope Science Institute and the National Astronomy Consortium for the opportunity to develop this project. In particular, we acknowledge funding and support from the 2023 version of the Space Astronomy Summer Program (SASP) at STScI that made it possible for the authors to work together on this project.

# References
