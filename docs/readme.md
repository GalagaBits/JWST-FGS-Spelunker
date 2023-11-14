---
title: 'spelunker: A guidestar pipeline for the JWST's fine guidence sensor ('
tags:
  - Python
  - astronomy
  - jwst
authors:
  - name: Derod Deal
    orcid: 0009-0006-6758-4751
    affiliation: 1
  - name: Néstor Espinoza
    orcid: 0000-0001-9513-1449
    affiliation: 2
affiliations:
 - name: University of Florida, Department of Astronomy, P.O. Box 112055, Gainesville, FL, USA
   index: 1
 - name: Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
   index: 2 
date: 6 Oct 2023
bibliography: paper.bib

aas-doi: 
aas-journal: The Astronomical Journal
---

# Summary

When exoplanets pass in front of their stars from our point of view on Earth, they imprint a transit signature on the stellar light curve which, to date, has been assumed to be symmetric in time, owing to the planet being modelled as a circular area occulting the stellar surface [see, e.g., @Mandel02; @Kreidberg15; @Luger19]. However this signature might be asymmetric due to several possible effects, one of which is the different temperature/pressure and/or chemical compositions the different terminator regions a transiting planet could have [see, e.g., @Powell19]. Being able to model these asymmetric signatures directly from transit light curves could give us an unprecedented glimpse into planetary 3-dimensional structure, helping constrain models of atmospheric evolution, structure and composition.

``catwoman`` is a Python package that models these asymmetric transit light curves, calculating light curves for any radially symmetric stellar limb darkening law and where planets are modelled as two semi-circles, of different radii, using the integration algorithm developed in [@Kreidberg15] and implemented in the ``batman`` library, from which ``catwoman`` builds upon. It is fast and efficient and open source with full documentation available to view at https://catwoman.readthedocs.io .
     
The light curves are modelled as follows: The decrease in flux, $\delta$, as a planet transits its star can be approximated by the sum 

\begin{eqnarray}
\label{eq:theproblem}
\delta = \sum_{i=1}^{N} I\left(x_m\right)\Delta A(x_m,R_{p,1},R_{p,2},\varphi,d),
\end{eqnarray}

splitting the semi-circles into iso-intensity bands centred on the star and for each intersectional segment (see \autoref{fig:strips}) you multiply its area, $\Delta A$, by the intensity of the star and then sum these strips to generate the full $\delta$ for a specific separation between the centre of the star and planet, $d$. The code then increments $d$ by a small pre-determined amount (based on the time array given by the user) and recalculates $\delta$.

<!-- ![Diagram of the geometric configuration during transit of two stacked semi-circles (one of radius $R_{p,1}$, and another of radius $R_{p,2}$) that model the different limbs of an exoplanet transiting in front of a star. The area of the star has been divided in different sections of radius $x_i$ (dashed circles) --- between each subsequent section, the star is assumed to have a radially symmetric intensity profile (e.g., blue band between $x_{i-1}$ and $x_i$ above). In order to obtain the light curve, the challenge is to calculate the sum of the intersectional areas between a given iso-intensity band and the semi-circles, $\Delta A$ (blue band with dashed grey lines). Note the stacked semi-circles are inclined by an angle $\varphi$ with respect to the planetary orbital motion.\label{fig:strips}](strips.png) -->

The width of the iso-intensity bands determines the truncation error of the model. The model is first initialised with parameters including a maximum truncation error either set by the user or taken as the pre-set value as 1ppm. As in ``batman``, ``catwoman`` first calculates many models, with varying widths and geometrically searches for a width that produces an error less than 1% away (and always less than) the specified level. The model then uses this width value to calculate the desired light curves. A lower specified error, and therefore thinner iso-intensity bands, produces more accurate light curves, however more steps are needed to calculate $\delta$ which takes more time.  

``catwoman`` also allows for $\varphi$, the angle of rotation of the semi-circles, to vary as a free parameter, which is something no other model has tried to implement, accounting for the possibility of spin-orbit misalignments of the planet. The two semi-circle radii, $R_{p,1}$ and $R_{p,2}$, and other orbital variables are also completely free parameters.

``catwoman`` was designed to be used by astronomical researchers. For a realistic light curve with 100 in-transit data points, ``catwoman`` takes around 340 seconds to produce 1 million quadratic-limb-darkened light curves on a single 1.3 GHz Intel Core i5 processor. It is used in Espinoza & Jones (in prep.).

# Acknowledgements
We would like to thank the Max Plank Institute of Astronomy, Heidelberg, for providing the funding for this project and hosting Kathryn Jones as a summer student at the Institute. 

# References