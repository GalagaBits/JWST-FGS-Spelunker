Spelunker
=========

.. image:: spelunker.png
    :width: 400pt
    :align: center
    

``spelunker`` is a package that assists and finds technical anomalies on the JWST and stellar properties from guidestars. This pipeline includes tools for
conveniently handling guidestar data products and provides methods for investigating its potential in scientific and technical observations.

This documentation showcases all features, tools, and current applications of `spelunker`. Some applications includes overlaying JWST engineering telemetry
with fitted spatial Gaussian parameters like full width half maximum (FWHM), spatial animations of guidestar frames, and optimizing photometry with pixel-level
decorrelation (PLD). 

The pipline is currently under development on `GitHub <https://github.com/GalagaBits/JWST-FGS-Spelunker>`_. For any bugs or requests, send us an email or `open an issue on GitHub <https://github.com/GalagaBits/JWST-FGS-Spelunker/issues>`_.

.. image:: https://img.shields.io/badge/GitHub-GalagaBits%2FJWST_FGS_Spelunker-blue
   :alt: Static Badge
   :target: https://github.com/GalagaBits/JWST-FGS-Spelunker

.. image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
    :target: https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/LICENSE

.. raw:: html
   <br>

User Guide
==========

For a quick start on how to use spelunker, go to the :doc:`Quickstart guide <simple_quickstart>`.

    * :doc:`Installation <installation>`
    * :doc:`API <api>`

.. toctree::
   :maxdepth: 3
   :hidden:

   Index <self>
   Installation <installation>
   Quickstart <simple_quickstart>
   Exploring Spelunker tools and features <user/fgs-spelunker_quickstart>
   Time Series Observations (TSO) <timeseries_observations>
   Engineering Telemetry <engineering_telemetry>
   API <api>


Tutorials
=========

    * `Pixel Centroid Mnemonics <pixel_centroid_mnemonics>``
    * Using ``spelunker`` to study :doc:`JWST Time Series Observations (TSOs) <user/fgs-spelunker-and-tsos>`.

License & Attribution
=====================

Copyright (c) 2023 Derod Deal & Nestor Espinoza under the `MIT License <https://github.com/GalagaBits/JWST-FGS-Spelunker/blob/main/LICENSE>`_.