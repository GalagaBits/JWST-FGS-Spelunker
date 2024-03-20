Installation
############

``spelunker`` can be installed using ``pip``.

.. code:: bash

   pip install spelunker

The following packages are required:

* numpy
* scipy
* matplotlib
* astropy
* jwst
* astroquery
* pandas
* ray (required for parallel processing of Gaussian fitting)
* astroplan (optional; for use with :func:`spelunker.load.guidestar_plot`)

For Python 3.12, ``ray`` fails to install with Spelunker when using ``pip install spelunker``. We recommend using Python 3.10 when installing ``ray`` and ``spelunker``.

Additonally, you can refer to the `GitHub repository <https://github.com/GalagaBits/JWST-FGS-Spelunker/>`_ for the source code.

JWST User Toolkit (`jwstuser`)
******************************

To use the `mnemonics` engineering telemetry features of Spelunker, the python toolkit ``jwstuser`` will need to be installed as well.
Follow the instructions on `https://github.com/spacetelescope/jwstuser <https://github.com/spacetelescope/jwstuser>`_ to install ``jwstuser``.