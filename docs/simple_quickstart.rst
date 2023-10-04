Quickstart
==========

The JWST telescope carries four different instruments: NIRCam, NIRSpec, MIRI and FGS/NIRISS --- the latter containing the 
Fine Guidance Sensor (FGS). FGS Spelunker is a package designed to conveniently analyze guidestar data from.

Installation
------------

To install ``spelunker``, use ``pip install``.

.. code:: bash

   pip install spelunker


For in-depth installation steps, visit :doc:`the installation documentation <installation>`.

Using the library
-----------------

Get started with ``spelunker`` with only two lines of code.

.. code:: python

   import spelunker

   spk = spelunker.load(pid=1534)

This will download guidestar data for Program ID 1534; the ``spk``
object itself can then be used to explore this guidestar data! For
example, let us make a plot of the guidestar time-series for the first
minutes of this PID:

.. code:: python

   import matplotlib.pyplot as plt

   # Convert times from MJD to minutes:
   plt.plot( ( spk.fg_time - spk.fg_time[0] ) * 24 * 60, spk.fg_flux )

   plt.xlim(0,10)
   plt.xlabel('Time from start (minutes)')
   plt.ylabel('Counts')

.. image:: simple_guidestar_files/timeseries.png
   :scale: 60%

.. raw:: html

   <p align="center">

.. raw:: html

   </p>

(See below on more information that can be extracted, including fitting
2D gaussians to each FGS integration!). We can even make a plot of the
tracked guidestars within this Program ID:

.. code:: python

   spk.guidestar_plot()


.. image:: simple_guidestar_files/guidestar_positions.png
   :scale: 60%

.. raw:: html

   <p align="center">

.. raw:: html

   </p>

Mnemonics from JWST technical events can be overplotted on any
timeseries, such as high-gain antenna (HGA) movement or to identify if
the FGS tracks a new guidestar if the `jwstuser package is also
installed <https://github.com/spacetelescope/jwstuser/>`_.

.. code:: python

   import matplotlib.pyplot as plt

   spk.mast_api_token = 'insert a token from auth.MAST here'

   fig, ax = plt.subplots(figsize=(12,4),dpi=200)

   ax = spk.mnemonics_local('GUIDESTAR')
   ax = spk.mnemonics('SA_ZHGAUPST', 60067.84, 60067.9) 
   ax.plot(spk.fg_time, spk.fg_flux)
   plt.legend(loc=3)
   plt.xlim(60067.84, 60067.9)
   plt.show()

.. image:: simple_guidestar_files/mnemonics.png
   :scale: 60%