Exploring Spelunker tools and features
========================

--------------

**Authors:** Derod Deal (dealderod@ufl.edu), Néstor Espinoza
(nespinoza@stsci.edu) \| **Last update:** July 24, 2023

**Program version:** ver. 1.1.9

The JWST telescope carries four different instruments: NIRCam, NIRSpec,
MIRI and FGS/NIRISS — the latter containing the Fine Guidance Sensor
(FGS). FGS Spelunker is a package designed to conveniently analyze
guidestar data. In this notebook, we cover the following main functions
of this package.

1. `Getting started <#getting-started>`__

   -  Installation
   -  Using ``Spelunker``

2. `Downloading data <#downloading-data>`__
3. Spatial fitting guide stars

   -  `Gaussian fitting <#gaussian-fitting>`__
   -  `Quick fitting <#quick-fitting>`__

4. `Plotting parameters <#plotting-parameters>`__
5. Periodograms

   -  `Creating a periodogram <#periodograms>`__

6. `Mnemonics <#mnemonics>`__
7. `Animations <#animations>`__
8. `Getting tables <#getting-tables>`__

Getting started
---------------

To get started with FGS Spelunker, first call ``spelunker.load`` into a
variable while setting a given Program ID.

.. code:: python

   import spelunker
   spk = spelunker.load(pid=1534)

Calling load without the pid parameter ``spelunker.load()`` will
initialize ``spelunker`` without downloading any of the files. This is
useful if you already have guidestar FITS available locally.

Downloading data
----------------

To load Spelunker with a given Program ID for JWST, simply call
``spelunker.load`` with the Program ID ``pid`` as a parameter. This will
create a directory called ``spelunker_results``, which is where the FITS
files from a selected Program ID and other data will be downloaded and
saved. You can define your own directory by using ``dir=``.

   The Program IDs that can be loaded are limited to programs without an
   exclusive access period or are available publicly
   (https://jwst-docs.stsci.edu/accessing-jwst-data/jwst-data-retrieval/data-access-policy#DataAccessPolicy-Exclusiveaccessperiod).

.. code:: python

    import spelunker
    spk = spelunker.load(dir='/Users/ddeal/JWST-Treasure-Chest/', pid='1534')


.. parsed-literal::

    Current working directory for spelunker: /Users/ddeal/JWST-Treasure-Chest/spelunker_outputs
    
    Connecting with astroquery...
    INFO: Found cached file ./mastDownload/JWST/jw01534002002_04101_00001_guider2/jw01534002002_gs-fg_2022338021919_cal.fits with expected size 10428480. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534001004_03101_00001_guider1/jw01534001004_gs-fg_2022340010755_cal.fits with expected size 8766720. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534002004_03101_00001_guider2/jw01534002004_gs-fg_2022338025056_cal.fits with expected size 8769600. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534001002_03101_00001_guider1/jw01534001002_gs-fg_2022340003651_cal.fits with expected size 8772480. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534002003_03101_00001_guider2/jw01534002003_gs-fg_2022338023521_cal.fits with expected size 8772480. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534001003_03101_00001_guider1/jw01534001003_gs-fg_2022340005224_cal.fits with expected size 8772480. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534002001_05101_00002_guider2/jw01534002001_gs-fg_2022338014704_cal.fits with expected size 10941120. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534002001_05101_00002_guider2/jw01534002001_gs-fg_2022338015941_cal.fits with expected size 7830720. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534001001_03101_00002_guider1/jw01534001001_gs-fg_2022340000825_cal.fits with expected size 9388800. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534001001_03101_00002_guider1/jw01534001001_gs-fg_2022340002102_cal.fits with expected size 7827840. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534004004_03101_00001_guider2/jw01534004004_gs-fg_2023123213436_cal.fits with expected size 8769600. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534004003_03101_00001_guider2/jw01534004003_gs-fg_2023123211905_cal.fits with expected size 8766720. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534004002_03101_00001_guider2/jw01534004002_gs-fg_2023123210335_cal.fits with expected size 8766720. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534004001_03101_00002_guider2/jw01534004001_gs-fg_2023123203053_cal.fits with expected size 12974400. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534004001_03101_00002_guider2/jw01534004001_gs-fg_2023123204330_cal.fits with expected size 7827840. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534003001_03101_00002_guider1/jw01534003001_gs-fg_2023125174543_cal.fits with expected size 9809280. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534003001_03101_00002_guider1/jw01534003001_gs-fg_2023125175812_cal.fits with expected size 7793280. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534003002_02101_00001_guider1/jw01534003002_gs-fg_2023125181351_cal.fits with expected size 8337600. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534003003_02101_00001_guider1/jw01534003003_gs-fg_2023125182911_cal.fits with expected size 8337600. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534003004_02101_00001_guider1/jw01534003004_gs-fg_2023125185519_cal.fits with expected size 8337600. [astroquery.query]


To download the data after initialization, use ``spk.download()`` with
given proposal id with the optional parameters observation number
``obs_num`` and visit number ``visit``. You can also set the calibration
level ``calib_level``. This information are required to use
``astroquery.mast`` to search and download the necessary files. The
download function will download the selected files in the given
directory and create a 2D array of the guidestar data as well as an
array of time and a flux timeseries. The same parameters work with
``spelunker.load``.

.. code:: python

    spk2 = spelunker.load(pid=1534, obs_num='2', visit='1', calib_level=2)
    spk2.download(1534, obs_num='2', visit='2', calib_level=2) # This overwrites the object data in spk2 with data from the input parameters


.. parsed-literal::

    Current working directory for spelunker: /Users/ddeal/JWST-Treasure-Chest/spelunker_outputs
    
    Connecting with astroquery...


.. parsed-literal::

    2023-08-02 21:11:34,101 - stpipe - INFO - Found cached file ./mastDownload/JWST/jw01534002001_05101_00002_guider2/jw01534002001_gs-fg_2022338014704_cal.fits with expected size 10941120.
    2023-08-02 21:11:34,195 - stpipe - INFO - Found cached file ./mastDownload/JWST/jw01534002001_05101_00002_guider2/jw01534002001_gs-fg_2022338015941_cal.fits with expected size 7830720.


.. parsed-literal::

    INFO: Found cached file ./mastDownload/JWST/jw01534002001_05101_00002_guider2/jw01534002001_gs-fg_2022338014704_cal.fits with expected size 10941120. [astroquery.query]
    INFO: Found cached file ./mastDownload/JWST/jw01534002001_05101_00002_guider2/jw01534002001_gs-fg_2022338015941_cal.fits with expected size 7830720. [astroquery.query]
    Connecting with astroquery...


.. parsed-literal::

    2023-08-02 21:11:41,186 - stpipe - INFO - Found cached file ./mastDownload/JWST/jw01534002002_04101_00001_guider2/jw01534002002_gs-fg_2022338021919_cal.fits with expected size 10428480.


.. parsed-literal::

    INFO: Found cached file ./mastDownload/JWST/jw01534002002_04101_00001_guider2/jw01534002002_gs-fg_2022338021919_cal.fits with expected size 10428480. [astroquery.query]


After we downloaded our data, we can access preprocessed spatial, time,
and flux arrays for all FITS files images under the specified Program
ID. Use the attributes ``spk.fg_array``, ``spk.fg_time``, and
``spk.fg_flux`` to access the arrays.

.. code:: python

    spk2.fg_array.shape, spk2.fg_time.shape, spk2.fg_flux.shape




.. parsed-literal::

    ((10240, 8, 8), (10240,), (10240,))



Previously downloaded FITS files in a given directory will not be
re-downloaded. If there are multiple files downloaded for the given
parameter, ``spk.download`` will automatically stitch the data from the
files into an array based on the date and time for each file, along with
the time and flux arrays.

FGS Spelunker can also handle guidestar FITS already stored locally
by using:

.. code:: python

   spk2 = spelunker.load()
   spk2.readfile(pid=1534)


Any files under the initialized directory and specified Program ID, observation number, and visit number will be loaded into `spk2`.

    `spelunker.load.readfile` now has access to the same attributes as `spelunker.load.download`. So, using `spk.object_properties` and `spk.fg_table`
    will work.

Spatial fitting guide stars
---------------------------

After downloading the data, we can perform spatial fitting gaussians to
each frame in a guidestar timeseries. This uses parallel processing
through ``ray`` to speed up the process. We can also perform quick fits
to speed through a given timeseries, though this method is a lot less
accurate in the fitting.

Gaussian fitting
~~~~~~~~~~~~~~~~

The downloaded data comes as a spatial timeseries of a selected
guidestar. To measure the centriods and PSF width of each frame, we need
to apply fitting. We will use Gaussian spatial fitting to measure x and
y pixel coordinates, x and y standard deviations, theta, and the
offset. To perform spatial gaussian fitting, use ``gauss2d_fit`` with guidestar arrays (the
timeseries needs to be in an 8 by 8 array, which should be the same for
all guidestar fine guidence products).

.. code:: python

   spk.gauss2d_fit() # ncpus sets the number of cpu cores your computer has. Defaults to 4 cores.

.. code:: python

    # We are going to limit the amount of frames that we input into gauss2d_fit and other methods
    # since the gauss2d_fit can take a few houts for very large arrays.
    spk.fg_array = spk.fg_array[0:10000]
    spk.fg_flux = spk.fg_flux[0:10000]
    spk.fg_time = spk.fg_time[0:10000]

.. code:: python

    table_gauss_fit = spk.gauss2d_fit(ncpus=6) 


.. parsed-literal::

    2023-08-02 21:12:50,384	INFO worker.py:1636 -- Started a local Ray instance.


The ``gauss2d_fit`` function outputs an astropy table, which can bee
accessed with the ``spk.gaussfit_results`` attribute. If ``gauss2d_fit``
fails to fit a frame, it will return ``nan`` for that frame.

.. code:: python

    spk.gaussfit_results




.. raw:: html

    <div><i>Table length=10000</i>
    <table id="table4415257088" class="table-striped table-bordered table-condensed">
    <thead><tr><th>amplitude</th><th>x_mean</th><th>y_mean</th><th>x_stddev</th><th>y_stddev</th><th>theta</th><th>offset</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>280706.15465765796</td><td>3.1774294356249997</td><td>2.7465302838135206</td><td>0.6350976070387301</td><td>0.614009020575321</td><td>-1.9103595130650228</td><td>3023.106318279726</td></tr>
    <tr><td>280706.15465765796</td><td>3.1774294356249997</td><td>2.7465302838135206</td><td>0.6350976070387301</td><td>0.614009020575321</td><td>-1.9103595130650228</td><td>3023.106318279726</td></tr>
    <tr><td>280963.5540504813</td><td>3.177604462333186</td><td>2.7483597462452547</td><td>0.6306454543965104</td><td>0.6193386849707871</td><td>-2.057972902746876</td><td>3149.3240730860866</td></tr>
    <tr><td>280963.5540504813</td><td>3.177604462333186</td><td>2.7483597462452547</td><td>0.6306454543965104</td><td>0.6193386849707871</td><td>-2.057972902746876</td><td>3149.3240730860866</td></tr>
    <tr><td>282706.5250312361</td><td>3.1764861837068716</td><td>2.749817871515913</td><td>0.6334273199822001</td><td>0.6145497343103167</td><td>-1.9504191092501943</td><td>3053.0948632606123</td></tr>
    <tr><td>282706.5250312361</td><td>3.1764861837068716</td><td>2.749817871515913</td><td>0.6334273199822001</td><td>0.6145497343103167</td><td>-1.9504191092501943</td><td>3053.0948632606123</td></tr>
    <tr><td>277126.33630266984</td><td>3.1748827601728564</td><td>2.7477495874396674</td><td>0.6189797899040209</td><td>0.6340116557887706</td><td>-3.48449959258196</td><td>3105.682301707251</td></tr>
    <tr><td>277126.33630266984</td><td>3.1748827601728564</td><td>2.7477495874396674</td><td>0.6189797899040209</td><td>0.6340116557887706</td><td>-3.48449959258196</td><td>3105.682301707251</td></tr>
    <tr><td>280742.3344982786</td><td>3.1719030737999923</td><td>2.756636337651271</td><td>0.6154040193075433</td><td>0.6363143600933248</td><td>-3.570644823307217</td><td>3017.796074602062</td></tr>
    <tr><td>280742.3344982786</td><td>3.1719030737999923</td><td>2.756636337651271</td><td>0.6154040193075433</td><td>0.6363143600933248</td><td>-3.570644823307217</td><td>3017.796074602062</td></tr>
    <tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
    <tr><td>288936.6587997144</td><td>3.1514848995974614</td><td>2.816421337922728</td><td>0.6078414078127158</td><td>0.6255153338398373</td><td>-0.724102219944298</td><td>3159.747623016102</td></tr>
    <tr><td>288936.6587997144</td><td>3.1514848995974614</td><td>2.816421337922728</td><td>0.6078414078127158</td><td>0.6255153338398373</td><td>-0.724102219944298</td><td>3159.747623016102</td></tr>
    <tr><td>287608.5204882826</td><td>3.148081209519121</td><td>2.8097574913336154</td><td>0.6092268378675755</td><td>0.6288855374510539</td><td>-0.6364418904422164</td><td>3098.4078599410695</td></tr>
    <tr><td>287608.5204882826</td><td>3.148081209519121</td><td>2.8097574913336154</td><td>0.6092268378675755</td><td>0.6288855374510539</td><td>-0.6364418904422164</td><td>3098.4078599410695</td></tr>
    <tr><td>286304.0727626729</td><td>3.1471623118694176</td><td>2.8102083208968813</td><td>0.6085355521172578</td><td>0.6298236704220975</td><td>-0.5591615297330863</td><td>3183.299010073181</td></tr>
    <tr><td>286304.0727626729</td><td>3.1471623118694176</td><td>2.8102083208968813</td><td>0.6085355521172578</td><td>0.6298236704220975</td><td>-0.5591615297330863</td><td>3183.299010073181</td></tr>
    <tr><td>284871.6486689821</td><td>3.1499465078006614</td><td>2.8072167275653706</td><td>0.6111915236092285</td><td>0.6277931861719188</td><td>-0.7047253049826113</td><td>3261.2487765038327</td></tr>
    <tr><td>284871.6486689821</td><td>3.1499465078006614</td><td>2.8072167275653706</td><td>0.6111915236092285</td><td>0.6277931861719188</td><td>-0.7047253049826113</td><td>3261.2487765038327</td></tr>
    <tr><td>288107.09702730743</td><td>3.14940434535617</td><td>2.807916552216667</td><td>0.6081505348286508</td><td>0.6295003348022744</td><td>-0.6030650650578055</td><td>3197.4098077599647</td></tr>
    <tr><td>288107.09702730743</td><td>3.14940434535617</td><td>2.807916552216667</td><td>0.6081505348286508</td><td>0.6295003348022744</td><td>-0.6030650650578055</td><td>3197.4098077599647</td></tr>
    </table></div>



----------


Quick fitting
~~~~~~~~~~~~~

There are some situations where you need to quickly obtain rough
statistics of changes in guidestar products overtime. Quick fitting fits
the x and y pixel locations and standard deviations as an astropy table
using centroid and variance calculations. To perform quick fitting, run
``quick_fit`` with an appropriate array.

.. code:: python

    table_quick_fit = spk.quick_fit()

.. code:: python

    spk.quickfit_results




.. raw:: html

    <div><i>Table length=10000</i>
    <table id="table4415251568" class="table-striped table-bordered table-condensed">
    <thead><tr><th>amplitude</th><th>x_mean</th><th>y_mean</th><th>x_stddev</th><th>y_stddev</th><th>theta</th><th>offset</th></tr></thead>
    <thead><tr><th>float32</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>int64</th><th>int64</th></tr></thead>
    <tr><td>254451.56</td><td>3.240314850861845</td><td>2.8033942297495758</td><td>1.74462175414244</td><td>1.8158228238188503</td><td>0</td><td>0</td></tr>
    <tr><td>254451.56</td><td>3.240314850861845</td><td>2.8033942297495758</td><td>1.74462175414244</td><td>1.8158228238188503</td><td>0</td><td>0</td></tr>
    <tr><td>255055.25</td><td>3.3206004778017384</td><td>2.8434574303565463</td><td>1.8543257785557397</td><td>1.8293394846671764</td><td>0</td><td>0</td></tr>
    <tr><td>255055.25</td><td>3.3206004778017384</td><td>2.8434574303565463</td><td>1.8543257785557397</td><td>1.8293394846671764</td><td>0</td><td>0</td></tr>
    <tr><td>256947.42</td><td>3.3505845162736376</td><td>2.925690858450849</td><td>1.8077292667969422</td><td>1.8943471255043283</td><td>0</td><td>0</td></tr>
    <tr><td>256947.42</td><td>3.3505845162736376</td><td>2.925690858450849</td><td>1.8077292667969422</td><td>1.8943471255043283</td><td>0</td><td>0</td></tr>
    <tr><td>251888.12</td><td>3.3039389301600726</td><td>2.886233231270987</td><td>1.854677926018813</td><td>1.8433178905598915</td><td>0</td><td>0</td></tr>
    <tr><td>251888.12</td><td>3.3039389301600726</td><td>2.886233231270987</td><td>1.854677926018813</td><td>1.8433178905598915</td><td>0</td><td>0</td></tr>
    <tr><td>257109.62</td><td>3.2835164773971806</td><td>2.774318082677534</td><td>1.837107063709473</td><td>1.7647732623026264</td><td>0</td><td>0</td></tr>
    <tr><td>257109.62</td><td>3.2835164773971806</td><td>2.774318082677534</td><td>1.837107063709473</td><td>1.7647732623026264</td><td>0</td><td>0</td></tr>
    <tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
    <tr><td>273886.84</td><td>3.307248070570433</td><td>2.9459581137888096</td><td>1.8638542966133642</td><td>1.8248573282234368</td><td>0</td><td>0</td></tr>
    <tr><td>273886.84</td><td>3.307248070570433</td><td>2.9459581137888096</td><td>1.8638542966133642</td><td>1.8248573282234368</td><td>0</td><td>0</td></tr>
    <tr><td>272548.8</td><td>3.303274024993382</td><td>2.888558147490168</td><td>1.8282836367085207</td><td>1.7580760556837993</td><td>0</td><td>0</td></tr>
    <tr><td>272548.8</td><td>3.303274024993382</td><td>2.888558147490168</td><td>1.8282836367085207</td><td>1.7580760556837993</td><td>0</td><td>0</td></tr>
    <tr><td>271490.1</td><td>3.228820447972362</td><td>3.055912219282716</td><td>1.8189049613644188</td><td>1.8755066513378191</td><td>0</td><td>0</td></tr>
    <tr><td>271490.1</td><td>3.228820447972362</td><td>3.055912219282716</td><td>1.8189049613644188</td><td>1.8755066513378191</td><td>0</td><td>0</td></tr>
    <tr><td>269606.9</td><td>3.328221486759065</td><td>2.963716959723631</td><td>1.8706223659386954</td><td>1.8586654374692335</td><td>0</td><td>0</td></tr>
    <tr><td>269606.9</td><td>3.328221486759065</td><td>2.963716959723631</td><td>1.8706223659386954</td><td>1.8586654374692335</td><td>0</td><td>0</td></tr>
    <tr><td>272629.9</td><td>3.304655431094987</td><td>2.9615702404863526</td><td>1.873261996709939</td><td>1.9288479581727678</td><td>0</td><td>0</td></tr>
    <tr><td>272629.9</td><td>3.304655431094987</td><td>2.9615702404863526</td><td>1.873261996709939</td><td>1.9288479581727678</td><td>0</td><td>0</td></tr>
    </table></div>



Plotting parameters
-------------------

We can plot a timeseries of a given parameter or flux from guidestars.
The method ``timeseries_binned_plot`` will generate a matplotlib axes
object of a given timeseries.

.. code:: python

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize = (12,4), dpi=200)
    
    ax = spk.timeseries_binned_plot()



.. image:: fgs-spelunker_quickstart_files/fgs-spelunker_quickstart_35_0.png
   :scale: 50%

Within guidestar data, changes in the PSF can impact the observed flux
of the star. Certain events might see changes in all fitted parameters.
In this case, subplots of each parameter will provide more information
to the user about the event, giving them the change of guidestar
position, brightness, and FWHM overtime.

.. code:: python

    ax = spk.timeseries_list_plot()



.. image:: fgs-spelunker_quickstart_files/fgs-spelunker_quickstart_37_0.png
   :scale: 50%

Periodograms
------------

FGS Spelunker comes with various tools to visualize and explore
guidestar data. Periodograms are useful for guidestar products to detect
periodicities not only within flux timeseries, but also within
centroids, FWHM, theta, and offset. From a selected fitting method, we
can use the table output to apply Lomb-Scargle periodograms to our
parameters.

``periodogram``
~~~~~~~~~~~~~~~

To obtain the power and frequencies of Lomb-Scargle periodograms for
each fitted parameter, use ``periodogram``. The periodograms for each
given parameter from a fit can be conveniently plotted in a single
figure with the same method.

.. code:: python

    ax = spk.periodogram()



.. image:: fgs-spelunker_quickstart_files/fgs-spelunker_quickstart_41_0.png
   :scale: 50%

To get the frequency and power for each fitted parameter, use
``spk.pgram_{parameter}``. > Available parameters: > -
``spk.pgram_amplitude`` > - ``spk.pgram_x_mean`` > -
``spk.pgram_y_mean`` > - ``spk.pgram_x_stddev`` > -
``spk.pgram_y_stddev`` > - ``spk.pgram_theta`` > - ``spk.pgram_offset``

.. code:: python

    freq = spk.pgram_x_mean[0] # periodogram frequency
    power = spk.pgram_x_mean[1] # periodogram power
    
    freq[0], power[0]




.. parsed-literal::

    (0.0003127661546504965, 0.005397779092056495)



Mnemonics
---------

When observing the timeseries of the guidestar, there might be technical
events from the JWST that causes changes in obtained data. For example,
high gain antenna or filter changes in NIRCAM can cause noticeable
changes in flux or other guidestar properties. We can overlay these
events onto fitted parameters using ``mnemonics`` and
``mnemonics_plot``. You will need a MAST API token to use ``mnemonics``,
as well as the ``jwstuser`` package: 
- https://auth.mast.stsci.edu/docs/
(MAST API TOKEN) 
- https://github.com/spacetelescope/jwstuser/tree/main
(jwstuser)

   Current supported mnemonics: 
   - *SA_ZHGAUPST* (high-gain antenna),
   - *INIS_FWMTRCURR* (NIRISS Filter Wheel Motor Current).

There are thousands of different mnemonics to explore on https://mast.stsci.edu/viz/api/v0.1/info/mnemonics. Use `spk.mnemonics` to try
the mnemonics you are interested in comparing with any JWST data, not just guidestars from the FGS.

.. code:: python

    spk2 = spelunker.load('/Users/ddeal/JWST-Treasure-Chest/', pid=1534)



.. code:: python

    spk2.mast_api_token = 'enter_mast_token_id_here' # input mast_api token here!
    
    fig, ax = plt.subplots(figsize=(12,4),dpi=200)
    
    ax = spk2.mnemonics_local('GUIDESTAR') # plots when the JWST tracks onto a new guidestars as a vertical line
    ax = spk2.mnemonics('SA_ZHGAUPST', 60067.84, 60067.9) # plots the start and end of high gain antenna movement
    
    ax.plot(spk2.fg_time, spk2.fg_flux)
    plt.legend(loc=3)
    
    plt.xlim(60067.84, 60067.9)






.. image:: fgs-spelunker_quickstart_files/fgs-spelunker_quickstart_47_1.png
   :scale: 50%

If you do have a MAST API token, you will have access to any program under that token.

If you do not have access to a MAST API token, you can only download and use publicly available Program IDs. 
However, with the readfile function, you can use fine guidance files you already have downloaded locally, 
and with the current version, with no drawbacks. 

Animations
----------

Spatial data of guidestar imaging can bring essential information about
how the point spread function changes overtime. Animations of the
spatial timeseries are convenient and helpful methods to analyze
guidestar data. To get a side by side comparison of the evolution of a
spatial timeseries and a parameter, use
``flux_spatial_timelapse_animation``.

   You may have to install ``ffmpeg`` on your computer to get ``mp4``
   formats.

.. code:: python

    plt.plot(spk2.fg_flux[2600:3100])




.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x1c16b7550>]




.. image:: fgs-spelunker_quickstart_files/fgs-spelunker_quickstart_50_1.png

.. code:: python

    spk2.flux_spatial_timelapse_animation(start=2600,stop=3100,) # to save an animation with a filename, use *filename=*. Defaults to movie.gif


.. parsed-literal::

    2023-08-02 21:19:50,803	INFO worker.py:1636 -- Started a local Ray instance.



.. image:: fgs-spelunker_quickstart_files/movie.gif
   :scale: 40%

Getting tables
--------------

After downloading a selected proposal id with ``download``, we can
easily output metadata about each downloaded file, including extracted
data from the filename including ``visit_group``,
``parallel_sequence_id``, and ``exposure_number``. The guide star used
in each file is also included, as well as filter magnitudes and other
stellar properties.

.. code:: python

    spk.fg_table # We can simply call this attribute after using spk.download() to obtain our table!




We can obtain a neat DataFrame of each tracked guidestar, which gives us
information such as the intergation start times and galactic
coordinates.

.. code:: python

    spk.object_properties




.. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>guidestar_catalog_id</th>
          <th>gaiadr1ID</th>
          <th>gaiadr2ID</th>
          <th>int_start</th>
          <th>int_stop</th>
          <th>ra</th>
          <th>dec</th>
          <th>Jmag</th>
          <th>Hmag</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>0</th>
          <td>S1HP079555</td>
          <td>4658077781377287680</td>
          <td>4658077781376437888</td>
          <td>59917.066396</td>
          <td>59917.074354</td>
          <td>80.837584</td>
          <td>-69.541124</td>
          <td>13.659</td>
          <td>12.898</td>
        </tr>
        <tr>
          <th>1</th>
          <td>S1HP080554</td>
          <td>4658077991763987712</td>
          <td>4658077991799023616</td>
          <td>59917.089163</td>
          <td>59917.096759</td>
          <td>80.806837</td>
          <td>-69.530972</td>
          <td>15.001</td>
          <td>14.282</td>
        </tr>
        <tr>
          <th>2</th>
          <td>S1HP078573</td>
          <td>4657983910572904320</td>
          <td>4657983910572904320</td>
          <td>59917.112547</td>
          <td>59917.118705</td>
          <td>80.807043</td>
          <td>-69.553474</td>
          <td>13.839</td>
          <td>13.078</td>
        </tr>
        <tr>
          <th>3</th>
          <td>S1HP079590</td>
          <td>4657986831103727872</td>
          <td>4657986835382982016</td>
          <td>59918.999015</td>
          <td>59919.005848</td>
          <td>80.510790</td>
          <td>-69.545479</td>
          <td>15.410</td>
          <td>14.839</td>
        </tr>
        <tr>
          <th>4</th>
          <td>S1HP079769</td>
          <td>4657986831078120832</td>
          <td>4657986835433225728</td>
          <td>59919.019436</td>
          <td>59919.025598</td>
          <td>80.518235</td>
          <td>-69.543415</td>
          <td>15.231</td>
          <td>14.341</td>
        </tr>
        <tr>
          <th>5</th>
          <td>S1HP078292</td>
          <td>4657986796681532672</td>
          <td>4657986801073794432</td>
          <td>59919.041018</td>
          <td>59919.047165</td>
          <td>80.519564</td>
          <td>-69.558464</td>
          <td>12.804</td>
          <td>11.883</td>
        </tr>
        <tr>
          <th>6</th>
          <td>S1HP077850</td>
          <td>4657986762384054144</td>
          <td>4657986766713867264</td>
          <td>60067.871344</td>
          <td>60067.877490</td>
          <td>80.573531</td>
          <td>-69.562862</td>
          <td>12.957</td>
          <td>12.227</td>
        </tr>
        <tr>
          <th>7</th>
          <td>S1HP197501</td>
          <td>4657986865463528832</td>
          <td>4657986869793061376</td>
          <td>60067.882117</td>
          <td>60067.888264</td>
          <td>80.571447</td>
          <td>-69.551750</td>
          <td>13.063</td>
          <td>12.168</td>
        </tr>
        <tr>
          <th>8</th>
          <td>S1HP773376</td>
          <td></td>
          <td>4658078124973829632</td>
          <td>60069.733171</td>
          <td>60069.740086</td>
          <td>80.794522</td>
          <td>-69.504084</td>
          <td>13.426</td>
          <td>12.654</td>
        </tr>
        <tr>
          <th>9</th>
          <td>S1HP081366</td>
          <td>4658078056254368128</td>
          <td>4658078056254368128</td>
          <td>60069.753592</td>
          <td>60069.759620</td>
          <td>80.758291</td>
          <td>-69.524143</td>
          <td>12.765</td>
          <td>11.899</td>
        </tr>
        <tr>
          <th>10</th>
          <td>S1HP082164</td>
          <td>4658077953064455552</td>
          <td>4658077957439332608</td>
          <td>60069.764246</td>
          <td>60069.770278</td>
          <td>80.865554</td>
          <td>-69.514107</td>
          <td>12.753</td>
          <td>11.871</td>
        </tr>
      </tbody>
    </table>
    </div>


