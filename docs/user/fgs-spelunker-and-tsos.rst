Using ``spelunker`` to study JWST Time Series Observations (TSOs)
=================================================================

--------------

**Authors:** Néstor Espinoza (nespinoza@stsci.edu), Derod Deal
(dealderod@ufl.edu) \| **Last update:** July 25, 2023

**Program version:** ver. 0.5.3

JWST Time Series Observations (TSOs) are multi-integration exposures
typically targetted at exploring time-varying phenomena: from
exoplanetary transits to accreeting material in distant objects.
Guidestar data such as the one ``spelunker`` can query can become very
helpful at exploring this data; this tutorial provides an introduction
on how to use the ``spelunker`` products to analyze it.

1. The case of HAT-P-14b NIRISS/SOSS observations
=================================================

The first dataset we will be analyzing below comes from an exoplanetary
transit obtained by `Program ID
1541 <https://www.stsci.edu/jwst/science-execution/program-information?id=1541>`__
(PI: Espinoza). This dataset was already introduced in `Albert et
al. (2023) <https://arxiv.org/abs/2306.04572>`__, and consisted of a
transit observation of the exoplanet HAT-P-14~b, which was used to test
the sensitivity and stability of the NIRISS/SOSS instrument during
commissioning.

To start the analysis, let’s load some libraries:

.. code:: ipython3

    import numpy as np
    import matplotlib.pyplot as plt
    
    import spelunker

1.1 Exploring the transit event of HAT-P-14b
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let’s first load and plot the NIRISS/SOSS HAT-P-14~b data:

.. code:: ipython3

    t, f, ferr = np.loadtxt('hp14_lightcurve.dat', unpack = True, usecols = (0,1,2))

.. code:: ipython3

    tstart = t[0]
    time_since_start = (t-tstart)*24
    
    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = '.', 
                             ms = 1, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.title('Exoplanet transit of HAT-P-14 b with NIRISS/SOSS', fontsize = 18)
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), np.max(time_since_start))




.. parsed-literal::

    (0.0, 6.09976451843977)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_4_1.png


All right. This time-series shows a nice transit event from about hour
2.5 since exposure start, all the way until about hour 5. Aside from the
transit event, however, there is an evident oscillation in the data,
which is evident if we do a zoom to the first three hours of data:

.. code:: ipython3

    tstart = t[0]
    time_since_start = (t-tstart)*24
    
    plt.figure(figsize=(10,4))
    
    plt.plot(time_since_start, (f-1)*1e6, color = 'black')
    plt.errorbar(time_since_start, (f-1)*1e6, ferr*1e6, fmt = '.', 
                             ms = 1, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.title('Exoplanet transit of HAT-P-14 b with NIRISS/SOSS (zoom)', fontsize = 18)
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux - 1 (ppm)', fontsize = 18)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), 3)
    plt.ylim(-500, 500)




.. parsed-literal::

    (-500.0, 500.0)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_6_1.png


The light curve has at least two oscillation patterns. One is a
long-term one, on which the light curve seems to rise at about hour 0.5
after start, then go down until about hour 1, and then oscillate up
again at about hour 3. The other is a short-frequency oscillation, with
a period of about ~5 minutes. The amplitude of those oscillations is
small — around 200 ppm.

Big question is: are those oscillations really happening on HAT-P-14
(the star)? Or is this an instrumental effect? Let’s now explore the
guidestar data to find this out.

1.2 Exploring the guidestar data of PID 1541
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let’s use ``spelunker`` to load the guidestar data for this Program
ID/observation/visit:

.. code:: ipython3

    spk = spelunker.load(pid=1541, obs_num='1', visit='1', save=True)


.. parsed-literal::

    Current working directory for spelunker: /Users/nespinoza/github/JWST-FGS-Spelunker/notebooks/spelunker_outputs
    
    Connecting with astroquery...


Let’s check the time-series of the guidestar data:

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    fg_time_since_start = (spk.fg_time + 2400000.5 - tstart) * 24
    
    plt.plot(fg_time_since_start, spk.fg_flux / np.nanmedian( spk.fg_flux ) , color = 'tomato')
    
    plt.title('Guidestar flux for PID 1541, visit 1, observation 1', fontsize = 18)
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    plt.ylim(0.5,1.5)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), np.max(time_since_start))




.. parsed-literal::

    (0.0, 6.09976451843977)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_10_1.png


The raw photometry from the guidestar varies by a lot more than our
target star. Let’s bin the FGS time-series to the same cadence as the
science time-series. To this end, let’s write a function that does this:

.. code:: ipython3

    def bin_fgs_to_science(tscience, tfgs, ffgs):
        """
        This function bins an FGS time-series defined by the times `tfgs` and fluxes `ffgs`, to times `tscience`. 
        The function assumes that (1) `tscience` are times obtained at pseudo-regular intervals (i.e., that times 
        on `tscience` next to each other are similar), and that (2) `tscience` is ordered in chronological order.
        """
    
        nscience = len(tscience)
        binned_fgs = np.zeros( nscience )
        binned_fgs_err = np.zeros( nscience )
        for i in range( nscience ):
    
            if i == 0:
    
                dt = tscience[1] - tscience[0] 
    
            elif i == nscience - 1:
    
                dt = tscience[-1] - tscience[-2]
    
            else:
    
                dt1 = tscience[i] - tscience[i-1]
                dt2 = tscience[i+1] - tscience[i]
                dt = ( dt1 + dt2 ) * 0.5
                
            idx = np.where( np.abs(tscience[i] - tfgs) < 0.5*dt )[0]
            binned_fgs[i] = np.mean( ffgs[idx] )
            binned_fgs_err[i] = np.sqrt( np.var( ffgs[idx] ) ) / np.sqrt( len(idx) )
    
        return binned_fgs, binned_fgs_err

.. code:: ipython3

    fbin, fbinerr = bin_fgs_to_science(time_since_start, fg_time_since_start, spk.fg_flux / np.nanmedian( spk.fg_flux ))

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, (fbin-1)*1e6, fbinerr*1e6, color = 'tomato')
    
    plt.title('(Binned) GS flux for PID 1541, visit 1, observation 1', fontsize = 18)
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux (ppm)', fontsize = 18)
    
    plt.ylim(1.0-0.05,1+0.05)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.ylim(-10000, 10000)
    plt.xlim(np.min(time_since_start), 3)




.. parsed-literal::

    (0.0, 3.0)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_14_1.png


This actually resembles the science time-series quite nicely, although
at a different amplitude. Let’s plot both on the same figure:

.. code:: ipython3

    tstart = t[0]
    time_since_start = (t-tstart)*24
    
    plt.figure(figsize=(10,4))
    
    plt.plot(time_since_start, (f-1)*1e6, color = 'black', label = 'NIRISS/SOSS TSO')
    
    plt.plot(time_since_start, (fbin-1)*1e6*0.05, color = 'tomato', label = r'FGS Guidestar Flux TSO $\times$ 0.05')
    
    plt.legend()
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux - 1 (ppm)', fontsize = 18)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), 3)
    plt.ylim(-500, 500)




.. parsed-literal::

    (-500.0, 500.0)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_16_1.png


Remarkable! The amplitude might need some tweaking, but it seems this
can, indeed, help track some lightcurve variations. Let’s look next at
other features that could be correlated with instrumental systematics.

1.3 More, more! Correlating PSF Guidestar properties to JWST TSOs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``spelunker`` can also fit Gaussians to each of the 2D FG frames, and
extract more precise parameters than the simple “crude” photometry
described above. This takes a while (a few minutes), so we need to be a
bit patient:

.. code:: ipython3

    spk.gauss2d_fit(ncpus=4)


.. parsed-literal::

    2023-07-28 15:06:29,069	ERROR services.py:1207 -- Failed to start the dashboard , return code 1
    2023-07-28 15:06:29,072	ERROR services.py:1232 -- Error should be written to 'dashboard.log' or 'dashboard.err'. We are printing the last 20 lines for you. See 'https://docs.ray.io/en/master/ray-observability/ray-logging.html#logging-directory-structure' to find where the log file is.
    2023-07-28 15:06:29,097	ERROR services.py:1276 -- 
    The last 20 lines of /tmp/ray/session_2023-07-28_15-06-25_792862_10261/logs/dashboard.log (it contains the error message from the dashboard): 
    2023-07-28 15:06:28,974	INFO head.py:226 -- Loading DashboardHeadModule: <class 'ray.dashboard.modules.usage_stats.usage_stats_head.UsageStatsHead'>
    2023-07-28 15:06:28,974	INFO head.py:239 -- Loaded 1 modules. [<ray.dashboard.modules.usage_stats.usage_stats_head.UsageStatsHead object at 0x106db63a0>]
    2023-07-28 15:06:28,974	INFO head.py:331 -- http server disabled.
    2023-07-28 15:06:28,974	ERROR head.py:346 -- Failed to connect to GCS. Please check `gcs_server.out` for more details.
    2023-07-28 15:06:28,987	ERROR dashboard.py:259 -- The dashboard on node rakisduam.stsci.edu failed with the following error:
    Traceback (most recent call last):
      File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/dashboard/dashboard.py", line 248, in <module>
        loop.run_until_complete(dashboard.run())
      File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/asyncio/base_events.py", line 647, in run_until_complete
        return future.result()
      File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/dashboard/dashboard.py", line 75, in run
        await self.dashboard_head.run()
      File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/dashboard/head.py", line 346, in run
        self.gcs_client.internal_kv_put(
      File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
      File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
      File "python/ray/_raylet.pyx", line 2221, in ray._raylet.GcsClient.internal_kv_put
      File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    ray.exceptions.RpcError: failed to connect to all addresses
    2023-07-28 15:06:50,551	INFO worker.py:1621 -- Started a local Ray instance.
    [2m[33m(raylet)[0m [2023-07-28 15:07:29,521 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26746404864; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:07:35,412 E 10444 36357357] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x000000010304c590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x000000010304c77c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x0000000102807df4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x00000001028ca850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x00000001028c7d4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x00000001005910d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m [2023-07-28 15:07:35,706 E 10446 36357362] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x0000000107ae8590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x0000000107ae877c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x00000001072a3df4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x0000000107366850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x0000000107363d4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x00000001051bd0d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:07:36,461 E 10445 36357360] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x000000010733c590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x000000010733c77c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x0000000106af7df4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x0000000106bba850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x0000000106bb7d4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x0000000102c050d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:07:36,811 E 10443 36357355] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x0000000105cc4590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x0000000105cc477c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x000000010547fdf4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x0000000105542850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x000000010553fd4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x0000000100a090d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m [2023-07-28 15:07:39,586 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26743943168; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:07:49,654 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26745270272; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:07:59,734 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26742882304; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:08:09,802 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26741587968; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:08:19,872 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26741252096; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:08:29,937 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26739695616; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:08:40,007 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26741256192; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:08:44,775 E 10406 36355151] (raylet) worker_pool.cc:544: Some workers of the worker process(10493) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:08:44,780 E 10406 36355151] (raylet) worker_pool.cc:544: Some workers of the worker process(10494) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:08:44,781 E 10406 36355151] (raylet) worker_pool.cc:544: Some workers of the worker process(10495) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:08:44,784 E 10406 36355151] (raylet) worker_pool.cc:544: Some workers of the worker process(10496) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:08:50,077 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26741645312; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:09:00,151 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26741493760; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:09:10,217 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26741489664; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:09:20,274 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26739683328; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:09:30,333 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26738573312; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:09:40,398 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26736566272; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:09:44,783 E 10406 36355151] (raylet) worker_pool.cc:544: Some workers of the worker process(10581) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:09:50,468 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26735161344; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[36m(ray_curve_fit pid=10640)[0m /Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/scipy/optimize/_minpack_py.py:881: OptimizeWarning: Covariance of the parameters could not be estimated
    [2m[36m(ray_curve_fit pid=10640)[0m   warnings.warn('Covariance of the parameters could not be estimated',
    [2m[33m(raylet)[0m [2023-07-28 15:10:00,518 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26729492480; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:10:10,566 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26731814912; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:10:20,621 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26731790336; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:10:30,684 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26730823680; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:10:40,732 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26754113536; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:10:50,787 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26752978944; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:11:00,838 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26749263872; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:11:10,892 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26749222912; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:11:20,939 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26748469248; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:11:30,989 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26749116416; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:11:41,044 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26748141568; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:11:51,096 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26747527168; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:12:01,156 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26744139776; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:12:11,200 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26744922112; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:12:21,245 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26742046720; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:12:31,294 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26744434688; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:12:41,351 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26741665792; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:12:51,400 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26740588544; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:13:01,455 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26739531776; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:13:11,508 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26738024448; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:13:21,548 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26735910912; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:13:31,586 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26735661056; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:13:41,615 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26734063616; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:13:51,643 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26736775168; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:14:01,675 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26733940736; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:14:11,700 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26732744704; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:14:21,745 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26731364352; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:14:31,780 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26731507712; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:14:41,796 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26729496576; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:14:51,834 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26730328064; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:15:01,881 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26726350848; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:15:11,903 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26726793216; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:15:21,940 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26726776832; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:15:31,970 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26722263040; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:15:42,000 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26722217984; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:15:52,033 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26724208640; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:16:02,067 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26718584832; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:16:12,098 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26722365440; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:16:22,121 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26721865728; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:16:32,145 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26718986240; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:16:42,167 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26717904896; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:16:52,206 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26717892608; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:17:02,228 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26714910720; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:17:12,260 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26716073984; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:17:22,296 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26713628672; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:17:32,332 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26712940544; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:17:42,356 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26711887872; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:17:52,388 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26712031232; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:18:02,432 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26710487040; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:18:12,446 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26709463040; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:18:22,480 E 10406 36357164] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-06-25_792862_10261 is over 95% full, available space: 26719281152; capacity: 994662584320. Object creation will fail if spilling is required.




.. raw:: html

    <div><i>Table length=401210</i>
    <table id="table11103717744" class="table-striped table-bordered table-condensed">
    <thead><tr><th>amplitude</th><th>x_mean</th><th>y_mean</th><th>x_stddev</th><th>y_stddev</th><th>theta</th><th>offset</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>96.13580818091697</td><td>2.9489472302249555</td><td>2.8459489641612317</td><td>0.6615126270671262</td><td>0.5982857173004258</td><td>0.5084229861100592</td><td>1.0716251998893254</td></tr>
    <tr><td>93.22999728230792</td><td>2.945971179314676</td><td>2.8321791067722404</td><td>0.6774501510971196</td><td>0.610760669468479</td><td>0.3324114335819944</td><td>0.7014456008012028</td></tr>
    <tr><td>94.83739704115126</td><td>2.94180146636494</td><td>2.82125715177069</td><td>0.6799938081616088</td><td>0.6012986779279296</td><td>0.40859458291440265</td><td>0.6536850854022012</td></tr>
    <tr><td>94.12867663693092</td><td>2.959077309637642</td><td>2.851779332482687</td><td>0.6862760348911927</td><td>0.5881373871809559</td><td>0.43462908400767114</td><td>1.193779403425883</td></tr>
    <tr><td>95.09568324476194</td><td>2.9776655111309482</td><td>2.843257410088901</td><td>0.6745053368568955</td><td>0.5948537560303414</td><td>0.4436336278798898</td><td>0.6849093750054519</td></tr>
    <tr><td>96.09149159005842</td><td>2.936042771168488</td><td>2.8397027851196017</td><td>0.6561865491418671</td><td>0.5921877323474933</td><td>0.31271186839202075</td><td>1.0385142236743439</td></tr>
    <tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
    <tr><td>87.74914558987756</td><td>2.6959405822169913</td><td>2.9659239812865574</td><td>0.7329552946836715</td><td>0.6397187407597072</td><td>0.4390018622237538</td><td>1.1684951454848664</td></tr>
    <tr><td>90.59196412930257</td><td>2.692654629778445</td><td>2.988546762625284</td><td>0.7198071514645737</td><td>0.6515556783383257</td><td>0.4521983178189343</td><td>0.6256298448573007</td></tr>
    <tr><td>88.12868979341526</td><td>2.690163957731549</td><td>3.0035667488552447</td><td>0.7117124119515923</td><td>0.6592995135570662</td><td>0.04656003733616875</td><td>1.127188290790881</td></tr>
    <tr><td>91.8260184010446</td><td>2.6757315498293726</td><td>2.97830644356315</td><td>0.7041102851009186</td><td>0.6372860805814633</td><td>0.24634331930260062</td><td>1.017394646015352</td></tr>
    <tr><td>91.13612563488523</td><td>2.715309467890296</td><td>2.988211023380153</td><td>0.6954974855383295</td><td>0.6308965496012459</td><td>0.35436411165274095</td><td>0.9773797864782711</td></tr>
    <tr><td>86.74058679676348</td><td>2.7005866029105525</td><td>2.9828734613586385</td><td>0.7116428778089073</td><td>0.664374229086644</td><td>0.1855252791897348</td><td>0.9938171960427725</td></tr>
    <tr><td>91.1500119503748</td><td>2.7045599931124356</td><td>2.97123181229028</td><td>0.7074935500971281</td><td>0.638799590023998</td><td>0.361961713119347</td><td>1.5019415061419104</td></tr>
    </table></div>



.. code:: ipython3

    results = {}
    
    for key in list(spk.gaussfit_results.keys()):
    
        results[key], _ = bin_fgs_to_science(time_since_start, 
                                             fg_time_since_start, 
                                             spk.gaussfit_results[key].value)

.. code:: ipython3

    print(results.keys())


.. parsed-literal::

    dict_keys(['amplitude', 'x_mean', 'y_mean', 'x_stddev', 'y_stddev', 'theta', 'offset'])


Let’s plot all of those parameters for the entire duration of the TSO:

.. code:: ipython3

    for key in list(results.keys()):
    
        plt.figure(figsize=(10,4))
    
        median = np.nanmedian(results[key])
        std = np.nanmedian(np.abs(results[key] - median)) * 1.4826
        plt.plot(time_since_start, results[key], color = 'tomato')
        
        plt.title(key+' for FGS 2D Gaussian Fit', fontsize = 18)
        plt.xlabel('Time since exposure start (hours)', fontsize = 18)
        plt.ylabel(key, fontsize = 18)
        
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        
        plt.xlim(np.min(time_since_start), np.max(time_since_start))
        plt.ylim(median-3*std,median+5*std)



.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_22_0.png



.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_22_1.png



.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_22_2.png



.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_22_3.png



.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_22_4.png



.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_22_5.png



.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_22_6.png


Neat! Many things to unpack.

First, note how the position in X (``x_mean``) and the standard
deviation in this direction (``x_stddev``), together with the standard
deviation on the y-direction (``y_stddev``) and the rotation angle of
the gaussian (``theta``) all oscillate in short frequency, in concert
with the science TSO. Let’s compare those time-series on top of the
science TSO for the first three hours. To do this, let’s create a helper
function that standarizes our regressors:

.. code:: ipython3

    def standarize(x):
    
        median = np.nanmedian(x)
        std = np.nanmedian(np.abs(x - median)) * 1.4826    
    
        return ( x - median ) / std

.. code:: ipython3

    variable = 'x_mean'
    
    plt.figure(figsize=(10,4))
    
    plt.plot(time_since_start, (f-1)*1e6, color = 'black', label = 'NIRISS/SOSS TSO')
    
    plt.plot(time_since_start, standarize(results[variable])*100*(-1), 
             color = 'tomato', label = r'FGS '+variable+r' $\times$ -100')
    
    plt.legend()
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux - 1 (ppm)', fontsize = 18)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), 3)
    plt.ylim(-500, 500)




.. parsed-literal::

    (-500.0, 500.0)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_25_1.png


Very nice correlation between variables! Also, note how the x-standard
deviation detects what appears to be a small tilt event:

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.plot(time_since_start, results['x_stddev'], color = 'tomato')
        
    plt.title('X-standard deviation for FGS 2D Gaussian Fit', fontsize = 18)
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('$\sigma_X$ (pix)', fontsize = 18)
        
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylim(0.7025,0.7115)
    plt.xlim(np.min(time_since_start), np.max(time_since_start))




.. parsed-literal::

    (0.0, 6.09976451843977)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_27_1.png


Seeing this from the actual TSO is quite difficult, because the tilt
event happened *just* during ingress:

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = '.', 
                             ms = 1, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    
    plt.xlim(np.min(time_since_start), np.max(time_since_start))
    plt.ylim(0.993, 1.001)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.twinx()
    variable = 'x_stddev'
    plt.plot(time_since_start, results[variable], 
             color = 'tomato', label = r'FGS '+variable)
    
    plt.ylim(0.7025,0.7115)
    
    plt.ylabel('$\sigma_X$ (pix)', fontsize = 18, color = 'tomato')
    
    plt.yticks(fontsize=16, color = 'tomato')




.. parsed-literal::

    (array([0.702, 0.704, 0.706, 0.708, 0.71 , 0.712]),
     [Text(1, 0.7020000000000001, '0.702'),
      Text(1, 0.7040000000000001, '0.704'),
      Text(1, 0.7060000000000001, '0.706'),
      Text(1, 0.7080000000000001, '0.708'),
      Text(1, 0.7100000000000001, '0.710'),
      Text(1, 0.7120000000000001, '0.712')])




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_29_1.png


To showcase further the beauty of ``spelunker`` as a TSO event detector,
let’s analyze one more dataset on which the tilt event is obvious: the
ERS observations of WASP-39 b with NIRSpec/G395H.

2. The case of WASP-39 b NIRSpec/G395H observations
===================================================

Let’s repeat the analysis for the transit WASP-39 b with NIRSpec/G395H.
Let’s study the NRS1 lightcurve presented in `Alderson et
al. (2023) <https://www.nature.com/articles/s41586-022-05591-3>`__:

.. code:: ipython3

    t, f, ferr = np.loadtxt('w39_lightcurve.dat', unpack = True, usecols = (0,1,2))

.. code:: ipython3

    tstart = t[0]
    time_since_start = (t-tstart)*24
    
    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = '.', 
                             ms = 1, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.title('Exoplanet transit of WASP-39b with NIRSpec/G395H', fontsize = 18)
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), np.max(time_since_start))




.. parsed-literal::

    (0.0, 8.256138402968645)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_32_1.png


Note that break in the transit light curve? That’s a tilt event. One
that ``spelunker`` can also detect! Let’s run the ``spelunker`` magic
for this program, which is `PID 1366 (PI: Batalha, co-PI: Bean,
Stevenson) <https://www.stsci.edu/jwst/science-execution/program-information?id=1366>`__.
This, in particular, is observation number 3, visit 1:

.. code:: ipython3

    spk = spelunker.load(pid=1366, obs_num='3', visit='1', save=True)


.. parsed-literal::

    Current working directory for spelunker: /Users/nespinoza/github/JWST-FGS-Spelunker/notebooks/spelunker_outputs
    
    Connecting with astroquery...


Let’s explore the guidestar (binned) photometry:

.. code:: ipython3

    fg_time_since_start = (spk.fg_time + 2400000.5 - tstart) * 24
    fbin, fbinerr = bin_fgs_to_science(time_since_start, fg_time_since_start, spk.fg_flux / np.nanmedian( spk.fg_flux ))

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, (fbin-1)*1e6, fbinerr*1e6, color = 'tomato')
    
    plt.title('(Binned) GS flux for PID 1541, visit 1, observation 1', fontsize = 18)
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux (ppm)', fontsize = 18)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), np.max(time_since_start))




.. parsed-literal::

    (0.0, 8.256138402968645)




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_37_1.png


Oh my. It’s not only one, but perhaps…two, three tilt events?:

.. code:: ipython3

    tstart = t[0]
    time_since_start = (t-tstart)*24
    
    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = '.', 
                             ms = 1, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.xlim(np.min(time_since_start), np.max(time_since_start))
    
    plt.twinx()
    
    plt.errorbar(time_since_start, (fbin-1)*1e6, fbinerr*1e6, color = 'tomato')
    
    plt.ylabel('Relative flux (FGS)', fontsize = 18, color = 'tomato')
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16, color = 'tomato')




.. parsed-literal::

    (array([-100000.,  -80000.,  -60000.,  -40000.,  -20000.,       0.,
              20000.,   40000.,   60000.,   80000.,  100000.]),
     [Text(1, -100000.0, '−100000'),
      Text(1, -80000.0, '−80000'),
      Text(1, -60000.0, '−60000'),
      Text(1, -40000.0, '−40000'),
      Text(1, -20000.0, '−20000'),
      Text(1, 0.0, '0'),
      Text(1, 20000.0, '20000'),
      Text(1, 40000.0, '40000'),
      Text(1, 60000.0, '60000'),
      Text(1, 80000.0, '80000'),
      Text(1, 100000.0, '100000')])




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_39_1.png


Very interesting. Let’s explore the gaussian fits to the data:

.. code:: ipython3

    spk.gauss2d_fit(ncpus=4)


.. parsed-literal::

    2023-07-28 15:30:48,400	ERROR node.py:605 -- Failed to connect to GCS. Please check `gcs_server.out` for more details.
    2023-07-28 15:30:56,510	INFO worker.py:1621 -- Started a local Ray instance.
    [2m[33m(raylet)[0m [2023-07-28 15:31:35,478 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26698432512; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:31:38,803 E 11745 36386431] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x0000000110968590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x000000011096877c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x0000000110123df4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x00000001101e6850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x00000001101e3d4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x0000000102b650d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:31:39,459 E 11744 36386429] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x000000010633c590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x000000010633c77c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x0000000105af7df4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x0000000105bba850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x0000000105bb7d4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x0000000102f990d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:31:40,156 E 11743 36386427] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x000000010333c590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x000000010333c77c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x0000000102af7df4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x0000000102bba850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x0000000102bb7d4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x000000010073d0d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/workers/default_worker.py", line 229, in <module>
    [2m[33m(raylet)[0m     ray._private.worker.connect(
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2282, in connect
    [2m[33m(raylet)[0m     tracing_hook_val = worker.gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    2023-07-28 15:31:40,958	WARNING worker.py:2006 -- A worker died or was killed while executing a task by an unexpected system error. To troubleshoot the problem, check the logs for the dead worker. RayTask ID: 71ae11604e3192d58a50af37db419c90501abe8f01000000 Worker ID: 370182201c86cff3ec027a5d3854a5ea6f8da94d220f5056be6a763d Node ID: fdc51f9019dd751dce1e007749df708720612e251c8d416e1984675e Worker IP address: 127.0.0.1 Worker port: 64730 Worker PID: 11742 Worker exit type: SYSTEM_ERROR Worker exit detail: Worker unexpectedly exits with a connection error code 2. End of file. There are some potential root causes. (1) The process is killed by SIGKILL by OOM killer due to high memory usage. (2) ray stop --force is called. (3) The worker is crashed unexpectedly due to SIGSEGV or other unexpected errors.
    [2m[33m(raylet)[0m [2023-07-28 15:31:40,956 E 11742 36386425] logging.cc:104: Stack trace: 
    [2m[33m(raylet)[0m  0   _raylet.so                          0x000000010aa68590 _ZN3raylsERNSt3__113basic_ostreamIcNS0_11char_traitsIcEEEERKNS_10StackTraceE + 84 ray::operator<<()
    [2m[33m(raylet)[0m 1   _raylet.so                          0x000000010aa6877c _ZN3ray16TerminateHandlerEv + 228 ray::TerminateHandler()
    [2m[33m(raylet)[0m 2   libc++abi.dylib                     0x00000001904c5ea4 _ZSt11__terminatePFvvE + 20 std::__terminate()
    [2m[33m(raylet)[0m 3   libc++abi.dylib                     0x00000001904c5e2c _ZSt9terminatev + 44 std::terminate()
    [2m[33m(raylet)[0m 4   _raylet.so                          0x000000010a223df4 _ZN3ray4core10CoreWorkerD2Ev + 1344 ray::core::CoreWorker::~CoreWorker()
    [2m[33m(raylet)[0m 5   _raylet.so                          0x000000010a2e6850 _ZN3ray4core21CoreWorkerProcessImplD2Ev + 664 ray::core::CoreWorkerProcessImpl::~CoreWorkerProcessImpl()
    [2m[33m(raylet)[0m 6   _raylet.so                          0x000000010a2e3d4c _ZN3ray4core17CoreWorkerProcess12HandleAtExitEv + 28 ray::core::CoreWorkerProcess::HandleAtExit()
    [2m[33m(raylet)[0m 7   libsystem_c.dylib                   0x00000001903f7de0 __cxa_finalize_ranges + 480 __cxa_finalize_ranges
    [2m[33m(raylet)[0m 8   libsystem_c.dylib                   0x00000001903f7b74 exit + 44 exit
    [2m[33m(raylet)[0m 9   libdyld.dylib                       0x0000000190518ec4 _ZNK5dyld416LibSystemHelpers6getenvEPKc + 0 dyld4::LibSystemHelpers::getenv()
    [2m[33m(raylet)[0m 10  dyld                                0x00000001055310d8 start + 596 start
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m [2023-07-28 15:31:45,486 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26702929920; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:31:55,487 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26702495744; capacity: 994662584320. Object creation will fail if spilling is required.


.. parsed-literal::

    [2m[1m[36m(autoscaler +28m44s)[0m Tip: use `ray status` to view detailed cluster status. To disable these messages, set RAY_SCHEDULER_EVENTS=0.
    [2m[1m[33m(autoscaler +28m44s)[0m Warning: The following resource request cannot be scheduled right now: {'CPU': 4.0}. This is likely due to all cluster resources being claimed by actors. Consider creating fewer actors or adding more nodes to this Ray cluster.


.. parsed-literal::

    [2m[33m(raylet)[0m [2023-07-28 15:32:05,491 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26702692352; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:32:15,508 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26704752640; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:32:25,516 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26704728064; capacity: 994662584320. Object creation will fail if spilling is required.


.. parsed-literal::

    [2m[1m[33m(autoscaler +29m19s)[0m Warning: The following resource request cannot be scheduled right now: {'CPU': 4.0}. This is likely due to all cluster resources being claimed by actors. Consider creating fewer actors or adding more nodes to this Ray cluster.


.. parsed-literal::

    [2m[33m(raylet)[0m [2023-07-28 15:32:35,522 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26705850368; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:32:40,621 E 11739 36386340] (raylet) worker_pool.cc:544: Some workers of the worker process(11906) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:32:40,623 E 11739 36386340] (raylet) worker_pool.cc:544: Some workers of the worker process(11907) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:32:40,623 E 11739 36386340] (raylet) worker_pool.cc:544: Some workers of the worker process(11908) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:32:40,960 E 11739 36386340] (raylet) worker_pool.cc:544: Some workers of the worker process(11910) have not registered within the timeout. The process is dead, probably it crashed during start.
    [2m[33m(raylet)[0m [2023-07-28 15:32:45,522 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26704105472; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m Failed to connect to GCS. Please check `gcs_server.out` for more details.
    [2m[33m(raylet)[0m Failed to publish error: Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m  [type version_mismatch]
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/worker.py", line 2113, in connect
    [2m[33m(raylet)[0m     node.check_version_info()
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/node.py", line 339, in check_version_info
    [2m[33m(raylet)[0m     cluster_metadata = ray_usage_lib.get_cluster_metadata(self.get_gcs_client())
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/usage/usage_lib.py", line 710, in get_cluster_metadata
    [2m[33m(raylet)[0m     gcs_client.internal_kv_get(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2132, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2120, in ray._raylet._auto_reconnect.wrapper
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2185, in ray._raylet.GcsClient.internal_kv_get
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 410, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.RpcError: failed to connect to all addresses
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m During handling of the above exception, another exception occurred:
    [2m[33m(raylet)[0m 
    [2m[33m(raylet)[0m Traceback (most recent call last):
    [2m[33m(raylet)[0m   File "/Users/nespinoza/opt/anaconda3/envs/newen/lib/python3.9/site-packages/ray/_private/utils.py", line 203, in publish_error_to_driver
    [2m[33m(raylet)[0m     gcs_publisher.publish_error(
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 2357, in ray._raylet.GcsPublisher.publish_error
    [2m[33m(raylet)[0m   File "python/ray/_raylet.pyx", line 402, in ray._raylet.check_status
    [2m[33m(raylet)[0m ray.exceptions.GetTimeoutError: Failed to publish after retries: failed to connect to all addresses
    [2m[33m(raylet)[0m [2023-07-28 15:32:55,523 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26709344256; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:33:05,613 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26709024768; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:33:15,704 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26711957504; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:33:25,795 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26710233088; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:33:35,890 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26707730432; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:33:45,979 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26712797184; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:33:56,067 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26709250048; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:34:06,158 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26707238912; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:34:16,247 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26705580032; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:34:26,338 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26704465920; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:34:36,429 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26704707584; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:34:46,518 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26703581184; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:34:56,604 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26702155776; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:35:06,693 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26701422592; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:35:16,783 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 26702127104; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:35:26,870 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30996672512; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:35:36,958 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30992228352; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:35:47,049 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30994817024; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:35:57,138 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30993162240; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:36:07,229 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30990024704; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:36:17,323 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30994264064; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:36:27,413 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30992416768; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:36:37,506 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30991323136; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:36:47,594 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30986346496; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:36:57,680 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30984134656; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:37:07,769 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30980599808; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:37:17,862 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30979162112; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:37:27,947 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30983151616; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:37:38,036 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30986510336; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:37:48,122 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30984798208; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:37:58,213 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30980804608; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:38:08,301 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30980399104; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:38:18,389 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30976733184; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:38:28,482 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30983811072; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:38:38,575 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30979407872; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:38:48,667 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30972837888; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:38:58,755 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30975447040; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:39:08,854 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30969364480; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:39:18,942 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30972436480; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:39:29,026 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30973825024; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:39:39,120 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30968504320; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:39:49,210 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30968496128; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:39:59,294 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30967545856; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:40:09,382 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30967238656; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:40:19,480 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30967599104; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:40:29,558 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30964912128; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:40:39,649 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30963363840; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:40:49,738 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30960181248; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:40:59,827 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 30959869952; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:41:09,916 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29886369792; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:41:20,002 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29884792832; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:41:30,087 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29884694528; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:41:40,172 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29889052672; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:41:50,244 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29904269312; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:42:00,264 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29905477632; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:42:10,352 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29908611072; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:42:20,441 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29905760256; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:42:30,525 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29905514496; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:42:40,601 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29905219584; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:42:50,678 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29903106048; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:43:00,762 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29902733312; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:43:10,843 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29898997760; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:43:20,930 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29896843264; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:43:31,017 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29896830976; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:43:41,105 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29900427264; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:43:51,191 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29897138176; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:44:01,274 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29895049216; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:44:11,364 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29895041024; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:44:21,449 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29895520256; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:44:31,531 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29894946816; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:44:41,620 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29895917568; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:44:51,704 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29896093696; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:45:01,785 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29893738496; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:45:11,872 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29902282752; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:45:21,954 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29898493952; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:45:32,045 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29890584576; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:45:42,143 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29885095936; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:45:52,233 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29883973632; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:46:02,315 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29880090624; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:46:12,401 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29882724352; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:46:22,485 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29881815040; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:46:32,565 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29877051392; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:46:42,650 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29878689792; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:46:52,742 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29875945472; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:47:02,829 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29876817920; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:47:12,915 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29875523584; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:47:22,994 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29874618368; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:47:33,077 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29872230400; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:47:43,159 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29871153152; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:47:53,238 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29871153152; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:48:03,321 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29868302336; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:48:13,407 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29871022080; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:48:23,497 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29869314048; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:48:33,581 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29865148416; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:48:43,675 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29863239680; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:48:53,768 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29864644608; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:49:03,856 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29860782080; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:49:13,936 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29862035456; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:49:24,030 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29875097600; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:49:34,111 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29873397760; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:49:44,192 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29872562176; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:49:54,278 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29872529408; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:50:04,363 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29869674496; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:50:14,454 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29867143168; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:50:24,541 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29867495424; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:50:34,627 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29864673280; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:50:44,716 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29868904448; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:50:54,805 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29857140736; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:51:04,892 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29845704704; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:51:14,977 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29836926976; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:51:25,063 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29835341824; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:51:35,146 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29833109504; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:51:45,247 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29837807616; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:51:55,335 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29834792960; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:52:05,421 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29833568256; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:52:15,510 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29833859072; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:52:25,597 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29828616192; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:52:35,689 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29827952640; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:52:45,772 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29827670016; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:52:55,860 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29823766528; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:53:05,869 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29829320704; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:53:15,957 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29823389696; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:53:26,043 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29819998208; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:53:36,131 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29819211776; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:53:46,218 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29819518976; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:53:56,307 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29844058112; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:54:06,390 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29842558976; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:54:16,475 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29840211968; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:54:26,565 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29838585856; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:54:36,650 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29836800000; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:54:46,742 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29835550720; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:54:56,823 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29834137600; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:55:06,916 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29829791744; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:55:16,995 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29828759552; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:55:27,084 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29827907584; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:55:37,167 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29827571712; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:55:47,256 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29826605056; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:55:57,347 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29825486848; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:56:07,434 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29821673472; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:56:17,526 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29817495552; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:56:27,614 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29818433536; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:56:37,702 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29813547008; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:56:47,792 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29810462720; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:56:57,879 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29810266112; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:57:07,966 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29824274432; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:57:18,051 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29826334720; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:57:28,136 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29822513152; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:57:38,224 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29823807488; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:57:48,314 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29822033920; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:57:58,402 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29820821504; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:58:08,488 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29816836096; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:58:18,581 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29825998848; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:58:28,665 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29820809216; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:58:38,750 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29824659456; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:58:48,835 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29822529536; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:58:58,916 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29820878848; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:59:09,004 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29817761792; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:59:19,094 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29814386688; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:59:29,184 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29812051968; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:59:39,272 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29814005760; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:59:49,354 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29811355648; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 15:59:59,448 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29809041408; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:00:09,538 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29814325248; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:00:19,627 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29811761152; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:00:29,707 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29810630656; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:00:39,791 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29814411264; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:00:49,876 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29812551680; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:00:59,898 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29806743552; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:01:09,995 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29805703168; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:01:20,080 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29804916736; capacity: 994662584320. Object creation will fail if spilling is required.
    [2m[33m(raylet)[0m [2023-07-28 16:01:30,164 E 11739 36386361] (raylet) file_system_monitor.cc:111: /tmp/ray/session_2023-07-28_15-30-41_512027_10261 is over 95% full, available space: 29805150208; capacity: 994662584320. Object creation will fail if spilling is required.




.. raw:: html

    <div><i>Table length=987136</i>
    <table id="table11109265360" class="table-striped table-bordered table-condensed">
    <thead><tr><th>amplitude</th><th>x_mean</th><th>y_mean</th><th>x_stddev</th><th>y_stddev</th><th>theta</th><th>offset</th></tr></thead>
    <thead><tr><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th></tr></thead>
    <tr><td>101.2130875542324</td><td>3.3337141869616698</td><td>3.1211539650536118</td><td>0.6455960581549347</td><td>0.6787948700644653</td><td>-0.9882345633774572</td><td>1.2600175942667253</td></tr>
    <tr><td>103.05907718346985</td><td>3.33371710766348</td><td>3.142691936508253</td><td>0.6215190964749095</td><td>0.6992885770745251</td><td>8.421751880874677</td><td>0.4315062171942471</td></tr>
    <tr><td>107.62423515840916</td><td>3.3329848169769485</td><td>3.120927859230632</td><td>0.6134283966398637</td><td>0.6706534574534939</td><td>2.1012697533767146</td><td>1.0340357839332834</td></tr>
    <tr><td>105.13550772812634</td><td>3.340446149057725</td><td>3.1465417724009326</td><td>0.6194982306101161</td><td>0.6743140206939372</td><td>8.33284963100652</td><td>0.3417050813233455</td></tr>
    <tr><td>105.81117632539961</td><td>3.3294497348834233</td><td>3.1567470988399315</td><td>0.687199451837556</td><td>0.619261893520507</td><td>-8.95975328933037</td><td>1.0360083315734216</td></tr>
    <tr><td>106.22326098792206</td><td>3.317701776351794</td><td>3.1573008366804722</td><td>0.6732302628935715</td><td>0.6200459189744063</td><td>-2.6905026977639435</td><td>1.393730020315556</td></tr>
    <tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
    <tr><td>92.73785019284357</td><td>3.1873002094970913</td><td>3.7376647173565347</td><td>0.6696140357514423</td><td>0.6149201200291508</td><td>0.5719390312820085</td><td>0.7572221044577753</td></tr>
    <tr><td>97.97532198761044</td><td>3.188108276531111</td><td>3.719662716727466</td><td>0.6720903777369073</td><td>0.6015700463433664</td><td>0.7071603924088367</td><td>1.2912940683375458</td></tr>
    <tr><td>97.95090537050895</td><td>3.1717576547844515</td><td>3.7113592931717414</td><td>0.6551480798503434</td><td>0.6090486591123077</td><td>0.39401270499008983</td><td>0.39912480068024053</td></tr>
    <tr><td>90.74370580762309</td><td>3.170361673138889</td><td>3.7116379447082983</td><td>0.6730234514605041</td><td>0.6319562125739976</td><td>0.6698834836744338</td><td>1.0697096845533383</td></tr>
    <tr><td>92.47649592912568</td><td>3.185794302482265</td><td>3.7194128606132377</td><td>0.6635212462728267</td><td>0.6253418596247133</td><td>0.43507819808824466</td><td>0.7636161336946441</td></tr>
    <tr><td>95.9522036162299</td><td>3.1873304341030315</td><td>3.7316706015647925</td><td>0.663425284985132</td><td>0.6000655892444249</td><td>0.6216247022228777</td><td>1.0357714545777559</td></tr>
    <tr><td>94.44828097982054</td><td>3.1809946513558147</td><td>3.732769547434784</td><td>0.6717574169327243</td><td>0.6213462199187172</td><td>0.7776036482787176</td><td>0.8617779884022134</td></tr>
    </table></div>



Let’s bin this to the science time-stamps:

.. code:: ipython3

    results = {}
    
    for key in list(spk.gaussfit_results.keys()):
    
        results[key], _ = bin_fgs_to_science(time_since_start, 
                                             fg_time_since_start, 
                                             spk.gaussfit_results[key].value)

Let’s go right away to the X standard deviation:

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = '.', 
                             ms = 1, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    
    plt.xlim(np.min(time_since_start), np.max(time_since_start))
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.twinx()
    variable = 'x_stddev'
    plt.plot(time_since_start, results[variable], 
             color = 'tomato', label = r'FGS '+variable)
    
    plt.ylabel('$\sigma_X$ (pix)', fontsize = 18, color = 'tomato')
    
    plt.yticks(fontsize=16, color = 'tomato')




.. parsed-literal::

    (array([0.662, 0.664, 0.666, 0.668, 0.67 , 0.672]),
     [Text(1, 0.662, '0.662'),
      Text(1, 0.664, '0.664'),
      Text(1, 0.666, '0.666'),
      Text(1, 0.668, '0.668'),
      Text(1, 0.67, '0.670'),
      Text(1, 0.672, '0.672')])




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_45_1.png


Very interesting! The “tilt” event is beautifully detected by the
guidestar data. Best of all, we can see the event at any resolution we
want thanks to it. Let’s write a function that can bin the data at any
temporal resolution so we can see this in action:

.. code:: ipython3

    def bin_data(x,y,n_bin):
        
        x_bins = []
        y_bins = []
        y_err_bins = []
        
        for i in range(0,len(x),n_bin):
            
            x_bins.append(np.median(x[i:i+n_bin-1]))
            y_bins.append(np.median(y[i:i+n_bin-1]))
            y_err_bins.append(np.sqrt(np.var(y[i:i+n_bin-1]))/np.sqrt(len(y[i:i+n_bin-1])))
            
        return np.array(x_bins),np.array(y_bins),np.array(y_err_bins)

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = 'o', 
                             ms = 5, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    
    plt.xlim(4.5,5.0)
    plt.ylim(0.973, 0.980)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.twinx()
    
    variable = 'x_stddev'
    plt.plot(time_since_start, results[variable], 
             color = 'tomato', label = r'FGS '+variable)
    
    tbin, ybin, _ = bin_data((spk.fg_time + 2400000.5 - tstart) * 24, spk.gaussfit_results[variable].value, n_bin = 300)
    plt.plot(tbin, ybin, '.-',
             color = 'red', alpha = 0.3)
    
    plt.plot(time_since_start, results[variable], 
             'o', color = 'tomato')
    
    plt.xlim(4.5,5.0)
    plt.ylim(0.66,0.675)
    
    plt.ylabel('$\sigma_X$ (pix)', fontsize = 18, color = 'tomato')
    
    plt.yticks(fontsize=16, color = 'tomato')




.. parsed-literal::

    (array([0.66 , 0.662, 0.664, 0.666, 0.668, 0.67 , 0.672, 0.674, 0.676]),
     [Text(1, 0.66, '0.660'),
      Text(1, 0.662, '0.662'),
      Text(1, 0.664, '0.664'),
      Text(1, 0.666, '0.666'),
      Text(1, 0.668, '0.668'),
      Text(1, 0.67, '0.670'),
      Text(1, 0.672, '0.672'),
      Text(1, 0.674, '0.674'),
      Text(1, 0.676, '0.676')])




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_48_1.png


Interestingly, in this case, the Y-standard deviation samples the event
even better:

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = 'o', 
                             ms = 5, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    
    plt.xlim(4.5,5.0)
    plt.ylim(0.973, 0.980)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.twinx()
    
    variable = 'y_stddev'
    plt.plot(time_since_start, results[variable], 
             color = 'cornflowerblue', label = r'FGS '+variable)
    
    tbin, ybin, _ = bin_data((spk.fg_time + 2400000.5 - tstart) * 24, spk.gaussfit_results[variable].value, n_bin = 300)
    plt.plot(tbin, ybin, '.-',
             color = 'cornflowerblue', alpha = 0.3)
    
    plt.plot(time_since_start, results[variable], 
             'o', color = 'cornflowerblue')
    
    plt.xlim(4.5,5.0)
    plt.ylim(0.58,0.63)
    
    plt.ylabel('$\sigma_Y$ (pix)', fontsize = 18, color = 'cornflowerblue')
    
    plt.yticks(fontsize=16, color = 'cornflowerblue')




.. parsed-literal::

    (array([0.58, 0.59, 0.6 , 0.61, 0.62, 0.63]),
     [Text(1, 0.5800000000000001, '0.58'),
      Text(1, 0.5900000000000001, '0.59'),
      Text(1, 0.6000000000000001, '0.60'),
      Text(1, 0.6100000000000001, '0.61'),
      Text(1, 0.6200000000000001, '0.62'),
      Text(1, 0.6300000000000001, '0.63')])




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_50_1.png


It is very interesting how while the change of the PSF in the guidestar
is almost instantaneous, the PSF change in the TSO is slightly smoother.
This smooth behavior is actually quite nicely tracked by the mean
positions:

.. code:: ipython3

    plt.figure(figsize=(10,4))
    
    plt.errorbar(time_since_start, f, ferr, fmt = 'o', 
                             ms = 5, mfc = 'black', mec = 'black', 
                             elinewidth = 1, ecolor = 'black')
    
    plt.xlabel('Time since exposure start (hours)', fontsize = 18)
    plt.ylabel('Relative flux', fontsize = 18)
    
    
    plt.xlim(4.5,5.0)
    plt.ylim(0.973, 0.980)
    
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    plt.twinx()
    
    variable = 'y_mean'
    plt.plot(time_since_start, results[variable], 
             color = 'cornflowerblue', label = r'FGS '+variable)
    
    tbin, ybin, _ = bin_data((spk.fg_time + 2400000.5 - tstart) * 24, spk.gaussfit_results[variable].value, n_bin = 300)
    plt.plot(tbin, ybin, '.-',
             color = 'cornflowerblue', alpha = 0.3)
    
    plt.plot(time_since_start, results[variable], 
             'o', color = 'cornflowerblue')
    
    plt.xlim(4.5,5.0)
    plt.ylim(3.71,3.74)
    
    plt.ylabel('$Y$ (pix)', fontsize = 18, color = 'cornflowerblue')
    
    plt.yticks(fontsize=16, color = 'cornflowerblue')




.. parsed-literal::

    (array([3.71 , 3.715, 3.72 , 3.725, 3.73 , 3.735, 3.74 ]),
     [Text(1, 3.71, '3.710'),
      Text(1, 3.715, '3.715'),
      Text(1, 3.72, '3.720'),
      Text(1, 3.725, '3.725'),
      Text(1, 3.73, '3.730'),
      Text(1, 3.735, '3.735'),
      Text(1, 3.74, '3.740')])




.. image:: fgs-spelunker-and-tsos_files/fgs-spelunker-and-tsos_52_1.png

