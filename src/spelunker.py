import numpy as np
import scipy.optimize as opt
import pandas as pd
import datetime

from jwst import datamodels

from astroquery.mast import Observations

from astropy.table import Table, hstack, vstack
from astropy.time import Time
from astropy import units as u
from astropy.timeseries import LombScargle

from astroplan.plots import plot_finder_image
from astropy import coordinates
from astropy.coordinates import SkyCoord
from astroquery.skyview import SkyView

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.animation as animation
import matplotlib.transforms as transforms

import ray
import glob
import os


try:
    from jwstuser import engdb
    jwstuser_installed = True

except ImportError:
    jwstuser_installed = False
    print('Warning: jwstuser is not installed. mnemonics() will not work.')


class load:

    """
    This class loads auxiliary FITS files of FG-GS guidestar data from MAST from a JWST Program ID. Using "load" will create attributes as np.array for the spatial data, time array,
    and flux array. After downloading the MAST files, other attributes become available, such as "object_properties" and "fg_table", which contains metadata and stellar properties 
    about the guidestar data. Example usage:

                >>> spk = spelunker.load(pid=1534,obs_num='2',token='abcdefg1234567')

    Parameters
    ----------
                
        :param dir: (optional, string)
            The parameter "dir" allows the user to change the directory of the script. Changing the directory dictates where the class loads, creates folders, and save data and results.
            If ignored, "load" will save the current directory as an attribute for use with other functions.

        :param pid: (optional, int or string)
            An Python ``int`` or ``string`` of the Program ID for a given JWST program (i.e.: 1296, '1296'). Program IDs that are publicly avaliable are able to be accessed. Program IDs 
            with exclusive rights will not be accessible unless the user has a MAST API token containing the rights to the program ID. The MAST API token must be entered with param "token"
            to use programs marked with exclusive rights. The program ID needs to be entered to use all of functions in "load". However, the pid can be ignored to use specific functions
            within the class, though functionality will be severely limited. You can also forgo inputting the program ID to create a 'blank' object, then later input data for a pid using
            the function ``download`` given ``pid`` and the following parameters. For example:
                    
                    >>> spk = spelunker.load()
                    >>> spk.download(pid=1534)

        :param obs_num: (optional, int or string)
            The integer or string of the observation number of a program ID. Inputting the observation number will filter the auxiliary FITS to only include the files with the given
            number before downloading all files. Note that when you specifify the observation number, you must also input the ``pid``. 
            
        :param visit: (optional, int or string)
            The visit parameter requires a Python ``int`` or ``string``. Inputting the visit number will fillter the auxiliary FITS to only include the files with the given parameter
            before downloading all files.

        :param visit_group: (optional, int or string)
            The visit group parameter requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``vist_group`` will filter the
            data to only include files with the given visit group before loading the data into the object.

        :param parallel_sequence_id: (optional, int or string)
            The parallel sequence ID parameter requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``parallel_sequence_id`` 
            will filter the data to only include files with the given parallel sequence ID before loading the data into the object.

        :param activity_number: (optional, int or string)
            The activity number requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``activity_number`` will filter the
            data to only include files with the given activity number before loading the data into the object.

        :param exposure_number: (optional, int or string)
            The exposure_number requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``exposure_number`` will filter the
            data to only include files with the given exposure number before loading the data into the object.

        :param dir_seg: (optional, int or string)
            The segment number requires a Python ``int`` or ``string``. Specifying the segment number will include guidestar FITS with only the input number in the header if the segmented
            files has 'seg' within the naming scheme for the filename. This parameter will filter associated program ID FITS after downloading has completed and before the data is loaded
            into the object.

        :param guider: (optional, int or string)
            Similarly to `dir_seg`, the guider number requires a Python ``int`` or ``string``. Specifying the guider number will include guidestar FITS with only the input number in the 
            header if the files has 'guider' within the naming scheme for the filename. This parameter will filter associated program ID FITS after downloading has completed and before 
            the data is loaded into the object.

        :param calib_level: (optional, int or string)
            The calibration level can be specified as a Python ``int`` or ``string``. This parameter determines the calibration level for the FG-GS FITS. The only inputs for this
            parameter are 1, and 2, referring to higher stages of calibration with the JWST pipeline. Inputting the observation number will filter the auxiliary FITS to only include the
            files with the given calibration level before downloading all files. If ``calib_level`` is ignored, the calibration level will default to 2.

        :param save: (optional, boolean)
            Python ``boolean`` that either saves Gaussian fit, periodogram results, and other useful data automatically or does not save results except for the downloaded data. If 
            ignored, ``save`` will default to false.


        :param token: (optional, string)
            This parameter takes a Python ``string`` of a given MAST API token. This is required for certain functions, such as accessing JWST mnemonics and technical information, as well
            as program IDs with exclusive rights. If ignored, a blank MAST API token attribute will be created and access to JWST mnemonics will be limited and any program ID marked with
            exclusive rights will not be loaded.
    """


    def __init__(self, dir='None', pid='None', obs_num='None',  visit='None', visit_group='None', parallel_sequence_id='None', 
                 activity_number='None', exposure_number='None', dir_seg='None', guider='None', calib_level=2, save=False, token = None):

        '''
        The initalization function for ``load``. This function reads the current directory and initalizes essential attributes.

        Parameters

            param dir: (optional, string)
                The parameter "dir" allows the user to change the directory of the script. Changing the directory dictates where the class loads, creates folders, and save data and results.
                If ignored, "load" will save the current directory as an attribute for use with other functions.

            :param pid: (optional, int or string)
                An Python ``int`` or ``string`` of the Program ID for a given JWST program (i.e.: 1296, '1296'). Program IDs that are publicly avaliable are able to be accessed. Program IDs 
                with exclusive rights will not be accessible unless the user has a MAST API token containing the rights to the program ID. The MAST API token must be entered with param "token"
                to use programs marked with exclusive rights. The program ID needs to be entered to use all of functions in "load". However, the pid can be ignored to use specific functions
                within the class, though functionality will be severely limited. You can also forgo inputting the program ID to create a 'blank' object, then later input data for a pid using
                the function ``download`` given ``pid`` and the following parameters. For example:
                        
                        >>> spk = spelunker.load()
                        >>> spk.download(pid=1534)

            :param obs_num: (optional, int or string)
                The integer or string of the observation number of a program ID. Inputting the observation number will filter the auxiliary FITS to only include the files with the given
                number before downloading all files. Note that when you specifify the observation number, you must also input the ``pid``. 
                
            :param visit: (optional, int or string)
                The visit parameter requires a Python ``int`` or ``string``. Inputting the visit number will fillter the auxiliary FITS to only include the files with the given parameter
                before downloading all files.

            :param visit_group: (optional, int or string)
                The visit group parameter requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``vist_group`` will filter the
                data to only include files with the given visit group before loading the data into the object.

            :param parallel_sequence_id: (optional, int or string)
                The parallel sequence ID parameter requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``parallel_sequence_id`` 
                will filter the data to only include files with the given parallel sequence ID before loading the data into the object.

            :param activity_number: (optional, int or string)
                The activity number requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``activity_number`` will filter the
                data to only include files with the given activity number before loading the data into the object.

            :param exposure_number: (optional, int or string)
                The exposure_number requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``exposure_number`` will filter the
                data to only include files with the given exposure number before loading the data into the object.

            :param dir_seg: (optional, int or string)
                The segment number requires a Python ``int`` or ``string``. Specifying the segment number will include guidestar FITS with only the input number in the header if the segmented
                files has 'seg' within the naming scheme for the filename. This parameter will filter associated program ID FITS after downloading has completed and before the data is loaded
                into the object.

            :param guider: (optional, int or string)
                Similarly to `dir_seg`, the guider number requires a Python ``int`` or ``string``. Specifying the guider number will include guidestar FITS with only the input number in the 
                header if the files has 'guider' within the naming scheme for the filename. This parameter will filter associated program ID FITS after downloading has completed and before 
                the data is loaded into the object.

            :param calib_level: (optional, int or string)
                The calibration level can be specified as a Python ``int`` or ``string``. This parameter determines the calibration level for the FG-GS FITS. The only inputs for this
                parameter are 1, and 2, referring to higher stages of calibration with the JWST pipeline. Inputting the observation number will filter the auxiliary FITS to only include the
                files with the given calibration level before downloading all files. If ``calib_level`` is ignored, the calibration level will default to 2.

            :param save: (optional, boolean)
                Python ``boolean`` that either saves Gaussian fit, periodogram results, and other useful data automatically or does not save results except for the downloaded data. If 
                ignored, ``save`` will default to false.


            :param token: (optional, string)
                This parameter takes a Python ``string`` of a given MAST API token. This is required for certain functions, such as accessing JWST mnemonics and technical information, as well
                as program IDs with exclusive rights. If ignored, a blank MAST API token attribute will be created and access to JWST mnemonics will be limited and any program ID marked with
                exclusive rights will not be loaded.
        '''
        
        self.init_dir = os.getcwd()

        created_dir = "spelunker_outputs"

        if dir != 'None':
            if not os.path.exists(dir+'/'+created_dir):  #https://www.geeksforgeeks.org/how-to-create-directory-if-it-does-not-exist-using-python/
                os.makedirs(created_dir)
        
            self.directory = dir+created_dir

        else:
            if not os.path.exists(self.init_dir+'/'+created_dir):
                os.makedirs(created_dir)
        
            self.directory = os.getcwd()+'/'+created_dir

        print('Current working directory for spelunker: '+self.directory+'\n')

        self.mast_api_token = token
        self.fg_table = None

        self.fg_array = None
        self.fg_time = None
        self.fg_flux = None

        self.gaussfit_results = None
        self.quickfit_results = None
        self.pgram_results = None

        self.object_properties = None

        self.fontsize = 14
        
        
        if pid != 'None':
            self.download(pid, obs_num=obs_num, visit=visit, visit_group=visit_group, parallel_sequence_id=parallel_sequence_id,
                          activity_number=activity_number, exposure_number=exposure_number, dir_seg=dir_seg, guider=guider, calib_level=calib_level, save=save, token = token)
            
            self.pid = pid

            if obs_num == 'None':
                self.obs_num = None
            else:
                self.obs_num = obs_num

            if visit == 'None':
                self.visit = None
            else:
                self.visit = visit
            

    def stitcher(self, fg_timeseries):
        '''
        Stitches multiple of JWST objects together and creates a data array from the guidestar frames. 

        Parameters
        ----------

            fg_timeseries : list
                (Mandatory) This list should include JWST datamodel objects in each element.

        
        Returns
        -------

            The function creates an np.array of guidestar frames.
            
        '''

        fg_cube = []
        for data in fg_timeseries:
            for i in data.data:
                fg_cube.append(i)

        fg_array = np.array(fg_cube)

        return fg_array


    '''
    PREPROCESSING METHODS
    ------------------------------------------------------------------------------------------------------
    '''

    def optimize_photometry(self, mask = None, min_pixel = 1, max_pixel = 30, n_bin = 100):
        '''
        This function optimizes the photometry in `self.fg_flux` via either a user-defined `mask` (of size [8,8]), which has values of 1 
        in the pixels that want to be included in the photometry, or automatically by testing masks of the brightest npixels from min_pixel to 
        max_pixel.

        Parameters
        ----------
            
            mask : np.array
                (Optional) mask array. Must have same dimensions as an FGS frame, ie, [8,8]. Values of 1 are pixels to be used, values of 0 won't be used.

            min_pixel : int
                (Optional) minimum range of pixels to explore to optimize a mask.

            max_pixel : int
                (Optional) maximum range of pixels to explore to optimize a mask.

            n_bin : int
                (Optional) number of datapoints to bin to compare the precision between pixel mask sizes.

        Returns
        -------

            The function rewrites the self.fg_flux attribute to have the optimized photometry. It also rewrites the self.photometry_mask containing the 
            optimized mask.

        '''

        if mask is None:

            # Find optimal mask. To do this, we first get the median frame:
            median_frame = np.nanmedian( np.array( self.fg_array ), axis = 0)
            pixels = []
            rms = []
            for i in range(min_pixel, max_pixel + 1):

                # Get current mask:
                mask = self.get_mask(median_frame, npixels = i)

                # Get photometry; bin it:
                photometry = np.sum( np.array(self.fg_array) * mask, axis = (1,2) )
                _, fbin, fbinerr = self.bin_data( self.fg_time, photometry, n_bin = n_bin )

                # Calculate rms; save it:
                rms.append( np.nanmedian( fbinerr / np.nanmedian(fbin) ) )
                pixels.append(i)

            # Find optimal mask:
            idx = np.where(np.array(rms) == np.min(rms))[0]

            # Get photometry with this mask:
            mask = self.get_mask(median_frame, npixels = pixels[idx[0]])

        self.photometry_mask = mask
        self.fg_flux = list(np.sum( np.array(self.fg_array) * mask, axis = (1,2) ))

    def negative_flux_preprocessing(self, fg_array, fg_time, fg_flux):
        '''
        The function ``negative_flux_preprocessing`` removes elements in the parameters that are marked with counts with negative values from the parameter
        ``fg_flux``. All parameters needs to have the same amount of elements in each array.

        Parameters
        ----------

            fg_array : np.array
                (Mandatory) The spatial array of guidestar data. This must have a shape of (X, 8, 8), with X being the number of elements in the array.

            fg_time : np.array
                (Mandatory) The time array of guidestar data. This must have a shape of (X,), with X being the number of elements in the array.

            fg_flux : np.array
                (Mandatory) The flux array of guidestar data. This must have a shape of (X,), with X being the number of elements in the array.

        Returns
        -------

            Any elements in fg_flux that have negative values will be removed for elements in the given parameters as well. This function will return the three
            processed arrays.

        '''

        bool_mask = fg_flux >= 0

        fg_array = fg_array[bool_mask]
        fg_time = fg_time[bool_mask]
        fg_flux = fg_flux[bool_mask]

        return fg_array, fg_time, fg_flux
    
    def normalization_flux_preprocessing(self,):
        '''
        Normalize star counts from guidestar datamodel objects. This function takes the data from the guidestar objects and divides it by the median of the
        guidestar data for each frame.

        Parameters
        ----------
        
            No parameters, but this function takes the self.fg_timeseries as an input.

        Returns
        -------
            
            This function returns two arrays: the normalized spatial array and the normalized flux array.


        '''

        norm_array = []
        for data in self.fg_timeseries:
            norm_array.append(data.data/np.nanmedian(data.data))

        norm_array = np.vstack(norm_array)
        norm_flux = np.nansum(norm_array, axis=(1,2))

        return norm_array, norm_flux

    def time_sort_preprocessing(self, fg_array, fg_time, fg_flux):
        '''
        Sorts the three arrays by ascending order of the values in ``fg_time``, regardless of filename or file order.

        Parameters
        ----------

            fg_array : np.array
                (Mandatory) The spatial array of guidestar data. This must have a shape of (X, 8, 8), with X being the number of elements in the array.

            fg_time : np.array
                (Mandatory) The time array of guidestar data. This must have a shape of (X,), with X being the number of elements in the array.

            fg_flux : np.array
                (Mandatory) The flux array of guidestar data. This must have a shape of (X,), with X being the number of elements in the array.

        Returns
        -------

            Three numpy arrays that are sorted by ``fg_time`` in ascending order.

        '''

        table = Table()
        table['fg_array'] = fg_array
        table['fg_time'] = fg_time
        table['fg_flux'] = fg_flux

        table.sort(['fg_time'])

        return np.array(table['fg_array']), np.array(table['fg_time']), np.array(table['fg_flux'])


    '''
    DOWNLOADING FROM ASTROQUERY AND READING FILES
    ------------------------------------------------------------------------------------------------------
    '''

    def duplicate_rm(self, products):
        '''
        Remove duplicate filenames in a table.

        Parameters
        ----------
            products : astropy table
                (Mandatory) The data products given from astroquery.

        Returns
        -------
            An astropy table that with all duplicate elements removed.


        '''

        # Used *not in* statement from https://www.dataquest.io/blog/how-to-remove-duplicates-from-a-python-list/

        duplicates_mask = []
        non_duplicates = []
        for idx, filename in enumerate(products['productFilename']):
            if filename not in non_duplicates:
                duplicates_mask.append(True)
                non_duplicates.append(filename)
            else:
                duplicates_mask.append(False)
        
        return products[duplicates_mask]


    def download(self, pid, obs_num='None',  visit='None', visit_group='None', parallel_sequence_id='None', 
                 activity_number='None', exposure_number='None', dir_seg='None', guider='None', calib_level=2, 
                 save=False, token=None):
        '''
        This function downloads FG-GS FITS for a given program ID, observtion number, visit number, and other parameters described as parameters in the ``load`` class. The
        function also creates a working directory within the current directory to download selected FITS files and save results and plots. Additionally, the guidestar is
        processed and packaged into objects including the timeseries of guidiesstars, time and spatial arrays, and guidestar properties.

        Parameters
        ----------

            :param dir: (optional, string)
                The parameter "dir" allows the user to change the directory of the script. Changing the directory dictates where the class loads, creates folders, and save data and results.
                If ignored, "load" will save the current directory as an attribute for use with other functions.

            :param pid: (optional, int or string)
                An Python ``int`` or ``string`` of the Program ID for a given JWST program (i.e.: 1296, '1296'). Program IDs that are publicly avaliable are able to be accessed. Program IDs 
                with exclusive rights will not be accessible unless the user has a MAST API token containing the rights to the program ID. The MAST API token must be entered with param "token"
                to use programs marked with exclusive rights. The program ID needs to be entered to use all of functions in "load". However, the pid can be ignored to use specific functions
                within the class, though functionality will be severely limited. You can also forgo inputting the program ID to create a 'blank' object, then later input data for a pid using
                the function ``download`` given ``pid`` and the following parameters. For example:
                        
                        >>> spk = spelunker.load()
                        >>> spk.download(pid=1534)

            :param obs_num: (optional, int or string)
                The integer or string of the observation number of a program ID. Inputting the observation number will filter the auxiliary FITS to only include the files with the given
                number before downloading all files. Note that when you specifify the observation number, you must also input the ``pid``. 
                
            :param visit: (optional, int or string)
                The visit parameter requires a Python ``int`` or ``string``. Inputting the visit number will fillter the auxiliary FITS to only include the files with the given parameter
                before downloading all files.

            :param visit_group: (optional, int or string)
                The visit group parameter requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``vist_group`` will filter the
                data to only include files with the given visit group before loading the data into the object.

            :param parallel_sequence_id: (optional, int or string)
                The parallel sequence ID parameter requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``parallel_sequence_id`` 
                will filter the data to only include files with the given parallel sequence ID before loading the data into the object.

            :param activity_number: (optional, int or string)
                The activity number requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``activity_number`` will filter the
                data to only include files with the given activity number before loading the data into the object.

            :param exposure_number: (optional, int or string)
                The exposure_number requires either a ``int`` or ``string``. After the FITS files downloads for a corresponding program ID, specifying ``exposure_number`` will filter the
                data to only include files with the given exposure number before loading the data into the object.

            :param dir_seg: (optional, int or string)
                The segment number requires a Python ``int`` or ``string``. Specifying the segment number will include guidestar FITS with only the input number in the header if the segmented
                files has 'seg' within the naming scheme for the filename. This parameter will filter associated program ID FITS after downloading has completed and before the data is loaded
                into the object.

            :param guider: (optional, int or string)
                Similarly to `dir_seg`, the guider number requires a Python ``int`` or ``string``. Specifying the guider number will include guidestar FITS with only the input number in the 
                header if the files has 'guider' within the naming scheme for the filename. This parameter will filter associated program ID FITS after downloading has completed and before 
                the data is loaded into the object.

            :param calib_level: (optional, int or string)
                The calibration level can be specified as a Python ``int`` or ``string``. This parameter determines the calibration level for the FG-GS FITS. The only inputs for this
                parameter are 1, and 2, referring to higher stages of calibration with the JWST pipeline. Inputting the observation number will filter the auxiliary FITS to only include the
                files with the given calibration level before downloading all files. If ``calib_level`` is ignored, the calibration level will default to 2.

            :param save: (optional, boolean)
                Python ``boolean`` that either saves Gaussian fit, periodogram results, and other useful data automatically or does not save results except for the downloaded data. If 
                ignored, ``save`` will default to false.


            :param token: (optional, string)
                This parameter takes a Python ``string`` of a given MAST API token. This is required for certain functions, such as accessing JWST mnemonics and technical information, as well
                as program IDs with exclusive rights. If ignored, a blank MAST API token attribute will be created and access to JWST mnemonics will be limited and any program ID marked with
                exclusive rights will not be loaded.        

        Returns
        -------

            The ``download`` function creates attributes for frame by frame flux timeseries of a guidestar, as well as a time array, and spatial array. Other data that are
            created as attributes includes a table of the metadata for a selected program ID and object properties, which includes stellar properties of guidestars present
            within the program. Other attributes are created and recorded using the processed guidestar data, such as ``self.datamodel`` that stores the objects for the JWST 
            datamodel for each opened FITS and ``self.photometry_mask``. Example usage:

                        >>> self.fg_array; self.fg_flux; self.fg_time; self.object_properties
        
        '''

        if token is not None:

            Observations.login(token=token)

        if token is None:

            if self.mast_api_token is not None:

                Observations.login(token=self.mast_api_token)

        self.pid = pid

        if obs_num == 'None' and visit != 'None':
            raise ValueError('When a visit is identified, the obs_num needs to be identified.')

        print("Connecting with astroquery...")
        

        matched_obs = Observations.query_criteria(
            obs_collection = 'JWST',
            proposal_id=str(pid),
        )
        
        data_products = [Observations.get_product_list(obs) for obs in matched_obs]
        data_products = vstack(data_products)

        products = Observations.filter_products(data_products,
                    productType=["AUXILIARY",],
                    extension="fits",
                    productSubGroupDescription = "GS-FG",
                    dataproduct_type='image'
        )

        if obs_num != 'None' and visit != 'None':

            obs_num_col = []
            visit_col = []

            for i in products['productFilename']:
                obs_num_col.append(int(i[7:10]))
                visit_col.append(int(i[10:13]))

            products['obs_num'], products['visit'] = obs_num_col, visit_col

            mask1 = products['calib_level'] == int(calib_level)
            productsx = products[mask1]

            mask2 = productsx['obs_num'] == int(obs_num)
            productsx = productsx[mask2]

            mask3 = productsx['visit'] == int(visit)
            productsx = productsx[mask3]

        elif obs_num != 'None':

            obs_num_col = []

            for i in products['productFilename']:
                obs_num_col.append(int(i[7:10]))

            products['obs_num'] = obs_num_col

            mask1 = products['calib_level'] == int(calib_level)
            productsx = products[mask1]

            mask2 = productsx['obs_num'] == int(obs_num)
            productsx = productsx[mask2]

        else:
        
            # Filter table to only include data with calib_level = 2
            mask = products['calib_level'] == int(calib_level)
            productsx = products[mask]

        productsx = self.duplicate_rm(productsx)

        manifest = Observations.download_products(productsx, download_dir=self.directory)

        lookup_directory = self.directory+'/mastDownload/JWST/'+'**/jw0'+str(pid)+'**_gs-fg_**_cal.fits'
        fg_raw = sorted(glob.glob(lookup_directory))

        if len(fg_raw) == 0:

            print('\t Could not find any files in '+lookup_directory)

            raise Exception('No files were downloaded for program '+str(pid)+'. Two common causes of this issue are: \n'+\
                           ' 1.- You do not have exclusive access rights to see the data. If you have a MAST API, ingest it via spelunker.load(..., token = "yourtoken").\n'+\
                           ' 2.- There is no guidestar data for your program as of yet in MAST.')

        fg = []
        sliced_directory = []

        for i in fg_raw:
                fg.append(i.rsplit('/')[-1])
                sliced_directory.append(i.split('/')[-2])

        fg_table = Table()

        fg_table['filenames'] = fg
        fg_table['sliced_directory'] = sliced_directory

        obs_num_col = []
        visit_col = []

        visit_group_col = []
        parallel_sequence_id_col = []
        activity_num = []
        exposure_number_col = []
        dir_segment_col = []
        guider_col = []

        for i in fg_table['filenames']:
                obs_num_col.append(int(i[7:10]))
                visit_col.append(int(i[10:13]))

        for i in fg_table['sliced_directory']:
            visit_group_col.append(int(i[14:16]))
            parallel_sequence_id_col.append(int(i[16]))
            activity_num.append(str(i[17:19]))
            exposure_number_col.append(int(i[20:25]))

            if 'seg' in i:
                dir_segment_col.append(int(i[29:32]))
            else:
                dir_segment_col.append(0)
            if 'guider' in i:
                guider_col.append(int(i[32]))
            else:
                guider_col.append(0)

        fg_table['visit_group'] = visit_group_col
        fg_table['parallel_sequence_id'] = parallel_sequence_id_col
        fg_table['activity_number'] = activity_num
        fg_table['exposure_number'] = exposure_number_col
        fg_table['dir_seg'] = dir_segment_col
        fg_table['guider'] = guider_col        

        fg_table['obs_num'], fg_table['visit'] = obs_num_col, visit_col

        if obs_num != 'None' and visit != 'None':
                mask2 = fg_table['obs_num'] == int(obs_num)
                fg_table = fg_table[mask2]

                mask3 = fg_table['visit'] == int(visit)
                fg_table = fg_table[mask3]

        elif obs_num != 'None':
                mask2 = fg_table['obs_num'] == int(obs_num)
                fg_table = fg_table[mask2]

        if visit_group != 'None':
                mask4 = fg_table['visit_group'] == int(visit_group)
                fg_table = fg_table[mask4]         
       
        if parallel_sequence_id != 'None':
                mask5 = fg_table['parallel_sequence_id'] == int(parallel_sequence_id)
                fg_table = fg_table[mask5]

        if activity_number != 'None':
                mask6 = fg_table['activity_number'] == str(activity_number)
                fg_table = fg_table[mask6]

        if exposure_number != 'None':
                mask7 = fg_table['exposure_number'] == int(exposure_number)
                fg_table = fg_table[mask7]
        
        if dir_seg != 'None':
                mask8 = fg_table['dir_seg'] == int(dir_seg)
                fg_table = fg_table[mask8]

        elif guider != 'None':
                mask8 = fg_table['guider'] == int(guider)
                fg_table = fg_table[mask8]    

        f_slash = []
        mastDownload_dir = []
        for i in range(len(fg_table['filenames'])):
            f_slash.append('/')
            mastDownload_dir.append(self.directory+'/mastDownload/JWST/')

        reformed_directory = np.char.add(np.char.add(fg_table['sliced_directory'], f_slash), fg_table['filenames'])
        fg_table['reformed_directory'] = np.char.add(mastDownload_dir, reformed_directory)

        gs_id = []
        guidestar_time = []
        object_fg = []

        self.fg_datamodel = datamodels.open(list(fg_table['reformed_directory']))

        for file in self.fg_datamodel:
            gs_id.append(file.meta.guidestar.gs_id)
            guidestar_time.append(file.meta.guidestar.data_start)
            object_fg.append(file)
            

        fg_table['gs_id'] = gs_id
        fg_table['guidestar_time'] = guidestar_time
        fg_table['object_fg'] = object_fg
        # Reads in all fits files from reformed_directory
        #fg1 = datamodels.open(list(reformed_directory))

        _ = fg_table.sort(['guidestar_time'])

        self.fg_datamodel = list(fg_table['object_fg'])
        self.fg_table = fg_table

        self.fg_timeseries = [fn for fn in self.fg_datamodel]

        all_times = []
        all_fluxes = []
        
        for i in self.fg_datamodel:
            all_times.append( np.linspace( i.meta.guidestar.data_start, 
                               i.meta.guidestar.data_end, 
                               i.data.shape[0])
                )
            all_fluxes.append( np.nansum( i.data, axis = (1,2) ) )

        ## Preprocessing

        # https://stackoverflow.com/a/3844833
        # If all elements in a list are the same, pass.
        if len(set(fg_table['gs_id'])) != 1:

            fg_time = self.stitcher(all_times)
            fg_flux = self.stitcher(all_fluxes)
            fg_array = self.stitcher(self.fg_timeseries)

            fg_array, fg_time, fg_flux = self.negative_flux_preprocessing(fg_array, fg_time, fg_flux)

        else:

            norm_array, norm_flux = self.normalization_flux_preprocessing()

            fg_time = self.stitcher(all_times)

            fg_array, fg_time, fg_flux = self.negative_flux_preprocessing(norm_array, fg_time, norm_flux)

        data_table = Table()
        data_table['spatial'] = fg_array
        data_table['time'] = fg_time
        data_table['flux'] = fg_flux

        _ = data_table.sort(['time'])

        
        # Saving numpy arrays

        if save:

            data_arrays_dir = 'data_arrays'

            if not os.path.exists(self.directory+'/'+data_arrays_dir):

                os.makedirs(self.directory+'/'+data_arrays_dir)

            np.save(self.directory+'/'+data_arrays_dir+'/'+'pid_'+str(pid)+'_fg_array', list(data_table['spatial']))
            np.save(self.directory+'/'+data_arrays_dir+'/'+'pid_'+str(pid)+'_fg_time', list(data_table['time']))
            np.save(self.directory+'/'+data_arrays_dir+'/'+'pid_'+str(pid)+'_fg_flux', list(data_table['flux']))

        self.fg_array = data_table['spatial']
        self.fg_time = data_table['time']
        self.fg_flux = data_table['flux']
        self.photometry_mask = np.ones([data_table['spatial'].shape[1], data_table['spatial'].shape[2]])

        self.fg_table = self.table()
        self.object_properties = self.object_properties_func()

        if token is not None:

            Observations.logout()
    
    def object_properties_func(self,):
        '''
        This function neatly packages useful stellar properties of guidestars in a specified program ID from the ``self.fg_table``, created from running the ``download``
        function.

        Parameters
        ----------
        
            No function parameters, but takes the `self.fg_table`, `self.fg_datamodel`, and `self.stitcher` attributes as inputs.

        Returns
        -------

            This function returns a table with the guidestar or object properties.

        '''

        object_table = pd.DataFrame(columns=['guidestar_catalog_id','GAIAdr3sourceID', 'int_start', 'int_stop', 'ra', 'dec', 'Jmag', 'Hmag'])

        for idx, gs_id in enumerate(self.fg_table['gs_id']):

            row = []

            if gs_id not in list(object_table['guidestar_catalog_id']):

                row.append(self.fg_table['gs_id'][idx])
                row.append(self.fg_table['GAIAdr3sourceID'][idx])

                all_times = []
                all_times.append( np.linspace( self.fg_datamodel[idx].meta.guidestar.data_start, 
                                    self.fg_datamodel[idx].meta.guidestar.data_end, 
                                    self.fg_datamodel[idx].data.shape[0])            
                )
                object_time = self.stitcher(all_times)

                row.append(object_time[0])
                row.append(object_time[-1])

                row.append(self.fg_table['ra'][idx])
                row.append(self.fg_table['dec'][idx])
                row.append(self.fg_table['TmassJmag'][idx])
                row.append(self.fg_table['TmassHmag'][idx])

                object_table.loc[idx] = (row)

        return object_table.reset_index(drop=True)
    
    def table(self,):
        '''
        This function takes the guidestar ID from ``self.fg_table`` and adds essential information to a new table such as stellar properties, GAIA IDs, and magnitudes.

        Parameters
        ----------
        
            No function parameters, but takes the `self.fg_table` attribute as an input.

        Returns
        -------

            This function returns a new table that includes the information from ``self.table`` along with the added stellar properties as a ``pandas`` table.

        '''

        guidestar_id = self.fg_table['gs_id'][0]
        data = pd.read_csv('https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?id='+guidestar_id+'&format=csv', skiprows=[0])
        names=[]
        for k in data.keys():
            names.append(k)
        meta_table_df = pd.DataFrame(columns = names)
        meta_table = Table()

        for gs_id in self.fg_table['gs_id']:
            try:
                data = pd.read_csv('https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?id='+gs_id+'&format=csv', skiprows=[0])
                for k in data.keys():
                    meta_table[k] = [data[k].values[0]]
                meta_table_df.loc[len(meta_table_df)] = (list(meta_table[0]))
            except:
                print('Could not search for guidestar ID ' + str(gs_id)+'. It is probably no longer in the following catalog: https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?id='+str(gs_id)+'&format=csv')
                empty_strings = [''] * len(names)
                empty_strings[0] = str(gs_id)
                meta_table_df.loc[len(meta_table_df)] = empty_strings


        master_meta_table = Table.from_pandas(meta_table_df)

        master_table = hstack([self.fg_table, master_meta_table])
        
        return master_table
    

    def readfile(self, pid, obs_num='None', visit='None',):
        '''
        This function opens a singular FG-GS guidestar FITS and writes attributes for ``self.fg_array``, ``self.fg_time``, and ``self.flux``. Note that if a program ID
        is previously loaded using ``load``, the attributes already created will be overwritten. Currently, ``readfile`` does not create ``fg_table`` and ``object_properties``
        attributes.

        Parameters
        ----------

            :param dir: (optional, string)
                The parameter "dir" allows the user to change the directory of the script. Changing the directory dictates where the class loads, creates folders, and save data and results.
                If ignored, "load" will save the current directory as an attribute for use with other functions.

            :param pid: (optional, int or string)
                An Python ``int`` or ``string`` of the Program ID for a given JWST program (i.e.: 1296, '1296'). Program IDs that are publicly avaliable are able to be accessed. Program IDs 
                with exclusive rights will not be accessible unless the user has a MAST API token containing the rights to the program ID. The MAST API token must be entered with param "token"
                to use programs marked with exclusive rights. The program ID needs to be entered to use all of functions in "load". However, the pid can be ignored to use specific functions
                within the class, though functionality will be severely limited. You can also forgo inputting the program ID to create a 'blank' object, then later input data for a pid using
                the function ``download`` given ``pid`` and the following parameters. For example:
                        
                        >>> spk = spelunker.load()
                        >>> spk.download(pid=1534)

            :param obs_num: (optional, int or string)
                The integer or string of the observation number of a program ID. Inputting the observation number will filter the auxiliary FITS to only include the files with the given
                number before downloading all files. Note that when you specifify the observation number, you must also input the ``pid``. 
                
            :param visit: (optional, int or string)
                The visit parameter requires a Python ``int`` or ``string``. Inputting the visit number will fillter the auxiliary FITS to only include the files with the given parameter
                before downloading all files.

        Returns
        -------

            The ``readfile`` function creates attributes for frame by frame flux timeseries of a guidestar, as well as a time array, and spatial array. Other data that are
            created as attributes includes a table of the metadata for a selected program ID and object properties, which includes stellar properties of guidestars present
            within the program. Other attributes are created and recorded using the processed guidestar data, such as ``self.datamodel`` that stores the objects for the JWST 
            datamodel for each opened FITS and ``self.photometry_mask``. Example usage:

                        >>> self.fg_array; self.fg_flux; self.fg_time; self.object_properties
        '''

        pid = str(pid)

        self.pid = pid

        if obs_num == 'None':
            obs_num = '***'
        else:
            obs_num = str(obs_num).zfill(3)

            self.obs_num = obs_num

        if visit == 'None':
            visit = '***'
        else:
            visit = str(visit).zfill(3)  

            self.visit = visit  
        
        fg_raw = sorted(glob.glob(self.directory+'/mastDownload/JWST/'+'**/jw0'+pid+obs_num+visit+'**_gs-fg_**cal.fits'))

        if len(fg_raw) == 0:
            print('\t Could not find any files with these parameters.')
            raise Exception('No files were read for program '+str(pid))

        print(fg_raw)
        fg = []
        sliced_directory = []

        for i in fg_raw:
                fg.append(i.rsplit('/')[-1])
                sliced_directory.append(i.split('/')[-2])

        fg_table = Table()

        fg_table['filenames'] = fg
        print(fg_table['filenames'])
        fg_table['sliced_directory'] = sliced_directory

        obs_num_col = []
        visit_col = []

        # visit_group_col = []
        # parallel_sequence_id_col = []
        # activity_num = []
        # exposure_number_col = []
        # dir_segment_col = []
        # guider_col = []

        for i in fg_table['filenames']:
                obs_num_col.append(int(i[7:10]))
                visit_col.append(int(i[10:13]))
    
        # for i in fg_table['sliced_directory']:
        #     visit_group_col.append(int(i[14:16]))
        #     parallel_sequence_id_col.append(int(i[16]))
        #     activity_num.append(int(i[17:19]))
        #     exposure_number_col.append(int(i[20:25]))

        #     if 'seg' in i:
        #         dir_segment_col.append(int(i[29:32]))
        #     else:
        #         dir_segment_col.append(0)
        #     if 'guider' in i:
        #         guider_col.append(int(i[32]))
        #     else:
        #         guider_col.append(0)

        # fg_table['visit_group'] = visit_group_col
        # fg_table['parallel_sequence_id'] = parallel_sequence_id_col
        # fg_table['activity_number'] = activity_num
        # fg_table['exposure_number'] = exposure_number_col
        # fg_table['dir_seg'] = dir_segment_col
        # fg_table['guider'] = guider_col        

        fg_table['obs_num'], fg_table['visit'] = obs_num_col, visit_col

        if obs_num != '***' and visit != '***':
                mask2 = fg_table['obs_num'] == int(obs_num)
                fg_table = fg_table[mask2]

                mask3 = fg_table['visit'] == int(visit)
                fg_table = fg_table[mask3]

        elif obs_num != '***':
                mask2 = fg_table['obs_num'] == int(obs_num)

        f_slash = []
        mastDownload_dir = []
        for i in range(len(fg_table['filenames'])):
            f_slash.append('/')
            mastDownload_dir.append(self.directory+'/mastDownload/JWST/')

        print(fg_table['sliced_directory'])
        print(f_slash)
        print(fg_table['filenames'])

        reformed_directory = np.char.add(np.char.add(fg_table['sliced_directory'], f_slash), fg_table['filenames'])
        fg_table['reformed_directory'] = np.char.add(mastDownload_dir, reformed_directory)

        gs_id = []
        guidestar_time = []
        object_fg = []

        self.fg_datamodel = datamodels.open(list(fg_table['reformed_directory']))

        for file in self.fg_datamodel:
            gs_id.append(file.meta.guidestar.gs_id)
            guidestar_time.append(file.meta.guidestar.data_start)
            object_fg.append(file)
            

        fg_table['gs_id'] = gs_id
        fg_table['guidestar_time'] = guidestar_time
        fg_table['object_fg'] = object_fg
        # Reads in all fits files from reformed_directory
        #fg1 = datamodels.open(list(reformed_directory))

        _ = fg_table.sort(['guidestar_time'])

        self.fg_datamodel = list(fg_table['object_fg'])
        self.fg_table = fg_table

        self.fg_timeseries = [fn for fn in self.fg_datamodel]

        all_times = []
        all_fluxes = []
        
        for i in self.fg_datamodel:
            all_times.append( np.linspace( i.meta.guidestar.data_start, 
                               i.meta.guidestar.data_end, 
                               i.data.shape[0])
                )
            all_fluxes.append( np.nansum( i.data, axis = (1,2) ) )

        ## Preprocessing

        # https://stackoverflow.com/a/3844833
        # If all elements in a list are the same, pass.
        if len(set(fg_table['gs_id'])) != 1:

            fg_time = self.stitcher(all_times)
            fg_flux = self.stitcher(all_fluxes)
            fg_array = self.stitcher(self.fg_timeseries)

            fg_array, fg_time, fg_flux = self.negative_flux_preprocessing(fg_array, fg_time, fg_flux)

        else:

            norm_array, norm_flux = self.normalization_flux_preprocessing()

            fg_time = self.stitcher(all_times)

            fg_array, fg_time, fg_flux = self.negative_flux_preprocessing(norm_array, fg_time, norm_flux)

        data_table = Table()
        data_table['spatial'] = fg_array
        data_table['time'] = fg_time
        data_table['flux'] = fg_flux

        _ = data_table.sort(['time'])

        self.fg_array = data_table['spatial']
        self.fg_time = data_table['time']
        self.fg_flux = data_table['flux']
        self.photometry_mask = np.ones([data_table['spatial'].shape[1], data_table['spatial'].shape[2]])

        self.fg_table = self.table()
        self.object_properties = self.object_properties_func()

    def time_to_sec(self, fg_time):
        '''
        This function converts the values in ``fg_time`` from mjd epoch time to elapsed guidestar time in seconds.

        Parameters
        ----------

            fg_time : np.array
                (Mandatory) The guidestar time array.


        Returns
        -------

            A new numpy array ``np.array`` is returned containing the converted time values in elapsed guidestar time in seconds.

        '''
        return (fg_time - fg_time[0]) * 24 * 3600.


    '''
    FITTING TOOLS
    ------------------------------------------------------------------------------------------------------
    '''


    def gauss2d_fit(self, fg_array='None', ncpus=4, save=False):
        '''
        This function applies a spatial gaussian fit to ``self.fg_array`` for guidestar data. The Gaussian parameters include amplitude, centriods, stddev, theta,
        and offset. Amplitude represents the brightness in counts measured by the Gaussian fit, centriods are the x and y pixels measured, stddev is the standard deviations of the 
        pixels in the x and y coordinates, and theta is the orientation of the Gaussian fit. Additionally, the function ``gauss2d_fit`` will automatically apply ``ray`` parallel processing. For each frame in ``self.fg_array``, a gaussian
        (which is defined with the function ``gauss2d_fit``), is fit using ``scipy.optimize.curve_fit``.

        Parameters
        ----------

            fg_array : np.array
                (Optional) The spatial array can be changed into a user-defined input guidestar array for the fitting process.

            ncpus : int
                (Optional) The number of cores used for ray parallel processing. If ignored, the number of cpu cores will default to 4 cores.
            
            save : boolean
                (Optional) This parameter will determine whether or not to save the output Gaussian fit table as a ``.dat`` data file. If ignored, ``gauss2d_fit`` will not
                save as a data file.

        Returns
        -------

            The results from the Gaussian fitting will be saved into the following attribute as an astropy table: ``self.gaussfit_results``. This function will also return
            an astropy table.


        '''
        if type(fg_array) == str:
            if fg_array == 'None':
                fg_array = np.array(self.fg_array)

        # creates a gaussian function, from stack overflow
        # https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m

        def gaussian_2d(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
            x, y = xy
            xo = xo.astype(float)
            yo = yo.astype(float)
            a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
            b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
            c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
            g = amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                    + c*((y-yo)**2))) + offset
            return g.ravel()
        
        ray.shutdown()

        @ray.remote(num_cpus=ncpus, max_retries=3)
        def ray_curve_fit(gaussian_2d, xx, yy, datar, initial_guess):
            try:
                popt, pcov = opt.curve_fit(gaussian_2d, (xx, yy), datar, p0=initial_guess, maxfev = 2000)
                return popt
            
            except RuntimeError:
                popt = [np.nan]*7
                return popt
    
        ray.init(ignore_reinit_error=True, num_cpus = ncpus)
        
        x = np.linspace(0, 7, 8)
        y = np.linspace(0, 7, 8)
        xx, yy = np.meshgrid(x, y)

        if len(fg_array.shape) == 3 and fg_array.shape[0] != 1:

            initial_guess_main_obj = []
            for data in fg_array[0:100]:
                datar = data.ravel()
                zodical_light = np.nanmedian(data[0:3,5:8])
                coords = np.where(data==datar.max())

                initial_guess = np.array([datar.max(), int(coords[1]), int(coords[0]), 1, 1, 0, zodical_light])

                popt = ray_curve_fit.remote(gaussian_2d, xx, yy, datar, initial_guess)
                initial_guess_main_obj.append(popt)

            initial_guess_main = []

            for i in initial_guess_main_obj:
                popt2 = ray.get(i)
                initial_guess_main.append([popt2[0], popt2[1],popt2[2],popt2[3],popt2[4],popt2[5],popt2[6]])

            guess = Table(rows=initial_guess_main, names=['amplitude','x_mean', 'y_mean', 'x_stddev', 'y_stddev', 'theta', 'offset'])

            initial_guess = np.array([np.nanmedian(guess['amplitude']),np.nanmedian(guess['x_mean']),np.nanmedian(guess['y_mean']),
                                        np.nanmedian(guess['x_stddev']),np.nanmedian(guess['y_stddev']),np.nanmedian(guess['theta']),np.nanmedian(guess['offset'])])
            
            rows5_obj = []
            for idx, data in enumerate(fg_array):
                datar = data.ravel()
                zodical_light = np.nanmedian(data[0:3,5:8])
                coords = np.where(data==datar.max())

                initial_guess[6] = zodical_light
                initial_guess[1], initial_guess[2] = int(coords[1]), int(coords[0])

                popt = ray_curve_fit.remote(gaussian_2d, xx, yy, datar, initial_guess)
                rows5_obj.append(popt)

            rows5 = []
            for idx, i in enumerate(rows5_obj):
                popt3 = ray.get(i)
                rows5.append([popt3[0], popt3[1],popt3[2],popt3[3],popt3[4],popt3[5],popt3[6]])                

        elif len(fg_array.shape) == 2 and (fg_array.shape[0] & fg_array.shape[1]) == 8:
            
            data = np.copy(fg_array)
            datar = data.ravel()
            zodical_light = np.nanmedian(data[0:3, 5:8])
            coords = np.where(data==datar.max())
            initial_guess = np.array([datar.max(), int(coords[1]), int(coords[0]), 1, 1, 0, zodical_light])
            popt = ray.get(ray_curve_fit.remote(gaussian_2d, xx, yy, datar, initial_guess))

            rows5 = []
            rows5.append([popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]])

        ray.shutdown()

        table = Table(rows=rows5, names=['amplitude','x_mean', 'y_mean', 'x_stddev', 'y_stddev', 'theta', 'offset'])

        self.gaussfit_results = table
        
        if save:
            gaussfits_array_dir = 'gaussfits'
            if not os.path.exists(self.directory+'/'+gaussfits_array_dir):
                os.makedirs(self.directory+'/'+gaussfits_array_dir)
            table.write(self.directory+'/'+gaussfits_array_dir+'/'+str(self.pid)+'_'+gaussfits_array_dir+'.dat', format='ascii', overwrite=True)

        return table
    
    def quick_fit(self, fg_array='None', save=False):
        '''
        The function ``quick_fit`` fits the x and y pixel locations and standard deviations of a guidestar spatial array using centroid and variance calculations. This 
        function is 'quicker' than ``gauss2d_fit``. However, it is less accuarate and does not fit for offset or theta. The amplitude is calculated by simply measuring
        the brightest pixel in each guidestar frame.

        Parameters
        ----------

            fg_array : np.array
                (Optional) The spatial array can be changed into a user-defined input guidestar array for the fitting process.

            save : boolean
                (Optional) This parameter will determine whether or not to save the output Gaussian fit table as a ``.dat`` data file. If ignored, ``quick_fit`` will not
                save as a data file.

                
        Returns
        -------

            The results from the quick fitting will be saved into the following attribute as an astropy table: ``self.quickfit_results``. This function will also return
            an astropy table of the results.

        '''
        if type(fg_array) == str:
            if fg_array == 'None':
                fg_array = self.fg_array

        def get_centroid_and_sigma(axis, counts):
            
            centroid = np.sum(axis*counts) / np.sum(counts)
            variance = np.sum( ((axis-centroid)**2) * counts) / np.sum(counts)
            
            return centroid, np.sqrt(variance)
        
        x = np.linspace(0,7,8)
        y = np.linspace(0,7,8)

        amp = []

        cen_x = []
        std_x = []

        cen_y= []
        std_y = []

        for data in fg_array:
            cen_d, std_d = get_centroid_and_sigma(x, data[0,:])
            cen_x.append(cen_d)
            std_x.append(std_d)

        for data in fg_array:
            cen_d, std_d = get_centroid_and_sigma(y, data[:,0])
            cen_y.append(cen_d)
            std_y.append(std_d)

        for data in fg_array:
            amp_d = np.nanmax(data)
            amp.append(amp_d)


        row=[]
        quick_fit_table = Table()

        quick_fit_table['amplitude'] = np.nan_to_num(amp)
        quick_fit_table['x_mean'] = np.nan_to_num(cen_x)
        quick_fit_table['y_mean'] = np.nan_to_num(cen_y)
        quick_fit_table['x_stddev'] = np.nan_to_num(std_x)
        quick_fit_table['y_stddev'] = np.nan_to_num(std_y)
        quick_fit_table['theta'] = 0
        quick_fit_table['offset'] = 0

        self.quickfit_results = quick_fit_table
        
        if save:
            quickfits_array_dir = 'quickfits'
            if not os.path.exists(self.directory+'/'+quickfits_array_dir):
                os.makedirs(self.directory+'/'+quickfits_array_dir)
            quick_fit_table.write(self.directory+'/'+quickfits_array_dir+'/'+str(self.pid)+'_'+quickfits_array_dir+'.dat', format='ascii', overwrite=True)

        return quick_fit_table


    '''
    POSTPROCESSING METHODS 
    ------------------------------------------------------------------------------------------------------
    '''


    def minmax_gaussian_postprocessing(self, fg_array,):
        '''
        Calculates the minimum and maximum counts of a median guidestar frame using Gaussian fitting. The offset parameter will be the minimum value and
        the sum of the amplitude and offset counts will be the maximum value.


        Parameters
        ----------

            fg_array : np.array
                (Mandatory) The spatial array can be changed into a user-defined input guidestar array for the fitting process.

        Returns
        -------
            This function returns minimum and maximum float values for a median frame.
        '''

        median_fg_image = np.median(fg_array, axis=0)
        gauss_2d_table = self.gauss2d_fit(median_fg_image, ncpus=1)

        max = gauss_2d_table['amplitude'][0] + gauss_2d_table['offset'][0]
        min = gauss_2d_table['offset'][0]

        self.fg_minmax = (min, max)

        return min, max


    '''
    ANALYTICAL TOOLS 
    ------------------------------------------------------------------------------------------------------
    '''
    def timescale(self, fg_time, fg_time_sec,):
        '''
        Creates scaling to convert time from seconds to minutes, hours, or mjd based on the elapsed time.

        Parameters
        ----------

            fg_time : np.array
                (Mandatory) The guidestar time array.

            fg_time_sec : np.array
                (Mandatory) The guidestar time array in seconds.


        Returns
        -------

            A string and converted time returned based on how much time has elapsed.

        '''
        if fg_time_sec[-1] > 172800:
            # in mjd
            #ax2.plot(short_fg_time, short_fg_flux)
            xlabel = 'Time (mjd)'
            time = fg_time
            return xlabel, time
        elif fg_time_sec[-1] > 10800:
            # in hours
            #ax2.plot(short_fg_time_sec /60 /60, short_fg_flux)
            xlabel = 'Time (hours)'
            time = fg_time_sec /60 /60
            return xlabel, time
        elif fg_time_sec[-1] > 300:
            # in minutes
            #ax2.plot(short_fg_time_sec /60, short_fg_flux)
            xlabel = 'Time (mins)'
            time = fg_time_sec /60
            return xlabel, time
        else:
            # in secs
            #ax2.plot(short_fg_time_sec, short_fg_flux)
            xlabel = 'Time (sec)'
            time = fg_time_sec
            return xlabel, time

    def get_mask(self, frame, npixels):
        '''
        A specified number of pixels with the highest counts in a given frame are gathered within a mask.

        Parameters
        ----------

            frame : np.array
                (Mandatory) A guidestar frame with a shape of (8,8) is required.

        Returns
        -------

            Returns a boolean mask as a ``np.array`` type.

        '''


        # cpixels = 0
        # flat_mi = frame.flatten()
        
        # values = np.array([])
        # while cpixels < npixels:
        #     max_val = np.max(flat_mi)
        #     idx = np.where(flat_mi == max_val)[0]
        #     values = np.append(values, flat_mi[idx])
        #     flat_mi = np.delete(flat_mi,idx)
        #     cpixels += len(idx)
            
        # # With these values, now fill the mask:
        # mask = np.zeros(frame.shape)
        # mask_idxs = np.array([])
        # mask_idys = np.array([])
        
        # for i in range(len(values)):
        #     idx = np.where(frame == values[i])
        #     mask_idxs = np.append(mask_idxs,idx[0][0])
        #     mask_idys = np.append(mask_idys,idx[1][0])
        #     mask[idx] = 1.
            
        # return mask

        mask = np.zeros_like(frame.flatten())
        good_inds = np.argsort(frame.flatten(), axis=None)[-npixels:]
        mask[good_inds] = 1

        return mask.reshape(frame.shape)

    def bin_data(self, time, y, n_bin):
        '''
        Bins a timeseries given a number of bins. 

        Parameters
        ----------

            time : np.array
                (Mandatory) The time arrary that will be binned.
            
            y : np.array
                (Mandatory) The y output that will be binned to the time (usually a flux array from a timeseries).

            n_bin : int
                (Mandatory) Number of bins to bin the timeseries.

        Returns
        -------

            A binned timeseries in a tuple with elements `time_bins`, `y_bins`, and `y_err_bins`.

        '''

        # Stolen from juliet: https://github.com/nespinoza/juliet/blob/master/juliet/utils.py
        time_bins = []
        y_bins = []
        y_err_bins = []
        for i in range(0, len(time), n_bin):
            time_bins.append(np.median(time[i:i + n_bin - 1]))
            y_bins.append(np.median(y[i:i + n_bin - 1]))
            y_err_bins.append(
                np.sqrt(np.var(y[i:i + n_bin - 1])) /
                np.sqrt(len(y[i:i + n_bin - 1])))
            
        return np.array(time_bins), np.array(y_bins), np.array(y_err_bins)
    
    def timeseries_binned_plot(self, fg_time='None', fg_flux='None', start_time='None', end_time='None', n_bins=96):
        '''
        Creates a matplotlib plot of the guidestar flux timeseries along with the binned flux timeseries.

        Parameters
        ----------

            fg_time : np.array
                (Optional) The guidestar time array with the same shape as the guidestar flux array. If ignored, the fg_time will be taken from ``self.fg_time``.

            fg_flux : np.array
                (Optional) The guidestar flux array with the same shape as the guidestar time array. If ignored, the fg_flux will be taken from ``self.fg_flux``.

            start_time : float
                (Optional) A specified left bound for the time component of the timeseries. The float must in the modified julian day (mjd) format. The ``start_time``
                represents when the timeseries will start in the created plot.
            
            end_time : float
                (Optional) A specified right bound for the time component of the timeseries. The float must in the modified julian day (mjd) format. The ``end_time``
                represents when the timeseries will end in the created plot.

            n_bins : int
                (Optional) Controls the number of bins for the plotted timeseries.
        Returns
        -------

            This function will output a plot of the binned flux timeseries. Additionally, the functions returns the timeseries as a matplotlib ``axes`` object.
        '''

        if type(fg_flux) == str and type(fg_time) == str:
            if fg_time and fg_flux == 'None':
                fg_time = self.fg_time
                fg_flux = self.fg_flux

        if start_time != 'None' and end_time != 'None':
            start_time_idx = np.abs(fg_time - start_time).argmin() # https://www.geeksforgeeks.org/find-the-nearest-value-and-the-index-of-numpy-array/
            end_time_idx = np.abs(fg_time - end_time).argmin()

            fg_time = fg_time[start_time_idx:end_time_idx]
            fg_flux = fg_flux[start_time_idx:end_time_idx]
            
        fg_time_sec = self.time_to_sec(fg_time)
        xlabel, time_scaled = self.timescale(fg_time, fg_time_sec)

        ax = plt.subplot()

        norm_flux = fg_flux/ np.nanmedian(fg_flux)

        tbin, ybin, ybinerr = self.bin_data(time_scaled, norm_flux, n_bin = n_bins)
        ax.plot(time_scaled, norm_flux, color = 'black', alpha = 0.1)
        ax.errorbar(tbin, ybin, ybinerr, fmt = 'o', 
             mfc = 'white', mec = 'black', ecolor = 'black', elinewidth = 1, alpha=0.8)
        
        ax.set_ylabel('Relative flux', fontsize = self.fontsize)
        ax.set_xlabel(xlabel)
        ax.set_ylim(np.nanmean(norm_flux) - 2*np.nanstd(norm_flux), np.nanmean(norm_flux) + 2*np.nanstd(norm_flux))

        return ax
    
    def timeseries_list_plot(self, table='None', fg_time='None', start_time='None', end_time='None'):
        '''
        After Gaussian fitting with ``gauss2d_fit``, this function will plot each of the Gaussian parameters as a timeseries. This function can also handle quick fit fitting
        with the following example:
        
                    >>> spk.timeseries_list_plot(table=spk.quickfit_results)

        Parameters
        ----------

            table : np.array
                (Optional) Here, the user can input either the Gaussian fitting results of the quick fitting results as an astropy table. Note that the table has to have exactly
                seven columns with ordered columns representing amplitude, x_mean, y_mean, x_stddev, y_stddev, theta, and offset. If ignored, the fg_time will be taken from 
                ``self.gaussfit_results``.

            fg_time : np.array
                (Optional) The guidestar time array. If ignored, the fg_time will be taken from ``self.fg_time``.

            start_time : float
                (Optional) A specified left bound for the time component of the timeseries plots. The float must in the modified julian day (mjd) format. The ``start_time``
                represents when the timeseries will start in the created plot.
            
            end_time : float
                (Optional) A specified right bound for the time component of the timeseries plots. The float must in the modified julian day (mjd) format. The ``end_time``
                represents when the timeseries will end in the created plot.

        Returns
        -------

            This function will output a figure of seven timeseries plots for each of the Gaussian fitting parameters. Additionally, the functions returns the plot as a matplotlib
            ``axes`` object. Each plot within the figure can be manipulated using the ``axes`` object. For instance:

                        >>> ax[0,0].set_title('Centroid_x')
                        >>> ax[1,0].set_title('stddev_x')
                        >>> ax[0,1].set_title('Centroid_y')
                        >>> ax[1,1].set_title('stddev_x')
                        >>> ax[2,0].set_title('amplitude')
                        >>> ax[2,1].set_title('theta')
                        >>> ax[3,0].set_title('offset')
        
        '''

        if type(fg_time) == str and type(table) == str:
            if table and fg_time == 'None':
                table = self.gaussfit_results
                fg_time = self.fg_time

        if start_time != 'None' and end_time != 'None':
            start_time_idx = np.abs(fg_time - start_time).argmin()
            end_time_idx = np.abs(fg_time - end_time).argmin()
            fg_time = fg_time[start_time_idx:end_time_idx]

        fg_time_sec = self.time_to_sec(fg_time)
        xlabel, time_scaled = self.timescale(fg_time, fg_time_sec)

        fig, ax = plt.subplots(4,2, figsize = (12,21), dpi = 200)

        #ax[].plot(fg_time, cen_x, color = 'black', alpha = .4)
        ax[0,0].set_title('Centroid_x')
        ax[0,0].plot(time_scaled, table['x_mean'])
        ax[0,0].set_xlabel(xlabel)
        ax[0,0].set_ylabel('Pixel')
        ax[0,0].set_ylim(np.nanmean(table['x_mean']) - 5*np.nanstd(table['x_mean']), 
                         np.nanmean(table['x_mean']) + 5*np.nanstd(table['x_mean']))

        ax[1,0].set_title('stddev_x')
        ax[1,0].plot(time_scaled, table['x_stddev'])
        ax[1,0].set_ylim(np.nanmean(table['x_stddev']) - 1*np.nanstd(table['x_stddev']), 
                         np.nanmean(table['x_stddev']) + 5*np.nanstd(table['x_stddev']))
        ax[1,0].set_xlabel(xlabel)
        ax[1,0].set_ylabel('Pixel')

        ax[0,1].set_title('Centroid_y')
        ax[0,1].plot(time_scaled,table['y_mean'], color='orange')
        ax[0,1].set_ylim(np.nanmean(table['y_mean']) - 5*np.nanstd(table['y_mean']), 
                         np.nanmean(table['y_mean']) + 5*np.nanstd(table['y_mean']))
        ax[0,1].set_xlabel(xlabel)
        ax[0,1].set_ylabel('Pixel')
        
        ax[1,1].set_title('stddev_y')
        ax[1,1].plot(time_scaled,table['y_stddev'], color='orange')
        ax[1,1].set_ylim(np.nanmean(table['y_stddev']) - 1*np.nanstd(table['y_stddev']), 
                         np.nanmean(table['y_stddev']) + 5*np.nanstd(table['y_stddev']))
        ax[1,1].set_xlabel(xlabel)
        ax[1,1].set_ylabel('Pixel')

        ax[2,0].set_title('amplitude')
        ax[2,0].plot(time_scaled,table['amplitude'], color='blue')
        ax[2,0].set_ylim(np.nanmean(table['amplitude']) - 3*np.nanstd(table['amplitude']), 
                         np.nanmean(table['amplitude']) + 3*np.nanstd(table['amplitude']))
        ax[2,0].set_xlabel(xlabel)
        ax[2,0].set_ylabel('Counts')

        ax[2,1].set_title('theta')
        ax[2,1].plot(time_scaled,table['theta'], color='red')
        ax[2,1].set_ylim(np.nanmean(table['theta']) - 3*np.nanstd(table['theta']), 
                         np.nanmean(table['theta']) + 3*np.nanstd(table['theta']))
        ax[2,1].set_xlabel(xlabel)
        ax[2,1].set_ylabel('Radians')

        ax[3,0].set_title('offset')
        ax[3,0].plot(time_scaled,table['offset'], color='lightblue')
        ax[3,0].set_ylim(np.nanmean(table['offset']) - 3*np.nanstd(table['offset']), 
                         np.nanmean(table['offset']) + 3*np.nanstd(table['offset']))
        ax[3,0].set_xlabel('Time (mjd)')
        ax[3,0].set_ylabel('Counts')

        ax[3,1].set_visible(False)

        return ax
   
    def guidestar_plot(self,):
        '''
        Using astroplan, the guidestars within a given program ID are overplotted on reference stars from the Digitized Sky Survey (DSS). During a program ID, if multiple guidestars
        are used, a line is created between the positions of the guidestar.

        Parameters
        ----------

            No function parameters, but this function takes in `self.object_properties` as an input.


        Returns
        -------

            This function returns a matplotlib ``axes`` object of the guidestar plot. Additionally, the function outputs guidestar track plot.
        '''
        filted_properties = self.object_properties[self.object_properties['ra'].astype(bool)]

        coords = SkyCoord(filted_properties['ra'], filted_properties['dec'], unit='deg')
        target = SkyCoord(np.mean(coords.ra),np.mean(coords.dec),unit='deg')

        distance = []
        for coord in coords:
            distance.append(np.sqrt(  (target.ra.value - coord.ra.value)**2
                                + (target.dec.value - coord.dec.value)**2  ))

        fov_radius = 0.5*np.mean(distance)*u.deg + 1*np.std(distance)*u.deg
        fov_radius = 3 * u.deg if fov_radius > 3 * u.deg else fov_radius

        if fov_radius.value == 0: fov_radius = 2*u.deg # from https://github.com/GalagaBits/JWST-FGS-Spelunker/issues/19

        fig, ax1 = plt.subplots(figsize=(6,6),dpi=200)
        ax, hdu = plot_finder_image(target, survey='DSS', fov_radius=fov_radius,)

        ax1.set_axis_off()

        ax.scatter(coords.ra, coords.dec,  color='darkorange',marker='$\u25EF$', s=160, linewidth=0.5, alpha=0.4, transform=ax.get_transform('fk5'), label='guidestars')
        ax.plot(coords.ra, coords.dec,  color='orange', linewidth=1.2, alpha = 0.3, transform=ax.get_transform('fk5'), label='gs track')
        ax.text(coords.ra[0].value, coords.dec[0].value, s='start   ', horizontalalignment='right' , verticalalignment='center', weight='bold', transform=ax.get_transform('fk5'),)

        ax.set_title("Guidestar positions — "+str(self.pid))
        ax.legend()

        return ax


    def mnemonics(self, mnemonic, start, end, plot=True):
        '''
        Accesses JWST MAST for JWST technical mnemonics. This function searches for engineering telemetry from the Engineering Database (EDB) for
        a given time range in mjd. Currently, there are two mnemonics that are supported:

                    >>> 'SA_ZHGAUPST' (high gain antenna movement)
                    >>> 'INIS_FWMTRCURR' (NIS filter wheel current)

        Thousands of mnemonics are avaliable in the EDB, however only two are tested for this function. The output for each mnemonic can vary. Additionally,
        a MAST API token is required to use `mnemonics`.

        Parameters
        ----------

            mnemonic : string
                (Mandatory) This parameter contains a string with the mnemonic code for a selected engineering telemetry.
            
            start : float
                (Mandatory) The start time in mjd, which is the start of time range to search for mnemonic events.

            end : float
                (Mandatory) The end time in mjd, which is the end of time range to search for mnemonic events.

            plot : boolean
                (Optional) Turn on or off the ability to plot.

        Returns
        -------

            The function returns a matplotlib ``axes`` object for the mnemonics. For the ``SA_ZHGAUPST`` and ``INIS_FWMTRCURR``, the events are stored as a
            ``tuple`` in the respective attributes:

                        >>> self.mnemonics_event_hga
                        >>> self.mnemonics_event_nisfil

            For other mnemonics, the events will be stored within the following attribute:

                        >>> self.mnemonics_event   
        '''

        if not jwstuser_installed:
            raise ImportError("jwstuser needs to be installed before using the mnemonics function. Visit https://github.com/spacetelescope/jwstuser for installation instructions.")
        if self.mast_api_token is None:
            raise ImportError('The mast_api_token is not set. To set the token, use self.mast_api_token = token.')

        self.EDB = engdb.EngineeringDatabase(mast_api_token = self.mast_api_token)

        start, end = Time(start, format='mjd',scale = 'utc'), Time(end, format='mjd',scale = 'utc')

        action = self.EDB.timeseries(mnemonic, start.iso, end.iso)

        hga_times = action.time
        hga_t = Time(hga_times, scale='utc')

        event_code = []
        event_time = []

        if mnemonic == 'SA_ZHGAUPST':        
            # https://stackoverflow.com/questions/67324974/list-of-index-where-corresponding-elements-of-two-lists-are-same
            for ix, (x,y) in enumerate(zip(hga_t.mjd,action.value)):
                if action.value[ix] == 'MOVING':
                    if action.value[ix-1] == 'FINISHED':
                        event_code.append('START')
                        event_time.append(hga_t.mjd[ix])

                if action.value[ix] == 'FINISHED':
                    if action.value[ix-1] == 'MOVING':
                        event_code.append('STOP')
                        event_time.append(hga_t.mjd[ix-1])

            event_t = Time(event_time, format='mjd', scale='utc')

            self.mnemonics_event_hga = (event_time, event_code)

            if plot:
                if event_time != []:
                    # https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list
                    def pairwise(iterable):
                        a = iter(iterable)
                        return zip(a, a)



                    ax = plt.subplot()

                    for x_time, y_time in pairwise(event_time):
                        ax.axvspan(x_time, y_time, alpha=0.3)
                        ax.axvline(x_time,color='g',)
                        ax.axvline(y_time, color='r',)

                    ax.axvline(x_time, color='g', label='hga_start')
                    ax.axvline(y_time, color='r', label='hga_stop')

                    return ax
                elif event_time == []:
                    print('No HGA events found. Returning a blank plot.')
                    ax = plt.subplot()
                    return ax
        
        elif mnemonic == 'INIS_FWMTRCURR':

            
            action = self.EDB.timeseries(mnemonic, start.iso, end.iso)
            self.mnemonics_event_nisfil = action

            action_times = action.time
            action_t = Time(action_times, scale='utc')
            
            if plot:
                ax = plt.subplot()
                ax.plot(action_t.mjd, action.value, color='goldenrod', label='NIRISS Filter Wheel Motor Current')

                return ax
        
        else:
            
            action = self.EDB.timeseries(mnemonic, start.iso, end.iso)
            self.mnemonics_event = action

            action_times = action.time
            action_t = Time(action_times, scale='utc')

            if plot:
                ax = plt.subplot()
                ax.plot(action_t.mjd, action.value, label=mnemonic,)

                return ax     
       
    def mnemonics_local(self, mnemonic, plot=True):
        '''
        There are 'local mnemonics' from the attribute ``fg_table`` that the function accesses. This function returns information such as the guidestar ID and
        filename as a mnemonic with a timestamp for which each event has occurred. A MAST API token is not required and will not impact functionality. Currently,
        there are only two local mnemonics that are supported:

                    >>> 'GUIDESTAR' (guidestar ID)
                    >>> 'FILENAME' (filename)

        Parameters
        ----------

            mnemonic : string
                (Mandatory) A string that includes one of the two supported local mnemonics.

            plot : boolean
                (Optional) Turn on or off the ability to plot.

        Returns
        -------

            Here, the function ``mnemonics_local`` returns matplotlib ``axes`` objects. The mnemonics are also stored as ``tuples`` in attributes depending
            on which local mnemonic is specified in the parameter:

                >>> self.guidestar_mnemonics
                >>> self.filename_mnemonics 

        '''

        event_code = []
        event_time = []

        if mnemonic == 'GUIDESTAR':
            gs_id_event = self.fg_table['gs_id'][0]
            event_code.append(self.fg_table['gs_id'][0])
            event_time.append(self.fg_table['guidestar_time'][0])
        
            for idx, (time, gs_id) in enumerate(zip(self.fg_table['guidestar_time'],self.fg_table['gs_id'])):
                if gs_id != gs_id_event:
                    gs_id_event = gs_id
                    event_code.append(gs_id)
                    event_time.append(time)

            self.guidestar_mnemonics = (event_time, event_code)

            if plot:
                ax = plt.subplot()
                trans = transforms.blended_transform_factory(ax.transData, ax.transAxes) #https://stackoverflow.com/a/63153806

                for (x, y) in zip(self.guidestar_mnemonics[0], self.guidestar_mnemonics[1]):
                    ax.axvline(x, color='lightgreen',)
                    ax.text(x, 0.9, s=y, transform=trans, clip_on=True)
                ax.axvline(self.guidestar_mnemonics[0][0], color='lightgreen', label='guidestar')
                return ax
        
        elif mnemonic == 'FILENAME':
            for (time, filename) in zip(self.fg_table['guidestar_time'],self.fg_table['filenames']):
                    event_code.append(filename)
                    event_time.append(time)

            self.filename_mnemonics = (event_time, event_code)

            if plot:
                ax = plt.subplot()
                trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

                for idx, (x, y) in enumerate(zip(self.filename_mnemonics[0], self.filename_mnemonics[1])):
                    ax.axvline(x, color='brown',)

                    ax.text(x, 0.75, s=str(' '+y[0:20]),transform=trans, clip_on=True, rotation=5, wrap=True)
                    ax.text(x, 0.7, s=str(' '+y[20:42]),transform=trans, clip_on=True, rotation=5, wrap=True)

                ax.axvline(self.filename_mnemonics[0][0], color='brown', label='filename')
                    
                return ax

    def periodogram(self, table='None', time='None', save=False, plot=True):
        '''
        Creates Lomb-Scargle periodogram plots for all fitted Gaussian or quick fit parameters.

        Parameters
        ----------

            table : np.array
                (Optional) Here, the user can input either the Gaussian fitting results of the quick fitting results as an astropy table. Note that the table has 
                to have exactly seven columns with ordered columns representing amplitude, x_mean, y_mean, x_stddev, y_stddev, theta, and offset. If ignored, the 
                fg_time will be taken from ``self.gaussfit_results``.

            time : np.array
                (Optional) The guidestar time array. If ignored, the fg_time will be taken from ``self.fg_time``.

            save : boolean
                (Optional) This parameter will determine whether or not to save the output periodogram as a png. 

            plot : boolean
                (Optional) Turn on or off the ability to plot.
            
        Returns
        -------

            This function returns a matplotlib ``axes`` object. Each periodogram plot for each of the fitted parameters are accessed through the ``axes`` object.
            For example:
                        >>> ax[0].text('amplitude')
                        >>> ax[1].text('x_mean')
                        >>> ax[2].text('y_mean')
                        >>> ax[3].text('x_stddev')
                        >>> ax[4].text('y_stddev')
                        >>> ax[5].text('theta')
                        >>> ax[6].text('offset')

            Additionally, the results from the periodogram analysis can be accessed through attributes as astropy tables:

                        >>> self.pgram_results
            
            The individual parameters can be accessed through attributes for usability. 

        '''

        if type(time) == str and type(table) == str:
            if table == 'None' and time == 'None':
                table = self.gaussfit_results
                time = self.fg_time

        table['time'] = time
        table = table[~np.isnan(table['amplitude'])]

        time_t, gs_xmean, gs_ymean = table['time'].value, table['x_mean'].value, table['y_mean'].value
        gs_xstddev, gs_ystddev = table['x_stddev'].value, table['y_stddev'].value
        gs_theta, gs_amplitude = table['theta'].value, table['amplitude'].value
        gs_offset = table['offset'].value

        time_t = self.time_to_sec(time_t)

        frequency_sub0, power_sub0 = LombScargle(time_t, gs_amplitude).autopower(samples_per_peak=5)
        frequency_sub1, power_sub1 = LombScargle(time_t, gs_xmean).autopower(samples_per_peak=5)
        frequency_sub2, power_sub2 = LombScargle(time_t, gs_ymean).autopower(samples_per_peak=5)
        frequency_sub3, power_sub3 = LombScargle(time_t, gs_xstddev).autopower(samples_per_peak=5)
        frequency_sub4, power_sub4 = LombScargle(time_t, gs_ystddev).autopower(samples_per_peak=5)
        frequency_sub5, power_sub5 = LombScargle(time_t, gs_theta).autopower(samples_per_peak=5)
        frequency_sub6, power_sub6 = LombScargle(time_t, gs_offset).autopower(samples_per_peak=5)

        pgram_table = Table()

        pgram_table['frequency_amplitude'], pgram_table['power_amplitude'] = frequency_sub0, power_sub0
        pgram_table['frequency_x_mean'], pgram_table['power_x_mean'] = frequency_sub1, power_sub1
        pgram_table['frequency_y_mean'], pgram_table['power_y_mean'] = frequency_sub2, power_sub2
        pgram_table['frequency_x_stddev'], pgram_table['power_x_stddev'] = frequency_sub3, power_sub3
        pgram_table['frequency_y_stddev'], pgram_table['power_y_stddev'] = frequency_sub4, power_sub4
        pgram_table['frequency_theta'], pgram_table['power_theta'] = frequency_sub5, power_sub5
        pgram_table['frequency_offset'], pgram_table['power_offset'] = frequency_sub6, power_sub6

        self.pgram_results = pgram_table

        self.pgram_amplitude = pgram_table['frequency_amplitude'], pgram_table['power_amplitude']
        self.pgram_x_mean = pgram_table['frequency_x_mean'], pgram_table['power_x_mean']
        self.pgram_y_mean = pgram_table['frequency_y_mean'], pgram_table['power_y_mean']
        self.pgram_x_stddev = pgram_table['frequency_x_stddev'], pgram_table['power_x_stddev']
        self.pgram_y_stddev = pgram_table['frequency_y_stddev'], pgram_table['power_y_stddev']
        self.pgram_theta = pgram_table['frequency_theta'], pgram_table['power_theta']
        self.pgram_offset = pgram_table['frequency_offset'], pgram_table['power_offset']

        if plot:

            fig, ax = plt.subplots(7,1, figsize=(12,16), dpi=200)

            fig.tight_layout()

            fig.suptitle('PID ' + str(self.pid) + ' — ' + str(self.obs_num) + ' — ' + str(self.visit) + ' — ' + 
                            str(self.object_properties['guidestar_catalog_id'][0]), fontsize = self.fontsize)

            ax[0].semilogx(1/frequency_sub0, power_sub0, linewidth=.2, color='black', alpha=0.6)
            ax[1].semilogx(1/frequency_sub1, power_sub1, linewidth=.2, color='black', alpha=0.6)
            ax[2].semilogx(1/frequency_sub2, power_sub2, linewidth=.2, color='black', alpha=0.6)
            ax[3].semilogx(1/frequency_sub3, power_sub3, linewidth=.2, color='black', alpha=0.6)
            ax[4].semilogx(1/frequency_sub4, power_sub4, linewidth=.2, color='black', alpha=0.6)
            ax[5].semilogx(1/frequency_sub5, power_sub5, linewidth=.2, color='black', alpha=0.6)
            ax[6].semilogx(1/frequency_sub6, power_sub6, linewidth=.2, color='black', alpha=0.6)

            ax[0].text(0.99,0.9,'amplitude',transform=ax[0].transAxes, horizontalalignment='right', fontweight='bold')
            ax[1].text(0.99,0.9,'x_mean',transform=ax[1].transAxes, horizontalalignment='right', fontweight='bold')
            ax[2].text(0.99,0.9,'y_mean',transform=ax[2].transAxes, horizontalalignment='right', fontweight='bold')
            ax[3].text(0.99,0.9,'x_stddev',transform=ax[3].transAxes, horizontalalignment='right', fontweight='bold')
            ax[4].text(0.99,0.9,'y_stddev',transform=ax[4].transAxes, horizontalalignment='right', fontweight='bold')
            ax[5].text(0.99,0.9,'theta',transform=ax[5].transAxes, horizontalalignment='right', fontweight='bold')
            ax[6].text(0.99,0.9,'offset',transform=ax[6].transAxes, horizontalalignment='right', fontweight='bold')
            
            ax[3].set_ylabel('Power', fontsize = self.fontsize)
            ax[6].set_xlabel('Period (s)', fontsize = self.fontsize)


            if save:

                base_fname = self.object_properties['guidestar_catalog_id'][0] +\
                        '_'+str(self.object_properties['int_start'][0])
                
                fig.savefig(self.directory+'/'+base_fname+'_pgram_results.png')
            
            return ax


    '''
    ANIMATIONS 
    ------------------------------------------------------------------------------------------------------
    '''


    def timelapse_animation(self,fg_array='None', start = 0, stop = -1, interval=100, filename='movie.gif'):
        '''
        Creates an animation of a timeseries for a fg_array.

        Parameters
        ----------
        
            fg_array : np.array
                (Optional) The guidestar spatial array. If ignored, the fg_array will be taken from ``self.fg_array``.

            start : int
                (Optional) The starting index for animation to start within ``fg_array``.

            end : int
                (Optional) The ending index for animation to end within ``fg_array``.

            interval : int 
                (Optional) The delay between frames in milliseconds, according to ``matplotlib.animation.ArtistAnimation``.

        Returns
        -------

            The animation saves as a specified file, either as ``.mp4`` or ``.gif``.
        '''

        if type(fg_array) == str:
            if fg_array == 'None':
                fg_array = self.fg_array

        fig, ax = plt.subplots(dpi=200)

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        short_fg_array = fg_array[start:stop]

        # From matplotlib https://matplotlib.org/stable/gallery/animation/dynamic_image.html
        ims = []

        min, max = self.minmax_gaussian_postprocessing(fg_array)

        for i in range(len(short_fg_array)):
            im = ax.imshow(short_fg_array[i], animated=True)
            im.set_clim(vmin=min, vmax=max)
            if i == 0:
                ax.imshow(short_fg_array[i])
            ims.append([im])

        fig.suptitle('Guidestar spatial timeseries animation')
        fig.colorbar(im, label='Counts', ax = ax)
        fig.tight_layout()

        ani = animation.ArtistAnimation(fig, ims, interval=interval, blit=True,
                                        repeat_delay=1000)
        ani.save(filename)

    def flux_spatial_timelapse_animation(self,fg_array='None', fg_time='None', fg_flux='None', start = 0, stop = 100, interval=100, filename='movie.gif',):
        '''
        This function animates the spatial timeseries from ``self.fg_array``. The function will generate two plots, one with the flux timeseries, the other with
        the spatial timeseries sharing the same time array. Thus, this function will animate events in both plots simultaneously.

        Parameters
        ----------

            fg_array : np.array
                (Optional) The guidestar spatial array with the same number of elements as ``fg_time`` and ``fg_flux``. If ignored, the fg_array will be taken
                from ``self.fg_array``.

            fg_time : np.array
                (Optional) The guidestar time array with the same shape as the guidestar flux array. If ignored, the fg_time will be taken from ``self.fg_time``.

            fg_flux : np.array
                (Optional) The guidestar flux array with the same shape as the guidestar time array. If ignored, the fg_flux will be taken from ``self.fg_flux``.
            
            start : int
                (Optional) The starting index for animation to start within ``fg_array``.

            end : int
                (Optional) The ending index for animation to end within ``fg_array``.

            interval : int 
                (Optional) The delay between frames in milliseconds, according to ``matplotlib.animation.ArtistAnimation``.

            filename : string
                (Optional) A filename to save the animation. To change the file type, the extension can be changed between ``.mp4`` and ``.gif``.

        Returns
        -------

            The animation saves as a specified file, either as ``.mp4`` or ``.gif``.

        '''
        if type(fg_flux) == str and type(fg_array) == str and type(fg_time) == str:
            if fg_array == 'None' and fg_flux == 'None':
                fg_array = self.fg_array
                fg_flux = self.fg_flux
                fg_time = self.fg_time

        fig = plt.figure(figsize=(14,6), dpi=200)

        ax1 = fig.add_subplot(1,2,2)
        ax2 = fig.add_subplot(1,2,1)

        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)

        ax2.set_ylabel('Counts')
        ax2.get_xaxis().set_visible(True)
        
        short_fg_array = fg_array[start:stop]
        short_fg_flux = fg_flux[start:stop]
        short_fg_time = fg_time[start:stop]

        short_fg_time_sec = self.time_to_sec(short_fg_time)

        xlabel, short_fg_time_animated = self.timescale(short_fg_time, short_fg_time_sec)

        ax2.set_xlabel(xlabel)
            
        ims = []

        min, max = self.minmax_gaussian_postprocessing(fg_array)

        im2, = ax2.plot(short_fg_time_animated, short_fg_flux, animated=True, color='black', alpha=0.5)
        for idx, i in enumerate(short_fg_time_animated):
            
            im = ax1.imshow(short_fg_array[idx], animated=True)
            im.set_clim(vmin=min, vmax=max)

            im3 = ax2.vlines(i, np.min(short_fg_flux), np.max(short_fg_flux),  animated=True, color='red')
            ims.append([im, im2, im3])

        fig.suptitle('Guidestar spatial timeseries animation')
        fig.colorbar(im, label='Counts', ax = ax1)
        fig.tight_layout()

        ani = animation.ArtistAnimation(fig, ims, interval=interval, blit=True,
                                        repeat_delay=1000)
        
        ani.save(filename,writer='ffmpeg')


    '''  
    UTILITY FUNCTIONS 
    ------------------------------------------------------------------------------------------------------
    '''

    def save(self, suffix = None):
        '''
        This function saves all the time-series information available from the guidestar in an easy-to-read format. Right now works only for 
        the case on which one has a single guidestar:

        Parameters
        ----------

            :param suffix: (optional, string)
            String that will be added in front of the output filename.

        Returns
        -------

            This function saves outputs to a .txt file.


        '''

        header = '# spelunker guidestar outputs\n'+\
                 '# ---------------------------\n#\n' 

        for i in range( len(self.object_properties['guidestar_catalog_id']) ):

            base_fname = self.object_properties['guidestar_catalog_id'][i] +\
                         '_'+str(self.object_properties['int_start'][0]) + '.txt'

            header += '# Guidestar ID: {0:} | RA: {1:.6f}, DEC: {2:.6f}\n'.format(self.object_properties['guidestar_catalog_id'][i], \
                                                                                 self.object_properties['ra'][i], \
                                                                                 self.object_properties['dec'][i])
            if suffix is not None:

                base_fname = suffix + '_' +  base_fname

            header += '# Column 0: Time (JD-2400000.5)\n'
            header += '# Column 1: Total flux (Counts)\n'

            if self.gaussfit_results is not None:

                header += '# Column 2: Gaussian Amplitude (counts) \n'
                header += '# Column 3: Gaussian X location (pix) \n'
                header += '# Column 4: Gaussian Y location (pix) \n'
                header += '# Column 5: Gaussian stdev X (pix) \n'
                header += '# Column 6: Gaussian stdev Y (pix) \n'
                header += '# Column 7: Gaussian theta (degs)\n'
                header += '# Column 8: Background counts \n'

            fout = open(self.directory+'/'+base_fname, 'w')
            fout.write(header)

            for j in range( len(self.fg_time) ):

                if self.gaussfit_results is None:

                    fout.write('{0:.12f} {1:.5f}\n'.format(self.fg_time[j], self.fg_flux[j]))

                else:

                    fout.write('{0:.12f} {1:.5f} {2:.12f} {3:.12f} {4:.12f} {5:.12f} {6:.12f} {7:.12f} {8:.12f} \n'.format(self.fg_time[j], self.fg_flux[j], 
                                                                                                                   self.gaussfit_results['amplitude'].value[j],
                                                                                                                   self.gaussfit_results['x_mean'].value[j],
                                                                                                                   self.gaussfit_results['y_mean'].value[j],
                                                                                                                   self.gaussfit_results['x_stddev'].value[j],
                                                                                                                   self.gaussfit_results['y_stddev'].value[j],
                                                                                                                   self.gaussfit_results['theta'].value[j],
                                                                                                                   self.gaussfit_results['offset'].value[j]
                                                                                                                   ))

            fout.close()

        if self.pgram_results is not None:

            header = '# spelunker guidestar periodogram results\n'+\
                    '# ---------------------------\n#\n' 

            for i in range( len(self.object_properties['guidestar_catalog_id']) ):

                base_fname = self.object_properties['guidestar_catalog_id'][i] +\
                            '_'+str(self.object_properties['int_start'][0]) + '_pgram' +'.txt'

                header += '# Guidestar ID: {0:} | RA: {1:.6f}, DEC: {2:.6f}\n'.format(self.object_properties['guidestar_catalog_id'][i], \
                                                                                    self.object_properties['ra'][i], \
                                                                                    self.object_properties['dec'][i])
                if suffix is not None:

                    base_fname = suffix + '_' +  base_fname

                header += '# Column 0: Gaussian Amplitude frequency \n'
                header += '# Column 1: Gaussian Amplitude power \n'
                header += '# Column 2: Gaussian X location frequency \n'
                header += '# Column 3: Gaussian X location power \n'
                header += '# Column 4: Gaussian Y location frequency \n'
                header += '# Column 5: Gaussian Y location power \n'
                header += '# Column 6: Gaussian stdev X frequency \n'
                header += '# Column 7: Gaussian stdev X power \n'
                header += '# Column 8: Gaussian stdev Y frequency \n'
                header += '# Column 9: Gaussian stdev Y power \n'
                header += '# Column 10: Gaussian theta frequency\n'
                header += '# Column 11: Gaussian theta power\n'
                header += '# Column 12: Background counts frequency \n'
                header += '# Column 13: Background counts power\n'

                fout = open(self.directory+'/'+base_fname, 'w')
                fout.write(header)

                for j in range( len(self.pgram_results) ):
                        fout.write('{0:.12f} {1:.12f} {2:.12f} {3:.12f} {4:.12f} {5:.12f} {6:.12f} {7:.12f} {8:.12f}\n'.format( 
                                                                                                                    self.pgram_results['frequency_amplitude'].value[j],
                                                                                                                    self.pgram_results['power_amplitude'].value[j],
                                                                                                                    self.pgram_results['frequency_x_mean'].value[j],
                                                                                                                    self.pgram_results['power_x_mean'].value[j],
                                                                                                                    self.pgram_results['frequency_y_mean'].value[j],
                                                                                                                    self.pgram_results['power_y_mean'].value[j],
                                                                                                                    self.pgram_results['frequency_x_stddev'].value[j],
                                                                                                                    self.pgram_results['power_x_stddev'].value[j],
                                                                                                                    self.pgram_results['frequency_y_stddev'].value[j],
                                                                                                                    self.pgram_results['power_y_stddev'].value[j],
                                                                                                                    self.pgram_results['frequency_theta'].value[j],
                                                                                                                    self.pgram_results['power_theta'].value[j],
                                                                                                                    self.pgram_results['frequency_offset'].value[j],
                                                                                                                    self.pgram_results['power_offset'].value[j],
                                                                                                                    ))

                fout.close()

        if self.quickfit_results is not None:

            header = '# spelunker guidestar quickfit results\n'+\
                    '# ---------------------------\n#\n' 

            for i in range( len(self.object_properties['guidestar_catalog_id']) ):

                base_fname = self.object_properties['guidestar_catalog_id'][i] +\
                            '_'+str(self.object_properties['int_start'][0]) + '_quickfit' +'.txt'

                header += '# Guidestar ID: {0:} | RA: {1:.6f}, DEC: {2:.6f}\n'.format(self.object_properties['guidestar_catalog_id'][i], \
                                                                                    self.object_properties['ra'][i], \
                                                                                    self.object_properties['dec'][i])
                if suffix is not None:

                    base_fname = suffix + '_' +  base_fname

                header += '# Column 2: Quickfit Amplitude (counts) \n'
                header += '# Column 3: Quickfit X location (pix) \n'
                header += '# Column 4: Quickfit Y location (pix) \n'
                header += '# Column 5: Quickfit stdev X (pix) \n'
                header += '# Column 6: Quickfit stdev Y (pix) \n'
                header += '# Column 7: Quickfit theta (degs)\n'
                header += '# Column 8: Background counts \n'

                fout = open(self.directory+'/'+base_fname, 'w')
                fout.write(header)

                for j in range( len(self.quickfit_results) ):
                    fout.write('{0:.12f} {1:.5f} {2:.12f} {3:.12f} {4:.12f} {5:.12f} {6:.12f} {7:.12f} {8:.12f} \n'.format(self.fg_time[j], self.fg_flux[j], 
                                                                                                                   self.quickfit_results['amplitude'].value[j],
                                                                                                                   self.quickfit_results['x_mean'].value[j],
                                                                                                                   self.quickfit_results['y_mean'].value[j],
                                                                                                                   self.quickfit_results['x_stddev'].value[j],
                                                                                                                   self.quickfit_results['y_stddev'].value[j],
                                                                                                                   self.quickfit_results['theta'].value[j],
                                                                                                                   self.quickfit_results['offset'].value[j]
                                                                                                                   ))

                fout.close()