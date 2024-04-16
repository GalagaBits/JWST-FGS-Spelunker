# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [1.1.10] - 2024-04-16 (Current release)

- Added ability to control number of bins in `timeseries_binned_plot.
- Changed Python requirements.

## [1.1.9] - 2024-03-19

- Changes and improvements to `spelunker.load.guidestar_plot`.

## [1.1.8] - 2024-03-19

- Added `if fov_radius.value == 0: fov_radius = 2*u.deg` to `spelunker.load.guidestar_plot` as recommended by @hdiamondlowe from issue #19.
- Added GitHub actions for `general_test.py`.

## [1.1.7] - 2024-03-19

- Added `jwstuser` to `setup.py`.

## [1.1.6] - 2024-03-18

- Fixed an issue where the searching the guidestar catalog with certain guidestars IDs will not be found, thus breaking `spelunker.load`. 

## [1.1.5] - 2024-03-18

- Updates to `object_properties_func`.

## [1.1.4] - 2024-03-18

- Updates to `object_properties_func`.

## [1.1.3] - 2024-03-17

- Updates to `object_properties_func`.

## [1.1.2] - 2024-03-17 

- Updates to `object_properties_func`.

## [1.1.1] - 2024-03-17 
- Changed `python_requires` to 3.10.
- Removed dashboard.py
- Small optimizations to the source code.
- Added jwstuser to dependencies list.
- Changes to object properties table: Instead of GAIAdr1sourceID, and GAIAdr2sourceID, only GAIAdr3sourceID will be listed. There were some changes in the source Spelunker uses to get the GAIA ID, that being https://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx... The parameters GAIAdr1sourceID and GAIAdr2sourceID no longer exist and has been replaced by GAIAdr3sourceID. 

- Updated installation guide to include `jwstuser` instructions.
- Changed the readthdocs theme to be more user friendly


## [1.0.3] - (Unreleased)
### Added
- Added option to define field-of-view for `guidestar_plot()`.
- Added test in `testing` for the `guidestar_plot()` function.
- Added the `traceback` library for general checks in `testing` which spits the full traceback.

### Fixed
- Fixed bug on `guidestar_plot()` which didn't work for single targets.

## [1.0.2] - 2023-09-01 
### Added
- Photometry can now be optimized with a mask. `spk.optimize_photometry()` does the trick via an automatic mask generation; 
  a manual mask can also be passed instead via `spk.optimize_photometry(mask = mask)`. TODO: handle optimal masks for 
  multiple guidestars.
- Added test function in `testing`, useful for testing before releases/updates. 

### Fixed
- Fixed bug introduced when storing `data_arrays` in the `spelunker_outputs` directory.
- Fixed bug introduced by directory changes.

## [1.0.1] - 2023-08-02 
- Error message when data is not downloaded either to not having propietary access period or no data.
- Possibility to download data with exclusive access period by specifying a `token` in `spelunker.load`. This is a string.
- Time is now sorted by guidestar epochs.
- Spelunker now does not use `os.chdir`, so the directory will not change unexpectedly while using Jupyter notebooks.

## [1.0.0] - 2023-07-23
### Added

- Added better file management systems. Now `spelunker` will create directories for data processing and downloads.
- Added a scalable x axis that changes timescales based on the input time array.
- You can specify a directory for `spelunker.load` with the `dir` parameter. Specifying a directory is now optional and you can call `spelunker.load` without any parameters.
- Added a `save` parameter for `spelunker.load` to toggle saving the `fg_array`, `fg_flux`, and `fg_time` arrays as a `.npy` file.
- Added saving methods for `gauss2d_fit`, `quick_fit`, and  with the `save` parameter.
- Added try and except statements to `gauss2d_fit` to prevent code crashing if a fit fails
- Added `save` method to save Gaussian fit results and object properties to a `txt` file.

### Fixed
- Fixed initial guesses of gaussian fits in `gauss2d_fit` to handle guidestar changes in position. Gaussian fits should run faster and be more reliable.


## [0.5.5] - 2023-07-20

- Changes to setup.py. Now `pip install spelunker` should run more properly.

## [0.5.4] - 2023-07-19

- Added changes to testpypi

## [0.5.3] - 2023-07-19

### Added

- `spelunker` is added to PyPI!
