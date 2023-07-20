import re
import numpy
from setuptools import setup, Extension, find_packages

VERSIONFILE='src/_version.py'
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(name='spelunker',
      version=verstr,
      description='spelunker: a library to extract guidestar data and observe technical and stellar events',
      url='https://github.com/GalagaBits/JWST-FGS-Spelunker',
      author='Derod Deal',
      author_email='dealderod@ufl.edu',
      license='MIT',
      packages=['spelunker'],
      package_dir={'spelunker': 'src'},
      requires=['numpy','scipy', 'pandas', 'jwst', 'astroquery', 'astropy','astroplan','matplotlib', 'ray'],
      install_requires=['numpy','scipy', 'pandas', 'jwst', 'astroquery', 'astropy','astroplan','matplotlib', 'ray', 'wheel'],
      keywords='guidestars',
      python_requires='>=3.0',
      zip_safe=False)