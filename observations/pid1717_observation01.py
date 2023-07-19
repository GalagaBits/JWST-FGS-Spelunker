import sys
sys.path.append('/Users/ddeal/JWST-FGS-Spelunker/JWST-FGS-Spelunker-main/')
import Spelunker
from astropy.io import fits

spk = Spelunker.load('/Users/ddeal/JWST-Treasure-Chest/',pid=1717, obs_num=1)
table_gauss2d_fit = spk.gauss2d_fit(spk.fg_array, ncpus=6) # ncpus sets the number of cpu cores your computer has. Defaults to 4 cores.
table_gauss2d_fit.write('/Users/ddeal/JWST-Treasure-Chest/pid1717_observation01.dat', format='ascii', overwrite=True)