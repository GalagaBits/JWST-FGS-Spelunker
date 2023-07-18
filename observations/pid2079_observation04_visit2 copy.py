import sys
sys.path.append('/Users/ddeal/JWST-Treasure-Chest-2023/JWST-FGS-Spelunker-Repos/JWST-FGS-Spelunker/')
import Spelunker
from astropy.io import fits

spk = Spelunker.load('/Users/ddeal/JWST-Treasure-Chest-2023/',pid=1717)
table_gauss2d_fit = spk.gauss2d_fit(spk.fg_array, ncpus=6) # ncpus sets the number of cpu cores your computer has. Defaults to 4 cores.
table_gauss2d_fit.write('/Users/ddeal/JWST-Treasure-Chest-2023/pid2079_observation04_visit2.dat', format='ascii', overwrite=True)