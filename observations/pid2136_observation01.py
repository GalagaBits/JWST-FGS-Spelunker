import spelunker.spelunker as spelunker
from astropy.io import fits

import sys
sys.path.append('/Users/ddeal/JWST-FGS-Spelunker/JWST-FGS-Spelunker-main/')

spk = spelunker.load('/Users/ddeal/JWST-Treasure-Chest/',pid=2136,obs_num=4,visit=2)

table_gauss2d_fit = spk.gauss2d_fit(spk.fg_array, ncpus=6) # ncpus sets the number of cpu cores your computer has. Defaults to 4 cores.
table_gauss2d_fit.write('/Users/ddeal/JWST-Treasure-Chest/pid2136_observation04_visit2.dat', format='ascii', overwrite=True)