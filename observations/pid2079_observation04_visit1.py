from astropy.io import fits
import fgs_spelunker_ver_0405
import numpy as np


spk = fgs_spelunker_ver_0405.fgs_splelunker('/Users/ddeal/JWST-Treasure-Chest-2023/', 2079)

pid2079_array, pid2079_time, pid2079_flux = spk.download('2079', obs_num=4, exposure_number=1, visit=1)

table_gauss2d_fit = spk.gauss2d_fit(pid2079_array, ncpus=6) # ncpus sets the number of cpu cores your computer has. Defaults to 4 cores.
table_gauss2d_fit.write('/Users/ddeal/JWST-Treasure-Chest-2023/pid2079_observation04_visit1.dat', format='ascii', overwrite=True)

np.save("pid2079_visit1_array",pid2079_array)
np.save("pid2079_visit1_time", pid2079_time)
np.save("pid2079__visit1_flux", pid2079_flux)