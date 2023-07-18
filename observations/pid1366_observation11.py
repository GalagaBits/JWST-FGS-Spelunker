from astropy.io import fits
import fgs_spelunker_ver_0405

spk = fgs_spelunker_ver_0405.fgs_splelunker('/Users/ddeal/JWST-Treasure-Chest-2023/', 1366)

pid1366_array, pid1366_time, pid1366_flux = spk.download('1366', obs_num=11,)

table_gauss2d_fit = spk.gauss2d_fit(pid1366_array, ncpus=6) # ncpus sets the number of cpu cores your computer has. Defaults to 4 cores.
table_gauss2d_fit.write('/Users/ddeal/JWST-Treasure-Chest-2023/pid1366_observation11_01.dat', format='ascii', overwrite=True)