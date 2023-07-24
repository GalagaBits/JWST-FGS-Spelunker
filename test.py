import os
import sys
os.chdir('/Users/ddeal/JWST-Treasure-Chest/')
sys.path.append('/Users/ddeal/JWST-FGS-Spelunker/JWST-FGS-spk-main/src/')

import matplotlib.pyplot as plt
import spelunker

spk = spelunker.load('/Users/ddeal/JWST-Treasure-Chest/', pid=1534)

fig, ax = plt.subplots(figsize=(6,6))
ax.plot(spk.fg_time, spk.fg_flux)
plt.show()