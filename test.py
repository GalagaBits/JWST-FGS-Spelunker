import os
import sys

import spelunker

spk = spelunker.load('/Users/ddeal/JWST-Treasure-Chest/', pid=1534)

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6,6))
ax = spk.timeseries_binned_plot()

plt.show()