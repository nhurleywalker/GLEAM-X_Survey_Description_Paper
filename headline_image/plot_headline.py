#!/usr/bin/env python

__author__ = "Natasha Hurley-Walker"
__date__ = "15/10/2018"

import os
import sys
import shutil
import glob

import matplotlib
matplotlib.use('Agg') # So does not use display -- only good if just making plots
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8})

cm = 1/2.54

from astropy.io import fits
from astropy import wcs
#from astropy.table import Table, Column
from astropy import units as u

from astropy.coordinates import SkyCoord
#from astropy.visualization import PercentileInterval
from astropy.visualization import AsinhStretch

from astropy.nddata import Cutout2D

import numpy as np

def normalize(arr, vmin, vmax):
    nor = (arr - vmin) / (vmax - vmin)
    nor[np.where(nor<0.0)] = 0.0
    nor[np.where(nor>1.0)] = 1.0
    return nor

stretch = AsinhStretch(a=0.1)

freqs = ["072-103", "103-134", "139-170", "170-231"]
fits_files = ["headline_{0}MHz.fits".format(x) for x in freqs]
hdus = [fits.open(ff) for ff in fits_files]
#ras = ["11:30", "08:30", "05:30"]
ras = ["05:30", "08:30", "11:30"]
#ras = ["10:45", "06:15"]

# Two-column portrait figure, whole page
fig = plt.figure(figsize=(19*cm,27*cm))

# Top plot first, work our way down
# Start at bottom, work your way up
ystart = 0.1
for ra in ras:
    pos = SkyCoord(ra, "-26:42:00", unit=(u.hour, u.deg), frame="fk5")
    cutouts = [Cutout2D(hdu[0].data, pos, (12*u.deg, 47*u.deg), wcs=wcs.WCS(hdu[0].header)) for hdu in hdus]
    rgb = [normalize(cutout.data, -0.009, 0.2) for cutout in cutouts[:-1]]
    d_rgb = np.dstack(rgb)
    d_w = stretch(normalize(cutouts[-1].data, -0.003, 0.2))
#    ax_rgb = fig.add_subplot(6,1,n, projection=cutouts[0].wcs)
#    ax_w = fig.add_subplot(6,1,n+1, projection=cutouts[-1].wcs, sharex=ax_rgb)
    xstart, xsize, ysize = 0.1, 0.9, 0.1
    ax_w = fig.add_axes([xstart, ystart, xsize, ysize], projection = cutouts[-1].wcs)
    ax_rgb = fig.add_axes([xstart, ystart+0.1, xsize, ysize], projection = cutouts[0].wcs)
    img_rgb = ax_rgb.imshow(stretch(d_rgb), origin="lower")
    img_w = ax_w.imshow(stretch(d_w), origin="lower", cmap="gray")
    lat = ax_rgb.coords["ra"]
    lat.set_ticklabel_visible(False)
    lat.set_ticks_position('b')
    lat.set_axislabel("")
    lat = ax_w.coords["ra"]
    lat.set_major_formatter('hh')
    if ystart == 0.1:
        lat.set_axislabel("Right Ascension")
    else:
        lat.set_axislabel(" ")
    for ax in ax_w, ax_rgb:
        lon = ax.coords['dec']
        lon.set_axislabel("Declination")
        lon.set_major_formatter('dd')
    ystart += 0.23
#    n += 2

#fig.tight_layout()
fig.savefig("headline_gleamx.pdf", bbox_inches="tight", dpi=900)
