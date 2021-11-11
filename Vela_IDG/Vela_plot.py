#!/usr/bin/python

import numpy
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['serif']})

#tables and votables
from astropy.io.votable import parse_single_table
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

mwa_bmaj=0.0208 # degrees
mwa_bmin=0.0165 # degrees
mwa_bpa=-31.88 # degrees

gleam="Vela_GLEAM.fits"
gleamx="Vela_GLEAM-X.fits"
idg="Vela_IDG.fits"

vmin_gleam = -0.1
vmax_gleam = 1.0
vmin_gleamx = -0.05
vmax_gleamx = 0.3
vmin_idg = -0.01
vmax_idg = 0.2

cmap="cubehelix"

# Unused colorbar stuff
#ax_gleam = fig.add_axes([0.025,0.1,0.525,0.85], projection = cutout_gleam.wcs)
#cax_gleam = fig.add_axes([0.515, 0.1, 0.015, 0.85])
#cb_gleam = plt.colorbar(im_gleam, cax = cax_gleam)
# Set up the giant wrapper fits figure
fig = plt.figure(figsize=(15,5))

# Left panel: GLEAM
hdu_gleam = fits.open(gleam)
ax_gleam = fig.add_subplot(131, projection = wcs.WCS(hdu_gleam[0].header, naxis=2))
im_gleam = ax_gleam.imshow(hdu_gleam[0].data, origin="lower", vmin = vmin_gleam, vmax = vmax_gleam, cmap = cmap)

# Middle panel: GLEAM-X
hdu_gleamx = fits.open(gleamx)
ax_gleamx = fig.add_subplot(132, projection = wcs.WCS(hdu_gleamx[0].header))
im_gleamx = ax_gleamx.imshow(hdu_gleamx[0].data, origin="lower", vmin = vmin_gleamx, vmax = vmax_gleamx, cmap = cmap)

# Right panel: IDG
hdu_idg = fits.open(idg)
ax_idg = fig.add_subplot(133, projection = wcs.WCS(hdu_idg[0].header, naxis=2))
im_idg = ax_idg.imshow(hdu_idg[0].data, origin="lower", vmin = vmin_idg, vmax = vmax_idg, cmap = cmap)


lon = ax_gleam.coords['dec']
lon.set_axislabel("Declination")
for ax in ax_gleamx, ax_idg:
    lon = ax.coords['dec']
    lon.set_axislabel("")
#    lon.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)
#    lon.set_major_formatter('dd')
for ax in ax_gleam, ax_gleamx, ax_idg:
    lat = ax.coords['ra']
    lat.set_axislabel("Right Ascension")

plt.subplots_adjust(wspace=0, hspace=0)

fig.savefig("Vela.pdf",bbox_inches='tight', dpi=300)

