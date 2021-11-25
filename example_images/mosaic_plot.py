#!/usr/bin/python

import numpy
#tables and votables
from astropy.io.votable import parse_single_table
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy import units as u

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8})

cm = 1/2.54

# Not used at the moment
mwa_bmaj=0.0208 # degrees
mwa_bmin=0.0165 # degrees
mwa_bpa=-31.88 # degrees

mosaic="XG.fits"
bkg="XG_bkg.fits"
rms="XG_rms.fits"
glm="GLEAM_rg.fits"

vmin_mos = -0.003e3
vmax_mos = 0.01e3
vmin_glm = -0.02e3
vmax_glm = 0.06e3
vmin_bkg = -0.0002e3
vmax_bkg = 0.0002e3
vmin_rms = 0.0009e3
vmax_rms = 0.0013e3

c = SkyCoord("10:30:00 -27:30:00", frame="fk5", unit=(u.hour, u.deg))
framesize = 3.5*u.deg

#cmap="cubehelix"
#cmap="Greys_r"

#table = parse_single_table("XG_mosaic.vot")
#my_data = table.array

# Set up the giant wrapper fits figure
fig = plt.figure(figsize=(19*cm,8.5*cm))

# Left-most panel: GLEAM
hdu_glm = fits.open(glm)
w_glm = wcs.WCS(hdu_glm[0].header)
cutout_glm = Cutout2D(1000*hdu_glm[0].data, c, framesize, w_glm)
ax_glm = fig.add_axes([0.085, 0.085, 0.3, 0.85], projection = cutout_glm.wcs)
im_glm = ax_glm.imshow(cutout_glm.data, origin="lower", vmin = vmin_glm, vmax = vmax_glm)
cax_glm = fig.add_axes([0.095, 0.86, 0.3-0.015, 0.015])
cb_glm = plt.colorbar(im_glm, cax = cax_glm, orientation="horizontal")
cb_glm.ax.xaxis.set_ticks_position('top')
cb_glm.ax.xaxis.set_label_position('top')

# Next panel: attractive slice of GLEAM-X
hdu_mos = fits.open(mosaic)
w_mos = wcs.WCS(hdu_mos[0].header)
cutout_mos = Cutout2D(1000*hdu_mos[0].data, c, framesize, w_mos)
ax_mos = fig.add_axes([0.395,0.085,0.3,0.85], projection = cutout_mos.wcs)
im_mos = ax_mos.imshow(cutout_mos.data, origin="lower", vmin = vmin_mos, vmax = vmax_mos)
cax_mos = fig.add_axes([0.4, 0.86, 0.3-0.015, 0.015])
cb_mos = plt.colorbar(im_mos, cax = cax_mos, orientation="horizontal")
cb_mos.ax.xaxis.set_ticks_position('top')
cb_mos.ax.xaxis.set_label_position('top')

## Top right panel: background
hdu_bkg = fits.open(bkg)
w_bkg = wcs.WCS(hdu_bkg[0].header)
cutout_bkg = Cutout2D(1000*hdu_bkg[0].data, c, framesize, w_bkg)
ax_bkg = fig.add_axes([0.74,0.6,0.2,0.3], projection = cutout_bkg.wcs)
im_bkg = ax_bkg.imshow(cutout_bkg.data, origin="lower", vmin = vmin_bkg, vmax = vmax_bkg)
cax_bkg = fig.add_axes([0.92, 0.6, 0.008, 0.3])
cb_bkg = plt.colorbar(im_bkg, cax = cax_bkg)
#
## Bottom right panel: RMS
hdu_rms = fits.open(rms)
w_rms = wcs.WCS(hdu_rms[0].header)
cutout_rms = Cutout2D(1000*hdu_rms[0].data, c, framesize, w_rms)
ax_rms = fig.add_axes([0.74,0.15,0.2,0.3], projection = cutout_rms.wcs)
im_rms = ax_rms.imshow(cutout_rms.data, origin="lower", vmin = vmin_rms, vmax = vmax_rms)
cax_rms = fig.add_axes([0.92, 0.15, 0.008, 0.3])
cb_rms = plt.colorbar(im_rms, cax = cax_rms)

#
lon = ax_glm.coords['dec']
lon.set_axislabel("Declination")
for ax in ax_glm, ax_mos:
    lat = ax.coords['ra']
    lat.set_axislabel("Right Ascension")
lon = ax_mos.coords['dec']
lon.set_axislabel(" ")
#lon.set_ticks_visible(False)
lon.set_ticklabel_visible(False)

for ax in ax_bkg, ax_rms:
    lon = ax.coords['dec']
    lon.set_axislabel("Dec")
#    lon.set_major_formatter('dd')
    lat = ax.coords['ra']
    lat.set_axislabel("RA")
# RA ticks overlap if left alone
    lat.set_ticks(number=3)


for cb in cb_glm, cb_mos, cb_bkg, cb_rms:
    cb.set_label("Brightness / mJy beam$^{-1}$")#, labelpad=-2)
#for cb in cb_bkg, cb_rms:
#    cb.set_label("Flux density / mJy beam$^{-1}$")

fig.savefig("XG_mosaic.pdf", dpi=900, bbox_inches='tight')
fig.savefig("XG_mosaic.png", dpi=900, bbox_inches='tight')

