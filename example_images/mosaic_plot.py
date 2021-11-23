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

vmin_mos = -0.003e3
vmax_mos = 0.01e3
vmin_bkg = -0.0002e3
vmax_bkg = 0.0002e3
vmin_rms = 0.0009e3
vmax_rms = 0.0013e3

c = SkyCoord("10:30:00 -27:30:00", frame="fk5", unit=(u.hour, u.deg))
framesize = 3.5*u.deg

cmap="cubehelix"
cmap="Greys_r"

#table = parse_single_table("XG_mosaic.vot")
#my_data = table.array

# Set up the giant wrapper fits figure
fig = plt.figure(figsize=(17*cm,8.5*cm))

# Left panel: attractive slice of mosaic
hdu_mos = fits.open(mosaic)
w_mos = wcs.WCS(hdu_mos[0].header)
cutout_mos = Cutout2D(1000*hdu_mos[0].data, c, framesize, w_mos)
ax_mos = fig.add_axes([0.025,0.1,0.525,0.85], projection = cutout_mos.wcs)
im_mos = ax_mos.imshow(cutout_mos.data, origin="lower", vmin = vmin_mos, vmax = vmax_mos)
cax_mos = fig.add_axes([0.515, 0.1, 0.015, 0.85])
cb_mos = plt.colorbar(im_mos, cax = cax_mos)

# Top right panel: background
hdu_bkg = fits.open(bkg)
w_bkg = wcs.WCS(hdu_bkg[0].header)
cutout_bkg = Cutout2D(1000*hdu_bkg[0].data, c, framesize, w_bkg)
ax_bkg = fig.add_axes([0.63,0.6,0.24,0.35], projection = cutout_bkg.wcs)
im_bkg = ax_bkg.imshow(cutout_bkg.data, origin="lower", vmin = vmin_bkg, vmax = vmax_bkg)
cax_bkg = fig.add_axes([0.85, 0.6, 0.015, 0.35])
cb_bkg = plt.colorbar(im_bkg, cax = cax_bkg)

# Bottom right panel: RMS
hdu_rms = fits.open(rms)
w_rms = wcs.WCS(hdu_rms[0].header)
cutout_rms = Cutout2D(1000*hdu_rms[0].data, c, framesize, w_rms)
ax_rms = fig.add_axes([0.63,0.1,0.24,0.35], projection = cutout_rms.wcs)
im_rms = ax_rms.imshow(cutout_rms.data, origin="lower", vmin = vmin_rms, vmax = vmax_rms)
cax_rms = fig.add_axes([0.85, 0.1, 0.015, 0.35])
cb_rms = plt.colorbar(im_rms, cax = cax_rms)

lon = ax_mos.coords['dec']
lon.set_axislabel("Declination")
lat = ax_mos.coords['ra']
lat.set_axislabel("Right Ascension")
for ax in ax_bkg, ax_rms:
    lon = ax.coords['dec']
    lon.set_axislabel("Dec")
#    lon.set_major_formatter('dd')
    lat = ax.coords['ra']
    lat.set_axislabel("RA")

cb_mos.set_label("Flux density / mJy beam$^{-1}$", labelpad=-2)
for cb in cb_bkg, cb_rms:
    cb.set_label("Flux density / mJy beam$^{-1}$")

fig.savefig("XG_mosaic.pdf",bbox_inches='tight')

