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

opt = "optional_filter.fits"
mand = "mandatory_filter_v3.fits"
vmin_opt = -0.01e3
vmax_opt = 0.01e3
vmin_mand = -0.01e3
vmax_mand = 0.05e3

pos_not_removed = fits.open("XG_170-231MHz_comp_rescaled_filtered.fits")[1].data
neg_not_removed = fits.open("XG_170-231MHz_psf_comp_negative_dafix_comp_filtered.fits")[1].data
pos_single_removed = fits.open("XG_170-231MHz_psf_comp_positive_dafix_comp_removed.fits")[1].data
neg_single_removed = fits.open("XG_170-231MHz_psf_comp_negative_dafix_comp_removed.fits")[1].data
neg_double_removed = fits.open("XG_170-231MHz_psf_comp_negative_dafix_comp_removed_by_both_filters.fits")[1].data

# Sensible color / symbol scheme
# positive sources can be black or red
# negative sources can be white or red
# loads of markers to choose from
# crosses make sense for sources that have been removed ("cross out")
# maybe a plus for the optional filter
# so normal sources are just circles

# Showing the optional removal of negative sources next to faint positive sources
fig = plt.figure(figsize=(8.5*cm,6.5*cm))
hdu_opt = fits.open(opt)
w_opt = wcs.WCS(hdu_opt[0].header)
ax_opt = fig.add_axes([0.15, 0.15, 0.7, 0.7], projection = w_opt)
im_opt = ax_opt.imshow(1000*hdu_opt[0].data, origin="lower", vmin = vmin_opt, vmax = vmax_opt)
cax_opt = fig.add_axes([0.8, 0.15, 0.03, 0.7])
cb_opt = plt.colorbar(im_opt, cax = cax_opt, orientation="vertical")
cb_opt.set_label("Brightness / mJy beam$^{-1}$")#, labelpad=-2)
lon = ax_opt.coords['dec']
lon.set_axislabel("Declination")
lat = ax_opt.coords['ra']
lat.set_axislabel("Right Ascension")
ax_opt.autoscale(False)
ax_opt.scatter(pos_not_removed["ra"], pos_not_removed["dec"], transform=ax_opt.get_transform("world"), marker="o", color="k", lw=0.5, s=5, facecolor="none")
ax_opt.scatter(neg_single_removed["ra"], neg_single_removed["dec"], transform=ax_opt.get_transform("world"), marker="+", color="w", lw=0.5, s=5)
fig.savefig("Optional_filter_example.pdf", bbox_inches = "tight")

# Showing the removal of artefacts around bright sources
fig = plt.figure(figsize=(8.5*cm,6.5*cm))
hdu_mand = fits.open(mand)
w_mand = wcs.WCS(hdu_mand[0].header)
ax_mand = fig.add_axes([0.15, 0.15, 0.7, 0.7], projection = w_mand)
im_mand = ax_mand.imshow(1000*hdu_mand[0].data, origin="lower", vmin = vmin_mand, vmax = vmax_mand)
cax_mand = fig.add_axes([0.8, 0.15, 0.03, 0.7])
cb_mand = plt.colorbar(im_mand, cax = cax_mand, orientation="vertical")
cb_mand.set_label("Brightness / mJy beam$^{-1}$")#, labelpad=-2)
lon = ax_mand.coords['dec']
lon.set_axislabel("Declination")
lat = ax_mand.coords['ra']
lat.set_axislabel("Right Ascension")
lat.set_major_formatter("hh:mm")
ax_mand.autoscale(False)
ax_mand.scatter(pos_not_removed["ra"], pos_not_removed["dec"], transform=ax_mand.get_transform("world"), marker="o", color="k", lw=0.5, s=5, facecolor="none")
ax_mand.scatter(neg_not_removed["ra"], neg_not_removed["dec"], transform=ax_mand.get_transform("world"), marker="o", color="w", lw=0.5, s=5, facecolor="none")
ax_mand.scatter(pos_single_removed["ra"], pos_single_removed["dec"], transform=ax_mand.get_transform("world"), marker="x", color="k", lw=0.5, s=5)
#ax_mand.scatter(neg_single_removed["ra"], neg_single_removed["dec"], transform=ax_mand.get_transform("world"), marker="+", color="r", lw=0.5, s=5)
ax_mand.scatter(neg_double_removed["ra"], neg_double_removed["dec"], transform=ax_mand.get_transform("world"), marker="x", color="w", lw=0.5, s=5)
fig.savefig("Mandatory_filter_example.pdf", bbox_inches = "tight")

# Next panel: attractive slice of GLEAM-X
#hdu_mos = fits.open(mosaic)
#w_mos = wcs.WCS(hdu_mos[0].header)
#cutout_mos = Cutout2D(1000*hdu_mos[0].data, c, framesize, w_mos)
#ax_mos = fig.add_axes([0.395,0.085,0.3,0.85], projection = cutout_mos.wcs)
#im_mos = ax_mos.imshow(cutout_mos.data, origin="lower", vmin = vmin_mos, vmax = vmax_mos)
#cax_mos = fig.add_axes([0.4, 0.86, 0.3-0.015, 0.015])
# Next panel: attractive slice of GLEAM-X
#hdu_mos = fits.open(mosaic)
#w_mos = wcs.WCS(hdu_mos[0].header)
#cutout_mos = Cutout2D(1000*hdu_mos[0].data, c, framesize, w_mos)
#ax_mos = fig.add_axes([0.395,0.085,0.3,0.85], projection = cutout_mos.wcs)
#im_mos = ax_mos.imshow(cutout_mos.data, origin="lower", vmin = vmin_mos, vmax = vmax_mos)
#cax_mos = fig.add_axes([0.4, 0.86, 0.3-0.015, 0.015])
#cb_mos = plt.colorbar(im_mos, cax = cax_mos, orientation="horizontal")
#cb_mos.ax.xaxis.set_ticks_position('top')
#cb_mos.ax.xaxis.set_label_position('top')
#
### Top right panel: background
#hdu_bkg = fits.open(bkg)
#w_bkg = wcs.WCS(hdu_bkg[0].header)
#cutout_bkg = Cutout2D(1000*hdu_bkg[0].data, c, framesize, w_bkg)
#ax_bkg = fig.add_axes([0.74,0.6,0.2,0.3], projection = cutout_bkg.wcs)
#im_bkg = ax_bkg.imshow(cutout_bkg.data, origin="lower", vmin = vmin_bkg, vmax = vmax_bkg)
#cax_bkg = fig.add_axes([0.92, 0.6, 0.008, 0.3])
#cb_bkg = plt.colorbar(im_bkg, cax = cax_bkg)
##
### Bottom right panel: RMS
#hdu_rms = fits.open(rms)
#w_rms = wcs.WCS(hdu_rms[0].header)
#cutout_rms = Cutout2D(1000*hdu_rms[0].data, c, framesize, w_rms)
#ax_rms = fig.add_axes([0.74,0.15,0.2,0.3], projection = cutout_rms.wcs)
#im_rms = ax_rms.imshow(cutout_rms.data, origin="lower", vmin = vmin_rms, vmax = vmax_rms)
#cax_rms = fig.add_axes([0.92, 0.15, 0.008, 0.3])
#cb_rms = plt.colorbar(im_rms, cax = cax_rms)
#
##
#
#for ax in ax_bkg, ax_rms:
#    lon = ax.coords['dec']
#    lon.set_axislabel("Dec")
##    lon.set_major_formatter('dd')
#    lat = ax.coords['ra']
#    lat.set_axislabel("RA")
## RA ticks overlap if left alone
#    lat.set_ticks(number=3)
#
#
#for cb in cb_opt, cb_mos, cb_bkg, cb_rms:
#    cb.set_label("Brightness / mJy beam$^{-1}$")#, labelpad=-2)
##for cb in cb_bkg, cb_rms:
##    cb.set_label("Flux density / mJy beam$^{-1}$")
#
#fig.savefig("XG_mosaic.pdf", dpi=900, bbox_inches='tight')
#fig.savefig("XG_mosaic.png", dpi=900, bbox_inches='tight')

