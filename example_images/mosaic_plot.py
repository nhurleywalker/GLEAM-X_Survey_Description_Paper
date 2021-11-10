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
fig = plt.figure(figsize=(10,5))

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
ax_bkg = fig.add_axes([0.6,0.6,0.24,0.35], projection = cutout_bkg.wcs)
im_bkg = ax_bkg.imshow(cutout_bkg.data, origin="lower", vmin = vmin_bkg, vmax = vmax_bkg)
cax_bkg = fig.add_axes([0.82, 0.6, 0.015, 0.35])
cb_bkg = plt.colorbar(im_bkg, cax = cax_bkg)

# Bottom right panel: RMS
hdu_rms = fits.open(rms)
w_rms = wcs.WCS(hdu_rms[0].header)
cutout_rms = Cutout2D(1000*hdu_rms[0].data, c, framesize, w_rms)
ax_rms = fig.add_axes([0.6,0.1,0.24,0.35], projection = cutout_rms.wcs)
im_rms = ax_rms.imshow(cutout_rms.data, origin="lower", vmin = vmin_rms, vmax = vmax_rms)
cax_rms = fig.add_axes([0.82, 0.1, 0.015, 0.35])
cb_rms = plt.colorbar(im_rms, cax = cax_rms)

for ax in ax_mos, ax_bkg, ax_rms:
    lon = ax.coords['dec']
    lon.set_axislabel("Declination")
#    lon.set_major_formatter('dd')
    lat = ax.coords['ra']
    lat.set_axislabel("Right Ascension")

for cb in cb_mos, cb_bkg, cb_rms:
    cb.set_label("Flux density / mJy beam$^{-1}$")



#f1.show_colorscale(cmap=cmap,vmin=mos_vmin,vmax=mos_vmax)
#f1.set_title("Wideband image")
#
#f2.show_colorscale(cmap=cmap,vmin=bkg_vmin,vmax=bkg_vmax)
#f2.set_title("Background")
#
#f3 = aplpy.FITSFigure(rms,figure=bigfig,subplot=[0.7,0.1,0.24,0.35])
#f3.show_colorscale(cmap=cmap,vmin=rms_vmin,vmax=rms_vmax)
#f3.set_title("RMS")
#
## Workaround for not being allowed colorbars in subplot mode
## mosaic
#af1 = bigfig.add_axes([0.56,0.1,0.015,0.85])
#nf1 = mpl.colors.Normalize(vmin=1000*mos_vmin, vmax=1000*mos_vmax)
#cf1 = mpl.colorbar.ColorbarBase(af1, cmap=cmap,orientation="vertical",norm=nf1)
#cf1.set_label('Flux density (mJy\,beam$^{-1}$)')#,size=35)
##cf1.ax.tick_params(labelsize='35')
#
## background
#af2 = bigfig.add_axes([0.95,0.6,0.015,0.35])
#nf2 = mpl.colors.Normalize(vmin=1000*bkg_vmin, vmax=1000*bkg_vmax)
#cf2 = mpl.colorbar.ColorbarBase(af2, cmap=cmap,orientation="vertical",norm=nf2)
#cf2.set_label('Flux density (mJy\,beam$^{-1}$)')#,size=35)
##cf2.ax.tick_params(labelsize='35')
#
## rms
#af3 = bigfig.add_axes([0.95,0.1,0.015,0.35])
#nf3 = mpl.colors.Normalize(vmin=1000*rms_vmin, vmax=1000*rms_vmax)
#cf3 = mpl.colorbar.ColorbarBase(af3, cmap=cmap,orientation="vertical",norm=nf3)
#cf3.set_label('Flux density (mJy\,beam$^{-1}$)')#,size=35)
##cf3.ax.tick_params(labelsize='35')
#
#f1.add_beam()
#f1.beam.show(zorder=100)
#f1.beam.set_major(mwa_bmaj)  # degrees
#f1.beam.set_minor(mwa_bmin)  # degrees
#f1.beam.set_angle(mwa_bpa)  # degrees
#f1.beam.set_corner('bottom left')
#f1.beam.set_frame(True)
#f1.beam.set_hatch('/')
#f1.beam.set_color('black')
#f1.beam.set_edgecolor('black')
#f1.beam.set_facecolor('black')
#
## Overplot sources
##f1.show_ellipses(my_data['RAJ2000'], my_data['DECJ2000'], my_data['a_deep']/3600.0, my_data['b_deep']/3600.0, my_data['pa_deep'], edgecolor='white')#, linewidth=2,marker='+',s=800)
#    
#f1.show_markers(my_data['RAJ2000'], my_data['DECJ2000'],marker='D',edgecolor="cyan",linewidth=0.5,s=70)
#    
#
#for fig in f1,f2,f3:
#    fig.recenter(RA,Dec,width=width,height=height)
#    fig.tick_labels.set_yformat("dd")
#    fig.tick_labels.set_xformat("hh:mm")
#
##    fig.axis_labels.set_font(size='large', weight='medium', \
##                         stretch='normal', family='serif', \
##                         style='normal', variant='normal')
##    fig.tick_labels.set_font(size='large', weight='medium', \
##                         stretch='normal', family='serif', \
##                         style='normal', variant='normal')
###    fig.ticks.set_xspacing(75)
###    fig.ticks.set_yspacing(30)
#
##fig.set_xaxis_coord_type("longitude")
##fig.tick_labels.set_yformat("mm.m")
##fig.tick_labels.set_font(size='xx-large')
##fig.axis_labels.set_font(size='xx-large')

fig.savefig("XG_mosaic.pdf",bbox_inches='tight')

