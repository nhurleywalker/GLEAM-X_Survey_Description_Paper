import os
import numpy as np
from matplotlib import pyplot as plt
#import pltab as plt
#from matplotlib.pyplot import ion
from astropy.io import fits
from scipy import *
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Distance,  SkyCoord

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8}
)

cm = 1/2.54

###
# RUN USING
# aegean \
# --seedclip = 4 \
# --maxsummits = 5 \
# --cores 1 \
# --autoload \
# --progress \
# --nopositive \
# --negative \
# --region /astro/mwasci/tvernstrom/GLEAMX/negative/prime_plus_edge.mim \
# --psf = /astro/mwasci/tgalvin/Mosaic_GLEAMX_Version2/Coadd_Mosaic/IDR_v1/DeepAegeanFix/XG_170-231MHz_projpsf_psf.fits \
# --table = /astro/mwasci/tvernstrom/GLEAMX/negative/XG_170-231MHz_psf_comp_negative_dafix1.fits \
# /astro/mwasci/tgalvin/Mosaic_GLEAMX_Version2/Coadd_Mosaic/IDR_v1/DeepAegeanFix/XG_170-231MHz.fits

cat_neg = fits.getdata('XG_170-231MHz_psf_comp_negative_dafix_comp.fits')

cat_full = fits.getdata('XG_170-231MHz_comp_rescaled.fits', 1)

# Don't download the enormous mosaic just to get the bmaj and bmin; hardcoding them here

#hd = fits.getheader('XG_170-231MHz.fits')
#bmaj = hd["BMAJ"]
#bmin = hd["BMIN"]
bmaj = 2.07931250334E-02
bmin = 1.63785628974E-02

nsamps = ((((13.-4.)*15*np.cos(np.radians(-26.7)))*12)/(bmaj*bmin*1.13))*2.

## FIND THOSE WITHIN THE CENTRAL REGION
##USING MIN RA OF 60 AND MAX 195 AND MIN DEC -32.7 MAX DEC -20.7
keepneg = (cat_neg['ra']>= 60)&(cat_neg['ra']<= 195)&(cat_neg['dec']>= -32.7)&(cat_neg['dec']<= -20.7)
keeppos = (cat_full['ra']>= 60)&(cat_full['ra']<= 195)&(cat_full['dec']>= -32.7)&(cat_full['dec']<= -20.7)


cat_neg = cat_neg[keepneg]
cat_full = cat_full[keeppos]


#CUT OUT ONES BELOW SNR OF 5
cat_neg = cat_neg[abs(cat_neg['int_flux']/cat_neg['local_rms'])>= 5]
cat_full = cat_full[abs(cat_full['int_flux']/cat_full['local_rms'])>= 5]


nneg = cat_neg.shape[0]
npos = cat_full.shape[0]

noise1 = np.median( cat_full['local_rms'])


## MAKE A REGION FILE FOR INSPECTION ##
with open('negatives.reg', 'w') as outf1:
    outf1.write("fk5\n")
    for i in range(cat_neg.size):
        outf1.write(f"point({cat_neg['ra'][i]},{cat_neg['dec'][i]}) # point = cross color = blue\n")

with open('positives.reg', 'w') as outf1:
    outf1.write("fk5\n")
    for i in range(cat_full.size):
        outf1.write(f"point({cat_full['ra'][i]},{cat_full['dec'][i]}) # point = X color = red\n")

## LOOK FOR NEGATIVES SOURCES NEAR POSTIVE ONES
coords_neg  =  SkyCoord(cat_neg['ra']*u.deg,  cat_neg['dec']*u.deg,  frame = 'fk5')
coords_pos = SkyCoord(cat_full['ra']*u.deg, cat_full['dec']*u.deg, frame = 'fk5')
search = coords_neg.search_around_sky(coords_pos, 15*u.arcmin)

## FIND THE BRIGHTEST POSTIVE SOURCE WITHIN 15 ARCMINUTE OF NEGATIVE AND THE CLOSEST POSTIVE
brightid = np.zeros(nneg)
sep_bright = np.zeros(nneg)

closeid = np.zeros(nneg)
sep_close = np.zeros(nneg)
for i in range(nneg):
    jj = np.argwhere(search[1] == i)[:, 0]
    brightid[i] = search[0][jj[np.argmax(cat_full['peak_flux'][search[0][jj]])]]
    sep_bright[i] = search[2][jj[np.argmax(cat_full['peak_flux'][search[0][jj]])]].value*60.
    sep_close[i] = (search[2][jj].value*60).min()
    closeid[i] = search[0][jj][np.argmin(search[2][jj])]

brightid = brightid.astype(int)
closeid = closeid.astype(int)
peak_flux_bright = cat_full['peak_flux'][brightid]
drange = abs(peak_flux_bright/cat_neg['peak_flux'])

## SET CRITERIA 

criteria1 = (sep_bright<= 5)&(drange>= 350)&(peak_flux_bright>= 2.)
criteria2 = (sep_bright<= 12)&(sep_bright>5)&(drange>= 650)&(peak_flux_bright>= 6.)

og_criteria1=(sep_bright<=10)&(drange>=30)&(peak_flux_bright>=0.1)
og_criteria1=(sep_bright<= 5)&(drange>= 350)&(peak_flux_bright>= 2.)
og_criteria2=sep_close<=2

good_neg = np.argwhere((criteria1 | criteria2) == False)[:,0]
og_good_neg = np.argwhere((og_criteria1 | og_criteria2) == False)[:,0]
comb_good_neg = np.argwhere((criteria1 | criteria2 | og_criteria1 | og_criteria2) == False)[:,0]

## DO SAME MATCHING AND FILTERING FOR POSTIVE CATALOGUE
search_pos = coords_pos.search_around_sky(coords_pos, 15*u.arcmin)

brightid_pos = np.zeros(npos).astype(int)
sep_bright_pos = np.zeros(npos).astype(int)

closeid_pos = np.zeros(npos)
sep_close_pos = np.zeros(npos)
for i in range(npos):
    jj = np.argwhere(search_pos[1] == i)[:, 0]
    kj = np.argwhere(search_pos[0][jj]!= i)[:, 0]
    jj = jj[kj]
    if len(jj) != 0:
        brightid_pos[i] = search_pos[0][jj[np.argmax(cat_full['peak_flux'][search_pos[0][jj]])]]
        sep_bright_pos[i] = search_pos[2][jj[np.argmax(cat_full['peak_flux'][search_pos[0][jj]])]].value*60.
        sep_close_pos[i] = (search_pos[2][jj].value*60).min()
        closeid_pos[i] = search_pos[0][jj][np.argmin(search_pos[2][jj])]

peak_flux_bright_pos = cat_full['peak_flux'][brightid_pos]
drange_pos = abs(peak_flux_bright_pos/cat_full['peak_flux'])

criteria1_pos = (sep_bright_pos<= 5)&(drange_pos>= 350)&(peak_flux_bright_pos>= 2.)
criteria2_pos = (sep_bright_pos<= 12)&(sep_bright_pos>5)&(drange_pos>= 650)&(peak_flux_bright_pos>= 6.)

# I don't agree with Tim in including this -- we never want to filter the positive sources with these "OG" criteria
#og_criteria1_pos = (sep_bright_pos<=10)&(drange_pos>=30)&(peak_flux_bright_pos>=0.1)
#og_criteria1_pos = (sep_bright_pos<= 5)&(drange_pos>= 350)&(peak_flux_bright_pos>= 2.)
#og_criteria2_pos = sep_close_pos<=2

# Positive sources that remain after the mandatory filter
good_pos = np.argwhere((criteria1_pos | criteria2_pos) == False )[:,0]

# The sources that are removed -- we want to record these somewhere
# Negative sources
bad_neg = np.argwhere((criteria1 | criteria2) == True)[:,0]
# Positive sources
bad_pos = np.argwhere((criteria1_pos | criteria2_pos) == True )[:,0]

# I don't agree with Tim in including this -- we never want to filter the positive sources with these "OG" criteria
#og_good_pos = np.argwhere((og_criteria1_pos | og_criteria2_pos) == False )[:,0]

# Region files
# Positive sources that remain after the mandatory filter
with open('positives_not_near_bright.reg', 'w') as outf1:
    outf1.write("fk5\n")
    for i in range(good_pos.size):
        outf1.write(f"point({cat_full['ra'][good_pos][i]},{cat_full['dec'][good_pos][i]}) # point = X color = green\n")

# Positive sources that are removed by the mandatory filter
with open('positives_near_bright.reg', 'w') as outf1:
    outf1.write("fk5\n")
    for i in range(bad_pos.size):
        outf1.write(f"point({cat_full['ra'][bad_pos][i]},{cat_full['dec'][bad_pos][i]}) # point = cross color = red\n")

# Negative sources that remain after the mandatory filter
with open('negatives_not_near_bright.reg', 'w') as outf1:
    outf1.write("fk5\n")
    for i in range(good_neg.size):
        outf1.write(f"point({cat_neg['ra'][good_neg][i]},{cat_neg['dec'][good_neg][i]}) # point = X color = magenta\n")

# Negative sources that are removed by the mandatory filter
with open('negatives_near_bright.reg', 'w') as outf1:
    outf1.write("fk5\n")
    for i in range(bad_neg.size):
        outf1.write(f"point({cat_neg['ra'][bad_neg][i]},{cat_neg['dec'][bad_neg][i]}) # point = cross color = yellow\n")

# Save these as FITS tables as well so we can load them into other code

cat_neghd = fits.getheader('XG_170-231MHz_psf_comp_negative_dafix_comp.fits', 1)
cat_fullhd = fits.getheader('XG_170-231MHz_comp_rescaled.fits', 1)

# Positive sources that remain after the mandatory filter
if not os.path.exists("XG_170-231MHz_comp_rescaled_filtered.fits"):
    fits.writeto('XG_170-231MHz_comp_rescaled_filtered.fits', cat_full[good_pos], header = cat_fullhd)

# Positive sources that are removed by the mandatory filter
if not os.path.exists("XG_170-231MHz_psf_comp_positive_dafix_comp_removed.fits"):
    fits.writeto('XG_170-231MHz_psf_comp_positive_dafix_comp_removed.fits', cat_full[bad_pos], header = cat_fullhd)

# Negative sources that remain after the mandatory filter
if not os.path.exists("XG_170-231MHz_psf_comp_negative_dafix_comp_filtered.fits"):
    fits.writeto('XG_170-231MHz_psf_comp_negative_dafix_comp_filtered.fits', cat_neg[good_neg], header = cat_neghd)

# Negative sources that are removed by the mandatory filter
if not os.path.exists("XG_170-231MHz_psf_comp_negative_dafix_comp_removed.fits"):
    fits.writeto('XG_170-231MHz_psf_comp_negative_dafix_comp_removed.fits', cat_neg[bad_neg], header = cat_neghd)


with open('Filtered_Cat_GOOD_IDS.dat', 'w') as outf1:
    for gpi in good_pos:
       print(good_pos, file=outf1)
    
### MAKE SOME PLOTS

fig = plt.figure(figsize = (8*cm, 8*cm))
ax = fig.add_subplot(111)
hxy = ax.hist(cat_neg['int_flux']/cat_neg['local_rms'], bins = 30, alpha = .5, color = 'blue', label = 'No Filtering')
hxy2 = ax.hist(cat_neg['int_flux'][good_neg]/cat_neg['local_rms'][good_neg], bins = hxy[1], alpha = 0.8, color = 'r', label = 'Filtered')
plt.axis([-10, -4, 0., 70.])
ax.legend(loc = 2)#, prop = {'size':13})

ax.set_ylabel('No.of Detections')
ax.set_xlabel('Signal to Noise Ratio')
fig.tight_layout()
fig.savefig('GLEAMX_neg.png')

#hsnrfull = np.histogram(cat_full['peak_flux'][good_pos]/cat_full['local_rms'][good_pos], bins = np.linspace(5, 8, 10))
#hsnrneg = np.histogram((cat_neg['peak_flux'][good_neg]*-1)/cat_neg['local_rms'][good_neg], bins = np.linspace(5, 8, 10))
#hsnrhx = (hsnrneg[1][1:]+hsnrneg[1][:-1])/2

#fig = plt.figure(figsize = (8*cm, 8*cm))
#ax = fig.add_subplot(111)
#ax.plot(hsnrhx, (hsnrneg[0]/hsnrfull[0].astype('float64'))*100, 'r-')
#ax.set_ylabel('% of Fake Sources')
#ax.set_xlabel('ABS(Peak SNR)')
#fig.tight_layout()
#fig.savefig('GLEAMX_neg-ratio.png')

bin_edge = np.linspace(5, 8, 10)
bin_cen = (bin_edge[:-1] + bin_edge[1:]) / 2

# Don't need these -- we ALWAYS do the filter on the bad sources near bright sources
#no_filt_hsnrfull = np.histogram(cat_full['int_flux']/cat_full['local_rms'], bins = bin_edge)
#no_filt_hsnrneg = np.histogram((cat_neg['int_flux']*-1)/cat_neg['local_rms'], bins = bin_edge)

# This is a bit random -- this is the filter ONLY removing the negative sources right next to positive sources, not very useful by itself
# And then also removing the positive sources nearby -- no, this is OTT, we agreed on that.
#hsnrfull = np.histogram(cat_full['int_flux'][og_good_pos]/cat_full['local_rms'][og_good_pos], bins = bin_edge)
#hsnrneg = np.histogram((cat_neg['int_flux'][og_good_neg]*-1)/cat_neg['local_rms'][og_good_neg], bins = bin_edge)
#hsnrhx = (hsnrneg[1][1:]+hsnrneg[1][:-1])/2

# What we actually want:
# Just filtering the spurious components near bright sources:
# For positive sources, that's the only filter we ever apply:
hsnr_pos_single = np.histogram(cat_full['int_flux'][good_pos]/cat_full['local_rms'][good_pos], bins = bin_edge)
# For negative sources, we can apply that filter:
hsnr_neg_single = np.histogram(cat_neg['int_flux'][good_neg]*-1/cat_neg['local_rms'][good_neg], bins = bin_edge)
# and also the filter of removing the ones that lie near other sources
hsnr_neg_double = np.histogram(cat_neg['int_flux'][comb_good_neg]*-1/cat_neg['local_rms'][comb_good_neg], bins = bin_edge)

# This just calcuates the bin centers for easy plotting
hsnrhx = (hsnr_neg_single[1][1:]+hsnr_neg_single[1][:-1])/2

fig = plt.figure(figsize = (8*cm, 8*cm))
ax = fig.add_subplot(111)
ax.plot(hsnrhx, 100*(1-(hsnr_neg_double[0]/hsnr_pos_single[0].astype('float64'))), 'r-', label='Both filters')
ax.plot(hsnrhx, 100*(1-(hsnr_neg_single[0]/hsnr_pos_single[0].astype('float64'))), 'b-', label='Single filter')
ax.set_ylabel('Reliability / \%')
ax.set_xlabel('Signal-to-noise ($S_\mathrm{int} / \sigma$)')
ax.legend(loc='lower right')
fig.tight_layout()
fig.savefig('reliability.pdf')
fig.savefig('reliability.png', transparent=False)
