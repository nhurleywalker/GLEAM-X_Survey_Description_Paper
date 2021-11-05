#!/usr/bin/env python

# I used AeRes to mask the sources in the image
# Now this script makes a plot of the noise distribution

from glob import glob

import copy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['serif']})
labelfontsize=12

from astropy.io import fits
import numpy as np
from scipy.stats import gaussian_kde
from astropy.modeling import models, fitting

import sys

def hist_norm_height(n,bins,const):
    ''' Function to normalise bin height by a constant. 
        Needs n and bins from np.histogram.'''

    n = np.repeat(n,2)
    n = np.float32(n) / const
    new_bins = [bins[0]]
    new_bins.extend(np.repeat(bins[1:],2))
    
    return n,new_bins[:-1]

# TODO sort out argparse

bins = 75
nx = np.empty(bins)
stdev = 7.5
nsig = 6
xmin = -nsig*stdev
xmax = nsig*stdev

# -rw-rw-r-- 1 tash tash 6.0M Nov  5 14:23 XG_residual_bkgsubtracted_sigma.fits
#-rw-rw-r-- 1 tash tash 6.0M Nov  5 14:23 XG_masked_bkgsubtracted_sigma.fits
#-rw-rw-r-- 1 tash tash 6.0M Nov  5 14:22 XG_bkgsubtracted_sigma.fits
#-rw-rw-r-- 1 tash tash 6.0M Nov  5 14:22 XG_sigma.fits

# Background-subtracted
tempd = fits.getdata("XG_bkgsubtracted_sigma.fits")
nss, binsss = np.histogram(tempd, bins = bins, range=(-nsig,nsig))
nss = np.array(nss,dtype="float32")/len(nss)
binsss_centers = (binsss+((binsss[1]-binsss[0])/2))[:-1]
g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=1.0)
fit_g = fitting.LevMarLSQFitter()
gss = fit_g(g_init, binsss_centers, nss)
print("Stdev of background-subtracted image is {0}".format(gss.stddev.value))


# Masked, background-subtracted
tempd = fits.getdata("XG_masked_bkgsubtracted_sigma.fits")
nms, binsms = np.histogram(tempd,bins = bins,range=(-nsig,nsig))
nms = np.array(nms,dtype="float32")/len(nms)
binsms_centers = (binsms+((binsms[1]-binsms[0])/2))[:-1]
g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=1.0)
fit_g = fitting.LevMarLSQFitter()
gms = fit_g(g_init, binsms_centers, nms)
print("Stdev of masked, background-subtracted image is {0}".format(gms.stddev.value))

# Residuals, background-subtracted
tempd = fits.getdata("XG_residual_bkgsubtracted_sigma.fits")
nrs, binsrs = np.histogram(tempd,bins = bins,range=(-nsig,nsig))
nrs = np.array(nrs,dtype="float32")/len(nrs)
binsrs_centers = (binsrs+((binsrs[1]-binsrs[0])/2))[:-1]
g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=1.0)
fit_g = fitting.LevMarLSQFitter()
grs = fit_g(g_init, binsrs_centers, nrs)
print("Stdev of residual, background-subtracted image is {0}".format(grs.stddev.value))

# unmasked, unbackgrounded
img = "XG.fits"
tempd = 1000*fits.getdata(img)
nu, binsu = np.histogram(tempd,bins = np.linspace(xmin,xmax,bins))
nu = np.array(nu,dtype="float32")/len(nu)
binsu_centers = (binsu+((binsu[1]-binsu[0])/2))[:-1]
g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
fit_g = fitting.LevMarLSQFitter()
gu = fit_g(g_init, binsu_centers, nu)
print("Stdev of original image is {0}Jy/beam (meaningless, here for completeness)".format(gu.stddev.value))

gmod = copy.deepcopy(gu)
tempr = 1000*fits.getdata("XG_rms.fits")
mean_rms = np.nanmean(tempr)
gmod.stddev.value = mean_rms
print("BANE detects {0} mJy/beam noise in this region.".format(mean_rms))

gmod1 = copy.deepcopy(gu)
gmod1.stddev.value = 1.0
gmod1.mean.value = 0.0
gmod1.amplitude.value = max(nss)
gmod2 = copy.deepcopy(gu)
gmod2.stddev.value = 1.0
gmod2.mean.value = 0.0
gmod2.amplitude.value = max(nrs)

fig = plt.figure(figsize=(12, 4))
# Background-subtracted
ax1 = fig.add_subplot(131)
ax1.bar(binsss[:-1], nss, color = 'lightgrey', edgecolor = "none", width=(binsss[1]-binsss[0])) # Histogram
ax1.plot(np.linspace(-nsig,nsig,100),gmod1(np.linspace(-nsig,nsig,100)), color='k',lw=2, label='BANE')
ax1.set_ylim([1.0,1.2*max(nss)])
ax1.set_title("Background-subtracted S/N image")
# Masked, background-subtracted
ax2=fig.add_subplot(132)
ax2.bar(binsms[:-1], nms, color = 'lightgrey', edgecolor = "none", width=(binsms[1]-binsms[0])) # Histogram
ax2.plot(np.linspace(-nsig,nsig,100),gmod1(np.linspace(-nsig,nsig,100)), color='k',lw=2, label='BANE')
ax2.set_ylim([1.0,1.2*max(nms)])
ax2.set_title("$>5\sigma$ Sources masked; background subtracted")
# Residuals, background-subtracted
ax3 = fig.add_subplot(133)
ax3.bar(binsrs[:-1], nrs, color = 'lightgrey', edgecolor = "none", width=(binsrs[1]-binsrs[0])) # Histogram
ax3.plot(np.linspace(-nsig,nsig,100),gmod2(np.linspace(-nsig,nsig,100)), color='k',lw=2, label='BANE')
ax3.set_ylim([1.0,1.2*max(nrs)])
ax3.set_title("$>5\sigma$ Sources and background subtracted")

for ax in ax1, ax2, ax3:
    ax.axvline(x=0.0, lw=1, color='k', linestyle='-')
    ax.axvline(x=-1, lw=1, color='k', linestyle='--')
    ax.axvline(x=-2, lw=1, color='k', linestyle='-.')
    ax.axvline(x=-5, lw=1, color='k', linestyle=':')
    ax.axvline(x=1, lw=1, color='k', linestyle='--')
    ax.axvline(x=2, lw=1, color='k', linestyle='-.')
    ax.axvline(x=5, lw=1, color='k', linestyle=':')
    ax.set_yscale("log", nonposy='clip')
    ax.set_xlim([-nsig,nsig])
    ax.tick_params(labelsize=labelfontsize)
    ax.set_xlabel("$\sigma$",fontsize=labelfontsize)

fig.tight_layout()
fig.savefig("noise_distribution.pdf",pad_inches=0.1,bbox_inches="tight")

# What area of sky is this?
#unnorm_nx, unnorm_binsx = np.histogram(tempd,bins = bins,range=(xmin,xmax))
#totalpixels = np.sum(unnorm_nx)
#header = fits.getheader("image.fits")
#cdelt = header['CD1_1'] # deg per pixel

# TODO Get from the PSF map
#bmaj = 0.0208277
#bmin = 0.016561
#pix_per_beam = (np.pi*bmaj*bmin/(4*np.log(2)))/(cdelt*cdelt)
#totalbeams = totalpixels/pix_per_beam

#fig = plt.figure(figsize=(15, 5))
#ax1 = fig.add_subplot(131)
#ax1.bar(binsu[:-1], nu, color = 'lightgrey', edgecolor = "none", width=(binsu[1]-binsu[0])) # Histogram
#ax1.plot(np.linspace(-nsig,nsig,100),gmod1(np.linspace(-nsig,nsig,100)), color='k',lw=2, label='BANE')
#ax1.set_ylim([1.0,1.2*max(nbs)])
#ax1.set_yscale("log", nonposy='clip')
#ax1.set_title("Unmodified image")
#fig.savefig("first_test.png")



# unbackgrounded
#img = XG_residual_masked.fits
#tempd = 1000*fits.getdata("image_subtracted.fits")
#nm, binsm = np.histogram(tempd,bins = bins,range=(xmin,xmax))
#nm = np.array(nm,dtype="float32")/len(nm)
#binsm_centers = (binsm+((binsm[1]-binsm[0])/2))[:-1]
#g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
#fit_g = fitting.LevMarLSQFitter()
#gm = fit_g(g_init, binsm_centers, nm)

# subtracted, backgrounded
#tempd=1000*fits.getdata("blanked_subtracted_fc1.fits")
#nx, binsx = np.histogram(tempd,bins = bins,range=(xmin,xmax))
##nx=np.array(nx,dtype="float32")/float(nx.max())
#nx=np.array(nx,dtype="float32")/len(nx)
#binsx_centers = (binsx+((binsx[1]-binsx[0])/2))[:-1]
#
#g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
#fit_g = fitting.LevMarLSQFitter()
#gx = fit_g(g_init, binsx_centers, nx)

# Residuals, backgrounded
#tempd=1000*fits.getdata("residual_subtracted.fits")
#nr, binsr = np.histogram(tempd,bins = bins,range=(xmin,xmax))
##nr=np.array(nr,dtype="float32")/float(nr.max())
#nr=np.array(nr,dtype="float32")/len(nr)
#binsr_centers = (binsr+((binsr[1]-binsr[0])/2))[:-1]
#
#g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
#fit_g = fitting.LevMarLSQFitter()
#gr = fit_g(g_init, binsr_centers, nr)

# sigma plots

#raw
#tempd=fits.getdata("sigma.fits")
#nss, binsss = np.histogram(tempd,bins = bins,range=(-nsig,nsig))
#nss=np.array(nss,dtype="float32")/len(nss)
#binsss_centers = (binsss+((binsss[1]-binsss[0])/2))[:-1]
#
#g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=1.0)
#fit_g = fitting.LevMarLSQFitter()
#gss = fit_g(g_init, binsss_centers, nss)
#print gss.stddev.value



# Make a modified Gaussian based on what BANE fits

# Another deep copy, this time just plotting 1-sigma
# matplotlib bars are defined by their *bottom left* corner
# so don't use bin_centers, use original bins

#ax1.bar(binsu[:-1], nu, color = 'lightgrey', edgecolor = "none", width=(binsu[1]-binsu[0])) # Histogram
#ax1.plot(np.linspace(xmin,xmax,100),gu(np.linspace(xmin,xmax,100)),color='k',lw=1, label='Gaussian')
#ax1.plot([gu.mean.value,gu.mean.value], [0,gu.amplitude.value], 'k-', lw=1)
#ax1.set_title("Raw image")
#ax1.set_ylabel("N$_\mathrm{pix}$",fontsize=labelfontsize)
#ax1.set_ylim([1.0,1.2*max(nu)])
#ax1.set_yscale("log", nonposy='clip')
#ax1.set_xlim([xmin,xmax])
#ax1.tick_params(labelsize=labelfontsize)
#ax1.set_xlabel("S / mJy\,beam$^{-1}$",fontsize=labelfontsize)

#ax2 = fig.add_subplot(132)
#ax2.bar(binsm[:-1], nm, color = 'lightgrey', edgecolor = "none", width=(binsm[1]-binsm[0])) # Histogram
##ax2.fill_between(np.linspace(xmin,xmax,100),gm(np.linspace(xmin,xmax,100)), facecolor="black", alpha=0.3, label='Gaussian')
#ax2.plot(np.linspace(xmin,xmax,100),gm(np.linspace(xmin,xmax,100)), color='k', lw=1, label='Gaussian')
#ax2.plot([gm.mean.value,gm.mean.value], [0,gm.amplitude.value], 'k-', lw=1)
#ax2.set_title("Background subtracted")

#ax3=fig.add_subplot(133)
#ax3.bar(binsx[:-1], nx, color = 'lightgrey', edgecolor = "none", width=(binsx[1]-binsx[0])) # Histogram
#ax3.plot(np.linspace(xmin,xmax,100),gx(np.linspace(xmin,xmax,100)), color='k',lw=1, label='Ensemble')
#ax3.plot(np.linspace(xmin,xmax,100),gmod(np.linspace(xmin,xmax,100)), color='k',lw=2, label='BANE')
#lg = ax3.legend(fontsize=12,loc=2) #'upper left'
#lg.draw_frame(False)
#ax3.plot(binsx_centers,gx(binsx_centers)-nx, color="c",lw=2, label='Gaussian')
#ax3.set_title("$5\sigma+$ Sources blanked; bkg subtracted")

# Residuals instead of blank

#ax4 = fig.add_subplot(134)
#ax4.bar(binsr[:-1], nr, color = 'lightgrey', edgecolor = "none", width=(binsr[1]-binsr[0])) # Histogram
#ax4.plot(np.linspace(xmin,xmax,100),gr(np.linspace(xmin,xmax,100)), color='k',lw=1, label='Ensemble')
#ax4.plot(np.linspace(xmin,xmax,100),gmod(np.linspace(xmin,xmax,100)), color='k',lw=2, label='BANE')
#lg = ax4.legend(fontsize=12,loc=2) #'upper left'
#lg.draw_frame(False)
#ax4.set_title("$5\sigma+$ Sources and background subtracted")





# Mean, stddev
#for ax in ax3, ax4:
#    ax.plot([gmod.mean.value,gmod.mean.value], [0,gmod.amplitude.value], 'k-', lw=1)
#    ax.plot([gmod.mean.value-gmod.stddev.value,gmod.mean.value-gmod.stddev.value], [0,gmod.amplitude.value], 'k--', lw=1)
#    ax.plot([gmod.mean.value-2*gmod.stddev.value,gmod.mean.value-2*gmod.stddev.value], [0,gmod.amplitude.value], 'k-.', lw=1)
#    ax.plot([gmod.mean.value-5*gmod.stddev.value,gmod.mean.value-5*gmod.stddev.value], [0,gmod.amplitude.value], 'k:', lw=1)
#    ax.plot([gmod.mean.value+gmod.stddev.value,gmod.mean.value+gmod.stddev.value], [0,gmod.amplitude.value], 'k--', lw=1)
#    ax.plot([gmod.mean.value+2*gmod.stddev.value,gmod.mean.value+2*gmod.stddev.value], [0,gmod.amplitude.value], 'k-.', lw=1)
#    ax.plot([gmod.mean.value+5*gmod.stddev.value,gmod.mean.value+5*gmod.stddev.value], [0,gmod.amplitude.value], 'k:', lw=1)


# Difference
#fig2 = plt.figure(figsize=(5, 5))
#ax4=fig2.add_subplot(111)
#ax4.set_xlabel("$\sigma$",fontsize=labelfontsize)
#ax4.set_ylabel("Ratio of data to Gaussian fit",fontsize=labelfontsize)
#ax4.plot(np.array(binsx_centers)/gx.stddev.value,((nx/gx(np.array(binsx_centers)))),label="Gaussian")
#ax4.plot(np.array(binsx_centers)/gmod.stddev.value,((nx/gmod(np.array(binsx_centers)))),label="BANE")
#ax4.plot([-5,-5], [0,10000], 'k--', lw=1)
#ax4.plot([5,5], [0,10000], 'k--', lw=1)
#ax4.set_yscale("log", nonposy='clip')
#ax4.set_ylim([0,10000])
#ax4.legend()
#fig2.savefig("difference.png",pad_inches=0.1,bbox_inches="tight")


# Bin everything in sigma rather than absolute value
#sigma_binsx_centers = binsx_centers/gmod.stddev.value
#diff = nx-gmod(binsx_centers)
#
## Work out the negative source numbers
#excess_neg = np.sum(diff[np.where(sigma_binsx_centers<-5.0)])
#total = np.sum(nx)
#excess_neg_fraction=excess_neg/total
#excess_neg_beams = excess_neg_fraction*totalbeams
#
## area of sky = total pixels * cdelt * cdelt
#area_of_sky = cdelt * cdelt * totalpixels
#area_of_survey = 24402
#total_neg_sources = (area_of_survey/area_of_sky) * excess_neg_beams
#
## Work out the positive source numbers
#excess_pos = np.sum(diff[np.where(sigma_binsx_centers>5.0)])
#excess_pos_fraction = excess_pos/total
#excess_pos_beams = excess_pos_fraction*totalbeams
#
#total_pos_sources = (area_of_survey/area_of_sky) * excess_pos_beams
#
#print "This area contains ", int(excess_neg_beams)," unreal negative sources."
#print "This area contains ", int(excess_pos_beams)," extra positive sources."
#print "This area contains ", int(excess_pos_beams-excess_neg_beams)," real confused positive beams."

#print "Survey should contain ", int(total_neg_sources)," unreal negative sources, and therefore unreal positive sources."
#print "Survey should contain ", int(total_pos_sources)," extra positive sources."
#print "Survey should contain ", int(total_pos_sources-total_neg_sources)," real confused positive sources."




# Looking at smaller regions
# highrms
#
#tempd=1000*fits.getdata("masked_subsection_highrms.fits")
#mean_highrms = np.nanmean(tempd)
#tempd=1000*fits.getdata("masked_subtracted_highrms.fits")
#nx, binsx = np.histogram(tempd,bins = bins,range=(xmin,xmax))
#nx=np.array(nx,dtype="float32")/float(nx.max())
#binsx_centers = (binsx+((binsx[1]-binsx[0])/2))[:-1]
#
#g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
#fit_g = fitting.LevMarLSQFitter()
#gx = fit_g(g_init, binsx_centers, nx)
#
#gmod=copy.deepcopy(gx)
#gmod.stddev.value=mean_highrms
#print mean_highrms
#
#fig3 = plt.figure(figsize=(10, 10))
#ax5=fig3.add_subplot(221)
#ax5.bar(binsx[:-1], nx, color = 'lightgrey', edgecolor = "none", width=(binsx[1]-binsx[0])) # Histogram
#ax5.plot(np.linspace(xmin,xmax,100),gx(np.linspace(xmin,xmax,100)), color='k',lw=2, label='Gaussian')
#ax5.plot(np.linspace(xmin,xmax,100),gmod(np.linspace(xmin,xmax,100)), color='k',lw=1, label='Gaussian')
##ax3.plot(binsx_centers,gx(binsx_centers)-nx, color="c",lw=2, label='Gaussian')
#ax5.set_title("High RMS region")
#
#ax7 = fig3.add_subplot(223)
#ax7.set_xlabel("$\sigma$",fontsize=labelfontsize)
#ax7.set_ylabel("Ratio of data to Gaussian fit",fontsize=labelfontsize)
#ax7.plot(np.array(binsx_centers)/gx.stddev.value,((nx/gx(np.array(binsx_centers)))),label="Gaussian")
#ax7.plot(np.array(binsx_centers)/gmod.stddev.value,((nx/gmod(np.array(binsx_centers)))),label="BANE")
#ax7.plot([-5,-5], [0,10000], 'k--', lw=1)
#ax7.plot([5,5], [0,10000], 'k--', lw=1)
#ax7.set_yscale("log", nonposy='clip')
#ax7.set_ylim([0,10000])
#
#tempd=1000*fits.getdata("masked_subsection_lowrms.fits")
#mean_lowrms = np.nanmean(tempd)
#tempd=1000*fits.getdata("masked_subtracted_lowrms.fits")
#nx, binsx = np.histogram(tempd,bins = bins,range=(xmin,xmax))
#nx=np.array(nx,dtype="float32")/float(nx.max())
#binsx_centers = (binsx+((binsx[1]-binsx[0])/2))[:-1]
#
#g_init = models.Gaussian1D(amplitude=1., mean=0., stddev=stdev)
#fit_g = fitting.LevMarLSQFitter()
#gx = fit_g(g_init, binsx_centers, nx)
#
#gmod=copy.deepcopy(gx)
#gmod.stddev.value=mean_lowrms
#print mean_lowrms
#
#ax6=fig3.add_subplot(222)
#ax6.bar(binsx[:-1], nx, color = 'lightgrey', edgecolor = "none", width=(binsx[1]-binsx[0])) # Histogram
#ax6.plot(np.linspace(xmin,xmax,100),gx(np.linspace(xmin,xmax,100)), color='k',lw=2, label='Gaussian')
#ax6.plot(np.linspace(xmin,xmax,100),gmod(np.linspace(xmin,xmax,100)), color='k',lw=1, label='Gaussian')
##ax3.plot(binsx_centers,gx(binsx_centers)-nx, color="c",lw=2, label='Gaussian')
#ax6.set_title("Low RMS region")
#
#ax8 = fig3.add_subplot(224)
#ax8.set_xlabel("$\sigma$",fontsize=labelfontsize)
#ax8.set_ylabel("Ratio of data to Gaussian fit",fontsize=labelfontsize)
#ax8.plot(np.array(binsx_centers)/gx.stddev.value,((nx/gx(np.array(binsx_centers)))),label="Gaussian")
#ax8.plot(np.array(binsx_centers)/gmod.stddev.value,((nx/gmod(np.array(binsx_centers)))),label="BANE")
#ax8.plot([-5,-5], [0,10000], 'k--', lw=1)
#ax8.plot([5,5], [0,10000], 'k--', lw=1)
#ax8.set_yscale("log", nonposy='clip')
#ax8.set_ylim([0,10000])
#
#fig3.savefig("cfrms.png",pad_inches=0.1,bbox_inches="tight")
##
