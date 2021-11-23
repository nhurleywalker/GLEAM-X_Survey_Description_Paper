#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from optparse import OptionParser
#from scipy.optimize import curve_fit
from matplotlib import gridspec
#from scipy.stats import gaussian_kde
from scipy.stats import lognorm
from astropy.modeling import models, fitting
import os
import sys
#tables and votables
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib import rc

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8}
)

cm = 1/2.54

def iqr(x):
    q75, q25 = np.percentile(x, [75 ,25])
    iqr = q75 - q25
    return iqr

def hist_norm_height(n,bins,const):
   ''' Function to normalise bin height by a constant. 
       Needs n and bins from np.histogram.'''

   n = np.repeat(n,2)
   n = np.float32(n) / const
   new_bins = [bins[0]]
   new_bins.extend(np.repeat(bins[1:],2))
   
   return n,new_bins[:-1]

hdu = fits.open("IDR_v1.1_joined_rescaled_cat_seds_subset.fits")
data = hdu[1].data
w = data['int_flux']/data['local_rms']
alpha = data['sp_alpha']

good_alpha = np.where((np.abs(alpha) < 4) & (alpha != np.nan))
S_200 = data['sp_norm'][good_alpha]
alpha_plot = alpha[good_alpha]

source_dim = np.where(S_200 <= 0.01)
alpha_dim = alpha_plot[source_dim]

source_bright = np.where((S_200 > 0.01) & (S_200 <= 0.05))
alpha_bright = alpha_plot[source_bright]

source_bright = np.where((S_200 > 0.05) & (S_200 <= 0.2))
alpha_bright_05 = alpha_plot[source_bright]

source_bright = np.where((S_200 > 0.2))
alpha_bright_1 = alpha_plot[source_bright]

nbins = 100

n_dim, bins_dim = np.histogram(alpha_dim, bins = nbins)
n, bins = np.histogram(alpha_bright, bins = nbins)
n_05, bins_05 = np.histogram(alpha_bright_05, bins = nbins)
n_1, bins_1 = np.histogram(alpha_bright_1, bins = nbins)
n_dim, bins_dim = hist_norm_height(n_dim,bins_dim,max(n_dim))
n, bins = hist_norm_height(n,bins,max(n))
n_05, bins_05 = hist_norm_height(n_05,bins_05,max(n_05))
n_1, bins_1 = hist_norm_height(n_1,bins_1,max(n_1))
# g = fit_g(g_init, n, bins)
fig = plt.figure(figsize=(8*cm,8*cm))
ax = fig.add_subplot(111)
patch_dim, = ax.step(bins_dim, n_dim, color = 'c',linewidth = 3)
patch_br, = ax.step(bins, n, color = 'k',linewidth = 2)
patch_05, = ax.step(bins_05, n_05, color = 'blue',linewidth = 1)
patch_1, = ax.step(bins_1, n_1, color = 'red',linewidth = 0.5)
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel('Normalised count')
ax.set_xlim([-2.5,1.0])
ax.set_ylim([0,1.1])
ax.tick_params(axis="both", which="both")
#ax.set_title('Distribution of spectral indices in GLEAM')
ax.axvline(np.median(alpha_dim), color='c', linestyle='--', lw=0.5)
ax.axvline(np.median(alpha_bright), color='k', linestyle='--', lw=0.5)
ax.axvline(np.median(alpha_bright_05), color='blue', linestyle='--', lw=0.5)
ax.axvline(np.median(alpha_bright_1), color='red', linestyle='--', lw=0.5)

lg = ax.legend([patch_dim,patch_br,patch_05,patch_1], ["$S_\mathrm{200 MHz} \leq 10$\,mJy","$10\,\mathrm{mJy} \leq S_\mathrm{200 MHz}<50$\,mJy","$50\,\mathrm{mJy} \leq S_\mathrm{200 MHz}<200$\,mJy","$S_\mathrm{200 MHz} \geq 200$\,mJy"], bbox_to_anchor=(0.1,1.3), loc="upper left")
#lg = ax.legend([patch_br,patch_05,patch_1], ["$0.16\,\mathrm{Jy} \leq S_\mathrm{200 MHz}<0.5$\,Jy","$0.5\,\mathrm{Jy} \leq S_\mathrm{200 MHz}<1.0$\,Jy","$S_\mathrm{200 MHz} \geq 1.0$\,Jy"])
lg.draw_frame(False)

fig.savefig('alpha_distribution.pdf', bbox_inches="tight")
print(alpha_dim.shape[0], np.median(alpha_dim), iqr(alpha_dim)/2.)
print(alpha_bright.shape[0], np.median(alpha_bright), iqr(alpha_bright)/2.)
print(alpha_bright_05.shape[0], np.median(alpha_bright_05), iqr(alpha_bright_05)/2.)
print(alpha_bright_1.shape[0], np.median(alpha_bright_1), iqr(alpha_bright_1)/2.)

