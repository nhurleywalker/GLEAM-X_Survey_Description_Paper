import matplotlib.pyplot as plt
import matplotlib.cm as cmap
from matplotlib.ticker import FuncFormatter
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
#import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

cm = 1/2.54

def S(nu, nu0, S0, alpha):
   return S0 * (nu/nu0)**alpha

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8}
)

# May need to adjust that last digit
viridis = cmap.get_cmap('viridis', 1000)

hdu_gx = fits.open("../alpha_distribution/IDR_v1.1_joined_rescaled_cat_seds_subset.fits")
gx = hdu_gx[1].data

hdu_gl = fits.open("GLEAM_downselect.fits")
gl = hdu_gl[1].data

ind = np.where(np.logical_not(np.isnan(gx["sp_alpha"])))
gx = gx[ind]

ind = np.where(np.logical_not(np.isnan(gl["alpha"])))
gl = gl[ind]

gx_srcs = SkyCoord(gx["ref_ra"], gx["ref_dec"], unit = (u.deg, u.deg), frame="fk5")
gl_srcs = SkyCoord(gl["RAJ2000"], gl["DEJ2000"], unit = (u.deg, u.deg), frame="fk5")

# Match the sources within 15"
# https://docs.astropy.org/en/stable/coordinates/matchsep.html
max_sep = 15. * u.arcsec
idx, d2d, d3d = gx_srcs.match_to_catalog_sky(gl_srcs)
sep_constraint = d2d <= max_sep
gxm = gx[sep_constraint]
glm = gl[idx[sep_constraint]]


makeS = True
# It takes four minutes to generate this plot because of the error bar loop, so make it optional
if makeS is True:
    xmin = 0.02
    xmax = 25.
    # Dummy x for f(x)
    x = np.arange(xmin, xmax+10)

    fig = plt.figure(figsize=(8*cm,8*cm))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    ax.set_aspect("equal")
    ax.set_xlabel("GLEAM $S_\mathrm{200MHz,fitted}$ / Jy")
    ax.set_ylabel("GLEAM-X $S_\mathrm{200MHz,fitted}$ / Jy")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([xmin, xmax])
    ax.plot(x, x, zorder = 20, color="k", lw=0.5, ls="-")
    x = glm["int_flux_fit_200"]
    y = gxm["sp_norm"]

    xerr = glm["err_int_flux_fit_200"]
    yerr = gxm["sp_norm_err"]

    xy = np.log10(np.vstack([x,y]))
    z = gaussian_kde(xy)(xy)
    # Normalise
    z = z / np.max(z)

    # Have to do this in a loop because of reasons...
    # https://github.com/matplotlib/matplotlib/issues/16317
    for xp, yp, xerrp, yerrp, zp in zip(x, y, xerr, yerr, z):
        ax.errorbar(xp, yp, xerr = xerrp, yerr = yerrp, fmt="none", elinewidth = 0.5, lw=0, zorder=zp, ecolor = viridis(zp))
        ax.scatter(xp, yp, marker=".", color = viridis(zp), s = 2, zorder=zp)

    ax.xaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))

    fig.savefig("GLEAM-X_GLEAM_S200_comparison.pdf", bbox_inches="tight")

makeAlpha = True
if makeAlpha is True:
    xmin = -2.5
    xmax = 2.5
    # Dummy alpha
    x = np.arange(xmin, xmax+0.1)

    fig = plt.figure(figsize=(8*cm,8*cm))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    ax.set_aspect("equal")
    ax.set_xlabel("GLEAM $\\alpha_\mathrm{fitted}$")
    ax.set_ylabel("GLEAM-X $\\alpha_\mathrm{fitted}$")
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([xmin, xmax])
    ax.plot(x, x, zorder = 20, color="k", lw=0.5, ls="-")
    x = glm["alpha"]
    y = gxm["sp_alpha"]

    xerr = glm["err_alpha"]
    yerr = gxm["sp_alpha_err"]

    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    # Normalise
    z = z / np.max(z)

    # Have to do this in a loop because of reasons...
    # https://github.com/matplotlib/matplotlib/issues/16317
    for xp, yp, xerrp, yerrp, zp in zip(x, y, xerr, yerr, z):
        ax.errorbar(xp, yp, xerr = xerrp, yerr = yerrp, fmt="none", elinewidth = 0.5, lw=0, zorder=zp, ecolor = viridis(zp))
        ax.scatter(xp, yp, marker=".", color = viridis(zp), s = 2, zorder=zp)

#    ax.xaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))
#    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:g}'.format(y)))

    fig.savefig("GLEAM-X_GLEAM_alpha_comparison.pdf", bbox_inches="tight")