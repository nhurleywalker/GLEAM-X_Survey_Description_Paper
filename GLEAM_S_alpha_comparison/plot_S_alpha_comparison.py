import matplotlib.pyplot as plt
import matplotlib.cm as cmap
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FormatStrFormatter
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
#import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

cm = 1/2.54

def S(nu, nu0, S0, alpha):
   return S0 * (nu/nu0)**alpha

def trueS(S0, q, r):
   return S0 * (0.5 + 0.5*np.sqrt(1 - ((4*q + 4) / (r**2))))

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

makeSr = True
# Plot the ratio of GLEAM to GLEAM-X 200MHz wideband detection flux densities as a function of S/N
# Overplot the curve of the ratio of sigma based on Stefan's deconvolution tests
# Overplot the GLEAM completeness curve
if makeSr is True:
#    data_gleam = np.genfromtxt("../Completeness/GLEAM_completeness.txt", usecols=[0, 4, 5]).T
#    noise_ratios = [4, 8, 16]
#    snapshot_sigmas = [3, 5, 10, 20, 1000]
#    snapshot_ratios = []
#    for s in snapshot_sigmas:
#        csv = np.loadtxt(f"../flux_density_recovery/c169/{s}sigma_multi_r0.5.csv", delimiter=",", comments="#", usecols=(5,8), skiprows=2, max_rows=1)
#        snapshot_ratios.append(csv[1]/csv[0])
#        print(f"../flux_density_recovery/c169/{s}sigma_multi_r0.5.csv")
#        print(s, csv[1], csv[0], csv[1]/csv[0])
    # Convert to arrays to avoid issues later
#    snapshot_ratios = np.array(snapshot_ratios, dtype="float32")
# Scale the signal-to-noise by the ratio of the snapshots to the mosaic noise
# Typical snapshot RMS at 170--200MHz is 4.5mJy/beam
# Typical mosaic RMS in source-finding image (where we care) is 1.2 mJy/beam
    xmin = 12
    xmax = 1.e4
    ymin = 0.3
    ymax = 2.1
    ratio_error = np.sqrt((glm["err_int_flux_wide"]/glm["int_flux_wide"])**2 + (gxm["err_int_flux"]/gxm["int_flux"])**2)

    fig = plt.figure(figsize=(8*cm,8*cm))
    ax = fig.add_axes([0.2,0.2,0.79,0.79])

    #ax.set_ylabel("$\frac{S_\mathrm{200MHz,fitted,GLEAM}}{S_\mathrm{200MHz,fitted,GLEAM-X}}$")
    ax.set_ylabel(r"$\frac{S_\mathrm{200MHz,GLEAM-X}}{S_\mathrm{200MHz,GLEAM}}$")
    ax.set_xlabel(r"$\frac{S_\mathrm{200MHz,GLEAM-X}}{\sigma_\mathrm{200MHz,GLEAM-X}}$")
    #ax.set_xlabel(r"$S_\mathrm{200MHz,fitted,GLEAM-X}$")

# Only plot this if the x-axis is S/N
#    s_gleam = data_gleam[0]
#    c_gleam = 1.05*data_gleam[1]/100
#    cerr_gleam = data_gleam[2]/100
#    ax.errorbar(
#       s_gleam,
#       c_gleam,
#       ms=1.0,
#       color="red",
#       marker=',',
#       linewidth=0.7,
#       yerr = cerr_gleam
#    )

#    ax.set_xscale("log")
#    ax.set_yscale("log")
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    #y = gxm["sp_norm"]/glm["int_flux_fit_200"]
    glm_corr = trueS(glm["int_flux_wide"], 1.54, glm["int_flux_wide"]/glm["local_rms_wide"])
#    glm_corr = trueS(glm["int_flux_151"], 1.54, glm["int_flux_151"]/glm["local_rms_151"])
    #y = gxm["int_flux"]/glm_corr #glm["int_flux_wide"]
    y = gxm["int_flux"]/glm_corr #glm["int_flux_wide"]
    #y = gxm["int_flux_N_147_154MHz"]/glm_corr #glm["int_flux_wide"]
    #x = gxm["sp_norm"]#/gxm["local_rms"]
    x = gxm["int_flux"]/gxm["local_rms"]
    #x = gxm["int_flux_N_147_154MHz"]#/gxm["local_rms"]
    #ax.scatter(x, y, zorder = 1, color="k", marker=".", alpha=0.5)#, ls="-")
    ax.hexbin(x, y, xscale="log", yscale="log", bins="log", zorder = 1, lw=0.2, cmap="binary")#, color="k", marker=".", alpha=0.5)#, ls="-")
#    ax.hexbin(x, ratio_error, xscale="log", yscale="log", bins="log", zorder = 2, lw=0.2, cmap="viridis")#, color="k", marker=".", alpha=0.5)#, ls="-")
#    for n in noise_ratios:
        #noise_ratio = (4.5/1.2)
#        noise_ratio = n
#        mosaic_sigmas = noise_ratio*np.array(snapshot_sigmas, dtype="float32")
#        ax.plot(mosaic_sigmas, snapshot_ratios, zorder=20, ls="-", lw=1, label=f"{n}")
    ax.axhline(1.0, color="k", lw=0.5, ls="-")
    ax.axhline(1.05, color="k", lw=0.5, ls="--")
#    ax.axvline(3*noise_ratio, color="k", lw=0.5, ls="--")
#    ax.axvline(.060, color="k", lw=0.5, ls="--")
    ax.axvline(100, color="k", lw=0.5, ls="--")
    ax.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
    ax.yaxis.set_minor_formatter(FormatStrFormatter('%3.1f'))
#    ax.legend()
    fig.savefig("GLEAM-X_GLEAMcorr_nofit_ratio_hexbins.pdf",bbox_inches="tight",dpi=1000)

makeS = False
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

makeAlpha = False
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
