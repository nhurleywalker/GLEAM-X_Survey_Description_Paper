#! /usr/bin/env python


import sys

import numpy as np
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc


rc('font', **{'family':'serif', 'serif':['Times'], 'weight':'medium'})
# rc("font", **{"family": "sans-serif", "sans-serif":["Helvetica"], "weight":"medium"})
rc('text', usetex=True)
params = {"text.latex.preamble": [r"\usepackage{siunitx}", \
          r"\sisetup{detect-family = true}"]}
plt.rcParams.update(params)

mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'


def extract_profile(image, radius):

    with fits.open(image) as f:
        data = np.squeeze(f[0].data)
        coords = (f[0].header["CRPIX1"], f[0].header["CRPIX2"])
        bin_size = 1  # 2 pixels
        radius /= f[0].header["CDELT2"]
        print(radius)
        x, y = np.indices(data.shape)
        x, y = y.flatten(), x.flatten()
        r = np.sqrt((x-coords[0])**2 + (y-coords[1])**2)
        cond = r < radius
        indices = np.where(cond)[0]
        print(len(indices))
        x_i, y_i, radii = x[indices], y[indices], r[indices]
        print(y_i)
        values = np.array([data[y_i[i], x_i[i]] for i in range(len(y_i))])
        
        nbins = radius // bin_size
        bin_values = np.array([
            np.nanmean(values[np.where((radii <= (i+1)*bin_size) & (radii > i*bin_size))[0]]) for i in range(int(nbins))
        ])
        bin_radii = np.array([
            ((i+1)*bin_size - 0.5*bin_size)*f[0].header["CDELT2"] for i in range(int(nbins))
        ])

    with open(image.replace(".fits", "_radialprof.txt"), "w+") as f:
        f.write("radius,value\n")
        for i in range(len(bin_values)):
            f.write("{},{}\n".format(bin_radii[i], bin_values[i]))

    return bin_values, bin_radii


def main():

    try:
        cmap = sys.argv[1]
    except IndexError:
        cmap = "viridis"

    weights = ["uniform", r"$r = 0.0$", 
        r"$r = +0.5$",
        r"$r = +1.0$",
        r"$r = +2.0$",
        
        "natural"]
    weightnames = ["uniform", "r0.0", "r0.5", "r1.0", "r2.0",  "natural"]
    font_labels = 8
    font_ticks = 8

    phaseII = "1202216384"
    phaseI = "1102697480"

    cm = 1/2.54
    axes = [0.18, 0.1, 0.82, 0.9]
    fig = plt.figure(figsize=(8*cm, 5*cm))
    ax = fig.add_axes(axes)

    cmap = plt.get_cmap(cmap)
    colors = [cmap(i*0.16) for i in range(6)][::-1]


    for i in range(6):
        if i == 2:
            ls = "--"
        else:
            ls = "-"
        v, r = extract_profile("{}.ms.{}-psf.fits".format(phaseII, weightnames[i]), 10./60.)
        ax.plot(r*60., v, ls=ls, lw=1., color=colors[i], zorder=10+i, label=weights[i])
    v, r = extract_profile("{}.ms.{}-psf.fits".format(phaseI, "r-1.0"), 10./60.)
    ax.plot(r*60., v, ls=":", lw=1., color="black", zorder=10+i, label=r"GLEAM, $r = -1.0$")

    ax.legend(ncol=2, loc="upper right", fontsize=font_ticks-2., frameon=False)
    ax.set_xlabel("Radius / arcmin", fontsize=font_labels)
    ax.set_ylabel("Azimuthally-averaged PSF response", fontsize=font_labels)
    ax.tick_params(which="major", axis="both", labelsize=font_ticks)

    ax.set_xlim([0, 5.])
    ax.axhline(0., ls="-.", color="grey", lw=0.8, zorder=-1)

    fig.savefig("psf_profiles.pdf", transparent=True, dpi=72, bbox_inches="tight")







if __name__ == "__main__":
    main()
