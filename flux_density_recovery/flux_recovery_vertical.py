#! /usr/bin/env python

####
# for r in r0.0 r0.5 r2.0 r1.0 uniform natural; do cd multi_${r}; s=20; ~/Dropbox/scripts/flux_recovery.py 154 1202216384_zero_SNR${s} $(seq 60 20 600) --radius-unit 3600 --snr $s --sigma 0.0 -A 2.0 -a --plot --name multi-${r}-SNR${s} --recovery --model 1202216384_zero_SNR${s}_models.txt ; cd .. ; done
####

import sys

from re import M
import numpy as np
from argparse import ArgumentParser
from subprocess import Popen

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import rgb2hex, LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

from matplotlib import rc

channel = sys.argv[1]

cmap = "viridis"

cm = 1/2.54

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8}
)

mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'

weights = ["uniform", "robust0", "robust0.5", "robust1", "robust2", "natural"]
WEIGHTS = ["uniform", "r0.0", "r0.5", "r1.0", "r2.0", "natural"]
WEIGHTS = ["uniform", "r0.0", "r0.5", "r1.0", "natural"]

me = ""  # legacy

d = "c{}/20sigma_".format(channel)

files_multi = [d + me + "multi_{}.csv".format(i) for i in WEIGHTS]
data_multi1 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]
files_multi = [d + me + "nomulti_{}.csv".format(i) for i in WEIGHTS]
data_nomulti1 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]

d = "c{}/1000sigma_".format(channel)
files_multi = [d + me +  "multi_{}.csv".format(i) for i in WEIGHTS]
data_multi2 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]
files_multi = [d + me +  "nomulti_{}.csv".format(i) for i in WEIGHTS]
data_nomulti2 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]


d = "c{}/10sigma_".format(channel)
files_multi = [d + me +  "multi_{}.csv".format(i) for i in WEIGHTS]
data_multi3 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]
files_multi = [d + me +  "nomulti_{}.csv".format(i) for i in WEIGHTS]
data_nomulti3 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]


d = "c{}/5sigma_".format(channel)
files_multi = [d + me +  "multi_{}.csv".format(i) for i in WEIGHTS]
data_multi4 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]
files_multi = [d + me +  "nomulti_{}.csv".format(i) for i in WEIGHTS]
data_nomulti4 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]


d = "c{}/3sigma_".format(channel)
files_multi = [d + me +  "multi_{}.csv".format(i) for i in WEIGHTS]
data_multi5 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]
files_multi = [d + me +  "nomulti_{}.csv".format(i) for i in WEIGHTS]
data_nomulti5 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]



font_labels = 8
font_ticks = 8

axes1 = [0.18, 0.1, 0.45, 0.45]
axes2 = [0.18, 0.57, 0.45, 0.45]
axes3 = [0.18+0.47, 0.1, 0.45, 0.45]
axes4 = [0.18+0.47, 0.57, 0.45, 0.45]
axes5 = [0.18, 0.1+0.47+0.47, 0.45, 0.45]

ax1len = axes1[0] + axes1[3]
tlen = ax1len + 0.39 + (0.55-ax1len) 

fig = plt.figure(figsize=(8*cm,8*cm))
ax1 = fig.add_axes(axes1)
ax2 = fig.add_axes(axes2)
ax3 = fig.add_axes(axes3)
ax4 = fig.add_axes(axes4)
ax5 = fig.add_axes(axes5)

cmap = plt.get_cmap(cmap)
colors = [cmap(i*0.2) for i in range(len(WEIGHTS))][::-1]
beams = []

# this is weird but bear with me ("for legacy reasons"): 
# middle right, top left, bottom left, middle left, bottom right
all_axes = [ax4, ax5, ax1, ax2, ax3]
all_data = [data_multi3, data_multi2, data_multi4, data_multi1, data_multi5]
all_nodata = [data_nomulti3, data_nomulti2, data_nomulti4, data_nomulti1, data_nomulti5]
all_names = ["$10\\sigma$", "$1000\\sigma$",  "$5\\sigma$",  "$20\\sigma$", "$3\\sigma$"]


for n in range(len(all_axes)):
    data = all_data[n]
    nodata = all_nodata[n]
    for i in range(len(WEIGHTS)):
        
        bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
        beams.append(bmaj)
        all_axes[n].plot(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["aegean"]/data[i]["model"], 
            ls="-",
            marker="+", 
            ms=4,
            markeredgewidth=0.5, 
            color=colors[i], 
            lw=0.5,
            zorder=10+i)

        all_axes[n].plot(np.sqrt(nodata[i]["radius"]**2 - bmaj**2), nodata[i]["aegean"]/nodata[i]["model"], 
            ls=":",
            marker="v", 
            ms=2,
            markeredgewidth=0, 
            mec="none",
            color=colors[i], 
            lw=0.5,
            zorder=i,
            alpha=0.7)


for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.set_ylim([0.35, 1.03])
    ax.set_xlim([1., 10.])

for ax in [ax1, ax2, ax3, ax4, ax5]:
    for spacing in [0.4, 0.6, 0.8, 1.0]:
        ax.axhline(spacing, ls="--", lw=1, color="grey", zorder=-1)

handles = [
    Line2D([0], [0], linestyle=":", marker="v", color=colors[0], lw=0.5, markeredgewidth=0, ms=3, alpha=0.75, label="uniform, no multi-scale"),
    Line2D([0], [0], linestyle=":", marker="v", color=colors[1], lw=0.5, markeredgewidth=0, ms=3,alpha=0.75, label=r"$r=0.0$, no multi-scale"),
    Line2D([0], [0], linestyle=":", marker="v",color=colors[2], lw=0.5, markeredgewidth=0, ms=3,alpha=0.75, label=r"$r=+0.5$, no multi-scale"),
    Line2D([0], [0], linestyle=":", marker="v",color=colors[3], lw=0.5, markeredgewidth=0, ms=3,alpha=0.75, label=r"$r=+1.0$, no multi-scale"),
    Line2D([0], [0], linestyle=":", marker="v",color=colors[4], lw=0.5, markeredgewidth=0, ms=3,alpha=0.75, label="natural, no multi-scale"),
    Line2D([0], [0], linestyle="-", marker="+", color=colors[0], lw=0.5, markeredgewidth=0.5, label=r"uniform, multi-scale"),
    Line2D([0], [0], linestyle="-", marker="+",color=colors[1], lw=0.5, markeredgewidth=0.5, label=r"$r=0.0$, multi-scale"),
    Line2D([0], [0], linestyle="-", marker="+",color=colors[2], lw=0.5, markeredgewidth=0.5, label=r"$r=+0.5$, multi-scale"),
    Line2D([0], [0], linestyle="-", marker="+",color=colors[3], lw=0.5, markeredgewidth=0.5, label=r"$r=+1.0$, multi-scale"),
    Line2D([0], [0], linestyle="-", marker="+",color=colors[4], lw=0.5, markeredgewidth=0.5, label=r"natural, multi-scale")
]

ax2.legend(ncol=1, bbox_to_anchor=(1.0,2.04), loc="upper left", fontsize=font_ticks-1,
           handles=handles[::-1], handlelength=3., frameon=False)

fig.text(0.04, 0.78, r"$S_\mathrm{measured} / S_\mathrm{model}$ ", 
    fontsize=font_labels, 
    va="center", 
    rotation=90.)
fig.text(0.38, 0, "Model Gaussian FWHM [arcmin]", 
    fontsize=font_labels, 
    va="center")

for ax in ax1, ax2, ax3, ax4, ax5:
    ax.set_yticks([0.4, 0.6, 0.8, 1.0])
for ax in ax1, ax3, ax5:
    ax.set_xticks([2.,4.,6.,8.,10])
for ax in ax3, ax4:
    ax.axes.get_yaxis().set_ticks([])
for ax in ax2, ax4, ax5:
    ax.axes.get_xaxis().set_ticks([])

for i in range(len(all_axes)):
    all_axes[i].text(0.03, 0.1, "{} {}".format(all_names[i], me), 
                     fontsize=font_labels, 
                     ha="left", 
                     va="bottom", 
                     transform=all_axes[i].transAxes)


fig.savefig("flux_recovery_3-5-10-20-1000_vertical_c{}.pdf".format(channel), 
    transparent=True, 
    dpi=300, 
    bbox_inches="tight")
