#! /usr/bin/env python

import numpy as np
from argparse import ArgumentParser
from subprocess import Popen


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import rgb2hex, LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

#import cmasher as cmr

from matplotlib import rc

import sys

try:
    cmap = sys.argv[1]
except IndexError:
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

d = "./20sigma_"
files_multi = [d + "multi_{}.csv".format(i) for i in WEIGHTS]
files = [d + "nomulti_{}.csv".format(i) for i in WEIGHTS]
data1 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files]
data_multi1 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]

d = "./1000sigma_"
files_multi = [d + "multi_{}.csv".format(i) for i in WEIGHTS]
files = [d + "nomulti_{}.csv".format(i) for i in WEIGHTS]
data2 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files]
data_multi2 = [np.genfromtxt(f, delimiter=",", names="radius,flux,eflux,peak,rms,model,modelpeak,arcmin,aegean") for f in files_multi]

font_labels = 8
font_ticks = 8

axes1 = [0.18, 0.1, 0.45, 0.45]
axes2 = [0.18, 0.57, 0.45, 0.45]
axes3 = [0.18+0.47, 0.1, 0.45, 0.45]
axes4 = [0.18+0.47, 0.57, 0.45, 0.45]

ax1len = axes1[0] + axes1[3]
tlen = ax1len + 0.39 + (0.55-ax1len) 

fig = plt.figure(figsize=(8*cm,8*cm))
ax1 = fig.add_axes(axes1)
ax2 = fig.add_axes(axes2)
ax3 = fig.add_axes(axes3)
ax4 = fig.add_axes(axes4)

cmap = plt.get_cmap(cmap)
colors = [cmap(i*0.16) for i in range(6)][::-1]

beams = []

for i in range(6):
    data = data1
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
#    ax1.errorbar(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["flux"]/data[i]["model"], 
#        yerr=data[i]["eflux"]/data[i]["model"], 
#        ls="-", 
#        marker="+",
#        markeredgewidth=0.25,
#        color=colors[i], 
#        lw=0.5,
#        zorder=10+i)
    ax1.plot(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["aegean"]/data[i]["model"], 
        ls="-.", 
        color=colors[i], 
        lw=0.5,
        zorder=10+i)

for i in range(6):
    data_multi = data_multi1
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
#    ax2.errorbar(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["flux"]/data_multi[i]["model"], 
#        yerr=data_multi[i]["eflux"]/data_multi[i]["model"], 
#        ls="-", 
#        lw=0.5,
#        marker="+",
#        markeredgewidth=0.25,
#        color=colors[i], 
#        zorder=10+i)
    ax2.plot(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["aegean"]/data_multi[i]["model"], 
        ls="-.", 
        lw=0.5,
        color=colors[i], 
        zorder=10+i)

for i in range(6):
    data = data2
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
#    ax3.errorbar(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["flux"]/data[i]["model"], 
#        yerr=data[i]["eflux"]/data[i]["model"], 
#        ls="-", 
#        lw=0.5,
#        marker="+",
#        markeredgewidth=0.25,
#        color=colors[i], 
#        zorder=10+i)
    ax3.plot(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["aegean"]/data[i]["model"], 
        ls="-.", 
        lw=0.5,
        color=colors[i], 
        zorder=10+i)

for i in range(6):
    data_multi = data_multi2
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
#    ax4.errorbar(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["flux"]/data_multi[i]["model"], 
#        yerr=data_multi[i]["eflux"]/data_multi[i]["model"], 
#        ls="-", 
#        lw=0.5,
#        marker="+",
#        markeredgewidth=0.25,
#        color=colors[i], 
#        zorder=10+i)
    ax4.plot(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["aegean"]/data_multi[i]["model"], 
        ls="-.", 
        lw=0.5,
        color=colors[i], 
        zorder=10+i)

for ax in ax1, ax2, ax3, ax4:
    ax.set_ylim([0.5, 1.05])
    ax.set_xlim([1., 10.])

for ax in [ax1, ax2, ax3, ax4]:
    for spacing in [0.6, 0.8, 1.0]:
        ax.axhline(spacing, ls="--", lw=1, color="grey")

#handles = [Line2D([0], [0], linestyle="-", marker="+", color=colors[0], lw=0.5, markeredgewidth=0.25, label=r"uniform"),
#           Line2D([0], [0], linestyle="-", marker="+",color=colors[1], lw=0.5, markeredgewidth=0.25, label=r"$r=0.0$"),
#           Line2D([0], [0], linestyle="-", marker="+",color=colors[2], lw=0.5, markeredgewidth=0.25, label=r"$r=+0.5$"),
#           Line2D([0], [0], linestyle="-", marker="+",color=colors[3], lw=0.5, markeredgewidth=0.25, label=r"$r=+1.0$"),
#           Line2D([0], [0], linestyle="-", marker="+",color=colors[4], lw=0.5, markeredgewidth=0.25, label=r"$r=+2.0$"),
#           Line2D([0], [0], linestyle="-", marker="+",color=colors[5], lw=0.5, markeredgewidth=0.25, label=r"natural"),
#           Line2D([0], [0], linestyle="-.", color=colors[0], lw=0.5, label=r"\texttt{aegean}, uniform"),
#           Line2D([0], [0], linestyle="-.", color=colors[1], lw=0.5, label=r"\texttt{aegean}, $r=0.0$"),
#           Line2D([0], [0], linestyle="-.", color=colors[2], lw=0.5, label=r"\texttt{aegean}, $r=+0.5$"),
#           Line2D([0], [0], linestyle="-.", color=colors[3], lw=0.5, label=r"\texttt{aegean}, $r=+1.0$"),
#           Line2D([0], [0], linestyle="-.", color=colors[4], lw=0.5, label=r"\texttt{aegean}, $r=+2.0$"),
#           Line2D([0], [0], linestyle="-.", color=colors[5], lw=0.5, label=r"\texttt{aegean}, natural"),]

handles = [Line2D([0], [0], linestyle="-.", color=colors[0], lw=0.5, label=r"uniform"),
           Line2D([0], [0], linestyle="-.", color=colors[1], lw=0.5, label=r"$r=0.0$"),
           Line2D([0], [0], linestyle="-.", color=colors[2], lw=0.5, label=r"$r=+0.5$"),
           Line2D([0], [0], linestyle="-.", color=colors[3], lw=0.5, label=r"$r=+1.0$"),
           Line2D([0], [0], linestyle="-.", color=colors[4], lw=0.5, label=r"$r=+2.0$"),
           Line2D([0], [0], linestyle="-.", color=colors[5], lw=0.5, label=r"natural")]
ax2.legend(ncol=2, bbox_to_anchor=(0.4,1.32), loc="upper left", fontsize=font_ticks-2.,
           handles=handles[::-1], handlelength=3., frameon=False)

fig.text(0.04, 0.5, r"$S_\mathrm{measured} / S_\mathrm{model}$ ", fontsize=font_labels, va="center", rotation=90.)
fig.text(0.35, 0.0, "Model Gaussian FWHM [arcmin]", fontsize=font_labels, va="center")


for ax in ax1, ax2, ax3, ax4:
    ax.set_yticks([0.6, 0.8, 1.0])
for ax in ax1, ax3:
    ax.set_xticks([2.,4.,6.,8.,10])
for ax in ax3, ax4:
    ax.axes.get_yaxis().set_ticks([])
for ax in ax2, ax4:
    ax.axes.get_xaxis().set_ticks([])

extra_name = "$20\\sigma$"
ax1.text(0.5, 0.03, "{}".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax1.transAxes)
ax2.text(0.5, 0.03, "{} multi-scale".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax2.transAxes)

extra_name = "$1000\\sigma$"
ax3.text(0.5, 0.03, "{}".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax3.transAxes)
ax4.text(0.5, 0.03, "{} multi-scale".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax4.transAxes)

fig.savefig("flux_recovery_20and1000_alt.pdf", transparent=True, dpi=300, bbox_inches="tight")
