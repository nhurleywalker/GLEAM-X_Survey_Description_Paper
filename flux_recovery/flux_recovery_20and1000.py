#! /usr/bin/env python

import numpy as np
from argparse import ArgumentParser
from subprocess import Popen


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import rgb2hex, LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

import cmasher as cmr

from matplotlib import rc

import sys

try:
    cmap = sys.argv[1]
except IndexError:
    cmap = "gnuplot2"

rc('font', **{'family':'serif', 'serif':['Times'], 'weight':'medium'})
# rc("font", **{"family": "sans-serif", "sans-serif":["Helvetica"], "weight":"medium"})
rc('text', usetex=True)
params = {"text.latex.preamble": [r"\usepackage{siunitx}", \
          r"\sisetup{detect-family = true}"]}
plt.rcParams.update(params)

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

font_labels = 22
font_ticks = 20

axes1 = [0.15, 0.1, 0.39*(8/14), 0.39]
axes2 = [0.15, 0.51, 0.39*(8/14), 0.39]
axes3 = [0.15+(0.39*(8/14))+0.02*8/14, 0.1, 0.39*(8/14), 0.39]
axes4 = [0.15+(0.39*(8/14))+0.02*8/14, 0.51, 0.39*(8/14), 0.39]

ax1len = axes1[0] + axes1[3]
tlen = ax1len + 0.39 + (0.55-ax1len) 

fig = plt.figure(figsize=(14,8))
ax1 = plt.axes(axes1)
ax2 = plt.axes(axes2)
ax3 = plt.axes(axes3)
ax4 = plt.axes(axes4)

cmap = plt.get_cmap(cmap)
colors = [cmap(i*0.16) for i in range(6)][::-1]

beams = []

for i in range(6):
    data = data1
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
    ax1.errorbar(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["flux"]/data[i]["model"], 
        yerr=data[i]["eflux"]/data[i]["model"], 
        ls="-", 
        marker="+",
        color=colors[i], 
        zorder=10+i)
    ax1.plot(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["aegean"]/data[i]["model"], 
        ls="-.", 
        color=colors[i], 
        zorder=10+i)

for i in range(6):
    data_multi = data_multi1
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
    ax2.errorbar(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["flux"]/data_multi[i]["model"], 
        yerr=data_multi[i]["eflux"]/data_multi[i]["model"], 
        ls="-", 
        marker="+",
        color=colors[i], 
        zorder=10+i)
    ax2.plot(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["aegean"]/data_multi[i]["model"], 
        ls="-.", 
        color=colors[i], 
        zorder=10+i)

for i in range(6):
    data = data2
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
    ax3.errorbar(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["flux"]/data[i]["model"], 
        yerr=data[i]["eflux"]/data[i]["model"], 
        ls="-", 
        marker="+",
        color=colors[i], 
        zorder=10+i)
    ax3.plot(np.sqrt(data[i]["radius"]**2 - bmaj**2), data[i]["aegean"]/data[i]["model"], 
        ls="-.", 
        color=colors[i], 
        zorder=10+i)

for i in range(6):
    data_multi = data_multi2
    bmaj = np.sqrt(data[i]["radius"][0]**2 - 1)
    beams.append(bmaj)
    ax4.errorbar(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["flux"]/data_multi[i]["model"], 
        yerr=data_multi[i]["eflux"]/data_multi[i]["model"], 
        ls="-", 
        marker="+",
        color=colors[i], 
        zorder=10+i)
    ax4.plot(np.sqrt(data_multi[i]["radius"]**2 - bmaj**2), data_multi[i]["aegean"]/data_multi[i]["model"], 
        ls="-.", 
        color=colors[i], 
        zorder=10+i)



ax1.set_ylim([0.5, 1.05])
ax2.set_ylim([0.5, 1.05])
ax2.set_xlim([1., 10.])
ax1.set_xlim([1., 10.])
ax3.set_ylim([0.5, 1.05])
ax4.set_ylim([0.5, 1.05])
ax4.set_xlim([1., 10.])
ax3.set_xlim([1., 10.])

for ax in [ax1, ax2, ax3, ax4]:
    for spacing in [0.6, 0.8, 1.0]:
        ax.axhline(spacing, ls="--", color="grey")


handles = [Line2D([0], [0], linestyle="-", marker="+", color=colors[0], lw=2., label=r"uniform"),
           Line2D([0], [0], linestyle="-", marker="+",color=colors[1], lw=2., label=r"$r=0.0$"),
           Line2D([0], [0], linestyle="-", marker="+",color=colors[2], lw=2., label=r"$r=+0.5$"),
           Line2D([0], [0], linestyle="-", marker="+",color=colors[3], lw=2., label=r"$r=+1.0$"),
           Line2D([0], [0], linestyle="-", marker="+",color=colors[4], lw=2., label=r"$r=+2.0$"),
           Line2D([0], [0], linestyle="-", marker="+",color=colors[5], lw=2., label=r"natural"),
           Line2D([0], [0], linestyle="-.", color=colors[0], lw=2., label=r"\texttt{aegean}, uniform"),
           Line2D([0], [0], linestyle="-.", color=colors[1], lw=2., label=r"\texttt{aegean}, $r=0.0$"),
           Line2D([0], [0], linestyle="-.", color=colors[2], lw=2., label=r"\texttt{aegean}, $r=+0.5$"),
           Line2D([0], [0], linestyle="-.", color=colors[3], lw=2., label=r"\texttt{aegean}, $r=+1.0$"),
           Line2D([0], [0], linestyle="-.", color=colors[4], lw=2., label=r"\texttt{aegean}, $r=+2.0$"),
           Line2D([0], [0], linestyle="-.", color=colors[5], lw=2., label=r"\texttt{aegean}, natural"),]
ax2.legend(ncol=2, bbox_to_anchor=(0,1.8), loc="upper left", fontsize=font_ticks-2.,
           handles=handles[::-1], handlelength=3., frameon=False)

fig.text(0.06, 0.5, r"$S_\mathrm{measured} / S_\mathrm{model}$ ", fontsize=font_labels, va="center", rotation=90.)
fig.text(0.23, 0.0, "Model Gaussian FWHM [arcmin]", fontsize=font_labels, va="center")

ax1.tick_params(axis="y", which="both", labelsize=font_ticks, 
                         pad=7.)
ax1.tick_params(axis="y", which="major", length=10.)
ax1.tick_params(axis="y", which="minor", length=5.)
ax2.tick_params(axis="y", labelsize=font_ticks, pad=7.)
ax2.tick_params(which="major", labelsize=font_labels, length=8.)
ax2.tick_params(which="minor", labelsize=font_labels, length=5.)
ax1.tick_params(axis="x", which="major", labelsize=font_ticks, 
                            length=8., pad=7.)

ax2.set_yticks([0.6, 0.8, 1.0])
ax2.set_yticks([0.6, 0.8, 1.0])
ax1.set_xticks([2.,4.,6.,8.,10])
ax3.set_xticks([2.,4.,6.,8.,10.])
extra_name = "$20\\sigma$"
ax1.text(0.5, 0.03, "{}".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax1.transAxes)
ax2.text(0.5, 0.03, "{} multi-scale".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax2.transAxes)
ax2.axes.get_xaxis().set_ticks([])

ax4.tick_params(axis="y", which="both", labelsize=font_ticks, 
                         pad=7.)
ax4.tick_params(axis="y", which="major", length=10.)
ax4.tick_params(axis="y", which="minor", length=5.)
ax3.tick_params(axis="y", labelsize=font_ticks, pad=7.)
ax3.tick_params(which="major", labelsize=font_labels, length=8.)
ax3.tick_params(which="minor", labelsize=font_labels, length=5.)
ax4.tick_params(axis="x", which="major", labelsize=font_ticks, 
                            length=8., pad=7.)

ax4.set_yticks([0.6, 0.8, 1.0])
ax3.set_yticks([0.6, 0.8, 1.0])
extra_name = "$1000\\sigma$"
ax3.text(0.5, 0.03, "{}".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax3.transAxes)
ax4.text(0.5, 0.03, "{} multi-scale".format(extra_name), fontsize=font_labels, ha="center", va="bottom", transform=ax4.transAxes)
ax4.axes.get_xaxis().set_ticks([])
ax4.axes.get_yaxis().set_ticks([])
ax3.axes.get_yaxis().set_ticks([])

fig.savefig("flux_recovery_20and1000_alt.pdf", transparent=True, dpi=300, bbox_inches="tight")