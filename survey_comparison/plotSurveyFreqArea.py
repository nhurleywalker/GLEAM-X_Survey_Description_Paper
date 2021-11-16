import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.size": 12}
)

#plt.rcParams["patch.force_edgecolor"] = True
#plt.rcParams['grid.linewidth'] = 0.1
#plt.rcParams['axes.axisbelow'] = True

dataset = pd.read_csv("RadioContinuumSurveys_downselected.csv", delimiter=",", header=0, comment="#")
#dataset = np.loadtxt("RadioContinuumSurveys_downselected.csv", delimieter=",")
#dataset = dataset.values

# plt.xkcd()


fig = plt.figure(figsize=(5,5))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
cax = fig.add_axes([0.91,0.1,0.02,0.8])
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Central frequency / MHz")
ax.set_ylabel("RMS noise / mJy beam$^{-1}$")
c = ax.scatter(dataset["mean_freq"], dataset["S-min"]/5, s = 20*dataset["resolution"], c = dataset["area"]/1000, cmap="cool", alpha=0.7)
for i, txt in enumerate(dataset["Name"]):
    if txt == "GLEAM-X":
        txt = "\\textbf{GLEAM-X}"
    ax.annotate(txt, (0.7*dataset["mean_freq"][i], 1.2*dataset["S-min"][i]/5))
ax.errorbar(dataset["mean_freq"], dataset["S-min"]/5, xerr = dataset["bandwidth"]/2, fmt="", color="k", linestyle="", lw=0.5)
cb = plt.colorbar(c, cax = cax)
cb.set_label("Sky area ($\\times1000$ sq.deg.)")
# Removes weird striping in colorbar
# https://stackoverflow.com/questions/15003353/why-does-my-colorbar-have-lines-in-it
cb.solids.set_edgecolor("face")
fig.savefig("Surveys.png", bbox_inches="tight")
fig.savefig("Surveys.pdf", bbox_inches="tight")


#plt.figure(0)
#
#date = dataset[np.where(dataset[:,3] == "dish"),1]
#sources = dataset[np.where(dataset[:,3] == "dish"),-2]
#plt.scatter(date.astype(int), sources.astype(int), marker = "o", label = "Single Dish", color="red")
#
#date = dataset[np.where(dataset[:,3] == "array"),1]
#sources = dataset[np.where(dataset[:,3] == "array"),-2]
#plt.scatter(date.astype(int), sources.astype(int), marker = "s", label = "Interferometer Array", color="deepskyblue")
#
#date = dataset[np.where(dataset[:,3] == "synth"),1]
#sources = dataset[np.where(dataset[:,3] == "synth"),-2]
#plt.scatter(date.astype(int), sources.astype(int), marker = "P", label = "Single Pixel Synthesis Array", color = "royalblue")
#
#date = dataset[np.where(dataset[:,3] == "MPA"),1]
#sources = dataset[np.where(dataset[:,3] == "MPA"),-2]
#plt.scatter(date.astype(int), sources.astype(int), marker = "^", label = "Phased Array", color="springgreen")
#
#date = dataset[np.where(dataset[:,3] == "PAF"),1]
#sources = dataset[np.where(dataset[:,3] == "PAF"),-2]
#plt.scatter(date.astype(int), sources.astype(int), marker = "d", label = "PAF Synthesis Array", color="darkgreen")
#
#date = dataset[np.where(dataset[:,3] == "cylinder"),1]
#sources = dataset[np.where(dataset[:,3] == "cylinder"),-2]
#plt.scatter(date.astype(int), sources.astype(int), marker = "*", label = "Cylinder Synthesis Array", color="orange")
#
#date = dataset[np.where(dataset[:,3] == "other"),1]
#sources = dataset[np.where(dataset[:,3] == "other"),-2]
#plt.scatter(date.astype(int), sources.astype(int), marker = "X", label = "Other", color="purple")
#
#plt.text(2013,70000000, "EMU", fontsize=6, color="darkgreen")
#plt.text(2011,40000000, "LOFAR", fontsize=6, color="springgreen")
#plt.text(2012,17500000, "Apertif", fontsize=6, color="darkgreen")
#plt.text(2012,9000000, "VLASS", fontsize=6, color="royalblue")
#plt.text(2017,550000, "TGSS", fontsize=6, color="royalblue")
#plt.text(2017,300000, "GLEAM", fontsize=6, color="springgreen")
#plt.text(2019,200000, "MIGHTEE", fontsize=6, color="royalblue")
#plt.text(2018,100000, "ASKAP-ESP3", fontsize=6, color="darkgreen")
#plt.text(2012,150000, "MSSS", fontsize=6, color="springgreen")
#plt.text(1999,1550000, "NVSS", fontsize=6, color="royalblue")
#plt.text(1999,700000, "FIRST", fontsize=6, color="royalblue")
#plt.text(1991,200000, "WENSS", fontsize=6, color="royalblue")
#plt.text(2000,100000, "SUMSS", fontsize=6, color="deepskyblue")
#plt.text(1979,18000, "MRC", fontsize=6, color="orange")
#plt.text(1969,18000, "B2", fontsize=6, color="purple")
#plt.text(1962,4000, "4C", fontsize=6, color="deepskyblue")
#plt.text(1956,3500, "MSH", fontsize=6, color="royalblue")
#plt.text(1952,1700, "2C", fontsize=6, color="deepskyblue")
#plt.text(1963,300, "3CR", fontsize=6, color="deepskyblue")
#plt.text(1955,100, "BSS", fontsize=6, color="purple")
#plt.text(1948,70, "Mills", fontsize=6, color="purple")
#plt.text(1947,45, "1C", fontsize=6, color="deepskyblue")
#plt.text(1942,8, "Reber", fontsize=6, color="red")
#plt.text(1949,8, "Bolton", fontsize=6, color="purple")
#plt.text(1947,0.9, "Hey", fontsize=6, color="purple")
#plt.text(1938,1.4, "Reber", fontsize=6, color="red")
#
#plt.xlabel("Date of first publication [year]")
#plt.ylabel("Number of Sources [count]")
#ax = plt.gca()
#ax.set_yscale("log")
#plt.tight_layout()
#plt.xlim(1935,2027)
#plt.xticks(range(1940,2030,10))
#
#
#plt.legend(frameon=False, fontsize="small")
#plt.grid()
#plt.savefig("chap1_surveys_number_sources.pdf")
#
#
#plt.figure(1)
#dataset = np.delete(dataset, 0, axis=0) #Removing first Reber survey, which is uncomfortably in the bottom right, where we want the legend
#
#plt.plot([0.001,3.5],[0.018,110000], "--", color="k")
#
#sensitivity = dataset[np.where(dataset[:,3] == "dish"),6]
#area = dataset[np.where(dataset[:,3] == "dish"),-1]
#plt.scatter(sensitivity.astype(np.float32), area.astype(np.float32), marker = "o", label = "Single Dish", color="red")
#
#sensitivity = dataset[np.where(dataset[:,3] == "array"),6]
#area = dataset[np.where(dataset[:,3] == "array"),-1]
#plt.scatter(sensitivity.astype(np.float32), area.astype(np.float32), marker = "s", label = "Interferometer Array", color="deepskyblue")
#
#sensitivity = dataset[np.where(dataset[:,3] == "synth"),6]
#area = dataset[np.where(dataset[:,3] == "synth"),-1]
#plt.scatter(sensitivity.astype(np.float32), area.astype(np.float32), marker = "P", label = "Single Pixel Synthesis Array", color = "royalblue")
#
#sensitivity = dataset[np.where(dataset[:,3] == "MPA"),6]
#area = dataset[np.where(dataset[:,3] == "MPA"),-1]
#plt.scatter(sensitivity.astype(np.float32), area.astype(np.float32), marker = "^", label = "Phased Array", color="springgreen")
#
#sensitivity = dataset[np.where(dataset[:,3] == "PAF"),6]
#area = dataset[np.where(dataset[:,3] == "PAF"),-1]
#plt.scatter(sensitivity.astype(np.float32), area.astype(np.float32), marker = "d", label = "PAF Synthesis Array", color="darkgreen")
#
#sensitivity = dataset[np.where(dataset[:,3] == "cylinder"),6]
#area = dataset[np.where(dataset[:,3] == "cylinder"),-1]
#plt.scatter(sensitivity.astype(np.float32), area.astype(np.float32), marker = "*", label = "Cylinder Synthesis Array", color="orange")
#
## sensitivity = dataset[np.where(dataset[:,3] == "other"),6]
## area = dataset[np.where(dataset[:,3] == "other"),-1]
## plt.scatter(np.log10(sensitivity.astype(np.float32)), np.log10(area.astype(np.float32)), marker = "X", label = "Other", color="purple")
#
#
#plt.text(0.0275,27500, "EMU", fontsize=6, color="darkgreen")
#plt.text(0.033,9000, "Apertif", fontsize=6, color="darkgreen")
#plt.text(0.4,30000, "VLASS", fontsize=6, color="royalblue")
#plt.text(0.57,18000, "LOFAR", fontsize=6, color="springgreen")
#plt.text(2.9,30000, "NVSS", fontsize=6, color="royalblue")
#plt.text(14,33000, "TGSS", fontsize=6, color="royalblue")
#plt.text(2.7,9000, "SUMSS", fontsize=6, color="orange")
#plt.text(7.5,9000, "WENSS", fontsize=6, color="royalblue")
#plt.text(1.15,9500, "FIRST", fontsize=6, color="royalblue")
#plt.text(0.5,2000, "ASKAP-ESP3", fontsize=6, color="darkgreen")
#plt.text(1.2,475, "LoTSS-HETDEX", fontsize=6, color="springgreen")
#plt.text(2.15,140, "ASKAP-BETA", fontsize=6, color="darkgreen")
#plt.text(0.3,130, "Stripe82", fontsize=6, color="royalblue")
#plt.text(0.27,55, "ATLAS-SPT", fontsize=6, color="royalblue")
#plt.text(0.11,39, "GLASS", fontsize=6, color="royalblue")
#plt.text(0.0032,26, "MIGHTEE", fontsize=6, color="royalblue")
#plt.text(0.055,9, "ATLAS", fontsize=6, color="royalblue")
#plt.text(0.03,6, "LHW", fontsize=6, color="royalblue")
#plt.text(0.013,1.85, "COSMOS", fontsize=6, color="royalblue")
#plt.text(0.013,1.4, "LH-Ibar", fontsize=6, color="royalblue")
#plt.text(0.012,0.65, "FLS-WSRT", fontsize=6, color="royalblue")
#plt.text(0.005,0.35, "LH-Owen", fontsize=6, color="royalblue")
#plt.text(0.0011,0.25, "Chilesconpol", fontsize=6, color="royalblue")
#plt.text(0.0038,0.18, "eMERGE", fontsize=6, color="royalblue")
#plt.text(0.0031,0.034, "Frontier Fields", fontsize=6, color="royalblue")
#
#plt.legend(frameon=False, fontsize="small")
#plt.grid()
#plt.xlabel("Limiting Sensitivity [mJy]")
#plt.ylabel(r"Area [degrees$^2$]")
#plt.tight_layout()
#plt.xlim(0.001,180)
#plt.ylim(0.03, 100000)
#ax = plt.gca()
#ax.set_yscale("log", basey=10)
#ax.set_xscale("log", basex=10)
#
#
#plt.savefig("chap1_surveys_sensitivity.pdf")
