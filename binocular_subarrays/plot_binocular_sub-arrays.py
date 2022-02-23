import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import sys

def S(nu, nu0, S0, alpha):
   return S0 * (nu/nu0)**alpha

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 8}
)

cm = 1 / 2.54

# These antenna files are txt outputs from the metafits file containing the antenna layouts
# I have confirmed that the ordering in the metafits file matches the ANTENNA table in the
# measurement set for the same obsid, as displayed in casabrowser
if sys.argv[1]:
    intxt = sys.argv[1]
else:
    intxt = "1204310304_antennas.txt"

# To select these antennas I used the following criteria:
# ants_east: East > 700 (43 antennas)
# ants_west: East < -120 (44 antennas)
# ants_north: North > 720 (44 antennas)
# ants_south: North < -70 (43 antennas)

#dataset = pd.read_csv("MWAx.txt", delimiter=" ", header=0, comment="#")
#dataset = pd.read_csv("1286466120_antennas.txt", delimiter=" ", header=0, comment="#")
dataset = pd.read_csv(intxt, delimiter=" ", header=0, comment="#")
east = np.squeeze(np.where(dataset["East"] > 700))
west = np.squeeze(np.where(dataset["East"] < -120))
north = np.squeeze(np.where(dataset["North"] > 720))
south = np.squeeze(np.where(dataset["North"] < -70))

flags_east =  np.squeeze(np.where(dataset["East"] <= 700))
with open("flags_MWALB_east.txt", "w") as f:
    for e in flags_east:
        f.write(f"{e}\n")
flags_west =  np.squeeze(np.where(dataset["East"] >= -120))
with open("flags_MWALB_west.txt", "w") as f:
    for e in flags_west:
        f.write(f"{e}\n")
flags_north = np.squeeze(np.where(dataset["North"] <= 720))
with open("flags_MWALB_north.txt", "w") as f:
    for e in flags_north:
        f.write(f"{e}\n")
flags_south = np.squeeze(np.where(dataset["North"] >= -70))
with open("flags_MWALB_south.txt", "w") as f:
    for e in flags_south:
        f.write(f"{e}\n")

makeplot = True
if makeplot == True:
    fig = plt.figure(figsize=(8*cm,8*cm))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (km)")
    ax.set_ylabel("Northing (km)")
    ax.scatter(dataset["East"]/1.e3, dataset["North"]/1.e3, marker=".", color="k", s = 2)
    # Colors selected using https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=5
    ax.scatter(dataset["East"][east]/1.e3, dataset["North"][east]/1.e3, marker="o", edgecolor="#66a61e", facecolor="none", s = 10)
    ax.scatter(dataset["East"][west]/1.e3, dataset["North"][west]/1.e3, marker="o", edgecolor="#7570b3", facecolor="none", s = 10)
    ax.scatter(dataset["East"][north]/1.e3, dataset["North"][north]/1.e3, marker="s", edgecolor="#e7298a", facecolor="none", s = 30)
    ax.scatter(dataset["East"][south]/1.e3, dataset["North"][south]/1.e3, marker="s", edgecolor="#d95f02", facecolor="none", s = 30)

    fig.savefig("binocular_layout.pdf", bbox_inches="tight")
