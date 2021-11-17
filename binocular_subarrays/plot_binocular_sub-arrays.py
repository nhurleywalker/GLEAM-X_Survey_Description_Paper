import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np

def S(nu, nu0, S0, alpha):
   return S0 * (nu/nu0)**alpha

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12}
)

#plt.rcParams["patch.force_edgecolor"] = True
#plt.rcParams['grid.linewidth'] = 0.1
#plt.rcParams['axes.axisbelow'] = True

# To select these antennas I used the following criteria:
# ants_east: East > 700 (43 antennas)
# ants_west: East < -120 (44 antennas)
# ants_north: North > 720 (44 antennas)
# ants_south: North < -70 (43 antennas)

dataset = pd.read_csv("MWAx.txt", delimiter=" ", header=0, comment="#")
east = np.squeeze(np.where(dataset["East"] > 700))
west = np.squeeze(np.where(dataset["East"] < -120))
north = np.squeeze(np.where(dataset["North"] > 720))
south = np.squeeze(np.where(dataset["North"] < -70))

fig = plt.figure(figsize=(5,5))
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
