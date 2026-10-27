"""Fig. 5 of the paper from olabarrieta_consistent.npz (written by olabarrieta_consistent.py): runs of Fig. 8 of
Olabarrieta et al. (2018) after 150 years, salinity gradient and depth reduction along the estuary.
Usage: python3 -I make_fig5.py <folder with olabarrieta_consistent.npz> <output pdf>
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = np.load(sys.argv[1].rstrip("/") + "/olabarrieta_consistent.npz")
CM = 1 / 2.54
BLUE, ORANGE, GREY = "#1d4e89", "#c2410c", "0.55"
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "legend.fontsize": 8, "legend.frameon": False, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.6, "lines.linewidth": 1.4,
                     "pdf.fonttype": 42})

x, x_up = R["x"], float(R["x_up"])
fig, ax = plt.subplots(figsize=(8.5 * CM, 5.0 * CM))
ax.axvspan(0, x_up, color=GREY, alpha=0.2, lw=0)
l1, = ax.plot(x, R["Sx"] / R["Sx"].max(), color=BLUE)
l2, = ax.plot(x, R["dD"] / R["dD"].min(), color=ORANGE)
ax.set_ylabel("normalized")
ax.set_xlabel("distance from the head (km)")
ax.set_xlim(0, 150)
ax.set_ylim(0, 1.5)
ax.set_yticks([0, 0.5, 1])
ax.legend([l2, l1], ["depth reduction", "salinity gradient"], loc="upper left", handlelength=1.6,
          borderaxespad=0.2)
fig.tight_layout(pad=0.3)
fig.savefig(sys.argv[2], metadata={"CreationDate": None})
print("written", sys.argv[2])
