"""Sec. VI B and Fig. 5: diagnostics of the simulations of Olabarrieta et al. (2018, JGR Earth Surface 123,
2901-2924), from their published data set (Figshare, doi:10.6084/m9.figshare.5975164).

Runs used: the baroclinic and barotropic runs of their Fig. 8 (river discharge 100 m3/s, tidal range 1 m,
w_s = 0.6 mm/s, 41 N, after 150 years).
- depth: FIGURE8/h_mean.mat (cross-channel averaged water depth, baroclinic and barotropic);
- salinity of the same baroclinic run: FIGURE11/salt_thal.mat, column 0 (w_s = 0.6 mm/s, vertically averaged
  salinity of their Fig. 11d). Its depth (FIGURE11/h_sec.mat, column 0) is identical to that of FIGURE8.

Not used: the time-series files of FIGURE2, 3, 5 and 6. They do not reproduce the published reference
simulation of their Figs. 1-6 (Q = 150 m3/s, tidal range 1.5 m): their M2 amplitude at the mouth is 0.35 m
instead of 0.75 m, and the transport proxy of FIGURE6/q_thal_time.mat differs in shape and size from their
Fig. 6 and from the figure files in the same folder.

x is the distance from the estuary head, as in the data set.
Usage: python3 -I olabarrieta_consistent.py <DATASET folder> <output folder>
"""
import sys
import numpy as np
import scipy.io as sio
from scipy.ndimage import uniform_filter1d

DS, OUT = sys.argv[1].rstrip("/") + "/", sys.argv[2].rstrip("/") + "/"
L = lambda f: sio.loadmat(DS + f)

hm = L("FIGURE8/h_mean.mat")
x = hm["xx"][:, 36] / 1000.0
dx = np.diff(x).mean() * 1000.0
sm = lambda f, km=10.0: uniform_filter1d(f, int(round(km * 1000 / dx)), mode="nearest")
m = x <= 150.0                         # estuary proper (the ebb-tidal delta lies beyond)

D_bc, D_bt = hm["h_mean"][:, 0], hm["h_mean"][:, 1]
assert np.allclose(L("FIGURE11/h_sec.mat")["h_mean"][:, 0], D_bc), "salinity file is not from the Fig. 8 run"
S = L("FIGURE11/salt_thal.mat")["salt_thal"][:, 0]
dD = sm(D_bc - D_bt)
Sx = sm(np.gradient(S, x))             # psu/km, > 0 (salinity increases seaward)
x_up = x[m][np.argmax(Sx[m] > 0.02 * Sx[m].max())]
flat = m & (x < x_up)
k = np.argmax(np.abs(dD[m]))
G = Sx[m] / Sx[m].max()
a = np.sum(dD[m] * G) / np.sum(G * G)  # best local fit dD = a G, amplitude and sign free
R2 = 1 - np.sum((dD[m] - a * G) ** 2) / np.sum(dD[m] ** 2)

print("RUNS OF FIG. 8 (Q = 100 m3/s, TR = 1 m, 150 yr)")
print(f"  max salinity gradient at {x[m][np.argmax(Sx[m])]:.1f} km; max |dD| at {x[m][k]:.1f} km: "
      f"{dD[m][k]:.1f} m = {100 * dD[m][k] / D_bt[m][k]:.0f}% of the barotropic depth")
print(f"  dD < 0 over {100 * np.mean(dD[m] < 0):.0f}% of the estuary")
print(f"  gradient < 2% of its maximum within {x_up:.1f} km of the head = {100 * flat.sum() / m.sum():.0f}% of the estuary")
print(f"  there: mean dD = {dD[flat].mean():.1f} m = {100 * np.mean(dD[flat] / sm(D_bt)[flat]):.0f}% of the barotropic"
      f" depth; share of the integrated |dD| = {100 * np.abs(dD[flat]).sum() / np.abs(dD[m]).sum():.0f}%")
print(f"  best local fit a*G: a = {a:.1f} m, explains {100 * R2:.0f}% of sum(dD^2)")
np.savez(OUT + "olabarrieta_consistent.npz", x=x[m], Sx=Sx[m], dD=dD[m], x_up=x_up)
