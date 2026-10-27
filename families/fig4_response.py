"""Response to prescribed salinity distributions (Sec. V B and Fig. 4 of the paper).

Reference river: depth D_u = 10 m, velocity U_u = 1 m/s, roughness z0 = 1 mm; friction coefficient of
the depth-averaged logarithmic law with linear eddy viscosity, c_f = (kappa/(ln(D_u/z0) - 3/2))^2, and
slope s = c_f U_u^2/(g D_u) (Eq. 1), about 2.9e-5. Sea salinity 30 psu, haline contraction
beta = 7.7e-4 /psu, salt-intrusion length L_s = 40 km, hence eps_L = beta D_u S_sea/(2 s L_s) = 0.10.

Salinity distributions (x landward from the mouth):
  Van der Burgh family, Eq. (19): S = S_sea (1 - x/L_s)^(1/K) for x < L_s, with K = 0.6 and K = 1 (linear);
  front: S = S_sea [1 - tanh((x - x_f)/L_f)]/2, x_f = 20 km, L_f = 8 km.
First-order response, Eq. (14): local depth change Delta D = (5/11) gamma delta_0 D_u, with
delta_0 = beta D_u (dS/dx)/(2 s); set-up Delta h = (1 + 17 gamma/33)(beta D_u/2)[S(0) - S(x)].

Usage: python fig4_response.py [output.pdf]     (default: figures/fig4_families.pdf)
Expected (Sec. V B): max |delta_0| = 0.17, 0.10, 0.25; max shoaling per unit gamma 76, 46, 114 cm;
set-up 12.0 cm (gamma = 0.074) and 13.3 cm (gamma = 0.3), hydrostatic head 11.6 cm, about 1% less for the front.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
G, KAPPA, BETA = 9.81, 0.41, 7.7e-4
D, U, Z0, S_SEA, LS = 10.0, 1.0, 1e-3, 30.0, 40e3
CF = (KAPPA / (np.log(D / Z0) - 1.5)) ** 2
SLOPE = CF * U ** 2 / (G * D)
x = np.linspace(0.0, 60e3, 6001)


def van_der_burgh(K):
    return S_SEA * np.clip(1 - x / LS, 0, None) ** (1 / K)


def front(xf=20e3, lf=8e3):
    return S_SEA * (1 - np.tanh((x - xf) / lf)) / 2


def response(S, gamma):
    delta0 = BETA * D * np.gradient(S, x) / (2 * SLOPE)
    dD = 5 / 11 * gamma * delta0 * D
    dh = (1 + 17 / 33 * gamma) * BETA * D / 2 * (S[0] - S)
    return delta0, dD, dh


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "figures" / "fig4_families.pdf"
    cases = [("Van der Burgh, K = 0.6", van_der_burgh(0.6)), ("linear, K = 1", van_der_burgh(1.0)),
             ("front", front())]
    print(f"c_f = {CF:.5f}, s = {SLOPE:.2e}, eps_L = {BETA * D * S_SEA / (2 * SLOPE * LS):.3f}, "
          f"hydrostatic head = {100 * BETA * S_SEA * D / 2:.1f} cm")
    fig, ax = plt.subplots(3, 1, figsize=(3.4, 6.0), sharex=True)
    for name, S in cases:
        delta0, dD1, _ = response(S, 1.0)
        _, _, h074 = response(S, 0.074)
        _, _, h03 = response(S, 0.3)
        print(f"{name:24s} max|delta0| = {np.abs(delta0).max():.3f}   max shoaling per unit gamma = "
              f"{-100 * dD1.min():.0f} cm   set-up = {100 * h074[-1]:.1f} cm (gamma 0.074), "
              f"{100 * h03[-1]:.1f} cm (gamma 0.3)")
        ax[0].plot(x / 1e3, S, label=name)
        ax[1].plot(x / 1e3, 100 * dD1)
        ax[2].plot(x / 1e3, 100 * h074)
    ax[2].axhline(100 * BETA * S_SEA * D / 2, color="0.4", ls=":", lw=1)
    ax[0].set_ylabel("salinity (psu)")
    ax[1].set_ylabel(r"$\Delta D/\gamma$ (cm)")
    ax[2].set_ylabel(r"set-up (cm)")
    ax[2].set_xlabel("distance from the mouth (km)")
    ax[0].legend(frameon=False, fontsize=7)
    for a, lab in zip(ax, "abc"):
        a.set_title(f"({lab})", loc="left", fontsize=9)
    fig.tight_layout()
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out)
    print("written:", out)


if __name__ == "__main__":
    main()
