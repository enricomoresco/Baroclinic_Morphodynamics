"""Fig. 2 of the paper: friction factor F(delta) from Eq. (7), the 1DV k-epsilon model and ROMS.

Reads results/kepsilon_sweep.json (run_sweep.py) and, if present, results/roms_equilibrium.json
(roms/scripts/analyse.py) for the ROMS value of gamma. Writes figures/fig2_friction_factor.(pdf|png).
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
LAMBDA = np.log(10.0 / 1e-3)                     # D = 10 m, z0 = 1 mm
ALPHA = {"linear": (LAMBDA - 11 / 6) / (LAMBDA - 1.5), "parabolic": (LAMBDA - 1.5) / (LAMBDA - 1.0)}


def F_analytical(delta, alpha):
    return (1 + delta) ** 2 / (1 + alpha * delta) ** 2


def main():
    sw = json.loads((ROOT / "results" / "kepsilon_sweep.json").read_text())
    roms_file = ROOT / "results" / "roms_equilibrium.json"
    g_roms = json.loads(roms_file.read_text())["eps"]["0.1"]["gamma_local"] if roms_file.exists() else 0.074
    d_max = sw["runaway_threshold"]["last_steady"]

    cm = 1 / 2.54
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8, "legend.frameon": False, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 0.6, "lines.linewidth": 1.4,
                         "pdf.fonttype": 42})
    blue, orange, grey, dark = "#1d4e89", "#c2410c", "0.55", "0.15"
    fig, ax = plt.subplots(figsize=(8.5 * cm, 6.0 * cm))
    dd = np.linspace(d_max, 0, 200)
    ax.fill_between(dd, F_analytical(dd, ALPHA["linear"]), F_analytical(dd, ALPHA["parabolic"]),
                    color=grey, alpha=0.35, lw=0, label="analytical")
    for name, mk, c, lab in [("neutral", "o-", blue, r"$k$–$\varepsilon$, neutral"),
                             ("stratified", "s-", orange, r"$k$–$\varepsilon$, stratified")]:
        rows = [r for r in sw["cases"][name] if r["delta"] >= d_max]
        ax.plot([r["delta"] for r in rows], [r["F"] for r in rows], mk, color=c, ms=3, label=lab)
    dr = np.linspace(-0.3, 0, 10)
    ax.plot(dr, 1 + g_roms * dr, "--", color=dark, lw=1.0, label="ROMS")
    ax.set_xlim(-0.33, 0.005)
    ax.set_ylim(0.78, 1.01)
    ax.set_xlabel(r"$\delta$")
    ax.set_ylabel(r"$F=\tau_b/(\rho_0 c_{f}\bar u^2)$")
    ax.legend(loc="lower right", handlelength=1.8)
    fig.tight_layout(pad=0.3)
    (ROOT / "figures").mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(ROOT / "figures" / f"fig2_friction_factor.{ext}", dpi=300,
                    metadata={"CreationDate": None} if ext == "pdf" else {})
    print(f"written: figures/fig2_friction_factor.pdf/.png  (ROMS gamma = {g_roms:.3f})")


if __name__ == "__main__":
    main()
