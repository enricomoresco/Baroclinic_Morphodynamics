"""Fig. 2 of the paper.

(a) Friction factor F(delta) from Eq. (7), the 1DV k-epsilon model and ROMS. Reads
    results/kepsilon_sweep.json (run_sweep.py) and, if present, results/roms_equilibrium.json
    (roms/scripts/analyse.py) for the ROMS value of gamma.
(b, c) Response of the neutral k-epsilon model at delta = -0.2, relative to the barotropic run with
    the same bed stress: change of the turbulent length scale and of the local friction velocity
    (tau/rho0)^(1/2), and change of the velocity, also with the length scale held at its barotropic
    profile. Reads results/kepsilon_closures.json (closure_diagnostics.py).
Writes figures/fig2_friction_factor.(pdf|png).
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
    closures_file = ROOT / "results" / "kepsilon_closures.json"
    if not closures_file.exists():
        raise SystemExit("results/kepsilon_closures.json not found: run kepsilon_1dv/closure_diagnostics.py first")
    ke = json.loads(closures_file.read_text())["k_epsilon"]
    d_max = sw["runaway_threshold"]["last_steady"]

    # same style as Fig. 1 (verification/run_verification.py): 13.9 cm wide, 8-pt fonts
    cm = 1 / 2.54
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                         "lines.linewidth": 1.4, "pdf.fonttype": 42, "axes.titlelocation": "left",
                         "axes.titlepad": 3})
    blue, orange, grey, dark = "#1d4e89", "#c2410c", "0.55", "0.15"
    fig, ax = plt.subplots(1, 3, figsize=(13.9 * cm, 5.0 * cm), gridspec_kw={"width_ratios": [1.3, 1, 1]})

    # (a) friction factor
    a = ax[0]
    dd = np.linspace(d_max, 0, 200)
    a.fill_between(dd, F_analytical(dd, ALPHA["linear"]), F_analytical(dd, ALPHA["parabolic"]),
                   color=grey, alpha=0.35, lw=0, label="analytical")
    for name, mk, c, lab in [("neutral", "o-", blue, r"$k$–$\varepsilon$, neutral"),
                             ("stratified", "s-", orange, r"$k$–$\varepsilon$, stratified")]:
        rows = [r for r in sw["cases"][name] if r["delta"] >= d_max]
        a.plot([r["delta"] for r in rows], [r["F"] for r in rows], mk, color=c, ms=2.5, lw=1.2, label=lab)
    dr = np.linspace(-0.3, 0, 10)
    a.plot(dr, 1 + g_roms * dr, "--", color=dark, lw=1.0, label="ROMS")
    a.set_xlim(-0.33, 0.005)
    a.set_ylim(0.70, 1.01)
    a.set_xticks([-0.3, -0.2, -0.1, 0])
    a.set_xlabel(r"$\delta$")
    a.set_ylabel(r"$F=\tau_b/(\rho_0 c_{f}\bar u^2)$")
    a.set_title("(a) friction factor")
    a.legend(loc="lower right", handlelength=1.5, borderaxespad=0.1, labelspacing=0.25)

    # (b, c) neutral k-epsilon response at equal bed stress
    zeta = np.array(ke["profiles"]["zeta"])
    ustar_change = 100 * (np.sqrt(np.array(ke["profiles"]["stress_ratio"])) - 1)
    length_change = 100 * (np.array(ke["profiles"]["length_scale_ratio"]) - 1)
    du = np.array(ke["profiles"]["velocity_change_percent"])
    du_fixed = np.array(ke["profiles"]["velocity_change_fixed_length_scale_percent"])
    b, c = ax[1], ax[2]
    # curves labelled directly, in their colours (positions in data units, for delta = -0.2)
    b.plot(length_change, zeta, color=blue)
    b.plot(ustar_change, zeta, "--", color=dark, lw=1.0)
    b.text(4.2, 0.12, "length scale", color=blue, va="center")
    b.text(5.6, 0.86, r"$(\tau/\rho_0)^{1/2}$", color=dark, va="center")
    b.set_xlabel("change (%)")
    b.set_ylabel(r"$\zeta=z^*/D^*$")
    b.set_title("(b) turbulence")
    c.axvline(0, color=grey, lw=0.6)
    c.plot(du, zeta, color=blue)
    c.plot(du_fixed, zeta, "--", color=dark, lw=1.0)
    c.text(0.70, 0.86, r"$k$–$\varepsilon$", color=blue, va="center")
    c.text(0.55, 0.14, "fixed length\nscale", color=dark, va="center")
    c.set_xlabel(r"$\Delta u/\bar u$ (%)")
    c.set_title("(c) velocity")
    for p in (b, c):
        p.set_ylim(0, 1)
        p.set_yticks([0, 0.5, 1])
    fig.tight_layout(pad=0.3, w_pad=1.0)
    (ROOT / "figures").mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(ROOT / "figures" / f"fig2_friction_factor.{ext}", dpi=300,
                    metadata={"CreationDate": None} if ext == "pdf" else {})
    print(f"written: figures/fig2_friction_factor.pdf/.png  (ROMS gamma = {g_roms:.3f}; "
          f"depth-averaged velocity change {ke['depth_mean_change_percent']:+.3f} %, "
          f"{ke['depth_mean_change_fixed_length_scale_percent']:+.3f} % with fixed length scale)")


if __name__ == "__main__":
    main()
