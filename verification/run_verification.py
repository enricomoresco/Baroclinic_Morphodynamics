"""Verification of the first-order solution (Sec. III B and Fig. 1 of the paper).

Solves Eqs. (9)-(10) exactly for a hyperbolic-tangent density transition and compares
the result with the first-order solution, Eq. (13), for 0.005 <= eps <= 0.5 and both
eddy-viscosity shapes. Reference river: D = 10 m, z0 = 1 mm, so Lambda = ln(D/z0) = 9.21.

Output: results/verification.json and figures/fig1_verification.(pdf|png).
Run time: a few seconds.
"""
import json
from pathlib import Path

import numpy as np

from equilibrium_model import closure, first_order, solve_full

ROOT = Path(__file__).resolve().parents[1]
LAMBDA = np.log(10.0 / 1e-3)
EPS = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5]
EPS_FIG = 0.2

# density transition rho ~ 1 - tanh((x - 1)/0.25); G = normalized gradient, min(G) = -1
x = np.linspace(0.0, 3.0, 1501)
G = -1.0 / np.cosh((x - 1.0) / 0.25) ** 2


def main():
    out = dict(Lambda=LAMBDA, x_transition=1.0, width=0.25, shapes={})
    fields = {}
    for shape in ("linear", "parabolic"):
        c = closure(LAMBDA, shape)
        D1, h1, eta1 = first_order(x, G, c["gamma"])
        rows = []
        for e in EPS:
            D, h = solve_full(x, G, e, c["alpha"])
            errD = np.abs(D - 1 - e * D1).max()
            errh = np.abs(h - x - e * h1).max()
            erreta = np.abs((h - D) - (x - 1) - e * eta1).max()
            rows.append(dict(eps=e, err_D=errD, err_h=errh, err_eta=erreta,
                             rel_err_depth_correction=errD / np.abs(e * D1).max()))
            if shape == "linear" and e == EPS_FIG:
                fields = dict(D=D, h=h, D1=D1, eta1=eta1)
        small = np.log(EPS[:4])
        order = {k: float(np.polyfit(small, np.log([r[k] for r in rows[:4]]), 1)[0])
                 for k in ("err_D", "err_h", "err_eta")}
        out["shapes"][shape] = dict(alpha=c["alpha"], gamma=c["gamma"], cf=c["cf"],
                                    rows=rows, observed_order=order)
        print(f"\n{shape} eddy viscosity: alpha = {c['alpha']:.4f}, gamma = {c['gamma']:.4f}")
        print("  observed order of the error: " + ", ".join(f"{k[4:]} {v:.2f}" for k, v in order.items()))
        print("  eps     max err D   max err h   rel. err depth correction")
        for r in rows:
            print(f"  {r['eps']:<6}  {r['err_D']:.2e}    {r['err_h']:.2e}    {100 * r['rel_err_depth_correction']:5.2f} %")

    (ROOT / "results").mkdir(exist_ok=True)
    with open(ROOT / "results" / "verification.json", "w") as f:
        json.dump(out, f, indent=1)
    plot(out, fields)


def plot(out, fields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cm = 1 / 2.54
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                         "lines.linewidth": 1.4, "pdf.fonttype": 42, "axes.titlelocation": "left",
                         "axes.titlepad": 3})
    blue, orange, dark = "#1d4e89", "#c2410c", "0.15"
    rows = out["shapes"]["linear"]["rows"]
    E = np.array([r["eps"] for r in rows])
    eD = np.array([r["err_D"] for r in rows])
    eh = np.array([r["err_h"] for r in rows])
    D, h = fields["D"], fields["h"]

    fig, ax = plt.subplots(1, 3, figsize=(13.9 * cm, 5.0 * cm))
    ax[0].plot(x, ((h - D) - (x - 1)) / EPS_FIG, color=blue, lw=2.2)
    ax[0].plot(x, fields["eta1"], "--", color=dark, lw=1.0)
    ax[0].set_ylabel(r"$(\eta-\eta_0)/\epsilon$")
    ax[0].set_title("(a) bed elevation")
    ax[1].plot(x, (D - 1) / EPS_FIG, color=blue, lw=2.2)
    ax[1].plot(x, fields["D1"], "--", color=dark, lw=1.0)
    ax[1].set_ylabel(r"$(D-D_0)/\epsilon$")
    ax[1].set_title("(b) flow depth")
    for a in ax[:2]:
        a.set_xlabel(r"$x$")
        a.set_xticks([0, 1, 2, 3])
    ax[2].loglog(E, eD, "o-", color=blue, ms=3, label=r"$D$")
    ax[2].loglog(E, eh, "s-", color=orange, ms=3, label=r"$h$")
    ax[2].loglog(E, eh[0] * (E / E[0]) ** 2, ":", color=dark, lw=1.0, label="slope 2")
    ax[2].set_xlabel(r"$\epsilon$")
    ax[2].set_ylabel("maximum error")
    ax[2].set_title("(c) convergence")
    ax[2].set_ylim(5e-8, 1e-1)
    ax[2].legend(loc="upper left", handlelength=1.6, borderaxespad=0.2)
    fig.tight_layout(pad=0.3, w_pad=1.0)
    (ROOT / "figures").mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(ROOT / "figures" / f"fig1_verification.{ext}", dpi=300,
                    metadata={"CreationDate": None} if ext == "pdf" else {})
    print("\nwritten: results/verification.json, figures/fig1_verification.pdf/.png")


if __name__ == "__main__":
    main()
