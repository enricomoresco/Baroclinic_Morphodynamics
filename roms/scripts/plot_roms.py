"""Fig. 3 of the paper from results/roms_equilibrium.json (written by analyse.py).

Writes figures/fig3_roms_equilibrium.(pdf|png).
"""
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from roms_tools import REPO, RESULTS


def main():
    R = json.loads((RESULTS / "roms_equilibrium.json").read_text())
    cm = 1 / 2.54
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "legend.fontsize": 8, "legend.frameon": False,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
                         "lines.linewidth": 1.4, "pdf.fonttype": 42, "axes.titlelocation": "left",
                         "axes.titlepad": 3})
    blue, orange, grey = "#1d4e89", "#c2410c", "0.55"
    fig, ax = plt.subplots(1, 3, figsize=(13.9 * cm, 5.2 * cm), gridspec_kw=dict(width_ratios=[1.35, 1, 1]))

    x = np.array(R["x_km"])
    e1 = R["eps"]["0.1"]
    for nm, c, ls, lab in [("flat", grey, "-", "uniform bed"), ("k+0.0", blue, "--", "cumulative"),
                           ("k+1.0", orange, "-", "first order")]:
        ax[0].plot(x, 100 * np.array(e1["departure"][nm]), color=c, ls=ls, label=lab)
    ax[0].axhline(0, color="0.7", lw=0.6)
    ax[0].set_xlabel(r"$x^*$ (km)")
    ax[0].set_ylabel(r"$\tau_b/\tau_{b,\mathrm{ref}}-1$ (%)")
    ax[0].set_title(r"(a) stress mismatch, $\epsilon=0.1$")
    ax[0].set_ylim(-1.0, 2.1)
    ax[0].legend(loc="upper right", handlelength=1.6, borderaxespad=0.2)

    for key, mk, c in [("0.1", "o", blue), ("0.3", "s", orange)]:
        e = R["eps"][key]
        ax[1].plot(e["k"], 100 * np.array(e["c_sech2"]) / float(key), mk + "-", color=c, ms=3,
                   label=rf"$\epsilon={key}$")
    ax[1].axhline(0, color="0.7", lw=0.6)
    ax[1].axvline(1, color=grey, lw=0.6, ls=":")
    ax[1].set_xlabel(r"$D_1$ / theory")
    ax[1].set_ylabel(r"local mismatch / $\epsilon$ (%)")
    ax[1].set_title("(b) local part")
    ax[1].legend(loc="lower right", handlelength=1.6, borderaxespad=0.2, frameon=True,
                 facecolor="white", edgecolor="none", framealpha=1)

    s = R["setup_scan"]
    ax[2].plot(s["factor"], 100 * np.array(s["c_tanh"]), "o-", color=blue, ms=3)
    ax[2].axhline(0, color="0.7", lw=0.6)
    ax[2].axvline(1, color=grey, lw=0.6, ls=":")
    ax[2].set_xlabel("set-up / theory")
    ax[2].set_ylabel("cumulative mismatch (%)")
    ax[2].set_title("(c) cumulative part")
    fig.tight_layout(pad=0.3, w_pad=1.0)
    out = REPO / "figures"
    out.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig3_roms_equilibrium.{ext}", dpi=300,
                    metadata={"CreationDate": None} if ext == "pdf" else {})
    print("written: figures/fig3_roms_equilibrium.pdf/.png")


if __name__ == "__main__":
    main()
