"""Effect of a salinity relaxation (nudging) on the stratified friction factor (Appendix A of the paper).

The ROMS runs nudge the salinity on a time scale of 86 s, which limits the stratification. To mimic it,
a relaxation term -S'/tau_n is added to Eq. (A2) of the 1DV model (argument tau_n of run()), and
gamma_eff = (F - 1)/delta is computed at delta = -0.02 for the stratified case (c3 = -0.4, Pr_t = 1).

Output: results/kepsilon_nudging.json.  Run time: under a minute.
Expected (Appendix A): gamma_eff = 0.235 without relaxation and 0.218, 0.162, 0.059 for
tau_n = 3600, 600 and 86 s.
"""
import json
from pathlib import Path

from kepsilon_1dv import friction_coefficient, run

ROOT = Path(__file__).resolve().parents[1]
DELTA = -0.02
TAUS = [None, 3600.0, 600.0, 86.0]


def main():
    cf0 = friction_coefficient(run(t_end=5e4))
    rows = []
    for tau in TAUS:
        r = run(delta=DELTA, t_end=1.2e5, tau_n=tau)
        F = friction_coefficient(r) / cf0
        rows.append(dict(tau_n=tau, F=F, gamma_eff=(F - 1) / DELTA, dS_bottom_minus_top=r["dS"]))
        label = "none" if tau is None else f"{tau:.0f} s"
        print(f"tau_n = {label:>7}   F = {F:.5f}   gamma_eff = {(F - 1) / DELTA:.3f}   dS = {r['dS']:.4f} psu")
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "kepsilon_nudging.json").write_text(json.dumps(dict(delta=DELTA, cf0=cf0, runs=rows), indent=1))
    print("written: results/kepsilon_nudging.json")


if __name__ == "__main__":
    main()
