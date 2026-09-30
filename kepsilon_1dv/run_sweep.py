"""Friction factor F(delta) from the 1DV k-epsilon model (Fig. 2 and Table II of the paper).

Cases: neutral (no buoyancy), stratified with c3 = -0.4 and Pr_t = 1 (reference), and two
stratified variants (c3 = -0.6; Pr_t = 0.7). Each run is integrated for 0.8e5 and 1.2e5 s;
the difference between the two values of F measures the distance from steady state.
Also computes the threshold beyond which the stratified runs have no steady state.

Output: results/kepsilon_sweep.json.  Run time: about ten minutes.
"""
import json
from pathlib import Path

import numpy as np

from kepsilon_1dv import friction_coefficient, run

ROOT = Path(__file__).resolve().parents[1]
DELTAS = [-0.005, -0.01, -0.02, -0.05, -0.1, -0.15, -0.2, -0.3, -0.4]
CASES = {"neutral": dict(buoyancy=False),
         "stratified": dict(buoyancy=True),
         "stratified, c3=-0.6": dict(buoyancy=True, c3=-0.6),
         "stratified, Pr_t=0.7": dict(buoyancy=True, prt=0.7)}


def steady_F(delta, cf0, **kw):
    a = run(delta=delta, t_end=8e4, **kw)
    b = run(delta=delta, t_end=1.2e5, **kw)
    Fa, Fb = friction_coefficient(a) / cf0, friction_coefficient(b) / cf0
    return Fb, abs(Fb - Fa), b


def is_steady(delta, cf0, **kw):
    """A run is steady if F settles to a finite value (runaway stratification drives F to 0)."""
    F, change, _ = steady_F(delta, cf0, **kw)
    return F > 0.5 and change < 1e-3


def main():
    r0 = run(t_end=5e4)
    cf0 = friction_coefficient(r0)
    out = dict(cf0=cf0, ubar0=r0["ubar"], deltas=DELTAS, cases={})
    for name, kw in CASES.items():
        rows = []
        for d in DELTAS:
            F, change, r = steady_F(d, cf0, **kw)
            rows.append(dict(delta=d, F=F, dF_last_4e4s=change, dS_bottom_minus_top=r["dS"]))
        out["cases"][name] = rows
        print(f"\n{name}")
        for q in rows:
            print(f"  delta = {q['delta']:+.3f}   F = {q['F']:.4f}   (F-1)/delta = {(q['F'] - 1) / q['delta']:.3f}"
                  f"   dS = {q['dS_bottom_minus_top']:.3g} psu   |dF| = {q['dF_last_4e4s']:.1e}")

    # threshold of runaway stratification (reference stratified case), by bisection
    lo, hi = -0.30, -0.40                      # steady at lo, runaway at hi
    assert is_steady(lo, cf0) and not is_steady(hi, cf0)
    for _ in range(6):
        mid = 0.5 * (lo + hi)
        if is_steady(mid, cf0):
            lo = mid
        else:
            hi = mid
    out["runaway_threshold"] = dict(last_steady=lo, first_runaway=hi)
    print(f"\nstratified runs: steady for delta >= {lo:.4f}, runaway for delta <= {hi:.4f}")
    print(f"barotropic reference: ubar = {r0['ubar']:.3f} m/s, c_f0 = {cf0:.5f}; "
          f"Eq. (8) with this c_f: gamma = {2 / 3 * np.sqrt(cf0) / 0.41:.3f} (linear), "
          f"{np.sqrt(cf0) / 0.41:.3f} (parabolic)")
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "kepsilon_sweep.json").write_text(json.dumps(out, indent=1))
    print("written: results/kepsilon_sweep.json")


if __name__ == "__main__":
    main()
