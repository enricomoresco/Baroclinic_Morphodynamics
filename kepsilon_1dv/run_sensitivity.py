"""Sensitivity tests of the 1DV k-epsilon model (Appendix A of the paper) and the
amplitude of the exchange flow (Sec. IV A).

- F at delta = -0.2, neutral and stratified, with 50, 100 and 200 cells;
- gamma_eff = (F - 1)/delta at delta = -0.05 (stratified) with 50, 100 and 200 cells,
  with c3 = -0.6, with Pr_t = 0.7 and with surface mixing lengths z0s = 0.02 and 0.5 m;
- exchange flow at delta = -0.1: maximum departure of u/ubar from the barotropic profile.

Output: results/kepsilon_sensitivity.json.  Run time: a few minutes.
"""
import json
from pathlib import Path

import numpy as np

from kepsilon_1dv import friction_coefficient, run

ROOT = Path(__file__).resolve().parents[1]
T_END = 1.2e5


def main():
    out = dict(grid={}, gamma_eff_at_delta_m005={}, exchange_flow_at_delta_m01={})
    ref = {}
    for N in (50, 100, 200):
        ref[N] = run(N=N, t_end=5e4)
        cf0 = friction_coefficient(ref[N])
        Fn = friction_coefficient(run(delta=-0.2, N=N, buoyancy=False, t_end=T_END)) / cf0
        Fs = friction_coefficient(run(delta=-0.2, N=N, t_end=T_END)) / cf0
        g = (friction_coefficient(run(delta=-0.05, N=N, t_end=T_END)) / cf0 - 1) / -0.05
        out["grid"][N] = dict(F_neutral_delta_m02=Fn, F_stratified_delta_m02=Fs, gamma_eff_stratified=g)
        print(f"N = {N:3d}:  F(-0.2) neutral {Fn:.4f}, stratified {Fs:.4f};  gamma_eff(-0.05) stratified {g:.3f}")

    cf0 = friction_coefficient(ref[100])
    for name, kw in {"reference (c3=-0.4, Pr_t=1, z0s=0.1 m)": {}, "c3=-0.6": dict(c3=-0.6),
                     "Pr_t=0.7": dict(prt=0.7), "z0s=0.02 m": dict(z0s=0.02), "z0s=0.5 m": dict(z0s=0.5)}.items():
        g = (friction_coefficient(run(delta=-0.05, t_end=T_END, **kw)) / cf0 - 1) / -0.05
        out["gamma_eff_at_delta_m005"][name] = g
        print(f"gamma_eff(-0.05), {name}: {g:.3f}")

    u0 = ref[100]["u"] / ref[100]["ubar"]
    for name, kw in {"neutral": dict(buoyancy=False), "stratified": {}, "stratified, c3=-0.6": dict(c3=-0.6),
                     "stratified, Pr_t=0.7": dict(prt=0.7)}.items():
        r = run(delta=-0.1, t_end=T_END, **kw)
        dep = float(np.abs(r["u"] / r["ubar"] - u0).max())
        out["exchange_flow_at_delta_m01"][name] = dep
        print(f"exchange flow at delta = -0.1, {name}: max |u/ubar - u0/ubar0| = {100 * dep:.2f}%")

    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "kepsilon_sensitivity.json").write_text(json.dumps(out, indent=1))
    print("written: results/kepsilon_sensitivity.json")


if __name__ == "__main__":
    main()
