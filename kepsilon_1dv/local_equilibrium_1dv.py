"""Local equilibrium depth with the friction law of the stratified 1DV k-epsilon model (Sec. IV C of the paper).

Eqs. (9)-(10) of the paper are algebraic at each x:
    sediment:   D = F(delta)^(5/11)
    momentum:   h' (1 + delta) = D^(-10/3) F(delta)
    definition: delta = delta0 D / h',   delta0 = eps G (reference-state value)
Eliminating D and h' gives  delta = delta0 (1 + delta) F(delta)^(32/33),  solved here for each delta0.

F(delta) of the stratified run (c3 = -0.4, Pr_t = 1) is read from results/kepsilon_sweep.json (written by
run_sweep.py) and interpolated with a monotone cubic (PCHIP), with F(0) = 1. The first-order shoaling,
5 gamma |delta0| / 11 with gamma = 0.25, is given for comparison.

Usage: python kepsilon_1dv/local_equilibrium_1dv.py [path to kepsilon_sweep.json]
Output: results/local_equilibrium_1dv.json.  Run time: seconds.
Expected: shoaling 3.89% at delta0 = -0.3 (local delta -0.216) and 8.48% at delta0 = -0.5 (local delta -0.293).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "kepsilon_sweep.json"
    rows = [r for r in json.loads(src.read_text())["cases"]["stratified"] if r["delta"] >= -0.3]
    dl = np.array([0.0] + [r["delta"] for r in rows])[::-1]
    Fv = np.array([1.0] + [r["F"] for r in rows])[::-1]
    F = PchipInterpolator(dl, Fv)
    out = []
    print(" delta0 | local delta | F(local) | shoaling, 1DV law | shoaling, first order (gamma = 0.25)")
    for d0 in (-0.05, -0.1, -0.2, -0.3, -0.4, -0.5):
        dloc = brentq(lambda x: x - d0 * (1 + x) * F(x) ** (32 / 33), -0.3, 0.0)
        Floc = float(F(dloc))
        shoal, first = 1 - Floc ** (5 / 11), 5 / 11 * 0.25 * abs(d0)
        out.append(dict(delta0=d0, delta_local=dloc, F_local=Floc, shoaling=shoal, shoaling_first_order=first))
        print(f" {d0:+.2f}  |  {dloc:+.3f}     |  {Floc:.4f}  |  {100 * shoal:5.2f}%           |  {100 * first:5.2f}%")
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "local_equilibrium_1dv.json").write_text(json.dumps(out, indent=1))
    print("written: results/local_equilibrium_1dv.json")


if __name__ == "__main__":
    main()
