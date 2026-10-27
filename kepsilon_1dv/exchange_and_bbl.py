"""Two diagnostics of the 1DV runs quoted in Secs. IV A and V B of the paper.

1. Salt carried by the exchange flow (Sec. V B). At steady state the landward salt flux per unit
   width is F_ex = int (u - ubar) S' dz (negative = landward; x points along the river flow). It is
   expressed as an equivalent longitudinal dispersion coefficient K_ex = -F_ex / (D dS/dx) and compared
   with the coefficient K_req = ubar S_mean / (dS/dx) that balances the seaward salt flux of the river
   for a depth-averaged salinity S_mean = 15 psu (half the sea salinity of the cases of Sec. V B).
2. Height of the bottom boundary layer (Sec. IV A): lowest height above 0.05 D at which the flux
   Richardson number Rf = -B/P reaches 0.2, compared with the scaling h/D = (Rf/Ri_x)^(1/2), Rf = 0.2,
   of Stacey and Ralston (2005, J. Phys. Oceanogr. 35, 55-71), with Ri_x = 2|delta|/(1 + delta).

Output: results/kepsilon_exchange_bbl.json.  Run time: about two minutes.
Expected: at delta = -0.1 (stratified) K_ex = 5.0 m2/s against K_req = 2.3e4 m2/s;
h/D = 1.00, 0.92, 0.70, 0.49 at delta = -0.05, -0.1, -0.2, -0.3, against 1.00, 0.95, 0.63, 0.48.
"""
import json
from pathlib import Path

import numpy as np

from kepsilon_1dv import friction_coefficient, run

ROOT = Path(__file__).resolve().parents[1]
DELTAS = [-0.05, -0.1, -0.15, -0.2, -0.25, -0.3]
S_MEAN = 15.0          # psu
RF_C = 0.2             # critical flux Richardson number


def diagnostics(r, cf0):
    z, u, s, D = r["z"], r["u"], r["s"], r["D"]
    dz = z[1] - z[0]
    F_ex = float(np.sum((u - r["ubar"]) * s) * dz)
    K_ex = -F_ex / (D * r["Sx"])
    K_req = r["ubar"] * S_MEAN / r["Sx"]
    Rf = np.where(r["P"] > 0, -r["B"] / np.maximum(r["P"], 1e-30), 0.0)
    above = np.where((Rf >= RF_C) & (z > 0.05 * D))[0]
    h = float(z[above[0]] / D) if len(above) else 1.0
    Rix = 2 * abs(r["delta"]) / (1 + r["delta"])
    return dict(delta=r["delta"], F=friction_coefficient(r) / cf0, dS_bottom_minus_top=r["dS"],
                F_ex=F_ex, K_ex=K_ex, K_req=K_req, Ri_x=Rix, h_over_D=h,
                h_over_D_SR05=float(min(1.0, np.sqrt(RF_C / Rix))))


def main():
    cf0 = friction_coefficient(run(t_end=5e4))
    out = dict(S_mean=S_MEAN, Rf_c=RF_C, cases={})
    for name, kw in (("neutral", dict(buoyancy=False)), ("stratified", dict(buoyancy=True))):
        rows = [diagnostics(run(delta=d, t_end=1.2e5, **kw), cf0) for d in DELTAS]
        out["cases"][name] = rows
        print(f"\n{name}")
        for q in rows:
            line = (f"  delta = {q['delta']:+.2f}  Ri_x = {q['Ri_x']:.2f}  F = {q['F']:.4f}  dS = {q['dS_bottom_minus_top']:.3f} psu"
                    f"  K_ex = {q['K_ex']:5.1f} m2/s  K_req = {q['K_req']:.2e} m2/s")
            if name == "stratified":
                line += f"  h/D = {q['h_over_D']:.2f} (scaling {q['h_over_D_SR05']:.2f})"
            print(line)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "kepsilon_exchange_bbl.json").write_text(json.dumps(out, indent=1))
    print("written: results/kepsilon_exchange_bbl.json")


if __name__ == "__main__":
    main()
