"""End of the stratified branch of the 1DV model (Sec. IV A and Appendix A of the paper).

Near the threshold the approach to steady state is slow, so the runs are integrated for 5e5 and 6e5 s
and the friction factor F is compared between the two. A run is steady if F no longer changes; beyond
the threshold F keeps decreasing (runaway stratification).

Output: results/kepsilon_threshold.json.  Run time: about four minutes.
Expected (Appendix A): F = 0.502 at delta = -0.328 and 0.212 at delta = -0.333, both steady, with a
salinity difference between bed and surface of 7.1 and 46 psu; at delta = -0.334 F keeps decreasing.
"""
import json
from pathlib import Path

from kepsilon_1dv import friction_coefficient, run

ROOT = Path(__file__).resolve().parents[1]
DELTAS = [-0.328, -0.333, -0.334]


def main():
    cf0 = friction_coefficient(run(t_end=5e4))
    rows = []
    for d in DELTAS:
        a, b = run(delta=d, t_end=5e5), run(delta=d, t_end=6e5)
        Fa, Fb = friction_coefficient(a) / cf0, friction_coefficient(b) / cf0
        steady = bool(abs(Fb - Fa) < 1e-4)
        rows.append(dict(delta=d, F_5e5=float(Fa), F_6e5=float(Fb), steady=steady,
                         dS_bottom_minus_top=float(b["dS"])))
        print(f"delta = {d:+.3f}   F(5e5 s) = {Fa:.4f}   F(6e5 s) = {Fb:.4f}   dS = {b['dS']:.1f} psu   "
              f"{'steady' if steady else 'no steady state'}")
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "kepsilon_threshold.json").write_text(json.dumps(dict(cf0=cf0, runs=rows), indent=1))
    print("written: results/kepsilon_threshold.json")


if __name__ == "__main__":
    main()
