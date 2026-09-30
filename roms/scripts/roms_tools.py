"""Helpers shared by the ROMS scripts: paths, amplitudes, reading the output, fits."""
import os
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

from make_run import D0, ELL, SSEA, XS, ZOB, make_run, scoef_for

HERE = Path(__file__).resolve().parent.parent           # roms/
REPO = HERE.parent
RUNS = HERE / "runs"
RESULTS = REPO / "results"

P = 2.0 / (np.log(D0 / ZOB) - 1.0)     # p = -d ln c_f / d ln D for the log law, Eq. (16)
GAMMA_REF = 0.08                        # gamma used to size the imposed local depth change: median of
                                        # (F - 1)/delta on the fixed-bed twins (printed by analyse.py)
J = 2                                   # middle row of the three-row channel
WINDOW = (5e3, 95e3)                    # fitting window (m from the sea boundary)
KS = [-1.0, 0.0, 0.5, 1.0, 1.5, 2.0]    # local depth change imposed, in units of the theory
SETUP_FACTORS = [0.6, 0.8, 1.0, 1.2, 1.4]


def amplitudes(eps, gamma=GAMMA_REF):
    """Theoretical amplitudes of the first-order bed correction for bed-load (MPM) transport.

    A_D: peak local depth reduction, Eq. (16) with delta_max = eps  -> gamma eps D0 / (2 + p)
    A_h: half of the cumulative rise of bed and surface, Eq. (15) with 1 + gamma/(2 + p)
         in place of 1 + 17 gamma/33; beta S_sea D0 / 2 is the hydrostatic head.
    """
    beta = scoef_for(eps)
    A_D = gamma / (2.0 + P) * eps * D0
    A_h = (1.0 + gamma / (2.0 + P)) * beta * SSEA * D0 / 4.0
    return A_D, A_h


def tag(eps):
    return f"{int(round(eps * 100)):02d}"


def is_done(run, runs=RUNS):
    log = Path(runs) / run / "log.txt"
    return log.exists() and "ROMS: DONE" in log.read_text(errors="ignore")


def run_roms(name, eps, AD, Ah, exe, runs=RUNS, force=False):
    """Create the run directory and run ROMS (serial executable reading roms.in from stdin)."""
    if is_done(name, runs) and not force:
        return name
    run, _ = make_run(name, eps, AD, Ah, days=3.0, out_days=0.5, runs_dir=runs)
    with open(run / "roms.in") as fin, open(run / "log.txt", "w") as fout:
        subprocess.run([str(exe)], stdin=fin, stdout=fout, stderr=subprocess.STDOUT, cwd=run, check=False)
    if not is_done(name, runs):
        raise RuntimeError(f"ROMS run {name} did not finish, see {run / 'log.txt'}")
    return name


# ------------------------------------------------------------------ reading the output
def _nc(run, runs=RUNS):
    import netCDF4
    return netCDF4.Dataset(Path(runs) / run / "his.nc")


def last(run, var, runs=RUNS):
    with _nc(run, runs) as d:
        return np.asarray(d[var][-1])


def coords(runs=RUNS, ref="eq_bt"):
    with _nc(ref, runs) as d:
        return np.asarray(d["x_rho"][J, :]), np.asarray(d["x_u"][J, :])


def stress_departure(run, ref="eq_bt", runs=RUNS):
    """tau_b / tau_b,ref - 1 along the channel axis (u points), last output."""
    return last(run, "bustr", runs)[J] / last(ref, "bustr", runs)[J] - 1.0


def nonboussinesq(run, ref="eq_bt", runs=RUNS):
    """Stress departure due to the density weight of the surface-slope term (Appendix B).

    ROMS multiplies the surface-slope term by (1000 + rho'_s)/rho0. For the same stress the
    baroclinic run needs a slope larger by 1/f, f = ratio of the two weights; integrating from
    the sea gives a depth anomaly, and the log-law friction turns it into a stress anomaly
    -(2 + p) times its ratio to the depth. Returns this anomaly and 1/f - 1 at u points.
    """
    xr, xu = coords(runs, ref)
    rho_s = last(run, "rho", runs)[-1, J]              # surface level, density anomaly (kg/m3)
    rho_s_ref = last(ref, "rho", runs)[-1, J]
    f = (1000.0 + rho_s) / (1000.0 + rho_s_ref)
    zeta = last(ref, "zeta", runs)[J]
    depth = zeta + last(ref, "bath", runs)[J]
    dD = cumulative_trapezoid(np.gradient(zeta, xr) * (1.0 / f - 1.0), xr, initial=0.0)
    return np.interp(xu, xr, -(2.0 + P) * dD / depth), np.interp(xu, xr, 1.0 / f - 1.0)


def corrected_departure(run, ref="eq_bt", runs=RUNS):
    """Stress departure with the non-Boussinesq part removed."""
    return stress_departure(run, ref, runs) - nonboussinesq(run, ref, runs)[0]


def window_mask(xu):
    return (xu > WINDOW[0]) & (xu < WINDOW[1])


def fit_shapes(xu, r):
    """Least-squares fit r = c0 + c1 sech^2((x - x_s)/ell) + c2 tanh((x - x_s)/ell) in the window."""
    m = window_mask(xu)
    arg = (xu[m] - XS) / ELL
    X = np.vstack([np.ones(m.sum()), 1.0 / np.cosh(arg) ** 2, np.tanh(arg)]).T
    return np.linalg.lstsq(X, r[m], rcond=None)[0]


def rms_departure(xu, r):
    """rms of the departure in the window, after removing its mean."""
    m = window_mask(xu)
    a = r[m] - r[m].mean()
    return float(np.sqrt(np.mean(a ** 2)))


def roms_exe(arg=None):
    exe = Path(arg) if arg else Path(os.environ.get("ROMS_EXE", HERE / "romsS"))
    if not exe.exists():
        raise SystemExit(f"ROMS executable not found: {exe}  (build it with build_roms.sh or pass --exe)")
    return exe.resolve()
