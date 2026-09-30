"""Create a ROMS run directory for the river-dominated channel (equilibrium test).

usage:  python make_run.py NAME EPS [--AD A_D] [--Ah A_h] [--days 3] [--out-days 0.5]

NAME  run name; the directory is roms/runs/NAME
EPS   peak ratio of baroclinic to barotropic pressure gradient in the uniform flow,
      epsilon in Eq. (11) of the paper (0 gives the barotropic twin)
--AD  amplitude (m) of the local depth reduction imposed on the bed, sech^2 shape
--Ah  half of the cumulative rise (m) of bed and free surface, tanh shape

The ROMS input files are generated from the templates distributed with ROMS
(ROMS/External/roms_estuary_test.in and sediment_estuary_test.in), changing only
the entries listed below. The ROMS source directory is taken from --roms-src, or
from the environment variable ROMS_SRC, or from $ROMS_ROOT_DIR/roms.
"""
import argparse
import json
import math
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent          # .../roms

# ---- physical set-up (Sec. IV B and Appendix B of the paper)
D0, U = 10.0, 1.0                     # uniform-flow depth (m) and river velocity (m/s)
ZOB = 0.001                           # bed roughness length z0 (m)
KAPPA, G = 0.41, 9.81
CF = (KAPPA / (math.log(D0 / ZOB) - 1.0)) ** 2        # log-law friction coefficient
S0 = CF * U ** 2 / (G * D0)           # uniform-flow slope, Eq. (1)
SSEA, XS, ELL = 30.0, 50000.0, 15000.0                # salinity transition (psu, m, m)
R0 = 1027.0                           # ROMS reference density (kg/m3)


def scoef_for(eps):
    """Haline contraction coefficient that gives the requested epsilon.

    epsilon = D0 max|d rho/dx| / (2 R0 s0), with max|d rho/dx| = R0 SCOEF S_sea / (2 ell)
    for the tanh salinity profile.
    """
    return eps * 2.0 * R0 * S0 / D0 * 2.0 * ELL / (R0 * SSEA)


def setkey(text, key, val):
    """Replace the value of 'key == ...' (or 'key = ...') in a ROMS input file, keeping comments."""
    pat = re.compile(r"^([ \t]*" + re.escape(key) + r"[ \t]*==?[ \t]*)([^!\n]*?)([ \t]*(![^\n]*)?)$", re.M)
    new, n = pat.subn(lambda m: m.group(1) + val + ("   " + m.group(4) if m.group(4) else ""), text, count=1)
    if n != 1:
        raise KeyError(key)
    return new


def fortran(x):
    return ("%.6e" % x).replace("e", "d")


def roms_src(arg):
    if arg:
        return Path(arg)
    if os.environ.get("ROMS_SRC"):
        return Path(os.environ["ROMS_SRC"])
    if os.environ.get("ROMS_ROOT_DIR"):
        return Path(os.environ["ROMS_ROOT_DIR"]) / "roms"
    raise SystemExit("Set ROMS_SRC (ROMS source directory) or pass --roms-src.")


def make_run(name, eps, AD=0.0, Ah=0.0, days=3.0, out_days=0.5, dt=30.0, tnudg=0.001,
             roms_dir=None, runs_dir=None):
    src = roms_src(roms_dir) / "ROMS" / "External"
    run = Path(runs_dir or HERE / "runs") / name
    run.mkdir(parents=True, exist_ok=True)
    nt = int(round(days * 86400 / dt))
    nout = int(round(out_days * 86400 / dt))
    sc = scoef_for(eps)

    s = (src / "roms_estuary_test.in").read_text()
    for k, v in {
        "TITLE": "River-dominated channel, prescribed salinity, eps=%g" % eps,
        "VARNAME": str(src / "varinfo.yaml"),
        "NCS": "0", "NNS": "1",
        "LBC(isTvar)": "Cla     Clo     Cla     Clo \\",
        "NTIMES": str(nt), "DT": "%.1fd0" % dt, "NDTFAST": "20",
        "NRST": str(nt), "NHIS": str(nout), "NAVG": str(nout), "NDEFHIS": "0", "NDEFAVG": "0",
        "Zob": "%gd0" % ZOB,
        "TNUDG": "%gd0 %gd0" % (tnudg, tnudg),
        "TCOEF": "0.0d0", "SCOEF": fortran(sc),
        "LtracerCLM": "T T", "LnudgeTCLM": "T T",
        "NUSER": "8",
        # ana_m2obc spreads the discharge over four rows, the channel has three wet rows
        "USER": " ".join(fortran(v) for v in (D0, U * 4.0 / 3.0, SSEA, XS, ELL, S0, AD, Ah)),
        "Aout(idBath)": "T", "Aout(idUbms)": "T", "Aout(idVvis)": "T",
        "HISNAME": str(run / "his.nc"), "AVGNAME": str(run / "avg.nc"), "RSTNAME": str(run / "rst.nc"),
    }.items():
        s = setkey(s, k, v)
    # the tracer LBC has a continuation line for salinity: clamp it as well
    s = s.replace("RadNud  Clo     Cla     Clo         ! salinity",
                  "Cla     Clo     Cla     Clo         ! salinity", 1)
    s = re.sub(r"^(\s*SPARNAM\s*=\s*).*$", lambda m: m.group(1) + str(run / "sediment.in"), s, count=1, flags=re.M)
    (run / "roms.in").write_text(s)

    t = (src / "sediment_estuary_test.in").read_text()
    for k, v in {"SAND_SD50": "0.3d0", "SAND_SRHO": "2650.0d0", "SAND_WSED": "40.0d0",
                 "SAND_ERATE": "0.0d0", "SAND_TAU_CE": "0.2d0", "SAND_TAU_CD": "0.2d0",
                 "SAND_POROS": "0.4d0", "SAND_MORPH_FAC": "0.0d0",        # fixed bed
                 "BEDLOAD_COEFF": "1.0d0", "LBC(isTvar)": "Cla     Clo     Cla     Clo"}.items():
        t = setkey(t, k, v)
    (run / "sediment.in").write_text(t)

    prm = dict(eps=eps, scoef=sc, s0=S0, D0=D0, SSEA=SSEA, XS=XS, ELL=ELL, AD=AD, Ah=Ah)
    (run / "params.json").write_text(json.dumps(prm, indent=1))
    return run, prm


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("name")
    p.add_argument("eps", type=float)
    p.add_argument("--AD", type=float, default=0.0)
    p.add_argument("--Ah", type=float, default=0.0)
    p.add_argument("--days", type=float, default=3.0)
    p.add_argument("--out-days", type=float, default=0.5)
    p.add_argument("--roms-src", default=None)
    p.add_argument("--runs-dir", default=None)
    a = p.parse_args()
    run, prm = make_run(a.name, a.eps, a.AD, a.Ah, a.days, a.out_days, roms_dir=a.roms_src, runs_dir=a.runs_dir)
    print(run, "SCOEF=%.4e  s0=%.4e" % (prm["scoef"], prm["s0"]))
