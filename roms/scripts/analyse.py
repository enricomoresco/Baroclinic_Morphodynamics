"""Analysis of the ROMS equilibrium test: the ROMS numbers of Sec. IV B and Appendix B.

Reads roms/runs/*/his.nc (produced by run_equilibrium_test.py) and writes
results/roms_equilibrium.json, which plot_roms.py uses to draw Fig. 3.

Steps
1. gamma on the fixed-bed twins: median of (F - 1)/delta over 20-80 km, uniform bed, eps = 0.1.
2. Each baroclinic run is compared with the barotropic twin: tau_b/tau_b,ref - 1, after removing
   the non-Boussinesq contribution of the surface-slope term (roms_tools.nonboussinesq).
3. The departure is fitted with a constant, a sech^2 (local) and a tanh (cumulative) term.
   The sech^2 coefficient changes sign where the imposed local depth change is the equilibrium one
   (k_opt); on the bed without local change it gives gamma = -c_sech2/eps.
4. Set-up scan: zero crossing of the tanh coefficient.
5. Check of the non-Boussinesq correction on the full first-order bed (Appendix B).
"""
import json

import numpy as np

from roms_tools import (GAMMA_REF, J, KS, RESULTS, RUNS, SETUP_FACTORS, WINDOW, amplitudes, coords,
                        corrected_departure, fit_shapes, last, nonboussinesq, rms_departure,
                        stress_departure, tag, window_mask)


def gamma_fixed_bed(eps=0.1, runs=RUNS):
    """Median of (F - 1)/delta over 20-80 km in the fixed-bed twins on the uniform-flow bed."""
    prm = json.loads((runs / f"eq_{tag(eps)}_flat" / "params.json").read_text())
    xr, xu = coords(runs)

    def load(run):
        z, h = last(run, "zeta", runs)[J], last(run, "bath", runs)[J]
        S = last(run, "salt", runs)[:, J].mean(0)              # depth-averaged salinity
        return dict(ub=last(run, "ubar", runs)[J], bs=last(run, "bustr", runs)[J],
                    D=0.5 * (h[1:] + h[:-1] + z[1:] + z[:-1]), zx=np.diff(z) / np.diff(xr),
                    Sx=np.diff(S) / np.diff(xr))
    bt, bc = load("eq_bt"), load(f"eq_{tag(eps)}_flat")
    F = (bc["bs"] / bc["ub"] ** 2) / (bt["bs"] / bt["ub"] ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        delta = bc["D"] * prm["scoef"] * bc["Sx"] / (2.0 * bc["zx"])  # Eq. (4), rho_x = rho0 beta S_x
    m = (xu > 20e3) & (xu < 80e3) & (np.abs(delta) > 0.02)
    g = (F[m] - 1.0) / delta[m]
    return dict(median=float(np.median(g)), iqr=[float(np.percentile(g, 25)), float(np.percentile(g, 75))])


def local_fit(eps, runs=RUNS):
    """Fits of the k series for one eps: sech^2 and tanh coefficients, rms, k_opt, gamma."""
    xr, xu = coords(runs)
    cs, ct, rms = [], [], []
    for k in KS:
        r = corrected_departure(f"eq_{tag(eps)}_k{k:+.1f}", runs=runs)
        c = fit_shapes(xu, r)
        cs.append(float(c[1])); ct.append(float(c[2])); rms.append(rms_departure(xu, r))
    k_opt = float(np.interp(0.0, cs, KS))                   # sech^2 coefficient increases with k
    gamma = -cs[KS.index(0.0)] / eps                        # local stress deficit, no local depth change
    return dict(k=KS, c_sech2=cs, c_tanh=ct, rms_k=rms, k_opt=k_opt, gamma_local=gamma)


def main(runs=RUNS):
    xr, xu = coords(runs)
    m = window_mask(xu)
    out = dict(gamma_ref=GAMMA_REF, window_km=[w / 1e3 for w in WINDOW], x_km=(xu[m] / 1e3).tolist(),
               gamma_fixed_bed=gamma_fixed_bed(runs=runs), eps={})
    g = out["gamma_fixed_bed"]
    print(f"fixed-bed twins, eps = 0.1: gamma = {g['median']:.3f} (IQR {g['iqr'][0]:.3f}-{g['iqr'][1]:.3f})")

    for eps in (0.1, 0.3):
        A_D, A_h = amplitudes(eps)
        fit = local_fit(eps, runs)
        rms = {nm: rms_departure(xu, corrected_departure(f"eq_{tag(eps)}_{nm}", runs=runs))
               for nm in ("flat", "k+0.0", "k+1.0")}
        dep = {}
        for nm in ("flat", "k+0.0", "k+1.0"):
            r = corrected_departure(f"eq_{tag(eps)}_{nm}", runs=runs)[m]
            dep[nm] = (r - r.mean()).tolist()
        # non-Boussinesq check on the full first-order bed
        run = f"eq_{tag(eps)}_k+1.0"
        pred, finv = nonboussinesq(run, runs=runs)
        sim = stress_departure(run, runs=runs)
        i1, i2 = np.argmin(np.abs(xu - 47.5e3)), np.argmin(np.abs(xu - 92.5e3))
        nb = dict(one_over_f_minus_1=float(finv[i2]), predicted_drop=float(pred[i1] - pred[i2]),
                  simulated_drop=float(sim[i1] - sim[i2]))
        out["eps"][str(eps)] = dict(A_D_m=A_D, A_h_m=A_h, total_setup_m=2 * A_h, rms=rms, departure=dep,
                                    nonboussinesq_check=nb, **fit)
        print(f"\neps = {eps}: theory A_D = {100 * A_D:.2f} cm (local), cumulative rise = {100 * 2 * A_h:.2f} cm")
        print("  rms departure (%): uniform bed {:.3f}, cumulative part only {:.3f}, first-order bed {:.3f}"
              .format(*(100 * rms[k] for k in ("flat", "k+0.0", "k+1.0"))))
        print("  sech^2 coefficient (%) vs k:", " ".join(f"{100 * c:+.3f}" for c in fit["c_sech2"]))
        print(f"  k_opt = {fit['k_opt']:.3f}  ->  gamma = {fit['k_opt'] * GAMMA_REF:.4f};"
              f"  local deficit on the bed without local change -> gamma = {fit['gamma_local']:.4f}")
        print(f"  non-Boussinesq: 1/f - 1 = {100 * nb['one_over_f_minus_1']:.2f}% at 92.5 km; stress drop "
              f"47.5-92.5 km predicted {100 * nb['predicted_drop']:.2f}%, simulated {100 * nb['simulated_drop']:.2f}%")

    fit = local_fit(0.1, runs)
    ct, cs = [], []
    for f in SETUP_FACTORS:
        c = fit_shapes(xu, corrected_departure(f"setup_{f:.1f}", runs=runs))
        cs.append(float(c[1])); ct.append(float(c[2]))
    zero = float(np.interp(0.0, ct, SETUP_FACTORS))           # tanh coefficient increases with F
    out["setup_scan"] = dict(eps=0.1, k_local=fit["k_opt"], factor=SETUP_FACTORS, c_tanh=ct, c_sech2=cs,
                             zero_crossing=zero)
    print(f"\nset-up scan (eps = 0.1, local part at k = {fit['k_opt']:.3f}): tanh coefficient (%) "
          + " ".join(f"{100 * c:+.3f}" for c in ct) + f"  ->  zero at {zero:.3f} times the theory")

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "roms_equilibrium.json").write_text(json.dumps(out, indent=1))
    print("\nwritten: results/roms_equilibrium.json")


if __name__ == "__main__":
    main()
