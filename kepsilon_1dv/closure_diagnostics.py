"""Why the coefficient gamma of the friction law depends on the turbulence closure (Secs. II B and IV A,
Table II and Appendix A of the paper).

Without stratification, gamma in F = 1 + gamma*delta measures how the depth-averaged velocity at a
given bed stress responds to the baroclinic force, which increases the stress in the upper water
column relative to the bed stress. The script compares three ways in which the turbulence can
respond to this increase:

1. eddy viscosity of parabolic shape scaled by the bed friction velocity (Sec. II B): the viscosity
   ignores the increase, the shear grows in proportion to the stress and gamma = sqrt(c_f)/kappa;
2. mixing length l = kappa z (1 - z/D)^(1/2), which gives the same barotropic profile but lets the
   eddy viscosity l^2 |du/dz| follow the local stress: the shear grows as the square root of the
   stress, so gamma is halved, sqrt(c_f)/(2 kappa);
3. the 1DV k-epsilon model, in which the length scale adapts as well. The run at delta = -0.2 is
   compared with the barotropic run rescaled to the same bed stress (the neutral steady equations
   are invariant under u -> a u, k -> a^2 k, epsilon -> a^3 epsilon), and with the change that a
   length scale fixed at its barotropic profile would give, obtained by scaling the barotropic shear
   between cells with the square root of the stress ratio.
4. Sensitivity of the k-epsilon result to the constants of the epsilon equation: sigma_e = 1.1, 1.3
   and 1.5, with c1 set by kappa^2 = (c2 - c1) sigma_e c_mu^(1/2) so that the log layer keeps
   kappa = 0.41.

Cases 1 and 2 are computed by quadrature for the paper's depth and roughness (D = 10 m, z0 = 1 mm).
Output: results/kepsilon_closures.json.  Run time: under a minute.
Expected: gamma = 0.122 (prescribed viscosity) and 0.061 (mixing length) for delta -> 0; k-epsilon at
delta = -0.2: F = 0.9996; below a quarter of the depth the length scale grows more than the local
friction velocity, and the velocity decreases there (minimum -0.27 % of ubar) while it increases near
the surface (+0.62 %); depth-averaged velocity +0.02 % against +0.81 % with the barotropic length
scale; |gamma_eff| < 0.005 for the three values of sigma_e.
"""
import json
from contextlib import contextmanager
from pathlib import Path

import numpy as np

import kepsilon_1dv as ke

ROOT = Path(__file__).resolve().parents[1]
D, Z0 = 10.0, 1e-3
DELTA = -0.2
T_END = 1.2e5


def ubar_closure(delta, closure, n=40001):
    """kappa*ubar/u_* at a given bed stress, for 'prescribed' (parabolic viscosity scaled by the bed
    friction velocity) or 'mixing' (mixing length kappa z (1-z/D)^1/2); no slip at z0."""
    zeta = np.geomspace(Z0 / D, 1.0, n)
    one = np.maximum(1.0 - zeta, 1e-12)
    phi = (one + delta * one ** 2) / (1.0 + delta)          # tau / tau_b
    ratio = phi / one                                        # tau / (tau_b (1 - zeta))
    dudz = (ratio if closure == "prescribed" else np.sqrt(ratio)) / zeta
    u = np.concatenate([[0.0], np.cumsum(0.5 * (dudz[1:] + dudz[:-1]) * np.diff(zeta))])
    return np.trapezoid(u, zeta) / (zeta[-1] - zeta[0])


def gamma_closure(delta, closure):
    F = (ubar_closure(0.0, closure) / ubar_closure(delta, closure)) ** 2
    return (F - 1.0) / delta


@contextmanager
def constants(**values):
    """Temporarily change module constants of kepsilon_1dv (e.g. SIGE, C1)."""
    old = {name: getattr(ke, name) for name in values}
    for name, value in values.items():
        setattr(ke, name, value)
    try:
        yield
    finally:
        for name, value in old.items():
            setattr(ke, name, value)


def neutral_pair():
    bt = ke.run(delta=0.0, buoyancy=False, t_end=T_END)
    bc = ke.run(delta=DELTA, buoyancy=False, t_end=T_END)
    F = ke.friction_coefficient(bc) / ke.friction_coefficient(bt)
    return bt, bc, F


def compensation(bt, bc):
    """Profiles of the baroclinic run against the barotropic run rescaled to the same bed stress."""
    z = bc["z"]
    zeta = z / D
    a = np.sqrt(bc["taub"] / bt["taub"])                     # velocity scale factor
    u_bt, u_bc = a * bt["u"], bc["u"]
    ubar = u_bt.mean()
    stress = (1.0 + DELTA * (1.0 - zeta)) / (1.0 + DELTA)    # tau_bc / tau_bt at equal bed stress (exact)
    length = lambda r: ke.CMU ** 0.75 * r["k"] ** 1.5 / r["eps"]
    l_ratio = length(bc) / length(bt)
    nu_ratio = bc["nu"] / (a * bt["nu"])
    du = 100 * (u_bc - u_bt) / ubar                           # % of the depth-averaged velocity
    # length scale fixed at its barotropic profile: the shear between cells, as computed by the model,
    # scaled with the square root of the stress ratio at the cell faces (exact when the ratio is 1)
    z_face = 0.5 * (z[1:] + z[:-1])
    stress_face = (1.0 + DELTA * (1.0 - z_face / D)) / (1.0 + DELTA)
    u_fixed = u_bt[0] + np.concatenate([[0.0], np.cumsum(np.diff(u_bt) * np.sqrt(stress_face))])
    excess = l_ratio - np.sqrt(stress)                        # > 0 where l grows more than the local u_*
    i_cross = int(np.argmax(excess[1:] < 0)) + 1              # first cell above the wall cell where it turns negative
    i_min = int(np.argmin(du))
    i_zero = int(np.argmax((du > 0) & (np.arange(len(du)) > i_min)))
    summary = dict(
        stress_increase_top_percent=100 * (stress[-1] - 1),
        length_scale_exceeds_local_ustar_below_zeta=float(zeta[i_cross]),
        length_scale_increase_below_that_percent=[100 * float(l_ratio[1:i_cross].min() - 1),
                                                  100 * float(l_ratio[1:i_cross].max() - 1)],
        local_ustar_increase_below_that_percent=[100 * float(np.sqrt(stress[1:i_cross]).min() - 1),
                                                 100 * float(np.sqrt(stress[1:i_cross]).max() - 1)],
        length_scale_increase_max_percent=100 * float(l_ratio.max() - 1),
        eddy_viscosity_increase_max_percent=100 * float(nu_ratio.max() - 1),
        velocity_change_min_percent=float(du[i_min]), velocity_change_min_at_zeta=float(zeta[i_min]),
        velocity_change_sign_change_at_zeta=float(zeta[i_zero]),
        velocity_change_top_percent=float(du[-1]),
        depth_mean_change_percent=100 * (u_bc.mean() - ubar) / ubar,
        depth_mean_change_fixed_length_scale_percent=100 * (u_fixed.mean() - ubar) / ubar,
    )
    du_fixed = 100 * (u_fixed - u_bt) / ubar
    profiles = dict(zeta=zeta.tolist(), stress_ratio=stress.tolist(), length_scale_ratio=l_ratio.tolist(),
                    eddy_viscosity_ratio=nu_ratio.tolist(), velocity_change_percent=du.tolist(),
                    velocity_change_fixed_length_scale_percent=du_fixed.tolist())
    return summary, profiles


def main():
    lam = np.log(D / Z0)
    out = dict(depth=D, z0=Z0, analytic=dict(prescribed=1 / (lam - 1), mixing_length=1 / (2 * (lam - 1))))
    print(f"Lambda = ln(D/z0) = {lam:.3f}")
    print("1-2. friction-law coefficient by quadrature (analytic limit for delta -> 0 in brackets)")
    for closure, key in (("prescribed", "prescribed"), ("mixing", "mixing_length")):
        g = {f"{d:g}": gamma_closure(d, closure) for d in (-0.005, DELTA)}
        out[key] = g
        print(f"   {key:13s} gamma_eff = {g['-0.005']:.4f} at delta = -0.005, {g[f'{DELTA:g}']:.4f} at "
              f"delta = {DELTA}   [{out['analytic'][key]:.4f}]")

    print(f"3. k-epsilon, neutral, delta = {DELTA}, equal bed stress")
    bt, bc, F = neutral_pair()
    summary, profiles = compensation(bt, bc)
    out["k_epsilon"] = dict(delta=DELTA, F=F, gamma_eff=(F - 1) / DELTA, **summary, profiles=profiles)
    s = summary
    print(f"   F = {F:.5f}, gamma_eff = {(F - 1) / DELTA:.4f}")
    print(f"   stress near the surface +{s['stress_increase_top_percent']:.1f} %")
    print(f"   below zeta = {s['length_scale_exceeds_local_ustar_below_zeta']:.3f} the length scale grows by "
          f"{s['length_scale_increase_below_that_percent'][0]:.2f}-{s['length_scale_increase_below_that_percent'][1]:.2f} %,"
          f" the local friction velocity by {s['local_ustar_increase_below_that_percent'][0]:.2f}-"
          f"{s['local_ustar_increase_below_that_percent'][1]:.2f} %")
    print(f"   largest increase: length scale {s['length_scale_increase_max_percent']:.2f} %, eddy viscosity "
          f"{s['eddy_viscosity_increase_max_percent']:.2f} %")
    print(f"   velocity change: minimum {s['velocity_change_min_percent']:.3f} % of ubar at zeta = "
          f"{s['velocity_change_min_at_zeta']:.3f}, positive above zeta = {s['velocity_change_sign_change_at_zeta']:.3f},"
          f" {s['velocity_change_top_percent']:.3f} % at the surface")
    print(f"   depth-averaged velocity {s['depth_mean_change_percent']:+.3f} %, with the barotropic length scale "
          f"{s['depth_mean_change_fixed_length_scale_percent']:+.3f} %")

    print("4. k-epsilon, neutral: constants of the epsilon equation")
    rows = []
    for sige in (1.1, 1.3, 1.5):
        c1 = ke.C2 - ke.KAPPA ** 2 / (sige * np.sqrt(ke.CMU))
        with constants(SIGE=sige, C1=c1):
            _, _, Fs = neutral_pair()
        rows.append(dict(sigma_e=sige, c1=c1, F=Fs, gamma_eff=(Fs - 1) / DELTA))
        print(f"   sigma_e = {sige:.1f}, c1 = {c1:.3f}: F = {Fs:.5f}, gamma_eff = {(Fs - 1) / DELTA:+.4f}")
    out["epsilon_constants"] = rows

    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "kepsilon_closures.json").write_text(json.dumps(out, indent=1))
    print("written: results/kepsilon_closures.json")


if __name__ == "__main__":
    main()
