"""One-dimensional vertical (1DV) k-epsilon model of a river flow with a weak horizontal
density gradient (Sec. IV A and Appendix A of the paper).

Steady open-channel flow of depth D driven by a surface slope S and by a prescribed,
depth-uniform horizontal density gradient rho_x. Here x points in the direction of the
river flow and z upward from the bed, so rho_x > 0 (denser water seaward) opposes the flow,
more strongly near the bed. The model integrates in time to steady state

    du/dt  = g S - (g/rho0) rho_x (D - z) + d/dz(nu_t du/dz)                     (A1)
    dS'/dt = -(u - ubar) dS/dx + d/dz(K_s dS'/dz)                                (A2)

S' is the departure of salinity from its depth average: the first term of (A2) is the
straining of the horizontal gradient by the sheared flow, the second the mixing.
Turbulence follows the standard k-epsilon model (Rodi 1993) with shear production
P = nu_t (du/dz)^2 and buoyancy production B = -K_s N^2, N^2 = -g beta dS'/dz,
nu_t = c_mu k^2/epsilon, K_s = nu_t/Pr_t, and the buoyancy coefficient c3 of
Burchard & Baumert (1995): c3 = 1 for B > 0, c3 = -0.4 for B < 0.
In the 'neutral' runs B = 0 and the density enters only through the pressure gradient.

Boundary conditions: logarithmic wall function in the first cell (roughness z0);
zero stress and zero flux of k at the surface, where epsilon is set by the mixing length
kappa (dz/2 + z0s). Numerics: uniform cells, implicit diffusion, Patankar treatment of
the sink terms; the boundary values of k and epsilon are imposed inside the implicit
solve, so that the steady state does not depend on the time step.

The pressure-gradient ratio of the paper, Eq. (4), is delta = -D rho_x / (2 rho0 S) < 0.
The friction factor is F = (tau_b/ubar^2) / c_f0, with c_f0 from the barotropic run.
"""
import numpy as np
from scipy.linalg import solve_banded

G, RHO0, KAPPA, BETA = 9.81, 1000.0, 0.41, 7.7e-4
CMU, C1, C2, SIGK, SIGE = 0.09, 1.44, 1.92, 1.0, 1.3
C3_STABLE = -0.4          # c3 for stable stratification (B < 0)
PR_T = 1.0                # turbulent Prandtl (Schmidt) number


def _tridiag(a_lo, a_d, a_up, rhs):
    ab = np.zeros((3, len(rhs)))
    ab[0, 1:] = a_up[:-1]
    ab[1] = a_d
    ab[2, :-1] = a_lo[1:]
    return solve_banded((1, 1), ab, rhs)


def _diffusion_matrix(nu_f, dz, dt):
    """Implicit diffusion with zero flux through both ends (boundary fluxes are added separately)."""
    n = len(nu_f) - 1
    lo, up, d = np.zeros(n), np.zeros(n), np.ones(n)
    c = dt / dz ** 2
    lo[1:] = -c * nu_f[1:-1]
    up[:-1] = -c * nu_f[1:-1]
    d += c * (nu_f[:-1] + nu_f[1:])
    d[0] -= c * nu_f[0]
    d[-1] -= c * nu_f[-1]
    return lo, d, up


def _dirichlet(lo, d, up, rhs, i, value):
    """Replace row i of the tridiagonal system by the condition x_i = value."""
    lo, d, up, rhs = lo.copy(), d.copy(), up.copy(), rhs.copy()
    lo[i], d[i], up[i], rhs[i] = 0.0, 1.0, 0.0, value
    return lo, d, up, rhs


def run(S=2.5e-5, delta=0.0, D=10.0, z0=1e-3, N=100, buoyancy=True, c3=C3_STABLE,
        prt=PR_T, t_end=6e4, dt=5.0, z0s=0.1):
    """Integrate to t_end (s) and return the final state.

    S: surface slope; delta: pressure-gradient ratio (< 0 opposes the flow); D: depth (m);
    z0: bed roughness length (m); N: number of cells; buoyancy: False for the neutral runs;
    c3, prt: buoyancy coefficient and turbulent Prandtl number; z0s: surface mixing length (m).
    """
    dz = D / N
    z = (np.arange(N) + 0.5) * dz
    rho_x = -2 * RHO0 * S * delta / D              # delta < 0  ->  rho_x > 0
    Sx = rho_x / (RHO0 * BETA)                     # salinity gradient (psu/m)
    force = G * S - G / RHO0 * rho_x * (D - z)

    us0 = np.sqrt(G * D * S)
    u = us0 / KAPPA * np.log(np.maximum(z, 2 * z0) / z0)
    k = np.full(N, us0 ** 2 / np.sqrt(CMU)) * (1 - z / D) + 1e-6
    eps = CMU ** 0.75 * k ** 1.5 / (KAPPA * (z + z0) * np.sqrt(1 - z / D + 0.05))
    s = np.zeros(N)
    for _ in range(int(t_end / dt)):
        nu = np.maximum(CMU * k ** 2 / eps, 1e-6)
        nu_f = np.zeros(N + 1)
        nu_f[1:-1] = 0.5 * (nu[1:] + nu[:-1])
        # momentum, bottom stress implicit
        ustar = KAPPA * abs(u[0]) / np.log(z[0] / z0)
        lo, d, up = _diffusion_matrix(nu_f, dz, dt)
        d = d.copy()
        d[0] += dt / dz * ustar * KAPPA / np.log(z[0] / z0)
        u = _tridiag(lo, d, up, u + dt * force)
        ubar = u.mean()
        # salinity departure: straining and mixing
        lo, d, up = _diffusion_matrix(nu_f / prt, dz, dt)
        s = _tridiag(lo, d, up, s - dt * (u - ubar) * Sx)
        s -= s.mean()
        # shear and buoyancy frequency at cell centres
        du = np.zeros(N + 1)
        du[1:-1] = np.diff(u) / dz
        ds = np.zeros(N + 1)
        ds[1:-1] = np.diff(s) / dz
        M2 = 0.5 * (du[1:] ** 2 + du[:-1] ** 2)
        N2 = -G * BETA * 0.5 * (ds[1:] + ds[:-1]) if buoyancy else np.zeros(N)
        P = nu * M2
        B = -(nu / prt) * N2
        # wall function in the first cell, from the updated velocity
        ustar = KAPPA * abs(u[0]) / np.log(z[0] / z0)
        # k: wall value imposed in the first cell (Dirichlet row), no flux at the surface
        lo, d, up = _diffusion_matrix(nu_f / SIGK, dz, dt)
        d = d + dt * eps / k + dt * np.maximum(-B, 0) / k
        rhs = k + dt * (P + np.maximum(B, 0))
        lo, d, up, rhs = _dirichlet(lo, d, up, rhs, 0, ustar ** 2 / np.sqrt(CMU))
        k = np.maximum(_tridiag(lo, d, up, rhs), 1e-9)
        # epsilon: wall value in the first cell, mixing-length value in the surface cell
        c3v = np.where(B < 0, c3, 1.0)
        src = eps / k * (C1 * P + c3v * B)
        lo, d, up = _diffusion_matrix(nu_f / SIGE, dz, dt)
        d = d + dt * C2 * eps / k + dt * np.maximum(-src, 0) / eps
        rhs = eps + dt * np.maximum(src, 0)
        lo, d, up, rhs = _dirichlet(lo, d, up, rhs, 0, ustar ** 3 / (KAPPA * z[0]))
        lo, d, up, rhs = _dirichlet(lo, d, up, rhs, N - 1,
                                    CMU ** 0.75 * k[-1] ** 1.5 / (KAPPA * (dz / 2 + z0s)))
        eps = np.maximum(_tridiag(lo, d, up, rhs), 1e-12)
    ustar = KAPPA * abs(u[0]) / np.log(z[0] / z0)
    taub = ustar ** 2
    return dict(z=z, u=u, ubar=u.mean(), taub=taub, taub_exact=G * D * S * (1 + delta),
                s=s, dS=s[0] - s[-1], nu=CMU * k ** 2 / eps, D=D, S=S, delta=delta)


def friction_coefficient(r):
    """tau_b / ubar^2 (kinematic)."""
    return r["taub"] / r["ubar"] ** 2


if __name__ == "__main__":
    r = run()
    cf = friction_coefficient(r)
    print(f"barotropic run: ubar = {r['ubar']:.3f} m/s, c_f = {cf:.5f}, "
          f"tau_b/(g D S) = {r['taub'] / r['taub_exact']:.4f}, sqrt(c_f)/kappa = {np.sqrt(cf) / KAPPA:.4f}")
