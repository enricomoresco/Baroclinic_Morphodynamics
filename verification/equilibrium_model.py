"""Weakly baroclinic morphodynamic equilibrium of a river-dominated channel.

Nondimensional problem, Eqs. (9)-(10) of the paper:

    h'(1 + delta) = D^(-10/3) F(delta),     delta = eps D G(x) / h'      (momentum)
    D^(-11/2) F(delta)^(5/2) = 1                                          (sediment equilibrium)

D is the depth, h the free-surface elevation, h' = dh/dx, x the landward coordinate
(scaled with the backwater length), G(x) <= 0 the normalized density gradient and
eps its peak ratio to the barotropic pressure gradient. The friction law, Eq. (7), is

    F(delta) = (1 + delta)^2 / (1 + alpha delta)^2,   gamma = 2 (1 - alpha).

At each x the two equations are algebraic in D and h', so the exact solution needs no
expansion in eps: they are solved point by point and h follows by quadrature, h(0) = 0.
The first-order solution, Eq. (13), is D = 1 + eps D1, h = x + eps h1 with

    D1 = (5/11) gamma G,    h1 = -(1 + 17 gamma / 33) int_0^x G dx.
"""
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import fsolve

KAPPA = 0.41


def closure(Lambda, shape):
    """Profile coefficient alpha, friction-law coefficient gamma and c_f, Eqs. (6)-(8).

    Lambda = ln(D/z0); shape = 'linear' or 'parabolic' eddy viscosity.
    """
    if shape == "linear":
        alpha, sqcf = (Lambda - 11 / 6) / (Lambda - 1.5), KAPPA / (Lambda - 1.5)
    elif shape == "parabolic":
        alpha, sqcf = (Lambda - 1.5) / (Lambda - 1.0), KAPPA / (Lambda - 1.0)
    else:
        raise ValueError(shape)
    return dict(alpha=alpha, gamma=2 * (1 - alpha), cf=sqcf ** 2)


def friction_factor(delta, alpha):
    """F(delta) of Eq. (7): bed stress over rho0 c_f ubar^2."""
    return (1 + delta) ** 2 / (1 + alpha * delta) ** 2


def first_order(x, G, gamma):
    """First-order corrections D1, h1 and eta1 = h1 - D1, Eq. (13)."""
    R = cumulative_trapezoid(G, x, initial=0.0)
    D1 = 5 / 11 * gamma * G
    h1 = -(1 + 17 / 33 * gamma) * R
    return D1, h1, h1 - D1


def solve_full(x, G, eps, alpha):
    """Exact solution of Eqs. (9)-(10) for given eps: returns D(x) and h(x)."""
    D = np.empty_like(x)
    s = np.empty_like(x)                       # free-surface slope h'
    guess = np.array([1.0, 1.0])
    for i, g in enumerate(G):
        def residual(v):
            Di, si = v
            d = eps * Di * g / si
            F = friction_factor(d, alpha)
            return [si * (1 + d) - Di ** (-10 / 3) * F, Di ** (-5.5) * F ** 2.5 - 1]
        with warnings.catch_warnings():        # fsolve warns when it cannot improve on a
            warnings.simplefilter("ignore")     # converged root; the residual is checked below
            sol = fsolve(residual, guess, xtol=1e-14)
        if np.max(np.abs(residual(sol))) > 1e-11:
            raise RuntimeError(f"no convergence at x = {x[i]:.3f}, eps = {eps}")
        D[i], s[i] = sol
        guess = sol
    return D, cumulative_trapezoid(s, x, initial=0.0)
