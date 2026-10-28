# Density-driven morphodynamic equilibrium of river-dominated estuaries

Code accompanying the paper

> E. Moresco, L. Durante and N. Tambroni, *Effect of a weak longitudinal density gradient on the bed
> stress and equilibrium bed profile of river-dominated estuaries*, Physics of Fluids (in preparation).

The repository reproduces the numerical results and the figures of the paper:

| Folder | Paper | What it does | Run time |
|---|---|---|---|
| [`verification/`](verification) | Sec. III B, Fig. 1 | exact solution of Eqs. (9)–(10) compared with the first-order solution, Eq. (13) | seconds |
| [`kepsilon_1dv/`](kepsilon_1dv) | Secs. II B, IV A, IV C, V B, Table II, Appendix A, Fig. 2 | one-dimensional vertical k–ε model that tests the friction law, Eq. (7), and its diagnostics | ~25 min |
| [`roms/`](roms) | Sec. IV B, Appendix B, Fig. 3 | ROMS configuration and analysis of the equilibrium test | ~20 min (2 cores) |
| [`families/`](families) | Sec. V B, Fig. 4 | first-order response to prescribed salinity distributions, Eq. (14) | seconds |
| [`olabarrieta/`](olabarrieta) | Sec. VI B, Fig. 5 | analysis of the published simulations of Olabarrieta et al. (2018) | seconds |

The outputs of all scripts are stored in [`results/`](results) and [`figures/`](figures), so the figures
can be redrawn without running the models again.

## Quick start

Python 3.9 or later:

```bash
pip install -r requirements.txt
python verification/run_verification.py        # Fig. 1, results/verification.json
python kepsilon_1dv/run_sweep.py               # results/kepsilon_sweep.json
python kepsilon_1dv/run_sensitivity.py         # results/kepsilon_sensitivity.json
python kepsilon_1dv/run_nudging.py             # Appendix A, results/kepsilon_nudging.json
python kepsilon_1dv/run_threshold.py           # Sec. IV A and Appendix A, results/kepsilon_threshold.json
python kepsilon_1dv/exchange_and_bbl.py        # Secs. IV A and V B, results/kepsilon_exchange_bbl.json
python kepsilon_1dv/local_equilibrium_1dv.py   # Sec. IV C, results/local_equilibrium_1dv.json
python kepsilon_1dv/closure_diagnostics.py     # Secs. II B and IV A, results/kepsilon_closures.json
python kepsilon_1dv/plot_friction_factor.py    # Fig. 2, from results/kepsilon_sweep.json and kepsilon_closures.json
python roms/scripts/plot_roms.py               # Fig. 3 from results/roms_equilibrium.json
python families/fig4_response.py               # Fig. 4 and the numbers of Sec. V B
python olabarrieta/make_fig5.py results figures/fig5_olabarrieta.pdf   # Fig. 5 from results/olabarrieta_consistent.npz
```

The Python scripts run on Linux, macOS and Windows. Building and running ROMS needs Linux or macOS
(on Windows, WSL), see [below](#3-roms-equilibrium-test-roms).

## 1. Verification of the first-order solution (`verification/`)

In nondimensional form, with the depth $D$ and the free-surface elevation $h$ scaled with the
uniform-flow depth and the landward coordinate $x$ with the backwater length, the equilibrium problem
is (Eqs. 9–10 of the paper)

$$
\frac{\mathrm{d}h}{\mathrm{d}x}\,(1+\delta)=D^{-10/3}F(\delta),\qquad
D^{-11/2}F(\delta)^{5/2}=1,\qquad
\delta=\epsilon\,\frac{D\,G(x)}{\mathrm{d}h/\mathrm{d}x}.
$$

The first equation is the depth-integrated momentum balance with a Strickler friction coefficient; the
second is the equilibrium condition for Engelund–Hansen transport, a sediment flux uniform along the
channel. $\delta$ is the ratio of the baroclinic to the barotropic pressure gradient, $G\le 0$ the
normalized density gradient, $\epsilon$ its peak value and $F(\delta)=(1+\delta)^2/(1+\alpha\delta)^2$
the friction law of Eq. (7), with $\alpha$ given by Eq. (6) for a linear or parabolic eddy viscosity.

At each $x$ the two equations are algebraic in $D$ and $\mathrm{d}h/\mathrm{d}x$. They are solved
exactly, point by point, and $h$ follows by quadrature; no expansion in $\epsilon$ is made. The result
is compared with the first-order solution $D=1+\epsilon D_1$, $h=x+\epsilon h_1$ of Eq. (13), for a
hyperbolic-tangent density transition, $\ln(D/z_0)=9.21$ ($D=10$ m, $z_0=1$ mm) and
$0.005\le\epsilon\le 0.5$.

- `equilibrium_model.py`: friction law, first-order solution, exact solver;
- `run_verification.py`: error table, `results/verification.json`, `figures/fig1_verification.*`.

Expected output: the error decreases as $\epsilon^2$ (observed order 2.00); the relative error of the
depth correction is 1.2, 2.4 and 5.8 % at $\epsilon=0.1$, 0.2 and 0.5 with the linear eddy viscosity
(1.7, 3.4 and 7.9 % with the parabolic one).

![Fig. 1](figures/fig1_verification.png)

## 2. One-dimensional k–ε model (`kepsilon_1dv/`)

A vertical column 10 m deep with bed roughness $z_0=1$ mm, driven by a surface slope and by a
prescribed, depth-uniform horizontal density gradient, integrated to steady state (Appendix A of the
paper). Turbulence follows the standard k–ε model with the buoyancy coefficient of Burchard and
Baumert (1995). In the *neutral* runs the density enters only through the pressure gradient; in the
*stratified* runs the salinity stratification produced by straining is also computed and damps
turbulence. The friction factor $F$ is the ratio of $\tau_b/\bar u^2$ to its value in the barotropic
run with the same slope, and $\gamma_{\rm eff}=(F-1)/\delta$. The boundary values of $k$ and
$\varepsilon$ are imposed inside the implicit solve, so the steady state does not depend on the time
step (checked for $\Delta t=5$, 2.5 and 1.25 s).

- `kepsilon_1dv.py`: the model (equations and closure constants in the docstring);
- `run_sweep.py`: $F(\delta)$ for the neutral case, the stratified case ($c_3=-0.4$, $Pr_t=1$) and two
  stratified variants, plus the threshold beyond which the stratified runs have no steady state;
- `run_sensitivity.py`: grid (50, 100, 200 cells), $c_3$, $Pr_t$ and surface mixing length; amplitude of
  the exchange flow at $\delta=-0.1$;
- `plot_friction_factor.py`: Fig. 2, the friction factor of `run_sweep.py` in panel (a) and the k–ε
  response at $\delta=-0.2$ of `closure_diagnostics.py` in panels (b) and (c);
- `run_nudging.py`: relaxation of the salinity towards the depth average, $-S'/\tau_n$ (argument `tau_n`
  of `run()`), which mimics the nudging of the ROMS runs (Appendix A);
- `run_threshold.py`: end of the stratified branch, with integrations of $5\times10^5$ and $6\times10^5$ s
  (Sec. IV A and Appendix A);
- `exchange_and_bbl.py`: salt carried landward by the exchange flow, as an equivalent dispersion
  coefficient (Sec. V B), and height of the bottom boundary layer compared with the scaling of
  Stacey and Ralston (2005) (Sec. IV A);
- `local_equilibrium_1dv.py`: local equilibrium depth with the stratified $F(\delta)$ of `run_sweep.py`
  in place of Eq. (7) (Sec. IV C);
- `closure_diagnostics.py`: why, without stratification, $\gamma$ depends on the closure (Secs. II B and
  IV A, Table II, Appendix A). It compares an eddy viscosity of parabolic shape scaled by the bed friction
  velocity, a mixing length that lets the eddy viscosity follow the local stress, and the k–ε model, in
  which the length scale adapts as well, and tests the k–ε result against the constants of the ε equation.

Expected output of `run_nudging.py`, `run_threshold.py`, `exchange_and_bbl.py` and
`local_equilibrium_1dv.py`: $\gamma_{\rm eff}=0.235$ at $\delta=-0.02$ without relaxation
and 0.218, 0.162, 0.059 for $\tau_n=3600$, 600 and 86 s; the stratified runs are steady at
$\delta=-0.328$ ($F=0.50$) and $-0.333$ ($F=0.21$, bed–surface salinity difference 46 psu), while at
$\delta=-0.334$ $F$ keeps decreasing; at $\delta=-0.1$ the exchange flow carries salt
like a dispersion coefficient of 4.0 (neutral) and 5.0 m²/s (stratified), against $2.3\times10^4$ m²/s
needed to balance the river at a mean salinity of 15 psu; the height at which the flux Richardson number
reaches 0.2 is 0.92, 0.70 and 0.49 times the depth at $\delta=-0.1$, $-0.2$ and $-0.3$ (scaling: 0.95,
0.63, 0.48); the local shoaling with the stratified friction law is 3.89 % at $\delta_0=-0.3$ and
8.48 % at $\delta_0=-0.5$.

Expected output of the first scripts: the barotropic reference has $\bar u=1.03$ m/s and $c_f=2.31\times10^{-3}$.
Without stratification $|\gamma|<0.01$ (for the reason, see `closure_diagnostics.py` below). With
stratification $\gamma=0.23$–$0.25$ as $\delta\to0$ and $0.25$–$0.27$ at $\delta=-0.05$
($c_3=-0.4$ to $-0.6$), within 1 % on grids of 50–200 cells. The bisection of `run_sweep.py`, which
requires $F>0.5$ after $1.2\times10^5$ s, stops at $\delta=-0.325$; with the longer integrations of
`run_threshold.py` the steady branch extends to $\delta=-0.333$, as in the paper. At $\delta=-0.1$ the departure of the velocity profile from the barotropic one
(exchange flow) is 0.27 % of $\bar u$ without and 1.5–1.7 % with stratification.

Expected output of `closure_diagnostics.py`: $\gamma=0.122$ for the parabolic eddy viscosity and 0.061
for the mixing length ($\delta\to0$). In the k–ε model at $\delta=-0.2$ and equal bed stress, the length
scale grows more than the local friction velocity below a quarter of the depth, so the velocity decreases
there (by up to 0.27 % of $\bar u$) and increases near the surface (by 0.62 %); the depth-averaged
velocity changes by 0.02 %, against 0.81 % with the barotropic length scale. $|\gamma_{\rm eff}|<0.005$ for
$\sigma_\varepsilon=1.1$, 1.3 and 1.5, with $c_1$ adjusted to keep $\kappa=0.41$.

![Fig. 2](figures/fig2_friction_factor.png)

## 3. ROMS equilibrium test (`roms/`)

A straight channel 100 km long and 300 m wide, 20 terrain-following levels, k–ε closure of the generic
length-scale formulation with Kantha–Clayson stability functions, logarithmic bed drag with
$z_0=1$ mm, no tide, river velocity 1 m/s. Salinity is nudged, on a time scale of 86 s, to a
hyperbolic-tangent profile with 30 psu at the sea, centred at 50 km, half-width 15 km; the residual
bottom–surface salinity difference is below 0.02 psu. The haline contraction coefficient is set to give
$\epsilon=0.1$ and $0.3$, and to zero in the barotropic twin.

With bed-load (Meyer-Peter and Müller) transport, a bed profile is a morphodynamic equilibrium if and
only if the bed stress it produces is uniform. The test therefore imposes the bed predicted by the theory
on a fixed bed and measures the departure of the bed stress from that of the barotropic twin; the local
($\mathrm{sech}^2$) and cumulative ($\tanh$) parts of the correction are separated by a least-squares fit.
Before the fit, the contribution of a non-Boussinesq term of the ROMS pressure gradient, which the theory
does not contain, is computed from the simulated fields and removed (Appendix B).

Contents:

- `estuary_test.h`: CPP options; the case is built on the ROMS `ESTUARY_TEST` application;
- `Functionals/`: analytical functions that replace the ROMS defaults (grid and bed, initial state,
  boundary conditions, salinity climatology, sediment bed). The changes with respect to ROMS are
  listed in [`roms/CHANGES.md`](roms/CHANGES.md);
- `build_roms.sh`: the ROMS build script, set for a serial gfortran build;
- `scripts/make_run.py`: writes the input files of one run from the ROMS templates;
- `scripts/run_equilibrium_test.py`: the 20 runs of the test;
- `scripts/analyse.py`: all ROMS numbers of the paper, `results/roms_equilibrium.json`;
- `scripts/plot_roms.py`: Fig. 3.

### Build and run

Requirements: gfortran, the NetCDF-Fortran library, make, perl and the Python packages above.
The code is ROMS at revision `57aecf5` of <https://github.com/myroms/roms>; the configuration uses the
unmodified ROMS source.

```bash
git clone https://github.com/myroms/roms.git ~/src/roms
git -C ~/src/roms checkout 57aecf5
export ROMS_ROOT_DIR=~/src          # directory that contains roms/

cd roms
./build_roms.sh -j 2                # creates roms/romsS
cd scripts
python run_equilibrium_test.py --jobs 2
python analyse.py
python plot_roms.py
```

Each run simulates 3 days and takes one to two minutes on one core. The NetCDF output (about 20 MB
per run) is written to `roms/runs/` and is not part of the repository.

Expected output (Sec. IV B): the rms departure of the bed stress at $\epsilon=0.1$ is 0.73 % on the
uniform bed, 0.26 % with the cumulative part only and 0.023 % on the full first-order bed
(2.2, 0.76 and 0.10 % at $\epsilon=0.3$); the local mismatch changes sign at 0.94 and 0.91 times the
amplitude predicted with $\gamma=0.08$, which corresponds to $\gamma=0.074$; with the local part at
this amplitude, the cumulative mismatch vanishes at 1.01 times the theoretical set-up. The 15 runs of
the first stage are bit-for-bit identical to those used for the paper.

![Fig. 3](figures/fig3_roms_equilibrium.png)

## 4. Response to prescribed salinity distributions (`families/`)

First-order response, Eq. (14), of the reference river of Sec. V B (depth 10 m, velocity 1 m/s,
roughness 1 mm, sea salinity 30 psu, salt intrusion 40 km) to a Van der Burgh profile with $K=0.6$, a
linear profile and a hyperbolic-tangent front.

- `fig4_response.py`: Fig. 4 (`figures/fig4_families.*`) and the numbers quoted in Sec. V B.

Expected output: maximum $|\delta_0|$ of 0.17, 0.10 and 0.25; maximum shoaling per unit $\gamma$ of 76,
46 and 114 cm; set-up of 12.0 cm for $\gamma=0.074$ and 13.3 cm for $\gamma=0.3$ (11.9 and 13.2 cm for
the front), against a hydrostatic head of 11.6 cm.

![Fig. 4](figures/fig4_families.png)

## 5. Simulations of Olabarrieta et al. (2018) (`olabarrieta/`)

Depth reduction and salinity gradient in the baroclinic and barotropic runs of Fig. 8 of
Olabarrieta et al. (2018, *J. Geophys. Res. Earth Surf.* 123, 2901–2924), from their published data set
(Figshare, [doi:10.6084/m9.figshare.5975164](https://doi.org/10.6084/m9.figshare.5975164)). The data set
is not part of this repository: download and unzip it, then

```bash
python olabarrieta/olabarrieta_consistent.py <unzipped DATASET folder> results   # results/olabarrieta_consistent.npz
python olabarrieta/make_fig5.py results figures/fig5_olabarrieta.pdf               # Fig. 5
```

The second command alone redraws Fig. 5 from the stored `results/olabarrieta_consistent.npz`.
Expected output: depth change negative over 96 % of the estuary, maximum $-9.4$ m (59 % of the
barotropic depth) within 0.5 km of the maximum salinity gradient; within 36.5 km of the head (24 % of
the estuary) the salinity gradient is below 2 % of its maximum and the depth is reduced by 2.2 m on
average (34 %), 10 % of the integrated depth change; a local response fitted to the whole profile
explains 77 % of its squared magnitude.

![Fig. 5](figures/fig5_olabarrieta.png)

## Acknowledgements

Parts of the code and of this documentation were prepared with the help of an AI-based coding
assistant. All scripts were reviewed by the authors, and the results reported here were checked by
rerunning the full workflow.

## Citation

See [`CITATION.cff`](CITATION.cff). The paper reference will be updated on publication.

## License

MIT, see [`LICENSE`](LICENSE). The files derived from ROMS (`roms/Functionals/`, `roms/estuary_test.h`,
`roms/build_roms.sh`) remain under the ROMS license, reproduced in [`roms/LICENSE_ROMS.md`](roms/LICENSE_ROMS.md).
