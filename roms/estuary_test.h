/*
** River-dominated channel with a prescribed (nudged) longitudinal salinity
** field, used for the equilibrium test of Sec. IV B of the paper.
**
** The case is built on the ROMS ESTUARY_TEST application: the CPP name is the
** same, and the analytical functions in Functionals/ replace the default ones.
**
** Parameters passed through the USER array of the ROMS input file
** (set by scripts/make_run.py; x is measured from the sea boundary, west):
**   user(1) = D0     uniform-flow depth (m)
**   user(2) = river velocity (m/s) times 4/3, because ana_m2obc spreads the
**             discharge over four rows while the channel has three wet rows
**   user(3) = S_sea  salinity at the sea boundary
**   user(4) = x_s    centre of the salinity transition (m)
**   user(5) = ell    half-width of the salinity transition (m)
**   user(6) = s0     surface slope of the uniform flow, also the bed slope
**   user(7) = A_D    amplitude of the local depth reduction, sech^2 shape (m)
**   user(8) = A_h    half of the cumulative rise of bed and surface, tanh shape (m)
**
** Salinity: S(x) = S_sea/2 [1 - tanh((x - x_s)/ell)].
** Bed:      h(x) = D0 - s0 x - A_D sech^2((x - x_s)/ell)
**                  - A_h [tanh((x - x_s)/ell) - tanh(-x_s/ell)].
** The bed is fixed in the equilibrium test (SAND_MORPH_FAC = 0 in sediment.in).
*/

#define UV_ADV
#define UV_LOGDRAG
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define SALINITY
#define SOLVE3D
#define SEDIMENT
#ifdef SEDIMENT
# undef  SUSPLOAD
# define BEDLOAD_MPM
# define SED_MORPH
# define SED_UPWIND
# define SLOPE_LESSER
#endif
#define AVERAGES
#define GLS_MIXING
#undef  MY25_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# define N2S2_HORAVG
# define RI_SPLINES
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SEDIMENT
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_BPFLUX
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_TOBC
#define ANA_TCLIMA
#define ANA_NUDGCOEF
