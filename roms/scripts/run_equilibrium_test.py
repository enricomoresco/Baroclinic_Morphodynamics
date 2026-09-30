"""Run the ROMS equilibrium test of Sec. IV B (fixed bed, 3 days per run).

Stage 1, for eps = 0.1 and 0.3:
    eq_bt              barotropic twin on the uniform-flow bed (reference)
    eq_XX_flat         baroclinic run on the uniform-flow bed
    eq_XX_k+K.K        baroclinic run on the bed with the full cumulative part and K times the
                       local part of the first-order correction, K = -1, 0, 0.5, 1, 1.5, 2
Stage 2, eps = 0.1:
    setup_F.F          local part at its best-fit amplitude from stage 1, cumulative part
                       scaled by F = 0.6 ... 1.4

A bed profile is an equilibrium for bed-load transport if and only if the bed stress it
produces is uniform: analyse.py measures the departure of each run from the reference.

usage:  python run_equilibrium_test.py [--exe path/to/romsS] [--jobs 2]
Each run takes one to two minutes on one core; the 20 runs take about 20 minutes on 2 cores.
Completed runs are skipped, so the script can be restarted.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor

from roms_tools import KS, RUNS, SETUP_FACTORS, amplitudes, roms_exe, run_roms, tag


def run_all(jobs, exe, n_jobs):
    with ThreadPoolExecutor(n_jobs) as ex:
        for name in ex.map(lambda j: run_roms(*j, exe=exe), jobs):
            print("done:", name, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", default=None, help="ROMS executable (default roms/romsS or $ROMS_EXE)")
    p.add_argument("--jobs", type=int, default=2, help="runs executed in parallel")
    a = p.parse_args()
    exe = roms_exe(a.exe)

    # ---- stage 1
    jobs = [("eq_bt", 0.0, 0.0, 0.0)]
    for eps in (0.1, 0.3):
        A_D, A_h = amplitudes(eps)
        jobs.append((f"eq_{tag(eps)}_flat", eps, 0.0, 0.0))
        jobs += [(f"eq_{tag(eps)}_k{k:+.1f}", eps, k * A_D, A_h) for k in KS]
    print(f"stage 1: {len(jobs)} runs in {RUNS}", flush=True)
    run_all(jobs, exe, a.jobs)

    # ---- stage 2: set-up scan on the bed with the best-fit local part
    from analyse import local_fit
    k_opt = local_fit(0.1)["k_opt"]
    A_D, A_h = amplitudes(0.1)
    jobs = [(f"setup_{f:.1f}", 0.1, k_opt * A_D, f * A_h) for f in SETUP_FACTORS]
    print(f"stage 2: {len(jobs)} runs, local part at {k_opt:.3f} times the theory", flush=True)
    run_all(jobs, exe, a.jobs)
    print("all runs completed; now run analyse.py")


if __name__ == "__main__":
    main()
