#!/usr/bin/env python3
"""Check the results of a finished heated pipe flow run.

Reads only files written by the run (log.*, the last time directory,
postProcessing/); it does not run OpenFOAM and needs no Docker.

Usage:  python3 verify.py [case_directory]      (default: this directory)
Exit status: 0 if every check passes, 1 if a check fails, 2 if results are missing.
"""
import glob
import math
import os
import re
import sys

# ---- case data (must match 0/, constant/ and system/blockMeshDict) --------
T_IN = 300.0                 # K, inlet temperature        (0/T)
T_WALL = 350.0               # K, wall temperature         (0/T)
U_IN = 0.5                   # m/s, inlet velocity         (0/U)
P_REF = 101325.0             # Pa, absolute pressure       (0/p_rgh)
MOL_WEIGHT = 28.9            # kg/kmol                     (thermophysicalProperties)
R_UNIVERSAL = 8314.47        # J/(kmol K), OpenFOAM value
CP = 1005.0                  # J/(kg K)
MU = 1.84e-5                 # Pa s
PR = 0.71
HEIGHT = 0.1                 # m (y), blockMeshDict
DEPTH = 0.01                 # m (z), blockMeshDict
LENGTH = 2.0                 # m (x), blockMeshDict
N_INLET_FACES = 20
N_CELLS = 2000
END_TIME = 1000
RES_TARGET = 1e-4            # residualControl in fvSolution

results = []


def check(name, ok, measured, expected):
    results.append(ok)
    print("%s  %-30s measured=%s  expected=%s" % ("PASS" if ok else "FAIL", name, measured, expected))


def missing(msg):
    print("ERROR: %s\nRun ./run_docker.sh (or ./Allrun) first." % msg, file=sys.stderr)
    sys.exit(2)


def read_text(path):
    if not os.path.isfile(path):
        missing("missing file %s" % path)
    with open(path, errors="replace") as f:
        return f.read()


def latest_time(case):
    times = []
    for d in os.listdir(case):
        if os.path.isdir(os.path.join(case, d)):
            try:
                times.append((float(d), d))
            except ValueError:
                pass
    times = [t for t in times if t[0] > 0]
    if not times:
        missing("no time directory > 0 in %s" % case)
    return max(times)[1]


NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def scalar_list(block):
    """Values of 'nonuniform List<scalar> N ( ... )' or 'uniform X' in block."""
    m = re.search(r"nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)", block)
    if m:
        vals = [float(x) for x in m.group(2).split()]
        if len(vals) != int(m.group(1)):
            missing("malformed list in field file")
        return vals
    m = re.search(r"uniform\s+(%s)" % NUM, block)
    if m:
        return [float(m.group(1))]
    return None


def internal(text):
    i = text.index("internalField")
    j = text.index("boundaryField")
    return scalar_list(text[i:j])


def patch_values(text, patch, nfaces):
    b = text[text.index("boundaryField"):]
    m = re.search(r"\n\s*%s\s*\{(.*?)\n\s*\}" % patch, b, re.S)
    if not m:
        missing("patch %s not found in field file" % patch)
    mv = re.search(r"\bvalue\b(.*)", m.group(1), re.S)
    vals = scalar_list(mv.group(1)) if mv else None
    if vals is None:
        missing("no values for patch %s" % patch)
    return vals * nfaces if len(vals) == 1 else vals


def main():
    case = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))

    # ---- 1. mesh ---------------------------------------------------------
    chk = read_text(os.path.join(case, "log.checkMesh"))
    check("checkMesh reports Mesh OK", "Mesh OK" in chk and "Failed" not in chk,
          "Mesh OK" if "Mesh OK" in chk else "not OK", "Mesh OK")
    m = re.search(r"cells:\s+(\d+)", chk)
    n = int(m.group(1)) if m else -1
    check("cell count", n == N_CELLS, n, "%d (exact)" % N_CELLS)

    # ---- 2. solver run ---------------------------------------------------
    log = read_text(os.path.join(case, "log.buoyantSimpleFoam"))
    clean = ("FOAM FATAL" not in log and "Floating point exception" not in log
             and re.search(r"^End\s*$", log, re.M) is not None)
    check("solver finished cleanly", clean, "clean" if clean else "error/abort",
          "ends with 'End', no FATAL/FPE")

    m = re.search(r"SIMPLE solution converged in (\d+) iterations", log)
    it = int(m.group(1)) if m else None
    # Criterion: converges before endTime. 134 iterations was measured; no
    # reference iteration count exists, so only "converged at all" is tested.
    check("residualControl satisfied", it is not None and it < END_TIME,
          it if it is not None else "not converged", "< %d iterations" % END_TIME)

    last = {}
    for f, v in re.findall(r"Solving for (\w+), Initial residual = (%s)" % NUM, log):
        last[f] = float(v)
    worst = max(last.values()) if last else float("nan")
    # Same criterion as residualControl (1e-4); the final iteration must be below it.
    check("final initial residuals", bool(last) and worst <= RES_TARGET,
          "%.2e" % worst, "<= %.0e" % RES_TARGET)

    # ---- 3. fields in the last time directory ----------------------------
    t = latest_time(case)
    tdir = os.path.join(case, t)
    Ttxt = read_text(os.path.join(tdir, "T"))
    phitxt = read_text(os.path.join(tdir, "phi"))
    T = internal(Ttxt)
    if T is None or len(T) != N_CELLS:
        missing("could not read T internalField from %s" % tdir)

    # Maximum principle: no heat sources, so T stays between inlet and wall values.
    # 0.01 K covers the linear-solver tolerance (relTol 0.1, tol 1e-6 in h).
    tmin, tmax = min(T), max(T)
    tol = 0.01
    check("T within [T_in, T_wall]", tmin >= T_IN - tol and tmax <= T_WALL + tol,
          "%.3f .. %.3f K" % (tmin, tmax), "[%.0f, %.0f] K +- %.2f" % (T_IN, T_WALL, tol))

    # Mass flow. phi is negative on the inlet. rho_in from the ideal gas law.
    phi_in = patch_values(phitxt, "inlet", N_INLET_FACES)
    phi_out = patch_values(phitxt, "outlet", N_INLET_FACES)
    m_in = -sum(phi_in)
    rho_in = P_REF * MOL_WEIGHT / (R_UNIVERSAL * T_IN)
    m_exp = rho_in * U_IN * HEIGHT * DEPTH
    # 0.2 %: p at the inlet differs from 101325 Pa by ~0.2 Pa (2e-6); the rest is
    # rounding of R. Anything larger would indicate a wrong BC or property.
    err = abs(m_in - m_exp) / m_exp
    check("inlet mass flow", err <= 2e-3, "%.4e kg/s" % m_in, "%.4e +- 0.2%%" % m_exp)

    # Steady state: mass in = mass out. 0.1 % is >> the continuity error at the
    # 1e-4 residual level (log: sum local ~5e-6 of the throughput).
    m_out = sum(phi_out)
    err = abs(m_out - m_in) / m_in
    check("mass balance in/out", err <= 1e-3, "%.2e rel. error" % err, "<= 1e-3")

    # ---- 4. energy balance ----------------------------------------------
    dats = sorted(glob.glob(os.path.join(case, "postProcessing", "wallHeatFlux", "*", "*.dat")))
    if not dats:
        missing("postProcessing/wallHeatFlux not found (controlDict function object)")
    rows = {}
    for line in read_text(dats[-1]).splitlines():
        p = line.split()
        if len(p) == 5 and not line.startswith("#"):
            rows[p[1]] = (float(p[0]), float(p[4]))      # time, integral [W]
    if not {"topWall", "bottomWall"} <= set(rows):
        missing("wallHeatFlux.dat has no row for topWall/bottomWall")
    q_w = rows["topWall"][1] + rows["bottomWall"][1]

    # Bulk (mass-flux weighted) outlet temperature from the outlet face values.
    T_out = patch_values(Ttxt, "outlet", N_INLET_FACES)
    tb_out = sum(p * x for p, x in zip(phi_out, T_out)) / sum(phi_out)
    # Analytical: all wall heat goes into the stream, q_w = m_dot cp (Tb_out - T_in).
    # Neglected: change of kinetic energy (<1e-3 K here). The tolerance of 2 % of
    # the temperature rise absorbs the first-order upwind face interpolation of T
    # at the outlet and the 1e-4 residual level.
    dt_exp = q_w / (m_in * CP)
    dt_meas = tb_out - T_IN
    check("outlet bulk T rise (q_w / m cp)", abs(dt_meas - dt_exp) <= 0.02 * dt_exp,
          "%.2f K" % dt_meas, "%.2f K +- 2%%" % dt_exp)

    # ---- 5. heat transfer level (consistency, no validated reference) ------
    # Mean Nu from the log-mean temperature difference, against Dittus-Boelter for
    # fully developed flow. Entry-region flow has Nu >= fully developed value, and
    # even short ducts rarely exceed ~3x; this is a sanity band, not a validation.
    area = 2 * LENGTH * DEPTH
    q_flux = q_w / area
    d1, d2 = T_WALL - T_IN, T_WALL - tb_out
    lmtd = (d1 - d2) / math.log(d1 / d2)
    k = MU * CP / PR
    dh = 2 * HEIGHT
    nu_mean = q_flux / lmtd * dh / k
    nu_inlet_re = U_IN * dh * rho_in / MU
    nu_db = 0.023 * nu_inlet_re ** 0.8 * PR ** 0.4
    ratio = nu_mean / nu_db
    check("mean Nu / Dittus-Boelter", 1.0 <= ratio <= 3.0,
          "%.2f (Nu=%.1f, DB=%.1f)" % (ratio, nu_mean, nu_db), "1.0 .. 3.0 (sanity band)")

    n_fail = results.count(False)
    print("\n%d of %d checks passed (time directory %s)" % (len(results) - n_fail, len(results), t))
    print("RESULT: %s" % ("PASS" if n_fail == 0 else "FAIL"))
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
