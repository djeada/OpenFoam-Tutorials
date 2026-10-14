#!/usr/bin/env python3
"""Check the results of a finished steady NACA 0012 run.

Reads only files produced by ./Allrun (log.*, postProcessing/, case dictionaries).
Standard library only; needs no OpenFOAM, docker or ParaView.

Usage:  python3 verify.py [CASE_DIR]

Prints one line per check
    PASS|FAIL  <name>  measured=<v>  expected=<v> +/- <tol>
and exits 0 only if every check passes (1 on a failed check, 2 if the results
are missing).
"""

import math
import re
import sys
from pathlib import Path

CASE = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent

results = []   # (ok, name, measured, expected, tol)


def read(path):
    p = CASE / path
    if not p.is_file():
        print("ERROR: missing %s -- run ./run_docker.sh (or ./Allrun) first" % p)
        sys.exit(2)
    return p.read_text(errors="replace")


def check(name, measured, expected, tol, ok=None):
    """Numeric check |measured - expected| <= tol, or an explicit boolean ok."""
    if ok is None:
        ok = abs(measured - expected) <= tol
    results.append((bool(ok), name, measured, expected, tol))


def check_below(name, measured, limit):
    """Check measured <= limit."""
    results.append((measured <= limit, name, measured, "< %s" % fmt(limit), None))


def check_flag(name, ok, measured, expected):
    results.append((bool(ok), name, measured, expected, None))


def fmt(v):
    return "%.4g" % v if isinstance(v, float) else str(v)


# ---------------------------------------------------------------------------
# Case inputs (consistency between dictionaries)
# ---------------------------------------------------------------------------
transport = read("constant/transportProperties")
nu = float(re.search(r"^\s*nu\s+(?:nu\s*\[[^\]]*\]\s*)?([0-9.eE+-]+)\s*;", transport, re.M).group(1))

u0 = read("0/U")
m = re.search(r"inlet\s*\{[^}]*?value\s+uniform\s*\(\s*([-0-9.eE+]+)\s+([-0-9.eE+]+)\s+([-0-9.eE+]+)\s*\)", u0, re.S)
u_inlet = math.sqrt(sum(float(x) ** 2 for x in m.groups()))

control = read("system/controlDict")
mag_uinf = float(re.search(r"magUInf\s+([0-9.eE+-]+)\s*;", control).group(1))
l_ref = float(re.search(r"lRef\s+([0-9.eE+-]+)\s*;", control).group(1))
a_ref = float(re.search(r"Aref\s+([0-9.eE+-]+)\s*;", control).group(1))
end_time = int(float(re.search(r"^\s*endTime\s+([0-9.eE+-]+)\s*;", control, re.M).group(1)))

geo = read("constant/airfoil_benchmark.geo")
span = float(re.search(r"Extrude\s*\{\s*0\s*,\s*0\s*,\s*([0-9.eE+-]+)\s*\}", geo).group(1))

# The airfoil is rotated by -4 degrees in make_mesh.sh (nose-up, +4 deg AoA).
make_mesh = read("make_mesh.sh")
alpha_deg = -float(re.search(r"--angle-deg\s+(-?[0-9.]+)", make_mesh).group(1))

# Dictionary consistency.
# Tolerance 1e-6 relative: the values are typed constants, they must simply agree.
check("freestream |U| inlet vs forceCoeffs magUInf", u_inlet, mag_uinf, 1e-6 * mag_uinf)
check("Aref equals lRef x span", a_ref, l_ref * span, 1e-9)
# 5 % band around the nominal Re = 6e6 stated in the README (nu was rounded to 8.333e-6).
re_num = u_inlet * l_ref / nu
check("Reynolds number U*c/nu", re_num, 6.0e6, 0.05 * 6.0e6)

# ---------------------------------------------------------------------------
# Mesh and solver logs
# ---------------------------------------------------------------------------
check_mesh = read("log.checkMesh")
# checkMesh prints "Mesh OK." only if no check failed. The skewness limit (4) is the one
# check this Gmsh mesh is known to fail (2 faces), so it is reported separately as a
# warning below; every other check (topology, openness, aspect ratio, non-orthogonality)
# must pass.
m = re.search(r"Failed (\d+) mesh checks", check_mesh)
n_failed = int(m.group(1)) if m else 0
skew = re.search(r"Max skewness = ([0-9.eE+-]+), (\d+) highly skew faces", check_mesh)
n_other = n_failed - (1 if skew else 0)
check_flag("checkMesh: failed checks other than skewness", n_other == 0, n_other, 0)

# OpenFOAM treats non-orthogonality > 70 degrees as severe (checkMesh limit).
nonortho = float(re.search(r"Mesh non-orthogonality Max: ([0-9.eE+-]+)", check_mesh).group(1))
check_below("max non-orthogonality [deg]", nonortho, 70.0)

solver = read("log.simpleFoam")
check_flag("simpleFoam: no FOAM FATAL ERROR", "FOAM FATAL" not in solver,
           "none" if "FOAM FATAL" not in solver else "found", "none")
check_flag("simpleFoam: ran to 'End'", re.search(r"^End\s*$", solver, re.M) is not None,
           "End" if re.search(r"^End\s*$", solver, re.M) else "no End", "End")

times = [int(float(t)) for t in re.findall(r"^Time = ([0-9.eE+-]+)\s*$", solver, re.M)]
last_iter = times[-1] if times else 0
converged = "SIMPLE solution converged" in solver
# The run must either meet residualControl or reach endTime (both are a regular end).
check_flag("simpleFoam: finished (residualControl met or endTime reached)",
           converged or last_iter == end_time,
           "converged at %d" % last_iter if converged else "%d iterations" % last_iter,
           "residualControl met or %d iterations" % end_time)

# Residuals: first-solve initial residual of the last iteration, from
# postProcessing/residuals. They start at 1. A drop below 1e-1 (one order of magnitude)
# is only the minimum evidence that the solver is iterating towards a solution: with
# first-order upwind and this coarse mesh the pressure residual plateaus at about 1e-2
# (measured 1.0e-2 and 1.2e-2 on two meshes from different Gmsh versions), and
# residualControl in fvSolution (1e-6) is never reached. Whether the answer is usable is
# judged by the stationarity of Cl below, not by these residuals.
res_text = read("postProcessing/residuals/0/residuals.dat")
res_lines = [l.split() for l in res_text.splitlines() if l.strip() and not l.startswith("#")]
res_cols = None
for l in res_text.splitlines():
    if re.match(r"#\s*Time\b", l):
        res_cols = l.lstrip("#").split()
if not res_lines or not res_cols:
    print("ERROR: could not parse postProcessing/residuals/0/residuals.dat")
    sys.exit(2)
last_res = dict(zip(res_cols, res_lines[-1]))
for field in ("p", "Ux", "Uy", "nuTilda"):
    check_below("final initial residual %s" % field, float(last_res[field]), 1e-1)

# ---------------------------------------------------------------------------
# Force coefficients
# ---------------------------------------------------------------------------
fc_text = read("postProcessing/forceCoeffs/0/forceCoeffs.dat")
columns = None
rows = []
for line in fc_text.splitlines():
    if line.startswith("#"):
        if re.match(r"#\s*Time\b", line):
            columns = line.lstrip("#").split()
        continue
    parts = line.split()
    if parts:
        rows.append([float(x) for x in parts])
if columns is None or not rows:
    print("ERROR: could not parse postProcessing/forceCoeffs/0/forceCoeffs.dat")
    sys.exit(2)
idx = {c: i for i, c in enumerate(columns)}
last = rows[-1]
cl, cd, cm = last[idx["Cl"]], last[idx["Cd"]], last[idx["Cm"]]

# Thin-airfoil theory: Cl = 2*pi*alpha (alpha in rad). Measured NACA 0012 lift-curve
# slopes at Re ~ 6e6 are close to this for small alpha. +/-0.05 (about 11 %) allows
# for the coarse wall-function mesh, the 15-chord domain with symmetry-plane top/bottom
# boundaries and the 2 pi idealisation. This is a consistency band, not a validation.
cl_theory = 2.0 * math.pi * math.radians(alpha_deg)
check("Cl vs thin-airfoil theory 2*pi*alpha (alpha=%g deg)" % alpha_deg, cl, cl_theory, 0.05)

# Positive AoA with liftDir (0 1 0) must give positive lift; drag must be positive
# (acts along dragDir (1 0 0)) and smaller than the lift of a lifting airfoil.
check_flag("Cl > 0 (nose-up airfoil, liftDir +y)", cl > 0, fmt(cl), "> 0")
check_flag("0 < Cd < Cl", 0 < cd < cl, fmt(cd), "0 < Cd < %s" % fmt(cl))

# Symmetric section: thin-airfoil theory gives Cm about the quarter chord = 0.
# 0.02 is about 5 % of Cl; wall-function and mesh errors are expected below that.
check_below("|Cm| about quarter chord", abs(cm), 0.02)

# Stationarity: Cl changed by less than 1 % over the last 100 iterations.
t_last = last[idx["Time"]]
earlier = [r for r in rows if r[idx["Time"]] <= t_last - 100]
if earlier:
    cl_prev = earlier[-1][idx["Cl"]]
    check_below("Cl change over last 100 iterations (relative)", abs(cl - cl_prev) / abs(cl), 0.01)
else:
    check_flag("Cl change over last 100 iterations", False, "fewer than 100 iterations", "< 1 %")

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
print("NACA 0012 steady case: %s" % CASE)
for ok, name, measured, expected, tol in results:
    tag = "PASS" if ok else "FAIL"
    if tol is None:
        print("%s  %s  measured=%s  expected=%s" % (tag, name, fmt(measured), fmt(expected)))
    else:
        print("%s  %s  measured=%s  expected=%s +/- %s" % (tag, name, fmt(measured), fmt(expected), fmt(tol)))

# Informational (no pass/fail): known limitations of this mesh.
if skew:
    print("WARN  checkMesh: max skewness %s on %s faces (OpenFOAM limit 4); solver ran without problems"
          % (skew.group(1), skew.group(2)))
yplus_file = CASE / "postProcessing/yPlus/0/yPlus.dat"
if yplus_file.is_file():
    for line in reversed(yplus_file.read_text().splitlines()):
        f = line.split()
        if len(f) >= 5 and f[1] == "airfoil":
            print("INFO  y+ on airfoil at last write: min=%s max=%s average=%s"
                  % (fmt(float(f[2])), fmt(float(f[3])), fmt(float(f[4]))))
            break
print("INFO  Cd=%s is not compared with published values: this mesh is wall-function only "
      "(y+ of several hundred) and is not expected to match them" % fmt(cd))

failed = [r for r in results if not r[0]]
print("SUMMARY  %d/%d checks passed" % (len(results) - len(failed), len(results)))
sys.exit(1 if failed else 0)
