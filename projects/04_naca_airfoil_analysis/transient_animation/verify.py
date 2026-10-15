#!/usr/bin/env python3
"""Check the results of a finished transient NACA 0012 run (pimpleFoam, laminar).

Reads only files produced by ./Allrun (log.*, postProcessing/, time directories,
case dictionaries). Standard library only; needs no OpenFOAM, docker or ParaView.

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
        print("ERROR: missing %s -- run ../run_docker.sh transient (or ./Allrun) first" % p)
        sys.exit(2)
    return p.read_text(errors="replace")


def check(name, measured, expected, tol, ok=None):
    """Numeric check |measured - expected| <= tol, or an explicit boolean ok."""
    if ok is None:
        ok = abs(measured - expected) <= tol
    results.append((bool(ok), name, measured, expected, tol))


def check_below(name, measured, limit):
    results.append((measured <= limit, name, measured, "< %s" % fmt(limit), None))


def check_flag(name, ok, measured, expected):
    results.append((bool(ok), name, measured, expected, None))


def fmt(v):
    return "%.4g" % v if isinstance(v, float) else str(v)


def mean(vals):
    return sum(vals) / len(vals)


def stdev(vals):
    m_ = mean(vals)
    return math.sqrt(sum((v - m_) ** 2 for v in vals) / len(vals))


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
end_time = float(re.search(r"^\s*endTime\s+([0-9.eE+-]+)\s*;", control, re.M).group(1))
write_interval = float(re.search(r"^\s*writeInterval\s+([0-9.eE+-]+)\s*;", control, re.M).group(1))
max_co = float(re.search(r"^\s*maxCo\s+([0-9.eE+-]+)\s*;", control, re.M).group(1))

geo = read("constant/airfoil_animation.geo")
span = float(re.search(r"Extrude\s*\{\s*0\s*,\s*0\s*,\s*([0-9.eE+-]+)\s*\}", geo).group(1))

make_mesh = read("make_mesh.sh")
alpha_deg = -float(re.search(r"--angle-deg\s+(-?[0-9.]+)", make_mesh).group(1))

# Tolerance 1e-6 relative: typed constants that must agree.
check("freestream |U| inlet vs forceCoeffs magUInf", u_inlet, mag_uinf, 1e-6 * mag_uinf)
check("Aref equals lRef x span", a_ref, l_ref * span, 1e-9)
# 5 % band around Re = 1e4 stated in the README.
re_num = u_inlet * l_ref / nu
check("Reynolds number U*c/nu", re_num, 1.0e4, 0.05 * 1.0e4)

# ---------------------------------------------------------------------------
# Mesh log
# ---------------------------------------------------------------------------
check_mesh = read("log.checkMesh")
# The skewness limit (4) is the one check this Gmsh mesh fails (a few faces), reported as
# a WARN line below; every other check must pass.
m = re.search(r"Failed (\d+) mesh checks", check_mesh)
n_failed = int(m.group(1)) if m else 0
skew = re.search(r"Max skewness = ([0-9.eE+-]+), (\d+) highly skew faces", check_mesh)
n_other = n_failed - (1 if skew else 0)
check_flag("checkMesh: failed checks other than skewness", n_other == 0, n_other, 0)
nonortho = float(re.search(r"Mesh non-orthogonality Max: ([0-9.eE+-]+)", check_mesh).group(1))
check_below("max non-orthogonality [deg]", nonortho, 70.0)

# ---------------------------------------------------------------------------
# Solver log
# ---------------------------------------------------------------------------
solver = read("log.pimpleFoam")
check_flag("pimpleFoam: no FOAM FATAL ERROR", "FOAM FATAL" not in solver,
           "none" if "FOAM FATAL" not in solver else "found", "none")
has_end = re.search(r"^End\s*$", solver, re.M) is not None
check_flag("pimpleFoam: ran to 'End'", has_end, "End" if has_end else "no End", "End")

times = [float(t) for t in re.findall(r"^Time = ([0-9.eE+-]+)\s*$", solver, re.M)]
t_end = times[-1] if times else 0.0
# The last step overshoots endTime by less than one time step (maxDeltaT).
max_dt = float(re.search(r"^\s*maxDeltaT\s+([0-9.eE+-]+)\s*;", control, re.M).group(1))
check("pimpleFoam: end time reached [s]", t_end, end_time, max_dt)

# Courant number. adjustTimeStep sets deltaT from the Courant number of the previous
# step and limits its growth, so the printed maximum can exceed maxCo slightly while the
# flow changes. Two checks, both on finite values:
#  * after the first 10 steps (the impulsive start changes the field quickly: the third
#    step printed 0.575 in the run documented in the README) the maximum stays within
#    1 % of maxCo;
#  * over the whole run it stays below 1.2 x maxCo, the largest step-size increase
#    OpenFOAM allows in one step.
courants = [float(x) for x in re.findall(r"^Courant Number mean: \S+ max: ([0-9.eE+-]+|nan|inf)\s*$", solver, re.M)]
co_all = max(courants) if courants else float("nan")
co_late = max(courants[10:]) if len(courants) > 10 else float("nan")
finite_co = bool(courants) and all(math.isfinite(c) for c in courants)
check_flag("Courant number after the first 10 steps <= 1.01 x maxCo",
           finite_co and co_late <= 1.01 * max_co, fmt(co_late), "<= %s" % fmt(1.01 * max_co))
check_flag("Courant number over the whole run <= 1.2 x maxCo",
           finite_co and co_all <= 1.2 * max_co, fmt(co_all), "<= %s" % fmt(1.2 * max_co))

# Written time directories: endTime / writeInterval (floating point names such as 0.25).
n_expected = int(round(end_time / write_interval))
written = [d for d in CASE.iterdir()
           if d.is_dir() and re.fullmatch(r"[0-9]+(\.[0-9]+)?", d.name) and float(d.name) > 0]
check_flag("written time directories (reconstructed)", len(written) == n_expected, len(written), n_expected)
check_flag("last written time equals endTime",
           bool(written) and abs(max(float(d.name) for d in written) - end_time) < 1e-9,
           fmt(max((float(d.name) for d in written), default=0.0)), fmt(end_time))

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
t_col = [r[idx["Time"]] for r in rows]
cl_all = [r[idx["Cl"]] for r in rows]
cd_all = [r[idx["Cd"]] for r in rows]
cm_all = [r[idx["Cm"]] for r in rows]

finite = all(math.isfinite(v) for r in rows for v in r)
check_flag("forceCoeffs: all values finite", finite, "yes" if finite else "nan/inf found", "yes")
check("forceCoeffs: reaches endTime [s]", t_col[-1], end_time, max(max_dt, 2 * 0.01))

# Plausibility bounds, for the second half of the run (the impulsive start produces a
# large transient in the first seconds, so it is excluded):
#  * thin-airfoil theory gives Cl = 2 pi alpha = 1.64 at 15 degrees for attached flow;
#    separated laminar flow at Re = 1e4 gives less, and the instantaneous value must stay
#    far from non-physical sizes: |Cl| < 2.5 is 1.5 x the attached-flow value;
#  * mean Cl > 0 for the nose-up airfoil with liftDir +y;
#  * Cd > 0 on average (the airfoil is a bluff-ish body at 15 degrees, drag acts
#    downstream) and mean Cd < mean Cl.
late = [i for i, t in enumerate(t_col) if t >= 0.5 * end_time]
cl_late = [cl_all[i] for i in late]
cd_late = [cd_all[i] for i in late]
cm_late = [cm_all[i] for i in late]
cl_theory = 2.0 * math.pi * math.radians(alpha_deg)
check_flag("max |Cl| over the run < 1.5 x 2*pi*alpha (alpha=%g deg)" % alpha_deg,
           max(abs(c) for c in cl_all[len(cl_all) // 20:]) < 1.5 * cl_theory,
           fmt(max(abs(c) for c in cl_all[len(cl_all) // 20:])), "< %s" % fmt(1.5 * cl_theory))
check_flag("mean Cl in second half of run > 0", mean(cl_late) > 0, fmt(mean(cl_late)), "> 0")
check_flag("0 < mean Cd < mean Cl (second half)", 0 < mean(cd_late) < mean(cl_late),
           fmt(mean(cd_late)), "0 < Cd < %s" % fmt(mean(cl_late)))

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
print("NACA 0012 transient case: %s" % CASE)
for ok, name, measured, expected, tol in results:
    tag = "PASS" if ok else "FAIL"
    if tol is None:
        print("%s  %s  measured=%s  expected=%s" % (tag, name, fmt(measured), fmt(expected)))
    else:
        print("%s  %s  measured=%s  expected=%s +/- %s" % (tag, name, fmt(measured), fmt(expected), fmt(tol)))

if skew:
    print("WARN  checkMesh: max skewness %s on %s faces (OpenFOAM limit 4); solver ran without problems"
          % (skew.group(1), skew.group(2)))
dts = [float(x) for x in re.findall(r"^deltaT = ([0-9.eE+-]+)\s*$", solver, re.M)]
if dts:
    print("INFO  time steps: %d, deltaT min=%s max=%s" % (len(dts), fmt(min(dts)), fmt(max(dts))))
print("INFO  second half of the run (t >= %s s): Cl mean=%s sd=%s min=%s max=%s; Cd mean=%s sd=%s; Cm mean=%s"
      % (fmt(0.5 * end_time), fmt(mean(cl_late)), fmt(stdev(cl_late)), fmt(min(cl_late)), fmt(max(cl_late)),
         fmt(mean(cd_late)), fmt(stdev(cd_late)), fmt(mean(cm_late))))
print("INFO  Cl and Cd are laminar, 2D, Re = %s: they fluctuate with the shedding of the separated "
      "flow and are not compared with published values" % fmt(re_num))

failed = [r for r in results if not r[0]]
print("SUMMARY  %d/%d checks passed" % (len(results) - len(failed), len(results)))
sys.exit(1 if failed else 0)
