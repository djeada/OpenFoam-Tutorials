#!/usr/bin/env python3
"""Check the results of a finished dam-break run (read-only, standard library only).

Usage:  python3 verify.py [case_dir]      (default: directory of this file)
Run it after ./Allrun or ./run_docker.sh. Reads log.*, system/controlDict,
0/alpha.water and the written time directories (ASCII format).
Exit code 0 only if every check passes; 2 if results are missing.

There is no analytical reference solution for this case. The checks are
(a) solver/mesh health, (b) numerical limits taken from the case files and
(c) internal consistency (mass balance, bounds). Each tolerance is explained
next to its check.
"""
import os
import re
import sys

case = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__))
results = []   # (ok, name, measured, expected)


def read(path):
    with open(os.path.join(case, path)) as f:
        return f.read()


def check(ok, name, measured, expected):
    results.append((bool(ok), name, measured, expected))


def die(msg):
    print("ERROR: " + msg)
    print("Run ./Allrun (or ./run_docker.sh) first.")
    sys.exit(2)


def dict_value(text, key):
    m = re.search(r"^\s*%s\s+([^;\s]+);" % re.escape(key), text, re.M)
    return float(m.group(1)) if m else None


def read_internal_field(path):
    """Return list of floats of a nonuniform scalar internalField, or the uniform value."""
    text = read(path)
    m = re.search(r"internalField\s+uniform\s+([^;\s]+);", text)
    if m:
        return float(m.group(1))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not m:
        raise ValueError("cannot parse " + path)
    values = [float(x) for x in m.group(2).split()]
    if len(values) != int(m.group(1)):
        raise ValueError("bad list length in " + path)
    return values


# ---------------------------------------------------------------- inputs
for f in ("log.blockMesh", "log.interFoam", "log.checkMesh", "log.setFields",
          "system/controlDict", "0/alpha.water"):
    if not os.path.isfile(os.path.join(case, f)):
        die("missing " + f)

control = read("system/controlDict")
end_time = dict_value(control, "endTime")
write_interval = dict_value(control, "writeInterval")
max_co = dict_value(control, "maxCo")
max_alpha_co = dict_value(control, "maxAlphaCo")
max_dt = dict_value(control, "maxDeltaT")
fvsol = read("system/fvSolution")
m = re.search(r"p_rgh\s*\{[^}]*?tolerance\s+([^;\s]+);", fvsol)
p_tol = float(m.group(1))

log = read("log.interFoam")
logs_all = log + read("log.blockMesh") + read("log.setFields") + read("log.checkMesh")

# ---------------------------------------------------------------- mesh
m = re.search(r"nCells:\s*(\d+)", read("log.blockMesh"))
ncells = int(m.group(1)) if m else -1
# 46 x 46 x 1 cells are set in blockMeshDict: exact integer.
check(ncells == 46 * 46, "mesh cells", ncells, "2116 (exact)")

# checkMesh prints "Mesh OK." only if all geometry/topology checks pass.
check("Mesh OK." in read("log.checkMesh"), "checkMesh", "Mesh OK" if "Mesh OK." in read("log.checkMesh") else "failed", "Mesh OK")

# ---------------------------------------------------------------- run health
bad = re.findall(r"FOAM FATAL|FOAM Warning|Floating point exception|Segmentation", logs_all)
check(not bad, "fatal errors/warnings in logs", len(bad), "0 (exact)")

times = [float(x) for x in re.findall(r"^Time = ([0-9.eE+-]+)\s*$", log, re.M)]
check(len(times) > 0, "time steps found", len(times), "> 0")
if not times:
    die("no time steps in log.interFoam")
finished = log.rstrip().endswith("End")
check(finished, "solver finished (End)", "yes" if finished else "no", "yes")
# Final time is written exactly; timePrecision 6 -> 1e-6 is the print resolution.
check(abs(times[-1] - end_time) <= 1e-6, "end time reached [s]", times[-1], "%g +- 1e-6" % end_time)

# ---------------------------------------------------------------- Courant numbers
co = [float(x) for x in re.findall(r"^Courant Number mean: \S+ max: (\S+)", log, re.M)]
aco = [float(x) for x in re.findall(r"^Interface Courant Number mean: \S+ max: (\S+)", log, re.M)]
# The time step is adjusted after Co is evaluated and may grow by up to 20 % per
# step (interFoam setDeltaT), so Co can exceed maxCo briefly; 1.2*maxCo is that bound.
check(max(co) <= 1.2 * max_co, "max Courant number", "%.3f" % max(co), "<= %.2f (1.2 x maxCo)" % (1.2 * max_co))
check(max(aco) <= 1.2 * max_alpha_co, "max interface Courant", "%.3f" % max(aco),
      "<= %.2f (1.2 x maxAlphaCo)" % (1.2 * max_alpha_co))

# ---------------------------------------------------------------- time step
dts = [float(x) for x in re.findall(r"^deltaT = (\S+)", log, re.M)]
# Upper limit is maxDeltaT. Lower bound 1e-4 s corresponds to Co = 0.5 at
# |U| = 0.5*0.0127/1e-4 = 63 m/s, far above anything physical here (|U| ~ few m/s).
check(max(dts) <= max_dt + 1e-12 and min(dts) >= 1e-4, "time step range [s]",
      "%.2e .. %.2e" % (min(dts), max(dts)), "1e-4 .. %g" % max_dt)

# ---------------------------------------------------------------- alpha bounds
mins = [float(x) for x in re.findall(r"Min\(alpha\.water\) = (\S+)", log)]
maxs = [float(x) for x in re.findall(r"Max\(alpha\.water\) = (\S+)", log)]
# MULES with MULESCorr keeps alpha bounded to a small undershoot; 1e-3 is 0.1 % of
# the 0..1 range (a few times the usual undershoot) and changes density by < 1 kg/m3.
check(min(mins) >= -1e-3, "min(alpha.water)", "%.2e" % min(mins), ">= -1e-3")
check(max(maxs) <= 1 + 1e-3, "max(alpha.water)", "%.5f" % max(maxs), "<= 1.001")

# ---------------------------------------------------------------- pressure solver
# last p_rgh solve in each time step is p_rghFinal: relTol 0, so it must reach `tolerance`
per_step = {}
t = None
for line in log.splitlines():
    mt = re.match(r"^Time = (\S+)\s*$", line)
    if mt:
        t = mt.group(1)
        continue
    mp = re.search(r"Solving for p_rgh, Initial residual = \S+, Final residual = (\S+),", line)
    if mp and t is not None:
        per_step[t] = float(mp.group(1))
worst = max(per_step.values())
# Allow the printed value to sit at the tolerance (strictly <= tolerance if converged).
check(worst <= p_tol, "p_rgh final residual (max)", "%.3e" % worst, "<= %g (fvSolution tolerance)" % p_tol)

# ---------------------------------------------------------------- mass balance
fractions = re.findall(r"^Phase-1 volume fraction = (\S+)", log, re.M)
frac = [float(x) for x in fractions]
f0 = frac[0]
alpha0 = read_internal_field("0/alpha.water")
n_water = sum(1 for a in alpha0 if a > 0.5)
check(abs(f0 - n_water / float(ncells)) <= 1e-6, "initial water fraction = wet cells / cells",
      "%.6f" % f0, "%.6f +- 1e-6 (%d/%d cells)" % (n_water / float(ncells), n_water, ncells))

# setFields only selects whole cells, so the column differs from the exact
# geometry (0.1461 x 0.292 in 0.584 x 0.584). The error is at most about half a
# cell layer along the column perimeter: (0.146+0.292)*0.0127/2/0.341 = 0.008 -> 0.01.
m = re.search(r"box \(([^)]*)\) \(([^)]*)\)", read("system/setFieldsDict"))
lo = [float(x) for x in m.group(1).split()]
hi = [float(x) for x in m.group(2).split()]
exact = (hi[0] - lo[0]) * (hi[1] - lo[1]) / (0.584 * 0.584)
check(abs(f0 - exact) <= 0.01, "initial water fraction vs exact geometry",
      "%.4f" % f0, "%.4f +- 0.01" % exact)

# No water source exists: the fraction may only fall (water leaving through
# `atmosphere`). 1e-5 is the print precision of the log (6 significant digits).
check(max(frac) <= f0 + 1e-5, "water fraction never increases", "%.6f" % max(frac), "<= %.6f" % (f0 + 1e-5))

# Time directories in write order
tdirs = []
for d in os.listdir(case):
    p = os.path.join(case, d)
    if os.path.isdir(p) and os.path.isfile(os.path.join(p, "alpha.water")):
        try:
            tdirs.append((float(d), d))
        except ValueError:
            pass
tdirs.sort()
n_expected = int(round(end_time / write_interval)) + 1
check(len(tdirs) == n_expected, "written time directories (incl. 0)", len(tdirs), "%d (exact)" % n_expected)

NX = 46   # cells per row (blockMeshDict), cells numbered x fastest
first_wet_top = None
first_wall = None
for tv, d in tdirs:
    a = read_internal_field(d + "/alpha.water")
    if not isinstance(a, list):
        continue
    if first_wet_top is None and max(a[-NX:]) > 0.01:
        first_wet_top = tv
    if first_wall is None and a[NX - 1] > 0.5:   # bottom row, cell next to rightWall
        first_wall = tv

# Fraction is constant (to log precision) as long as the top cell row is dry.
if first_wet_top is None:
    dry_until = end_time
else:
    dry_until = max(tv for tv, _ in tdirs if tv < first_wet_top)
dev = 0.0
t = None
for line in log.splitlines():
    mt = re.match(r"^Time = (\S+)\s*$", line)
    if mt:
        t = float(mt.group(1))
        continue
    mf = re.match(r"^Phase-1 volume fraction = (\S+)", line)
    if mf and t is not None and t <= dry_until:
        dev = max(dev, abs(float(mf.group(1)) - f0))
check(dev <= 1e-5, "water conserved while top row dry (t <= %.2f s)" % dry_until,
      "%.1e" % dev, "<= 1e-5")

# Total loss: only the jet tip can leave through the open top, so a loss of more
# than a few percent would indicate a problem. 5 % is a plausibility bound
# (no reference exists); measured loss is about 2.5 %.
loss = (f0 - frac[-1]) / f0
check(-1e-4 <= loss <= 0.05, "water lost through atmosphere", "%.1f %%" % (100 * loss), "0 .. 5 %")

# ---------------------------------------------------------------- front arrival
# First written time at which the cell next to the right wall on the floor is wet.
# Lower bound: the frictionless dry-bed (Ritter) front speed 2*sqrt(g*H) = 3.38 m/s
# over 0.584-0.146 = 0.438 m gives 0.13 s; the real front cannot be faster.
# Upper bound 0.5 s: plausibility only (no experimental reference for this geometry;
# Martin & Moyce (1952) data are qualitative for this aspect ratio).
# Resolution is the write interval (0.05 s).
check(first_wall is not None and 0.13 <= first_wall <= 0.5, "front reaches right wall [s]",
      "%s" % first_wall, "0.13 .. 0.5 (resolution %g)" % write_interval)

# ---------------------------------------------------------------- report
width = max(len(r[1]) for r in results)
for ok, name, meas, exp in results:
    print("%s  %-*s  measured=%s  expected=%s" % ("PASS" if ok else "FAIL", width, name, meas, exp))
nfail = sum(1 for r in results if not r[0])
print("%d checks, %d passed, %d failed" % (len(results), len(results) - nfail, nfail))
sys.exit(1 if nfail else 0)
