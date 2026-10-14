#!/usr/bin/env python3
"""Check the results of a finished run of the wind turbine MRF case.

Usage:  python3 verify.py [case_directory]

Reads only files already written by ./Allrun (log.*, the latest time
directory, constant/polyMesh/cellZones, the case dictionaries and
README.md). It starts no solver. Standard library only.

Exit status: 0 if all checks pass, 1 if any check fails, 2 if results are
missing.

Scope: the case has no blades or walls inside the MRF zone, so there is no
torque or power to verify. The checks cover mesh, zone, convergence, mass
balance and the size and sign of the rotation-induced flow asymmetry.

The numerical bands (velocity range, pressure range) are measured values
from one reference run with OpenFOAM 6 plus an allowance (+-10 % velocity,
+-15 % pressure) for other versions, compilers and linear-solver details.
They are not derived from run-to-run scatter, which was not measured.
"""
import glob
import math
import os
import re
import sys

CASE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__))
results = []


def fail_missing(msg):
    print("ERROR: %s\nRun ./Allrun (or ./run_docker.sh) first." % msg, file=sys.stderr)
    sys.exit(2)


def path(*p):
    return os.path.join(CASE, *p)


def read(*p):
    f = path(*p)
    if not os.path.isfile(f):
        fail_missing("missing file %s" % os.path.relpath(f, CASE))
    with open(f, errors="replace") as fh:
        return fh.read()


def strip_comments(s):
    s = re.sub(r"/\*.*?\*/", "", s, flags=re.S)
    return re.sub(r"//[^\n]*", "", s)


def check(name, ok, measured, expected, tol=""):
    results.append(ok)
    print("%s  %-34s measured=%s  expected=%s%s"
          % ("PASS" if ok else "FAIL", name, measured, expected,
             (" +- " + tol) if tol else ""))


def within(name, val, exp, tol, fmt="%.4g", rel=False):
    t = tol * abs(exp) if rel else tol
    check(name, abs(val - exp) <= t, fmt % val, fmt % exp,
          (fmt % tol + (" (relative)" if rel else "")))


def parse_list(text, vec):
    """Parse 'N ( ... )' after a List<...> keyword; return list of floats/tuples."""
    m = re.search(r"List<\w+>\s*(\d+)\s*\(", text)
    if not m:
        return None
    n = int(m.group(1))
    body = text[m.end():]
    if vec:
        vals = re.findall(r"\(\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*\)", body)[:n]
        return [tuple(float(c) for c in v) for v in vals]
    vals = re.findall(r"^\s*([-+\d.eE]+)\s*$", body, re.M)[:n]
    return [float(v) for v in vals]


# ---------------------------------------------------------------- logs
logs = {n: read("log." + n) for n in ("blockMesh", "topoSet", "checkMesh", "simpleFoam")}

fatal = [n for n, t in logs.items() if "FOAM FATAL" in t]
check("no FATAL ERROR in logs", not fatal, "none" if not fatal else ",".join(fatal), "none")

# ---------------------------------------------------------------- mesh
m = re.search(r"cells:\s+(\d+)", logs["checkMesh"])
ncells = int(m.group(1)) if m else -1
bm = strip_comments(read("system", "blockMeshDict"))
nb = re.search(r"hex\s*\([^)]*\)\s*\((\d+)\s+(\d+)\s+(\d+)\)", bm)
nx, ny, nz = (int(v) for v in nb.groups())
check("mesh cell count (blockMeshDict)", ncells == nx * ny * nz, ncells, nx * ny * nz, "0")

check("checkMesh reports Mesh OK", "Mesh OK" in logs["checkMesh"],
      "Mesh OK" if "Mesh OK" in logs["checkMesh"] else "not found", "Mesh OK")

m = re.search(r"cell zones:\s+(\d+)", logs["checkMesh"])
check("number of cell zones", m is not None and int(m.group(1)) == 1,
      m.group(1) if m else "n/a", 1, "0")

# ---------------------------------------------------------------- zone
topo = strip_comments(read("system", "topoSetDict"))
radius = float(re.search(r"radius\s+([-+\d.eE]+);", topo).group(1))
m = re.search(r"cellZoneSet rotatingZone now size (\d+)", logs["topoSet"])
zone_n = int(m.group(1)) if m else -1
verts = re.findall(r"\(\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*\)", bm.split("blocks")[0])
lx = max(float(v[0]) for v in verts)
ly = max(float(v[1]) for v in verts)
dx, dy = lx / nx, ly / ny
expected_zone = math.pi * radius ** 2 / (dx * dy)
# A lattice-point count in a circle differs from pi r^2 / (dx dy) by far
# less than the perimeter cell count (2 pi r / dx ~ 94); 5 % (~35 cells)
# is a deliberately loose allowance.
within("cells in rotatingZone (topoSet)", zone_n, expected_zone, 0.05, "%.4g", rel=True)

zones = read("constant", "polyMesh", "cellZones")
zm = re.search(r"rotatingZone\s*\{.*?cellLabels\s+List<label>\s*(\d+)\s*\(([^)]*)\)", zones, re.S)
labels = [int(v) for v in zm.group(2).split()] if zm else []
check("cellZones file matches topoSet log", zm is not None and len(labels) == zone_n,
      len(labels), zone_n, "0")

check("MRF zone created by simpleFoam",
      re.search(r"creating MRF zone: MRF1", logs["simpleFoam"]) is not None,
      "MRF1" if "creating MRF zone: MRF1" in logs["simpleFoam"] else "not found",
      "MRF1")

# ---------------------------------------------------------------- omega / TSR vs README
mrf = strip_comments(read("constant", "MRFProperties"))
omega = float(re.search(r"\bomega\s+([-+\d.eE]+);", mrf).group(1))
origin = [float(v) for v in re.search(r"origin\s*\(([^)]*)\)", mrf).group(1).split()]
u0 = read("0", "U")
uin = float(re.search(r"inlet\s*\{[^}]*value\s+uniform\s*\(\s*([-+\d.eE]+)", u0).group(1))
readme = read("README.md")
rm = re.search(r"(\d+(?:\.\d+)?) RPM = ([\d.]+) rad/s", readme)
rpm_readme, omega_readme = float(rm.group(1)), float(rm.group(2))
# MRFProperties stores omega rounded to 5 significant digits, hence 1e-3.
within("omega vs README RPM [rad/s]", abs(omega), rpm_readme * 2 * math.pi / 60, 1e-3, "%.5g")
within("omega vs README rad/s value", abs(omega), omega_readme, 1e-3, "%.5g")
tsr = abs(omega) * radius / uin
tm = re.search(r"tip-speed\s+ratio\s+of\s+(\d+\.\d+)", readme)
# README quotes the TSR to two decimals, so +-0.005 is the rounding limit.
within("tip-speed ratio vs README", tsr, float(tm.group(1)), 0.005, "%.4f")
ox, oy = origin[0], origin[1]
cyl = re.search(r"p1\s*\(\s*([-+\d.eE]+)\s+([-+\d.eE]+)", topo)
within("MRF origin matches topoSet axis",
       math.hypot(ox - float(cyl.group(1)), oy - float(cyl.group(2))), 0.0, 1e-9, "%.3g")

# ---------------------------------------------------------------- convergence
cm = re.search(r"SIMPLE solution converged in (\d+) iterations", logs["simpleFoam"])
cd = strip_comments(read("system", "controlDict"))
end_time = float(re.search(r"\bendTime\s+([-+\d.eE]+);", cd).group(1))
it = int(cm.group(1)) if cm else -1
check("SIMPLE converged before endTime", cm is not None and it < end_time,
      it if cm else "not converged", "< %g" % end_time)

fs = strip_comments(read("system", "fvSolution"))
rc = re.search(r"residualControl\s*\{([^}]*)\}", fs).group(1)
tol_p = float(re.search(r"\bp\s+([-+\d.eE]+);", rc).group(1))
last = {}
for fld, val in re.findall(r"Solving for (\w+), Initial residual = ([-+\d.eE]+)", logs["simpleFoam"]):
    last[fld] = float(val)
worst = max(last[f] for f in ("Ux", "Uy", "p", "k", "epsilon"))
# SIMPLE stops when every initial residual is below residualControl (all
# equal to tol_p here), so the last-iteration values must be below it.
check("final initial residuals < residualControl", worst < tol_p, "%.3g" % worst,
      "< %g" % tol_p)

times = [d for d in os.listdir(CASE) if re.fullmatch(r"\d+(\.\d+)?", d) and d != "0"]
if not times:
    fail_missing("no time directory found")
latest = max(times, key=float)
check("latest time dir = convergence iteration", cm is not None and float(latest) == it,
      latest, it, "0")

# ---------------------------------------------------------------- fields
U = parse_list(read(latest, "U").split("internalField", 1)[1], True)
p = parse_list(read(latest, "p").split("internalField", 1)[1], False)
if U is None or p is None or len(U) != ncells or len(p) != ncells:
    fail_missing("could not read U/p fields of time %s" % latest)

speed = [math.sqrt(a * a + b * b + c * c) for a, b, c in U]
# Measured (OpenFOAM 6): min 7.93, max 11.80 m/s. Band +-10 %.
within("minimum |U| [m/s]", min(speed), 7.93, 0.10, "%.3g", rel=True)
within("maximum |U| [m/s]", max(speed), 11.80, 0.10, "%.3g", rel=True)
# Measured: p in [-26.6, +28.6] m^2/s^2. Band +-15 %.
within("minimum p [m^2/s^2]", min(p), -26.6, 0.15, "%.3g", rel=True)
within("maximum p [m^2/s^2]", max(p), 28.6, 0.15, "%.3g", rel=True)

# Mass balance from the face flux phi (kinematic: m^3/s for the 0.1 m depth).
phi = read(latest, "phi")


def patch_flux(patch):
    blk = re.search(r"\b%s\s*\{(.*?)\n    \}" % patch, phi.split("boundaryField", 1)[1], re.S).group(1)
    vals = parse_list(blk, False)
    if vals is None:
        v = re.search(r"value\s+uniform\s+([-+\d.eE]+)", blk)
        vals = [float(v.group(1))] * (int(ny))
    return sum(vals)


qin, qout = patch_flux("inlet"), patch_flux("outlet")
depth = max(float(v[2]) for v in verts)
within("inlet flux = -U*H*depth [m^3/s]", qin, -uin * ly * depth, 1e-6, "%.6g")
# Outflow equals inflow to the extent the pressure equation is converged;
# 0.5 % of the inflow is a loose limit for a 1e-4 residual target.
within("outlet flux + inlet flux [m^3/s]", qin + qout, 0.0, 0.005 * abs(qin), "%.3g")

# Rotation-induced asymmetry: mean p over the upper half of the zone minus
# the lower half (cells ordered i fastest, then j, for the single block).
up, lo = [], []
for c in labels:
    j = c // nx
    y = (j + 0.5) * dy
    (up if y > oy else lo).append(p[c])
dp = sum(up) / len(up) - sum(lo) / len(lo)
# Sign: measured upper - lower < 0 for omega > 0 (about +z); reversed when
# omega is negated (separate run, see README).
sign_ok = dp * omega < 0
scale = 2 * abs(omega) * uin * 2 * radius  # Coriolis estimate 2*omega*U*D
check("zone p asymmetry sign (upper-lower)", sign_ok, "%.3g" % dp,
      "opposite sign to omega")
# Order-of-magnitude check against the Coriolis scale 2*omega*U*D = 126
# m^2/s^2: measured |dp| is about 0.11 of it. Pass between 0.05 and 1.0.
check("zone p asymmetry magnitude", 0.05 * scale <= abs(dp) <= 1.0 * scale,
      "%.3g" % abs(dp), "[%.3g, %.3g]" % (0.05 * scale, scale))

n_pass = sum(results)
print("\n%d/%d checks passed" % (n_pass, len(results)))
print("Not verified: torque and power (the zone contains no blades or walls).")
sys.exit(0 if n_pass == len(results) else 1)
