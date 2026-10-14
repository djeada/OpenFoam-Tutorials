#!/usr/bin/env python3
"""Verify a finished lid-driven cavity run (read-only, standard library only).

Usage: python3 verify.py [case_dir]      (default: directory of this script)

Reads log.blockMesh, log.checkMesh, log.icoFoam, system/ and the time
directories written by ./Allrun. Prints one PASS/FAIL line per check and
exits 0 only if every check passes (2 if results are missing).
"""
import math
import os
import re
import sys

case = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__))
results = []


def fail_missing(msg):
    print("ERROR: %s\nRun ./Allrun (or ./run_docker.sh) first." % msg, file=sys.stderr)
    sys.exit(2)


def read(rel):
    path = os.path.join(case, rel)
    if not os.path.isfile(path):
        fail_missing("missing %s" % rel)
    with open(path, errors="replace") as f:
        return f.read()


def check(name, ok, measured, expected):
    results.append(ok)
    print("%s  %-34s measured=%s  expected=%s" % ("PASS" if ok else "FAIL", name, measured, expected))


def foam_value(text, key):
    m = re.search(r"^\s*%s\s+([^;\s]+)\s*;" % key, text, re.M)
    return float(m.group(1)) if m else None


def read_vector_field(path):
    """Return the internalField of a volVectorField as a list of (x, y)."""
    with open(path, errors="replace") as f:
        txt = f.read()
    body = txt.split("internalField", 1)[1].split("boundaryField", 1)[0]
    if body.lstrip().startswith("uniform"):
        fail_missing("%s has a uniform internalField" % path)
    vals = re.findall(r"\(\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*\)", body)
    return [(float(a), float(b)) for a, b, _ in vals]


def read_scalar_field(path):
    with open(path, errors="replace") as f:
        txt = f.read()
    body = txt.split("internalField", 1)[1].split("boundaryField", 1)[0]
    m = re.search(r"\n\s*\d+\s*\n\(\n(.*?)\n\)", body, re.S)
    if not m:
        fail_missing("cannot parse %s" % path)
    return [float(v) for v in m.group(1).split()]


# ---------------------------------------------------------------- inputs
log_mesh = read("log.blockMesh")
log_check = read("log.checkMesh")
log_ico = read("log.icoFoam")
control = read("system/controlDict")
block = read("system/blockMeshDict")

end_time = foam_value(control, "endTime")
dt = foam_value(control, "deltaT")
m = re.search(r"hex\s*\([\d\s]+\)\s*\(\s*(\d+)\s+(\d+)\s+1\s*\)", block)
nx, ny = (int(m.group(1)), int(m.group(2))) if m else (None, None)
if end_time is None or dt is None or nx is None:
    fail_missing("cannot read endTime/deltaT/cell counts from system/")
n_steps = round(end_time / dt)

times = []
for d in os.listdir(case):
    try:
        t = float(d)
    except ValueError:
        continue
    if os.path.isfile(os.path.join(case, d, "U")) and t > 0:
        times.append((t, d))
if not times:
    fail_missing("no time directories with a U field")
times.sort()
t_last, d_last = times[-1]

# ---------------------------------------------------------------- mesh
# Exact: blockMesh must create nx*ny*1 cells.
m = re.search(r"nCells:\s*(\d+)", log_mesh)
check("blockMesh cell count", m is not None and int(m.group(1)) == nx * ny,
      m.group(1) if m else "n/a", "%d (exact)" % (nx * ny))

# checkMesh verdict is exact.
check("checkMesh verdict", re.search(r"^Mesh OK\.", log_check, re.M) is not None,
      "Mesh OK" if "Mesh OK." in log_check else "not OK", "Mesh OK")

# The block is an axis-aligned box: non-orthogonality must be 0.
# Tolerance 1e-6 deg covers floating-point noise only.
m = re.search(r"non-orthogonality Max:\s*([-\d.eE+]+)", log_check)
v = float(m.group(1)) if m else float("nan")
check("max non-orthogonality [deg]", v <= 1e-6, "%g" % v, "0 +- 1e-6")

# ---------------------------------------------------------------- solver run
bad = re.search(r"FOAM FATAL|FOAM exiting|Floating point exception", log_ico)
finished = re.search(r"^End\s*$", log_ico, re.M) is not None
check("icoFoam finished, no fatal error", finished and not bad,
      "finished" if finished and not bad else "aborted", "End, no FOAM FATAL")

steps = len(re.findall(r"^Time = ", log_ico, re.M))
check("number of time steps", steps == n_steps, steps, "%d (endTime/deltaT, exact)" % n_steps)
check("last time directory", abs(t_last - end_time) < 1e-9, d_last, "%g" % end_time)

# Courant number: Co < 1 is the standard accuracy guideline (README, notes/08).
co = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log_ico)]
co_max = max(co) if co else float("nan")
check("max Courant number", co_max < 1.0, "%.3f" % co_max, "< 1.0")

# Continuity: the pressure equation enforces div(phi)=0 to the p tolerance
# (1e-6 relative to the initial residual). 1e-8 is two orders tighter than
# the tolerance because the residual is already tiny at steady state.
cont = re.findall(r"sum local = (\S+), global = (\S+), cumulative = (\S+?)\s*$", log_ico, re.M)
if cont:
    loc, glob, cum = (float(x.rstrip(",")) for x in cont[-1])
else:
    loc = glob = cum = float("nan")
check("continuity error (last step)", loc < 1e-8, "%.2e" % loc, "< 1e-8")
# A closed cavity must conserve mass globally; only round-off is allowed.
check("cumulative continuity error", abs(cum) < 1e-10, "%.2e" % abs(cum), "< 1e-10")

# Convergence: first-corrector initial residuals of the last time step.
# 1e-3 is the usual 'three orders of magnitude' criterion; the flow is only
# nearly steady at t = 0.5 s (see README), so a much tighter value is not expected.
res = {}
for var, val in re.findall(r"Solving for (Ux|Uy), Initial residual = (\S+),", log_ico):
    res[var] = float(val.rstrip(","))
for var in ("Ux", "Uy"):
    r = res.get(var, float("nan"))
    check("final initial residual %s" % var, r < 1e-3, "%.2e" % r, "< 1e-3")

# ---------------------------------------------------------------- fields
U = read_vector_field(os.path.join(case, d_last, "U"))
p = read_scalar_field(os.path.join(case, d_last, "p"))
check("field sizes", len(U) == nx * ny and len(p) == nx * ny, "%d/%d" % (len(U), len(p)), "%d" % (nx * ny))

# Maximum principle: |U| cannot exceed the lid speed 1 m/s (1e-6 round-off).
umax = max(math.hypot(a, b) for a, b in U)
check("max |U| (cell centres) [m/s]", all(math.isfinite(a + b) for a, b in U) and umax <= 1.0 + 1e-6,
      "%.3f" % umax, "<= 1.0")

# Pressure reference: pRefCell 0 / pRefValue 0 pins p of cell 0 to 0.
# Tolerance 1e-6: PISO re-applies the reference every corrector.
check("p in reference cell 0", abs(p[0]) < 1e-6 and all(math.isfinite(x) for x in p),
      "%.1e" % p[0], "0 +- 1e-6")

# Near-steady check: largest change of |U| over the last 0.1 s (one flow-through
# time L/U). The run is stopped at 0.5 s on purpose, where the field is still
# about 2.5 % of the lid speed away from the steady solution (README), so this is
# only a sanity check that the start-up transient has mostly decayed, not a
# convergence proof. Limit 0.05 m/s = 5 % of the lid speed; the measured value
# (about 0.018) shows the run is nearly, not fully, steady.
ref = [t for t in times if abs(t[0] - (t_last - 0.1)) < 1e-6]
if ref:
    U0 = read_vector_field(os.path.join(case, ref[0][1], "U"))
    dmax = max(math.hypot(a - c, b - d) for (a, b), (c, d) in zip(U, U0))
    check("max change of U over last 0.1 s", dmax < 0.05, "%.4f" % dmax, "< 0.05 m/s")
else:
    print("SKIP  max change of U over last 0.1 s  (time directory %.3f not written)" % (t_last - 0.1))

# ------------------------------------------------- Target: Ghia et al. (1982)
# Re = 100 vertical-centreline Ux (Table I of Ghia, Ghia and Shin).
GHIA = [(1.0000, 1.00000), (0.9766, 0.84123), (0.9688, 0.78871), (0.9609, 0.73722),
        (0.9531, 0.68717), (0.8516, 0.23151), (0.7344, 0.00332), (0.6172, -0.13641),
        (0.5000, -0.20581), (0.4531, -0.21090), (0.2813, -0.15662), (0.1719, -0.10150),
        (0.1016, -0.06434), (0.0703, -0.04775), (0.0625, -0.04192), (0.0547, -0.03717),
        (0.0000, 0.00000)]
# Ux on x/L = 0.5 is the mean of the two cell columns next to the centreline
# (exact linear interpolation on this uniform mesh). Wall values (0 and the
# lid speed 1) are added so that Ghia's points near the walls are bracketed.
if nx % 2 == 0:
    cols = (nx // 2 - 1, nx // 2)
    prof = [(0.0, 0.0)]
    for j in range(ny):
        prof.append(((j + 0.5) / ny, 0.5 * (U[cols[0] + nx * j][0] + U[cols[1] + nx * j][0])))
    prof.append((1.0, 1.0))
    worst, worst_y = 0.0, 0.0
    for y, ref_u in GHIA:
        for (y0, u0), (y1, u1) in zip(prof, prof[1:]):
            if y0 - 1e-12 <= y <= y1 + 1e-12:
                u = u0 + (y - y0) / (y1 - y0) * (u1 - u0)
                if abs(u - ref_u) >= worst:
                    worst, worst_y = abs(u - ref_u), y
                break
    # 0.03 m/s = 3 % of the lid speed: expected O(h^2) discretisation error of
    # a 20 x 20 mesh plus the linear interpolation between cell centres. The
    # measured value (about 0.02) is the case's own accuracy, not a fitted limit;
    # a 40 x 40 mesh (exercise 2) should be clearly lower.
    check("max |Ux - Ghia| on centreline", worst < 0.03, "%.4f (y/L=%.4f)" % (worst, worst_y), "< 0.03 m/s")

# Primary vortex centre = extremum of the stream function psi(x, y) built by
# integrating Ux upward through each cell column (psi = 0 on the bottom wall).
# Ghia: (0.6172, 0.7344). Tolerance 0.05 L = one cell width, the location
# resolution of a cell-centred field on this mesh.
psi_best = None
for i in range(nx):
    s = 0.0
    for j in range(ny):
        ux = U[i + nx * j][0]
        psi = s + 0.5 * ux / ny
        s += ux / ny
        if psi_best is None or psi < psi_best[0]:
            psi_best = (psi, (i + 0.5) / nx, (j + 0.5) / ny)
cx, cy = psi_best[1], psi_best[2]
err = max(abs(cx - 0.6172), abs(cy - 0.7344))
check("primary vortex centre (x/L, y/L)", err <= 0.05 + 1e-9, "(%.3f, %.3f)" % (cx, cy),
      "(0.617, 0.734) +- 0.05")

# ---------------------------------------------------------------- summary
n_fail = results.count(False)
print("\n%d checks, %d passed, %d failed -> %s" % (len(results), len(results) - n_fail, n_fail,
                                                   "PASS" if n_fail == 0 else "FAIL"))
sys.exit(0 if n_fail == 0 else 1)
