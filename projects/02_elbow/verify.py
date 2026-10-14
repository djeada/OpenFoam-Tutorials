#!/usr/bin/env python3
"""Check the results of a finished elbow run (log.*, time directories, mesh).

Only reads files already produced by ./Allrun (or ./run_docker.sh); it does not
run OpenFOAM. Prints one line per check and exits 0 only if all checks pass.

There is no analytical or published reference solution for this case. The
checks are therefore (a) run health, (b) conservation and consistency, and
(c) regression values measured with OpenFOAM 6 on this mesh. Each tolerance is
justified in a comment next to the check.

Usage:  python3 verify.py [case_directory]
"""

import math
import os
import re
import sys

case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
results = []


def path(*parts):
    return os.path.join(case, *parts)


def fail_missing(what):
    print(f"ERROR: {what} not found in {case}. Run ./run_docker.sh (or ./Allrun) first.")
    sys.exit(2)


def read(*parts):
    p = path(*parts)
    if not os.path.isfile(p):
        fail_missing(os.path.relpath(p, case))
    with open(p, errors="replace") as f:
        return f.read()


def check_num(name, measured, expected, tol):
    ok = abs(measured - expected) <= tol
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}  measured={measured:.6g}  expected={expected:.6g} ± {tol:.3g}")


def check_max(name, measured, limit):
    ok = measured < limit
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}  measured={measured:.6g}  expected=< {limit:.6g}")


def check_bool(name, measured):
    results.append(bool(measured))
    print(f"{'PASS' if measured else 'FAIL'}  {name}  measured={measured}  expected=True")


def dict_value(text, key):
    m = re.search(r"^\s*" + key + r"\s+([^;\s]+)\s*;", text, re.M)
    if not m:
        print(f"ERROR: entry '{key}' not found")
        sys.exit(2)
    return float(m.group(1))


# ── Case settings (used to derive expected values, not hard-coded) ─────
control = read("system", "controlDict")
end_time = dict_value(control, "endTime")
delta_t = dict_value(control, "deltaT")
write_interval = dict_value(control, "writeInterval")
n_steps = round(end_time / delta_t)
n_written = round(n_steps / write_interval) + 1  # includes time 0


def time_name(t):
    return f"{t:g}"


# ── Logs ────────────────────────────────────────────────────────────────
logs = {n: read("log." + n) for n in
        ("fluentMeshToFoam", "checkMesh", "icoFoam", "foamMeshToFluent", "foamDataToFluent")}

# A FOAM FATAL (IO) ERROR in any utility means the run is invalid. Zero is the
# only acceptable count, so no tolerance applies.
n_fatal = sum(len(re.findall(r"FOAM FATAL", t)) for t in logs.values())
check_num("fatal errors in all logs", n_fatal, 0, 0)

check_bool("icoFoam reached End", logs["icoFoam"].rstrip().endswith("End"))

# The number of time steps follows exactly from endTime / deltaT.
steps = len(re.findall(r"^Time = ", logs["icoFoam"], re.M))
check_num("time steps", steps, n_steps, 0)

time_dirs = [d for d in os.listdir(case)
             if os.path.isdir(path(d)) and re.fullmatch(r"[0-9.]+", d)]
check_num("written time directories", len(time_dirs), n_written, 0)

# ── Mesh quality (checkMesh) ────────────────────────────────────────────
cm = logs["checkMesh"]
check_bool("checkMesh: Mesh OK", "Mesh OK." in cm)
m = re.search(r"cells:\s+(\d+)", cm)
check_num("mesh cells", int(m.group(1)) if m else -1, 918, 0)  # elbow.msh has 918 cells
m = re.search(r"non-orthogonality Max:\s*([\d.eE+-]+)", cm)
# OpenFOAM flags severe non-orthogonality at 70 degrees.
check_max("max non-orthogonality (deg)", float(m.group(1)) if m else 1e9, 70.0)
m = re.search(r"Max skewness = ([\d.eE+-]+)", cm)
# checkMesh fails a mesh at skewness 4 (default).
check_max("max skewness", float(m.group(1)) if m else 1e9, 4.0)

# ── Solver stability ────────────────────────────────────────────────────
courant = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", logs["icoFoam"])]
# Rule of thumb for a PISO solver: Co < 1 (see notes/08_cfl_number.md).
check_max("max Courant number", max(courant) if courant else 1e9, 1.0)

cont = re.findall(r"cumulative = (\S+)", logs["icoFoam"])
# The cumulative continuity error is the sum of the time-step errors (mass
# imbalance of the corrected flux). A well-converged PISO run stays orders of
# magnitude below this limit (measured about 2e-7).
check_max("|cumulative continuity error|", abs(float(cont[-1])) if cont else 1e9, 1e-5)

# ── Helpers to read OpenFOAM ASCII fields ───────────────────────────────
def internal_field(t, name, vector):
    text = read(time_name(t), name)
    if "nan" in text.lower().split():
        return None
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*\n(\d+)\s*\n\(", text)
    if not m:
        print(f"ERROR: could not parse internalField of {name} at t={t}")
        sys.exit(2)
    n = int(m.group(1))
    body = text[m.end():]
    if vector:
        vals = re.findall(r"\(([^()\s]+) ([^()\s]+) ([^()\s]+)\)", body)[:n]
        return [tuple(float(c) for c in v) for v in vals]
    return [float(x) for x in body.split()[:n]]


def boundary_block(text, patch):
    b = text[text.index("boundaryField"):]
    i = b.index(patch)
    return b[i:b.index("}", i)]


def boundary_list(t, name, patch):
    block = boundary_block(read(time_name(t), name), patch)
    m = re.search(r"nonuniform List<\w+>\s*\d+\s*\(([^)]*)\)", block)
    if m:
        return [float(x) for x in m.group(1).split()]
    m = re.search(r"uniform\s+(\S+);", block)  # patch with uniform value
    return m.group(1) if m else None


t_end = end_time
U = internal_field(t_end, "U", True)
p = internal_field(t_end, "p", False)
check_bool("U and p contain no NaN", U is not None and p is not None)
umag = [math.sqrt(a * a + b * b + c * c) for a, b, c in U]

# ── Mass conservation from the boundary flux phi ────────────────────────
# Patch face count from the boundary file, to expand uniform phi values.
boundary = read("constant", "polyMesh", "boundary")


def n_faces(patch):
    i = boundary.index(patch)
    return int(re.search(r"nFaces\s+(\d+)", boundary[i:]).group(1))


def flux(patch):
    v = boundary_list(t_end, "phi", patch)
    if isinstance(v, list):
        return sum(v)
    return float(v) * n_faces(patch)


q5, q6, qo = flux("velocity-inlet-5"), flux("velocity-inlet-6"), flux("pressure-outlet-7")
# Inflow equals outflow to the pressure-solver tolerance. The relative limit
# 1e-4 is two orders above the final p tolerance (1e-6) scaled by the flow rate.
check_num("mass balance (outflow + inflows) / outflow", (qo + q5 + q6) / qo, 0.0, 1e-4)

# Inlet flux equals prescribed velocity times patch area (boundary flux of a
# fixedValue patch is exact up to round-off, so 1e-4 relative is generous).
def patch_area(patch):
    pts_txt = read("constant", "polyMesh", "points")
    fcs_txt = read("constant", "polyMesh", "faces")
    pts = [tuple(map(float, m)) for m in
           re.findall(r"\(([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)\)", pts_txt)]
    faces = [list(map(int, m.split())) for m in re.findall(r"\d+\(([\d ]+)\)", fcs_txt)]
    i = boundary.index(patch)
    start = int(re.search(r"startFace\s+(\d+)", boundary[i:]).group(1))
    area = 0.0
    for f in faces[start:start + n_faces(patch)]:
        # polygon area vector by fan triangulation
        o = pts[f[0]]
        sx = sy = sz = 0.0
        for a, b in zip(f[1:-1], f[2:]):
            u = [pts[a][k] - o[k] for k in range(3)]
            v = [pts[b][k] - o[k] for k in range(3)]
            sx += u[1] * v[2] - u[2] * v[1]
            sy += u[2] * v[0] - u[0] * v[2]
            sz += u[0] * v[1] - u[1] * v[0]
        area += 0.5 * math.sqrt(sx * sx + sy * sy + sz * sz)
    return area


u_text = read("0", "U")


def inlet_speed(patch):
    block = boundary_block(u_text, patch)
    m = re.search(r"value\s+uniform\s+\(([^)]*)\)", block)
    return math.sqrt(sum(float(x) ** 2 for x in m.group(1).split()))


for patch, q in (("velocity-inlet-5", q5), ("velocity-inlet-6", q6)):
    expected = inlet_speed(patch) * patch_area(patch)
    check_num(f"inflow {patch} / (U x area)", -q / expected, 1.0, 1e-4)

# ── Regression values (measured with OpenFOAM 6, this mesh and schemes) ─
# No reference solution exists. Tolerance 0.05 m/s is about 1 % of the peak:
# it allows round-off or compiler differences but detects a changed scheme,
# time step, boundary condition or mesh.
check_num("max |U| at end time (m/s)", max(umag), 4.067, 0.05)
check_num("min p at end time (m2/s2)", min(p), -6.165, 0.1)
check_num("max p at end time (m2/s2)", max(p), 1.236, 0.1)

# Plateau check: the peak speed changes by less than 0.25 % of its value
# between 15 s before the end and the end (measured change 0.007 m/s).
t_prev = end_time - 15
if os.path.isdir(path(time_name(t_prev))):
    U_prev = internal_field(t_prev, "U", True)
    peak_prev = max(math.sqrt(a * a + b * b + c * c) for a, b, c in U_prev)
    check_num(f"peak |U| change t={time_name(t_prev)} -> {time_name(t_end)} (m/s)",
              max(umag) - peak_prev, 0.0, 0.01)
else:
    check_bool(f"time directory {time_name(t_prev)} exists", False)

# ── Fluent export ───────────────────────────────────────────────────────
check_bool("fluentInterface/case.msh written", os.path.isfile(path("fluentInterface", "case.msh")))
check_bool(f"fluentInterface/case{time_name(t_end)}.dat written",
           os.path.isfile(path("fluentInterface", f"case{time_name(t_end)}.dat")))

# ── Summary ─────────────────────────────────────────────────────────────
n_pass = sum(results)
print(f"\n{n_pass}/{len(results)} checks passed")
if n_pass != len(results):
    print("RESULT: FAIL")
    sys.exit(1)
print("RESULT: PASS")
