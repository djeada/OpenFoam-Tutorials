#!/usr/bin/env python3
"""Check the results of a finished backward-facing-step run.

Usage:  python3 verify.py [case_dir]

Only reads files already written by ./Allrun (log.*, constant/polyMesh, the
latest time directory, system/*). No OpenFOAM, docker or third-party Python
packages are needed. Prints one line per check and exits 0 only if all pass.
"""
import os
import re
import sys

CASE = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

H = 0.0127          # step height [m]            (system/blockMeshDict, convertToMeters)
U_IN = 10.0         # inlet velocity [m/s]       (0/U)
DEPTH = 1.0 * H     # z thickness of the 2D mesh (blockMeshDict: z = -0.5 .. 0.5, in units of H)

results = []


def fail_missing(msg):
    print("ERROR: %s" % msg, file=sys.stderr)
    print("Run ./run_docker.sh (or ./Allrun inside OpenFOAM) first.", file=sys.stderr)
    sys.exit(2)


def read(rel):
    path = os.path.join(CASE, rel)
    if not os.path.isfile(path):
        fail_missing("missing file %s" % rel)
    with open(path, errors="replace") as f:
        return f.read()


def check(name, measured, expected, tol, ok, unit=""):
    results.append(ok)
    print("%s  %-34s measured=%s  expected=%s%s"
          % ("PASS" if ok else "FAIL", name, measured, expected,
             (" +- %s" % tol) if tol is not None else ""), flush=True)


def latest_time():
    times = []
    for d in os.listdir(CASE):
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(CASE, d)) and t > 0:
            times.append((t, d))
    if not times:
        fail_missing("no time directory > 0 found (solver did not write results)")
    return max(times)[1]


def boundary_values(text, patch):
    """Return the list of values of boundaryField/<patch>/value (floats or 3-tuples)."""
    m = re.search(r"^\s*%s\s*\{(.*?)\n\s*\}" % re.escape(patch), text, re.S | re.M)
    if not m:
        fail_missing("patch %s not found in field file" % patch)
    blk = m.group(1)
    m = re.search(r"value\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(\s*\n(.*?)\n\s*\)\s*;", blk, re.S)
    if not m:
        fail_missing("patch %s has no nonuniform value list" % patch)
    kind, n, body = m.group(1), int(m.group(2)), m.group(3).split("\n")
    vals = []
    for line in body:
        line = line.strip()
        if not line:
            continue
        if kind == "vector":
            vals.append(tuple(float(v) for v in line.strip("()").split()))
        else:
            vals.append(float(line))
    if len(vals) != n:
        fail_missing("patch %s: expected %d values, parsed %d" % (patch, n, len(vals)))
    return vals


def patch_face_range(patch):
    text = read("constant/polyMesh/boundary")
    m = re.search(r"^\s*%s\s*\{(.*?)\}" % re.escape(patch), text, re.S | re.M)
    n = int(re.search(r"nFaces\s+(\d+)", m.group(1)).group(1))
    s = int(re.search(r"startFace\s+(\d+)", m.group(1)).group(1))
    return s, n


def patch_face_centres_x(patch):
    """x of the centre of each face of a patch (mean of its vertex x), from polyMesh."""
    pts_txt = read("constant/polyMesh/points")
    i = pts_txt.index("\n(", pts_txt.index("FoamFile"))
    px = [float(m.group(1)) for m in
          re.finditer(r"^\([ \t]*(\S+)[ \t]+\S+[ \t]+\S+[ \t]*\)[ \t]*$", pts_txt[i:], re.M)]
    faces_txt = read("constant/polyMesh/faces")
    start, n = patch_face_range(patch)
    faces = re.findall(r"^\d+\(([\d ]+)\)\s*$", faces_txt, re.M)
    out = []
    for f in faces[start:start + n]:
        ids = [int(v) for v in f.split()]
        out.append(sum(px[k] for k in ids) / len(ids))
    return out


# ---------------------------------------------------------------- logs
log_mesh = read("log.blockMesh")
log_chk = read("log.checkMesh")
log_run = read("log.simpleFoam")

# 1. No fatal errors in any log
fatal = [n for n, t in (("blockMesh", log_mesh), ("checkMesh", log_chk), ("simpleFoam", log_run))
         if "FOAM FATAL" in t or "Floating point exception" in t or "Segmentation fault" in t]
check("no_fatal_error", "none" if not fatal else ",".join(fatal), "none", None, not fatal)

# 2. Mesh checks. 27200 = 80*40 + 2*300*40 from blockMeshDict.
m = re.search(r"cells:\s+(\d+)", log_chk)
ncells = int(m.group(1)) if m else -1
check("mesh_cell_count", ncells, 27200, "0", ncells == 27200)
check("checkMesh_Mesh_OK", "yes" if "Mesh OK." in log_chk else "no", "yes", None, "Mesh OK." in log_chk)

# 3. Convergence. The case stops on residualControl; if it ran to endTime it did not converge.
end_time = int(float(re.search(r"^endTime\s+([\d.e+]+);", read("system/controlDict"), re.M).group(1)))
m = re.search(r"SIMPLE solution converged in (\d+) iterations", log_run)
check("simple_converged", int(m.group(1)) if m else "not converged",
      "< %d (endTime)" % end_time, None, bool(m) and int(m.group(1)) < end_time)

# 4. Final initial residuals are below the residualControl values of system/fvSolution.
fvs = read("system/fvSolution")
rc = re.search(r"residualControl\s*\{(.*?)\}", fvs, re.S).group(1)
limits = {k: float(v) for k, v in re.findall(r"(\w+)\s+([\d.e+-]+);", rc)}
last = {}
for fld, res in re.findall(r"Solving for (\w+), Initial residual = ([\d.e+-]+)", log_run):
    last[fld] = float(res)
for fld, lim in limits.items():
    meas = max(last.get(fld + "x", 0), last.get(fld + "y", 0)) if fld == "U" else last.get(fld)
    if meas is None:
        fail_missing("no residual for %s in log.simpleFoam" % fld)
    check("final_residual_%s" % fld, "%.2e" % meas, "< %.0e" % lim, None, meas < lim)

# 5. y+ : wall-resolved mesh. Max y+ over all walls from the yPlus function object.
yp = {}
for patch, mx in re.findall(r"patch (\w+) y\+ : min = \S+, max = (\S+),", log_run):
    yp[patch] = float(mx)
ymax = max(yp.values()) if yp else -1
# Tolerance: y+ < 5 is the usual limit for a low-Re (viscous sublayer) near-wall treatment.
check("max_yplus", "%.2f" % ymax, "< 5", None, 0 < ymax < 5)

# ---------------------------------------------------------------- fields
t = latest_time()

# 6. Mass conservation from phi on inlet/outlet patches (m^3/s for the 1H-deep slab).
phi = read(os.path.join(t, "phi"))
q_in = sum(boundary_values(phi, "inlet"))
q_out = sum(boundary_values(phi, "outlet"))
q_ref = U_IN * H * DEPTH
# Inlet flux is fixed by the BC -> must equal U*H*depth to round-off (1e-6 relative).
check("inlet_flux", "%.6e" % -q_in, "%.6e" % q_ref, "1e-6 rel", abs(-q_in - q_ref) < 1e-6 * q_ref)
# Outlet flux must balance inlet; steady, incompressible, converged to p < 1e-4 -> 1e-3 relative.
imb = abs(q_in + q_out) / q_ref
check("mass_balance", "%.2e" % imb, "0", "1e-3 rel", imb < 1e-3)

# 7. Reattachment length from the sign change of tau_x on bottomWall.
# OpenFOAM's wallShearStress is the stress on the fluid, so tau_x > 0 where the near-wall
# flow is reversed and < 0 where it is attached. Reattachment = last + -> - crossing.
tau = [v[0] for v in boundary_values(read(os.path.join(t, "wallShearStress")), "bottomWall")]
xs = patch_face_centres_x("bottomWall")
pts = sorted(zip(xs, tau))
xr = None
for (x0, t0), (x1, t1) in zip(pts, pts[1:]):
    if t0 > 0 >= t1:
        xr = (x0 - t0 * (x1 - x0) / (t1 - t0)) / H
if xr is None:
    check("reattachment_found", "no sign change", "tau_x + -> -", None, False)
else:
    # Regression check: baseline measured with OpenFOAM 6 on this mesh (9.0 H).
    # +-0.5 H = about 7 cells at x ~ 9H, i.e. allows for different convergence history/compilers
    # but not for a changed model, mesh or inlet condition. This is NOT experimental validation.
    check("reattachment_regression", "%.2f" % xr, "9.0", "0.5", abs(xr - 9.0) <= 0.5)
    # Plausibility: turbulent BFS experiments with expansion ratio 1.1-2 give Xr/H of about
    # 6-8 (e.g. Driver & Seegmiller 6.26); RANS scatter of +-25% -> accept 5 .. 10.
    # No experiment matches this exact geometry/Re, so only a sanity band is possible.
    check("reattachment_plausible", "%.2f" % xr, "5 .. 10", None, 5.0 <= xr <= 10.0)

    # Independent cross-check: last sign change of Ux (reversed -> forward flow)
    # along the line sampled in the first cell row (system/sampleDict, written by Allrun).
    # Same physical event seen through a different field/discretisation; allowed difference
    # 0.2 H = about 3 cells at x ~ 9H (face-centre vs cell-centre location, sampling spacing 0.03 H).
    sample = os.path.join("postProcessing", "sampleDict", t, "nearWall_U.xy")
    ux = []
    for line in read(sample).splitlines():
        p = line.split()
        if len(p) >= 4:
            ux.append((float(p[0]), float(p[1])))
    cross = [(x0 - u0 * (x1 - x0) / (u1 - u0)) / H
             for (x0, u0), (x1, u1) in zip(ux, ux[1:]) if u0 < 0 <= u1]
    xu = cross[-1] if cross else None
    check("reattachment_nearwall_Ux", "%.2f" % xu if xu is not None else "none",
          "%.2f (from tau_x)" % xr, "0.2",
          xu is not None and abs(xu - xr) <= 0.2)

# ---------------------------------------------------------------- summary
npass = sum(results)
print("-" * 60)
print("%d of %d checks passed (latest time: %s)" % (npass, len(results), t))
print("RESULT: %s" % ("PASS" if npass == len(results) else "FAIL"))
sys.exit(0 if npass == len(results) else 1)
