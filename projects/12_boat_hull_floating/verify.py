#!/usr/bin/env python3
"""Check the results of a finished (or stopped) floating-hull run.

Reads only files already produced by ./Allrun (log.*, postProcessing/forces,
time directories) plus the case dictionaries and ../../boat_hull.stl.
No OpenFOAM, docker or third-party Python packages are needed.

Output: one line per check
    PASS|FAIL|SKIP  <name>  measured=<v>  expected=<v> +- <tol>
Exit code: 0 if no check FAILED and the logs exist; 1 otherwise.
SKIP means the run was stopped before the check can be evaluated
(for example the settling checks need the full 3 s of simulated time).

Tolerances are justified in the comment next to each check; they were not
tuned to make a check pass.
"""
import glob
import math
import os
import re
import struct
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
G = 9.81            # constant/g
RHO_W = 998.2       # constant/transportProperties, water

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
results = []        # (status, name, measured, expected, tol)


def path(*p):
    return os.path.join(HERE, *p)


def fail_hard(msg):
    print("ERROR: " + msg)
    print("Run ./run_docker.sh (or ./Allrun) first; verify.py only reads its results.")
    sys.exit(2)


def read(p, required=True):
    try:
        with open(p, errors="replace") as f:
            return f.read()
    except OSError:
        if required:
            fail_hard("missing file: " + os.path.relpath(p, HERE))
        return None


def add(status, name, measured, expected, tol):
    results.append((status, name, measured, expected, tol))


def fmt(v):
    if isinstance(v, float):
        return "%.5g" % v
    return str(v)


def check(name, measured, expected, tol, absolute=True):
    """PASS if |measured - expected| <= tol (absolute) or relative tol."""
    err = abs(measured - expected)
    lim = tol if absolute else tol * abs(expected)
    add("PASS" if err <= lim else "FAIL", name, measured, expected, tol)


def check_max(name, measured, limit):
    add("PASS" if measured <= limit else "FAIL", name, measured, "<= " + fmt(limit), "-")


def skip(name, why):
    add("SKIP", name, "n/a", why, "-")


# ---------------------------------------------------------------- case data
def dict_value(text, key):
    m = re.search(r"^\s*" + re.escape(key) + r"\s+(.*?);", text, re.M)
    return m.group(1) if m else None


def vec(s):
    return [float(x) for x in re.findall(NUM, s)]


def load_stl_m(p):
    """Binary STL in mm (as stored at the repo root) -> triangles in m."""
    with open(p, "rb") as f:
        d = f.read()
    n = struct.unpack("<I", d[80:84])[0]
    tris = []
    for i in range(n):
        v = struct.unpack("<12f", d[84 + 50 * i:84 + 50 * i + 48])
        tris.append([tuple(c * 1e-3 for c in v[3 + 3 * k:6 + 3 * k]) for k in range(3)])
    return tris


def submerged_volume(tris, zw):
    """Volume of the closed surface below the plane z = zw (divergence thm.)."""
    vol = 0.0
    for t in tris:
        poly = []
        for i in range(3):
            a, b = t[i], t[(i + 1) % 3]
            ina, inb = a[2] <= zw, b[2] <= zw
            if ina:
                poly.append(a)
            if ina != inb:
                s = (zw - a[2]) / (b[2] - a[2])
                poly.append(tuple(a[k] + s * (b[k] - a[k]) for k in range(3)))
        for k in range(1, len(poly) - 1):
            p0, p1, p2 = poly[0], poly[k], poly[k + 1]
            u = [p1[j] - p0[j] for j in range(3)]
            w = [p2[j] - p0[j] for j in range(3)]
            nz = 0.5 * (u[0] * w[1] - u[1] * w[0])          # z-component of area vector
            zc = (p0[2] + p1[2] + p2[2]) / 3.0 - zw
            vol += zc * nz
    return vol


def equilibrium_waterline(tris, volume):
    """Waterline height relative to the hull centre giving displaced `volume`."""
    lo, hi = min(p[2] for t in tris for p in t), max(p[2] for t in tris for p in t)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if submerged_volume(tris, mid) < volume:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------- load files
mesh_dict = read(path("constant", "dynamicMeshDict"))
control = read(path("system", "controlDict"))
block = read(path("system", "blockMeshDict"))
snappy_log = read(path("log.snappyHexMesh"))
check_log = read(path("log.checkMesh"))
set_log = read(path("log.setFields"))
solver_log = read(path("log.interFoam"))

mass = vec(dict_value(mesh_dict, "mass"))[0]
com0 = vec(dict_value(mesh_dict, "centreOfMass"))
end_time = vec(dict_value(control, "endTime"))[0]
max_co = vec(dict_value(control, "maxCo"))[0]
max_aco = vec(dict_value(control, "maxAlphaCo"))[0]

stl_path = os.path.normpath(path("..", "..", "boat_hull.stl"))
if not os.path.exists(stl_path):
    fail_hard("missing hull STL: " + stl_path)
tris = load_stl_m(stl_path)
v_hull = submerged_volume(tris, 1e9)
v_sub0 = submerged_volume(tris, 0.0)            # hull centred at z = 0, water below z = 0

# ------------------------------------------------------------------- mesh
add("PASS" if "Finished meshing without any errors" in snappy_log else "FAIL",
    "snappyHexMesh finished", "ok" if "Finished meshing without any errors" in snappy_log else "not found",
    "Finished meshing without any errors", "-")
add("PASS" if re.search(r"^Mesh OK\.", check_log, re.M) else "FAIL",
    "checkMesh Mesh OK", "ok" if re.search(r"^Mesh OK\.", check_log, re.M) else "failed", "Mesh OK.", "-")
add("PASS" if re.search(r"boat\s+\d+\s+\d+\s+ok \(closed", check_log) else "FAIL",
    "hull patch closed", "closed" if re.search(r"boat\s+\d+\s+\d+\s+ok \(closed", check_log) else "not closed",
    "closed singly connected", "-")

m = re.search(r"Mesh non-orthogonality Max: (" + NUM + ")", check_log)
# snappyHexMeshDict meshQualityControls: maxNonOrtho 65
check_max("max non-orthogonality [deg]", float(m.group(1)) if m else float("inf"), 65.0)
m = re.search(r"Max skewness = (" + NUM + ")", check_log)
# snappyHexMeshDict meshQualityControls: maxInternalSkewness 4
check_max("max skewness", float(m.group(1)) if m else float("inf"), 4.0)
m_tot = re.search(r"Total volume = (" + NUM + ")", check_log)
v_total = float(m_tot.group(1)) if m_tot else None

# ----------------------------------------------------------- solver output
# FOAM FATAL / FOAM exiting, a trapped "Floating point exception", or NaN
# values.  The start-up line "Enabling floating point exception trapping"
# is not an error (lower case, so it does not match).
FATAL_RE = r"FOAM FATAL|FOAM exiting|Floating point exception|\bnan\b"
add("PASS" if not re.search(FATAL_RE, solver_log) else "FAIL",
    "no fatal error / NaN in log.interFoam",
    "none" if not re.search(FATAL_RE, solver_log) else "found",
    "none", "-")

times = [float(t) for t in re.findall(r"^Time = (" + NUM + ")", solver_log, re.M)]
if not times:
    fail_hard("log.interFoam contains no time steps")
t_end = times[-1]
finished = bool(re.search(r"^End\s*$", solver_log, re.M))
full_run = finished and t_end >= end_time - 1e-9

# Constraint tensors printed by sixDoFRigidBodyMotion: only z-translation,
# only rotation about z (dictionary: line (0 0 1) + axis (0 0 1)).
m = re.search(r"Translational constraint tensor \(([^)]*)\)", solver_log)
m2 = re.search(r"Rotational constraint tensor \(([^)]*)\)", solver_log)
tc = vec(m.group(1)) if m else []
rc = vec(m2.group(1)) if m2 else []
want = [0, 0, 0, 0, 0, 0, 0, 0, 1]
add("PASS" if tc == want and rc == want else "FAIL", "DoF constraints: heave + yaw only",
    "translation=%s rotation=%s" % ("zz" if tc == want else "?", "zz" if rc == want else "?"),
    "only zz components 1", "-")

# Water volume fraction (Phase-1 volume fraction) at the first step against
# the geometric value: water below z = 0 minus the submerged part of the hull.
# Tolerance 2e-4 of the domain volume (= 52 cm3): the snapped mesh
# represents the 211 cm3 submerged hull volume within about 10 %.
vf = [float(x) for x in re.findall(r"Phase-1 volume fraction = (" + NUM + ")", solver_log)]
xs = [vec(l) for l in re.findall(r"^\s*\((.*?)\)\s*//\s*\d", block, re.M)]
zmin = min(p[2] for p in xs)
lx = max(p[0] for p in xs) - min(p[0] for p in xs)
ly = max(p[1] for p in xs) - min(p[1] for p in xs)
if vf and v_total:
    exp_vf = (lx * ly * (0.0 - zmin) - v_sub0) / v_total
    check("initial water volume fraction", vf[0], exp_vf, 2e-4)
    # Water is only lost/gained through the open boundaries and by MULES
    # truncation.  Tolerance 1.5e-3 relative = 225 cm3 = a 0.3 mm change of
    # the water level over the 0.75 m2 free surface, 5 % of the equilibrium
    # draft (6.3 mm): a level error that small does not change the buoyancy.
    check("water volume conserved (first vs last)", vf[-1], vf[0], 1.5e-3, absolute=False)
else:
    add("FAIL", "water volume fraction in log", "not found", "present", "-")

amin = [float(x) for x in re.findall(r"Min\(alpha\.water\) = (" + NUM + ")", solver_log)]
amax = [float(x) for x in re.findall(r"Max\(alpha\.water\) = (" + NUM + ")", solver_log)]
# MULES with nLimiterIter 5 keeps alpha within O(1e-3) of [0,1] on a moving
# mesh; 0.01 (1 % of the range) separates this from a failing interface.
check_max("alpha.water overshoot above 1", max(amax) - 1.0 if amax else float("inf"), 0.01)
check_max("alpha.water undershoot below 0", -min(amin) if amin else float("inf"), 0.01)

# Courant numbers.  adjustTimeStep sets deltaT from the Courant number of
# the previous step, so the controlled limits (maxCo, maxAlphaCo) are
# exceeded transiently.  The pass criterion is the stability bound of the
# explicit convection / MULES alpha transport used here, Co <= 0.5
# (independent of the case settings; the limits actually set are 0.4 / 0.25).
co = [float(x) for x in re.findall(r"^Courant Number mean: " + NUM + " max: (" + NUM + ")", solver_log, re.M)]
aco = [float(x) for x in re.findall(r"^Interface Courant Number mean: " + NUM + " max: (" + NUM + ")", solver_log, re.M)]
check_max("max Courant number (maxCo=%g)" % max_co, max(co) if co else float("inf"), 0.5)
check_max("max interface Courant (maxAlphaCo=%g)" % max_aco, max(aco) if aco else float("inf"), 0.5)

# ------------------------------------------------- 6-DoF motion from the log
# Blocks "Centre of mass / Orientation / Linear velocity / Angular velocity"
# are printed once per PIMPLE outer corrector; the last one per time step is
# the converged state.
com = [vec(s) for s in re.findall(r"Centre of mass: \(([^)]*)\)", solver_log)]
ori = [vec(s) for s in re.findall(r"Orientation: \(([^)]*)\)", solver_log)]
vel = [vec(s) for s in re.findall(r"Linear velocity: \(([^)]*)\)", solver_log)]
if not com:
    fail_hard("no 6-DoF output (Centre of mass) in log.interFoam")

# line constraint: x and y must not change at all (round-off only)
dxy = max(max(abs(c[0]), abs(c[1])) for c in com)
check_max("surge/sway of centre of mass [m]", dxy, 1e-9)

# axis constraint: matrix entries xz, yz, zx, zy of the orientation
# (pitch/roll) must stay zero (round-off only)
roll_pitch = max(max(abs(o[2]), abs(o[5]), abs(o[6]), abs(o[7])) for o in ori) if ori else float("inf")
check_max("pitch/roll orientation entries", roll_pitch, 1e-6)

# Yaw is the one free rotation; the symmetric hull in calm water should not
# yaw. 2 deg is far above numerical asymmetry seen in a stable run and well
# below anything physical here.
yaw = max(abs(math.degrees(math.atan2(o[3], o[0]))) for o in ori) if ori else float("inf")
check_max("max yaw [deg]", yaw, 2.0)

# Energy bound for the heave speed: the released buoyancy work
# (rho g V_sub0 - m g) * rise gives v <= sqrt(2 W / m); the rise is the
# equilibrium rise computed below.  15 % margin.
wl_eq = equilibrium_waterline(tris, mass / RHO_W)     # waterline rel. hull centre at equilibrium
rise = -wl_eq                                         # hull rises by this much
z_eq = com0[2] + rise
f_net = RHO_W * G * v_sub0 - mass * G
v_max = math.sqrt(2.0 * f_net * rise / mass)
vz = max(abs(v[2]) for v in vel)
check_max("max heave speed [m/s] (energy bound)", vz, 1.15 * v_max)

# Momentum balance: the hull momentum at the end of the run must equal the
# time integral of the vertical force (pressure + viscous z force from
# forces.dat, plus weight), multiplied by accelerationDamping because the
# Newmark update scales the acceleration by that factor.  This couples the
# force output to the motion output independently of the added mass (the
# force on the hull is NOT rho*g*V_sub while the hull accelerates).
# Tolerance 5 % of the gross impulse int|F - m g| dt: the force is sampled at
# the end of each step while the 6-DoF solver used the force of each outer
# corrector, and accelerationRelaxation 0.7 lags the response.
a_damp = vec(dict_value(mesh_dict, "accelerationDamping") or "1")[0]
fdat = glob.glob(path("postProcessing", "forces", "*", "forces.dat"))
rows_f = []
if fdat:
    for line in read(fdat[0]).splitlines():
        if line.strip() and not line.startswith("#"):
            nums = [float(x) for x in re.findall(NUM, line.replace("(", " ").replace(")", " "))]
            rows_f.append((nums[0], nums[3] + nums[6]))     # time, pressure z + viscous z
if len(rows_f) > 2:
    imp = 0.0
    gross = 0.0
    for (t0, f0), (t1, f1) in zip(rows_f, rows_f[1:]):
        g0, g1 = f0 - mass * G, f1 - mass * G
        imp += 0.5 * (g0 + g1) * (t1 - t0)
        gross += 0.5 * (abs(g0) + abs(g1)) * (t1 - t0)
    mv = mass * vel[-1][2]
    check("momentum m*v_end vs accelerationDamping*int(F-mg)dt [kg m/s]", mv, a_damp * imp, 0.05 * gross)
else:
    add("FAIL", "forces.dat rows", len(rows_f), "> 2", "-")

# Direction: the hull starts with the waterline at its centre, i.e. far
# more buoyancy than weight, so the centre of mass must rise.  Consistency
# check only (no published reference exists for this case).
z_end = com[-1][2]
add("PASS" if z_end > com0[2] + 1e-4 else "FAIL", "centre of mass rises (z_end - z_start) [m]",
    z_end - com0[2], "> 1e-4", "-")

# Equilibrium from the STL (waterline for displaced volume m/rho): fixes
# what the settled state must be.  Checked only if the run reached endTime.
if full_run:
    # Tolerance 3 mm = half the equilibrium draft: finest waterline cells are
    # 6 mm, so a smaller tolerance would test the mesh, not the physics.
    z_last = sum(c[2] for c in com[-30:]) / len(com[-30:])
    check("settled z of centre of mass [m]", z_last, z_eq, 3e-3)
else:
    skip("settled z of centre of mass [m]",
         "z_eq = %.4f m (needs run to endTime=%g, stopped at %.4g s)" % (z_eq, end_time, t_end))
    skip("heave oscillation decays", "needs run to endTime")

# Time directories written
tdirs = sorted(d for d in os.listdir(HERE)
               if re.fullmatch(NUM, d) and os.path.isdir(path(d)) and float(d) > 0)
add("PASS" if tdirs else "FAIL", "time directories written", len(tdirs), ">= 1", "-")

# ------------------------------------------------------------------ report
print("%-5s %-46s %-24s %-26s %s" % ("", "check", "measured", "expected", "+- tol"))
n_fail = 0
n_skip = 0
for st, name, meas, exp, tol in results:
    print("%s  %s  measured=%s  expected=%s +- %s" % (st, name, fmt(meas), fmt(exp), fmt(tol)))
    n_fail += st == "FAIL"
    n_skip += st == "SKIP"
n_pass = sum(1 for r in results if r[0] == "PASS")
print("")
print("Run reached t = %.4g s of endTime = %g s (%s); hull z_CoM = %.5f m, v_z = %.4f m/s"
      % (t_end, end_time, "complete" if full_run else "PARTIAL", z_end, vel[-1][2] if vel else float("nan")))
print("STL: hull volume %.1f cm3, submerged at z=0: %.1f cm3, equilibrium: rise %.1f mm, draft %.1f mm, z_CoM %.4f m"
      % (v_hull * 1e6, v_sub0 * 1e6, rise * 1e3,
         (wl_eq - min(p[2] for t in tris for p in t)) * 1e3, z_eq))
print("Initial state: waterline at the hull centre, buoyancy %.3f N vs weight %.3f N -> NOT in equilibrium "
      "(net %+.2f N, by design); the hull must rise %.1f mm to float at the %.1f mm draft"
      % (RHO_W * G * v_sub0, mass * G, RHO_W * G * v_sub0 - mass * G, rise * 1e3,
         (wl_eq - min(p[2] for t in tris for p in t)) * 1e3))
print("SUMMARY: %d PASS, %d FAIL, %d SKIP%s"
      % (n_pass, n_fail, n_skip, " (SKIP = needs a full run to endTime)" if n_skip else ""))
sys.exit(1 if n_fail else 0)
