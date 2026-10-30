#!/usr/bin/env python3
"""Check the results of a finished vertical-axis wind turbine run.

Usage:  python3 verify.py [case_dir]

Only reads files already written by ./Allrun (log.*, postProcessing/*,
system/*, constant/*). No OpenFOAM, docker or third-party Python packages are
needed. Prints one line per check and exits 0 only if all pass.
"""
import cmath
import glob
import math
import os
import re
import sys

CASE = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

U_INF = 8.0         # wind speed [m/s]                   (0/U)
OMEGA = 40.0        # rotor speed [rad/s]                (constant/dynamicMeshDict)
R = 0.5             # rotor radius [m]                   (generate_vawt_geo.py)
SPAN = 1.0          # mesh depth [m]                     (generate_vawt_geo.py)
RHO = 1.225         # air density used by the forces    (system/controlDict, rhoInf)
N_REV = 10          # revolutions simulated              (Allrun: endTime = 10 T)
N_AVG = 4           # revolutions averaged (7-10)
PERIOD = 2 * math.pi / OMEGA
BETZ = 16.0 / 27.0
# Mean power coefficient of revolutions 7-10 in the reference run (OpenFOAM 6,
# cfdengine/openfoam image, Gmsh 4.12.1, 8 cores); see the README.
CP_REF = 0.131

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


def check(name, measured, expected, tol, ok):
    results.append(ok)
    print("%s  %-30s measured=%s  expected=%s%s"
          % ("PASS" if ok else "FAIL", name, measured, expected,
             (" +- %s" % tol) if tol is not None else ""), flush=True)


def moment_z(name):
    """(time, M_z) of a forces function object, pressure + viscous, all restarts."""
    files = sorted(glob.glob(os.path.join(CASE, "postProcessing", name, "*", "forces.dat")),
                   key=lambda p: float(p.split(os.sep)[-2]))
    if not files:
        fail_missing("postProcessing/%s/<time>/forces.dat" % name)
    rows = {}
    for path in files:
        for line in open(path):
            if line.strip() and not line.startswith("#"):
                v = [float(x) for x in line.replace("(", " ").replace(")", " ").split()]
                rows[v[0]] = v[12] + v[15]  # moment: pressure z + viscous z
    return sorted(rows.items())


def mean_over(series, t0, t1):
    values = [m for t, m in series if t0 < t <= t1 + 1e-9]
    return sum(values) / len(values)


def cp(torque):
    return torque * OMEGA / (0.5 * RHO * U_INF ** 3 * 2 * R * SPAN)


# ---------------------------------------------------------------- logs
logs = {name: read("log." + name) for name in
        ("gmshToFoam", "changeDictionary", "checkMesh", "pimpleFoam.spinup", "pimpleFoam.record")}

# 1. No fatal errors in any log
fatal = [n for n, t in logs.items()
         if "FOAM FATAL" in t or "Floating point exception" in t or "Segmentation fault" in t]
check("no_fatal_error", "none" if not fatal else ",".join(fatal), "none", None, not fatal)

# 2. Mesh. The cell count depends a little on the Gmsh version, so it is a range.
m = re.search(r"cells:\s+(\d+)", logs["checkMesh"])
ncells = int(m.group(1)) if m else -1
check("mesh_cell_count", ncells, "100000 .. 130000", None, 100000 <= ncells <= 130000)
check("checkMesh_Mesh_OK", "yes" if "Mesh OK." in logs["checkMesh"] else "no", "yes", None,
      "Mesh OK." in logs["checkMesh"])
# The rotor and the stator are two separate regions: they touch only through AMI.
m = re.search(r"Number of regions: (\d+)", logs["checkMesh"])
check("two_mesh_regions", m.group(1) if m else "?", "2", None, bool(m) and m.group(1) == "2")

# 3. AMI: every face of each side is fully covered by the other side. sum(weights)
#    is 1 for a perfect overlap; the 10 micrometre gap between the circles and
#    the polygonal faces leave about 1e-5.
weights = [float(w) for pair in re.findall(
    r"sum\(weights\) min/max/average = ([\d.e+-]+), ([\d.e+-]+)",
    logs["pimpleFoam.spinup"] + logs["pimpleFoam.record"]) for w in pair]
worst = max(abs(w - 1) for w in weights) if weights else 1
check("AMI_weights", "max |sum - 1| = %.1e" % worst, "0", "1e-3", bool(weights) and worst < 1e-3)

# 4. Both stages reached their end time (8 and 10 revolutions).
for stage, revs in (("spinup", 8), ("record", N_REV)):
    times = [float(t) for t in re.findall(r"^Time = (\S+)", logs["pimpleFoam." + stage], re.M)]
    last = times[-1] if times else -1
    check("reached_end_%s" % stage, "%.4f s" % last, "%.4f s (%d rev)" % (revs * PERIOD, revs),
          "1e-4", abs(last - revs * PERIOD) < 1e-4 and "\nEnd" in logs["pimpleFoam." + stage])

# ---------------------------------------------------------------- torque and power
# The blades stall dynamically on every upwind pass and the shed vortices are never
# exactly the same twice, so the power changes by about +-5 % from one revolution to the
# next even once the wake has developed. Statistics therefore use the last four
# revolutions (7-10); forces are written every degree in both stages.
rotor = moment_z("rotor")
blades = [moment_z("blade%d" % k) for k in (1, 2, 3)]
t_end = N_REV * PERIOD
t_avg = t_end - N_AVG * PERIOD
cp_mean = cp(mean_over(rotor, t_avg, t_end))

# 5. The rotor extracts power, less than the Betz limit 16/27 of the wind's power.
check("power_coefficient", "%.3f" % cp_mean, "0 < Cp < %.3f (Betz)" % BETZ, None,
      0 < cp_mean < BETZ)
# Same result as the reference run. +-0.015: about twice the scatter of a 4-revolution
# mean, given the +-5 % scatter of single revolutions.
check("power_coefficient_regression", "%.3f" % cp_mean, "%.3f" % CP_REF, "0.015",
      abs(cp_mean - CP_REF) <= 0.015)

# 6. Settled: the mean power of revolutions 9-10 is within 10 % of revolutions 7-8.
#    (Early on it drifts by more: revolutions 3-6 give about 15 % more power than 7-10,
#    while the wake behind the rotor is still building up.)
cp_a = cp(mean_over(rotor, t_end - 4 * PERIOD, t_end - 2 * PERIOD))
cp_b = cp(mean_over(rotor, t_end - 2 * PERIOD, t_end))
change = abs(cp_b - cp_a) / abs(cp_b)
check("settled_last_revs", "Cp %.3f -> %.3f (%.1f %%)" % (cp_a, cp_b, 100 * change),
      "< 10 %", None, change < 0.10)

# 7. Three identical blades do the same work on average. Single revolutions scatter by
#    about +-15 % per blade, so over four revolutions the blades agree to ~15 %.
shares = [mean_over(b, t_avg, t_end) for b in blades]
spread = (max(shares) - min(shares)) / (sum(shares) / 3)
check("blades_share_equally", "spread %.1f %%" % (100 * spread), "< 15 %", None, spread < 0.15)

# 8. The rotor torque pulses three times per revolution (once per blade pass):
#    dominant frequency 3 * omega / (2 pi) = 19.1 Hz over the last four revolutions.
window = [(t, m) for t, m in rotor if t > t_avg]
mean = sum(m for _, m in window) / len(window)
n = len(window)
dt = (window[-1][0] - window[0][0]) / (n - 1)
spectrum = [(k / (n * dt), abs(sum((m - mean) * cmath.exp(-2j * math.pi * k * j / n)
                                   for j, (_, m) in enumerate(window))))
            for k in range(1, n // 4)]
f_peak = max(spectrum, key=lambda fp: fp[1])[0]
f_3p = 3 * OMEGA / (2 * math.pi)
check("torque_pulses_3_per_rev", "%.1f Hz" % f_peak, "%.1f Hz" % f_3p,
      "1 bin (%.1f Hz)" % (1 / (n * dt)), abs(f_peak - f_3p) <= 0.51 / (n * dt) + 1e-9)

# 9. Most of the work is done in the upwind half of the revolution. Blade k starts at
#    theta = 120 (k - 1) deg, where theta = 0 is the top of the rotor (moving into the
#    wind); the rotor turns anticlockwise, so 0 < theta < 180 deg is the upwind half.
up = down = 0.0
for k, blade in enumerate(blades):
    for t, m in blade:
        if t > t_avg:
            theta = (120 * k + math.degrees(OMEGA * t)) % 360
            if 0 < theta < 180:
                up += m
            else:
                down += m
check("upwind_half_dominates", "upwind %.0f %% of blade work" % (100 * up / (up + down)),
      "> 60 %", None, up / (up + down) > 0.6)

# ---------------------------------------------------------------- summary
npass = sum(results)
print("-" * 60)
print("Cp (revolutions 7-10) = %.3f, tip-speed ratio %.1f" % (cp_mean, OMEGA * R / U_INF))
print("%d of %d checks passed" % (npass, len(results)))
print("RESULT: %s" % ("PASS" if npass == len(results) else "FAIL"))
sys.exit(0 if npass == len(results) else 1)
