#!/usr/bin/env python3
"""Check the results of a finished backward-facing-step LES run.

Usage:  python3 verify.py [case_dir]

Only reads files already written by ./Allrun (log.*, postProcessing/*,
system/*). No OpenFOAM, docker or third-party Python packages are needed.
Prints one line per check and exits 0 only if all pass.
"""
import cmath
import glob
import math
import os
import re
import sys

CASE = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))

H = 0.01            # step height [m]                (system/blockMeshDict, convertToMeters)
U0 = 7.65           # free-stream speed at the inlet  (0/U, edge of the boundary layer)
T_AVG = 0.08        # start of the averaging window   (system/controlDict, fieldAverage)
XR_DNS = 6.28       # Le, Moin & Kim (1997), DNS, same geometry and Re_h = 5100
XR_EXP = 6.0        # Jovic & Driver (1994), experiment, +-0.15
XR_REF = 3.48       # this 2D LES (OpenFOAM 6, 8 cores): regression value

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


def series(pattern):
    """Rows of numbers from postProcessing files (all restart folders, in time order)."""
    files = sorted(glob.glob(os.path.join(CASE, pattern)),
                   key=lambda p: float(p.split(os.sep)[-2]))
    if not files:
        fail_missing("no file matches %s" % pattern)
    rows = {}
    for path in files:
        for line in open(path):
            if line.strip() and not line.startswith("#"):
                values = [float(v) for v in line.replace("(", " ").replace(")", " ").split()]
                rows[values[0]] = values  # a restart repeats its first time: keep one row
    return [rows[t] for t in sorted(rows)]


# ---------------------------------------------------------------- logs
logs = {name: read("log." + name) for name in
        ("blockMesh", "checkMesh", "pimpleFoam.startup", "pimpleFoam.statistics")}

# 1. No fatal errors in any log
fatal = [n for n, t in logs.items()
         if "FOAM FATAL" in t or "Floating point exception" in t or "Segmentation fault" in t]
check("no_fatal_error", "none" if not fatal else ",".join(fatal), "none", None, not fatal)

# 2. Mesh: 58,822 = 69*70 + 397*70 + 397*66 cells from blockMeshDict.
m = re.search(r"cells:\s+(\d+)", logs["checkMesh"])
ncells = int(m.group(1)) if m else -1
check("mesh_cell_count", ncells, 58822, "0", ncells == 58822)
check("checkMesh_Mesh_OK", "yes" if "Mesh OK." in logs["checkMesh"] else "no", "yes", None,
      "Mesh OK." in logs["checkMesh"])

# 3. Both stages reached their end time.
for stage, end in (("startup", 0.08), ("statistics", 0.4)):
    times = [float(t) for t in re.findall(r"^Time = (\S+)", logs["pimpleFoam." + stage], re.M)]
    last = times[-1] if times else -1
    check("reached_end_%s" % stage, "%.4f s" % last, "%.2f s" % end, "1e-6",
          abs(last - end) < 1e-6 and "\nEnd" in logs["pimpleFoam." + stage])

# 4. Courant number. The time step adapts to maxCo 0.8 (system/controlDict); the first
#    few steps after the impulsive start overshoot while it adapts, so they are skipped.
co = [float(v) for stage in ("startup", "statistics")
      for v in re.findall(r"Courant Number mean: \S+ max: (\S+)", logs["pimpleFoam." + stage])]
co_max = max(co[20:])
check("max_courant_after_start", "%.3f" % co_max, "<= 0.85", None, co_max <= 0.85)

# 5. Mass conservation: inlet and outlet flow rates (surfaceFieldValue, every 20 steps)
#    cancel. The inlet flux is fixed by the profile; any imbalance is solver error.
q_in = series("postProcessing/flowRates/*/surfaceFieldValue.dat")
q_out = series("postProcessing/outletFlow/*/surfaceFieldValue.dat")
out_at = {row[0]: row[1] for row in q_out}
imb = max(abs(row[1] + out_at[row[0]]) / abs(row[1]) for row in q_in if row[0] in out_at)
check("mass_balance_max", "%.1e" % imb, "0", "1e-4 rel", imb < 1e-4)
# The 1/7 power-law profile carries U0 * (5h - delta/8) per metre of span, times the
# 1h depth of the mesh; the table in 0/U is piecewise linear, so allow 0.5 %.
q_ref = U0 * (5 * H - 1.2 * H / 8) * H
check("inlet_flow_rate", "%.4e m3/s" % -q_in[-1][1], "%.4e" % q_ref, "0.5 %",
      abs(-q_in[-1][1] - q_ref) < 0.005 * q_ref)

# ---------------------------------------------------------------- mean flow
# 6. Mean reattachment length: where the time-averaged velocity in the first cell row
#    above the floor (UMean_x, system/sampleDict) turns from backwards to forwards for
#    the last time. Near the step a small corner eddy also has forward flow, so the
#    last - -> + crossing is the reattachment point of the main bubble.
latest = sorted(glob.glob(os.path.join(CASE, "postProcessing/sampleDict/*/floor_UMean.xy")),
                key=lambda p: float(p.split(os.sep)[-2]))
if not latest:
    fail_missing("postProcessing/sampleDict/<time>/floor_UMean.xy")
ux = [(float(p[0]), float(p[1])) for p in
      (line.split() for line in open(latest[-1])) if len(p) >= 4]
cross = [(x0 - u0 * (x1 - x0) / (u1 - u0)) / H
         for (x0, u0), (x1, u1) in zip(ux, ux[1:]) if u0 < 0 <= u1]
xr = cross[-1] if cross else None
if xr is None:
    check("reattachment_found", "no sign change", "UMean_x - -> +", None, False)
else:
    # Regression value of this 2D setup, NOT validation. The DNS of the same geometry
    # gives 6.28 h: without a third dimension the shed vortices stay large and
    # coherent and sweep the bubble out early, so a 2D run is expected to come out
    # short (see the README). +-10 %: run-to-run scatter of a 0.32 s average.
    check("reattachment_regression", "%.2f h" % xr, "%.2f h (DNS: %.2f h)" % (XR_REF, XR_DNS),
          "10 %", abs(xr - XR_REF) <= 0.1 * XR_REF)
    # A small counter-rotating corner eddy sits in the step corner: mean flow next to the
    # floor runs forwards there, then backwards under the main bubble. The main bubble
    # starts at the last + -> - crossing before Xr.
    starts = [(x0 - u0 * (x1 - x0) / (u1 - u0)) / H
              for (x0, u0), (x1, u1) in zip(ux, ux[1:]) if u0 > 0 >= u1 and x1 < xr * H]
    xc = starts[-1] if starts else None
    forward = xc is not None and max(u for x, u in ux if x < xc * H) > 0.01 * U0
    check("corner_eddy", "%.2f h" % xc if xc is not None else "none",
          "0.2 h .. Xr, forward flow before it", None, forward and 0.2 <= xc < xr)

# ---------------------------------------------------------------- unsteadiness
# 7. Vortex shedding: the vertical velocity at x = 4h on the step-lip line (probe 1)
#    must oscillate with a clear dominant frequency. Discrete Fourier transform of the
#    statistics window (t >= 0.08 s), resampled at 0.1 ms.
rows = [r for r in series("postProcessing/probes/*/U") if r[0] >= T_AVG]
if len(rows) < 100:
    fail_missing("too few probe samples after t = %g s" % T_AVG)
dt = 1e-4
t_end = rows[-1][0]
samples, i = [], 0
t = T_AVG
while t <= t_end:
    while rows[i + 1][0] < t:
        i += 1
    (t0, *a), (t1, *b) = rows[i], rows[i + 1]
    w = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
    samples.append(a[4] + w * (b[4] - a[4]))  # probe 1 = columns 3..5, Uy = index 4
    t += dt
mean = sum(samples) / len(samples)
samples = [s - mean for s in samples]
n = len(samples)
spectrum = []
for k in range(1, n // 2):
    f = k / (n * dt)
    if f > 500:
        break
    c = sum(s * cmath.exp(-2j * math.pi * k * j / n) for j, s in enumerate(samples[::1]))
    spectrum.append((f, abs(c) ** 2))
f_peak, p_peak = max(spectrum, key=lambda fp: fp[1])
background = sorted(p for _, p in spectrum)[len(spectrum) // 2]
st_h = f_peak * H / U0
check("shedding_peak_clear", "peak/median = %.0f" % (p_peak / background), "> 20", None,
      p_peak / background > 20)
# Strouhal number based on the step height, f h / U0: regression value of this
# setup (one vortex about every 23 ms).
check("shedding_strouhal_h", "%.3f (f = %.1f Hz)" % (st_h, f_peak), "0.056", "0.01",
      abs(st_h - 0.056) <= 0.01)
# ---------------------------------------------------------------- summary
npass = sum(results)
print("-" * 60)
print("%d of %d checks passed" % (npass, len(results)))
print("RESULT: %s" % ("PASS" if npass == len(results) else "FAIL"))
sys.exit(0 if npass == len(results) else 1)
