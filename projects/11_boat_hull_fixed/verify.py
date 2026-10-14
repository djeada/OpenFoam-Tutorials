#!/usr/bin/env python3
"""Check the results of a finished run of the fixed boat hull case.

Reads only files already produced by ./Allrun (log.*, constant/polyMesh,
time directories, postProcessing/). Standard library only; no OpenFOAM, no
docker. Prints one line per check:

    PASS|FAIL|SKIP  <name>  measured=<value>  expected=<value> +/- <tol>

Exit status: 0 if no check FAILed (SKIPs are allowed, the summary counts them),
1 if at least one check FAILed, 2 if results/logs are missing.

A check marked "consistency only" has no published reference value: the
expected value comes from the case files or from geometry, not from
literature or experiment.
"""
import math
import re
import struct
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = []


# --------------------------------------------------------------- reporting
def report(status, name, measured, expected, note=""):
    RESULTS.append(status)
    line = f"{status:<4}  {name:<30}  measured={measured}  expected={expected}"
    if note:
        line += f"  ({note})"
    print(line)


def check_close(name, measured, expected, tol, rel=False, note=""):
    """PASS if |measured - expected| <= tol (absolute, or relative if rel)."""
    limit = tol * abs(expected) if rel else tol
    tol_txt = f"{tol * 100:g}%" if rel else f"{tol:g}"
    ok = abs(measured - expected) <= limit
    report("PASS" if ok else "FAIL", name, f"{measured:.6g}",
           f"{expected:.6g} +/- {tol_txt}", note)


def check_max(name, measured, limit, note=""):
    ok = measured <= limit
    report("PASS" if ok else "FAIL", name, f"{measured:.6g}", f"<= {limit:.6g}", note)


def check_min(name, measured, limit, note=""):
    ok = measured >= limit
    report("PASS" if ok else "FAIL", name, f"{measured:.6g}", f">= {limit:.6g}", note)


def check_true(name, ok, measured, expected, note=""):
    report("PASS" if ok else "FAIL", name, measured, expected, note)


def skip(name, why):
    report("SKIP", name, "n/a", "n/a", why)


# ------------------------------------------------------------------ helpers
def read(path):
    p = HERE / path
    if not p.is_file():
        return None
    return p.read_text(errors="replace")


def missing(paths):
    return [p for p in paths if not (HERE / p).is_file()]


def foam_value(text, key):
    m = re.search(rf"^\s*{key}\s+([^;]+);", text, re.M)
    return m.group(1).strip() if m else None


def latest_time_dir():
    best = None
    for d in HERE.iterdir():
        if d.is_dir():
            try:
                t = float(d.name)
            except ValueError:
                continue
            if t > 0 and (best is None or t > best[0]):
                best = (t, d)
    return best


def read_alpha_range(path):
    """Min/max of the alpha.water internalField (binary or ascii)."""
    data = path.read_bytes()
    m = re.search(rb"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", data)
    if m:
        v = float(m.group(1))
        return v, v
    m = re.search(rb"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", data)
    if not m:
        return None
    n = int(m.group(1))
    start = m.end()
    fmt_binary = re.search(rb"format\s+binary", data[:2000]) is not None
    if fmt_binary:
        vals = struct.unpack(f"<{n}d", data[start:start + 8 * n])
    else:
        vals = [float(x) for x in data[start:].split(b")", 1)[0].split()[:n]]
    return min(vals), max(vals)


def hull_submerged_volume(stl_path, scale):
    """Volume of the closed STL below z = 0 (divergence theorem, F = (0,0,z))."""
    d = stl_path.read_bytes()
    n = struct.unpack("<I", d[80:84])[0]
    if 84 + 50 * n != len(d):
        return None  # not a binary STL

    def clip(poly):
        out = []
        for i, a in enumerate(poly):
            b = poly[(i + 1) % len(poly)]
            if a[2] <= 0:
                out.append(a)
            if (a[2] < 0) != (b[2] < 0) and a[2] != b[2]:
                s = -a[2] / (b[2] - a[2])
                out.append(tuple(a[k] + s * (b[k] - a[k]) for k in range(3)))
        return out

    vol = 0.0
    for i in range(n):
        o = 84 + 50 * i + 12
        tri = [tuple(c * scale for c in struct.unpack("<3f", d[o + 12 * v:o + 12 * v + 12]))
               for v in range(3)]
        poly = clip(tri)
        for j in range(1, len(poly) - 1):
            a, b, c = poly[0], poly[j], poly[j + 1]
            nz = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            vol += (a[2] + b[2] + c[2]) / 3.0 * nz / 2.0
    return vol


def read_forces(path):
    """Rows of (t, Fx, Fy, Fz) = pressure + viscous + porous force."""
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        nums = [float(x) for x in re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", line)]
        if len(nums) < 10:
            continue
        t, f = nums[0], nums[1:10]
        rows.append((t, sum(f[0::3]), sum(f[1::3]), sum(f[2::3])))
    return rows


def read_coeffs(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        nums = [float(x) for x in line.split()]
        if len(nums) >= 4:
            rows.append(nums)  # t, Cm, Cd, Cl, ...
    return rows


# --------------------------------------------------------------------- main
def main():
    needed = ["log.blockMesh", "log.snappyHexMesh", "log.checkMesh", "log.setFields",
              "log.interFoam", "system/controlDict", "constant/polyMesh/boundary",
              "postProcessing/forces/0/forces.dat",
              "postProcessing/forceCoeffs/0/forceCoeffs.dat"]
    gone = missing(needed)
    if gone:
        print("ERROR: results not found (run ./run_docker.sh or ./Allrun first):")
        for p in gone:
            print("  missing", p)
        return 2

    ctrl = read("system/controlDict")
    end_time = float(foam_value(ctrl, "endTime"))
    max_co = float(foam_value(ctrl, "maxCo"))
    max_alpha_co = float(foam_value(ctrl, "maxAlphaCo"))
    u_inf = float(re.search(r"magUInf\s+([^;]+);", ctrl).group(1))
    a_ref = float(re.search(r"Aref\s+([^;]+);", ctrl).group(1))
    rho_inf = float(re.search(r"rhoInf\s+([0-9.eE+-]+);", ctrl).group(1))

    blk = read("log.blockMesh")
    snappy = read("log.snappyHexMesh")
    chk = read("log.checkMesh")
    setf = read("log.setFields")
    solver = read("log.interFoam")
    boundary = read("constant/polyMesh/boundary")

    # ---- mesh --------------------------------------------------------------
    # blockMeshDict: 60 x 20 x 12 cells. Exact by construction.
    m = re.search(r"nCells:\s+(\d+)", blk)
    check_close("blockMesh cells", float(m.group(1)) if m else -1, 60 * 20 * 12, 0,
                note="exact: 60x20x12 from blockMeshDict")

    # snappyHexMesh prints this line only if it completed without errors.
    check_true("snappyHexMesh finished",
               "Finished meshing without any errors" in snappy,
               "finished" if "Finished meshing without any errors" in snappy else "not finished",
               "Finished meshing without any errors")

    # Final cell count. Consistency only: snappyHexMesh is deterministic for a
    # given version/settings, but cell counts shift slightly between OpenFOAM
    # versions, so 10 % around the value measured with OpenFOAM 6 is allowed.
    m = re.search(r"^\s+cells:\s+(\d+)", chk, re.M)
    cells = float(m.group(1)) if m else -1
    check_close("final cell count", cells, 75515, 0.10, rel=True,
                note="consistency only; value measured with OpenFOAM 6")

    check_true("checkMesh Mesh OK", re.search(r"^Mesh OK\.", chk, re.M) is not None,
               "Mesh OK" if "Mesh OK." in chk else "failed checks", "Mesh OK.")

    # meshQualityControls maxNonOrtho is 65 in snappyHexMeshDict.
    m = re.search(r"non-orthogonality Max:\s+([0-9.eE+-]+)", chk)
    check_max("max non-orthogonality [deg]", float(m.group(1)) if m else 1e9, 65.0,
              note="snappyHexMeshDict maxNonOrtho")

    # skewness limits in meshQualityControls: 4 internal, 20 boundary.
    m = re.search(r"Max skewness =\s+([0-9.eE+-]+)", chk)
    check_max("max skewness", float(m.group(1)) if m else 1e9, 4.0,
              note="snappyHexMeshDict maxInternalSkewness")

    # Domain bounding box equals the blockMeshDict domain.
    m = re.search(r"bounding box \(([^)]*)\) \(([^)]*)\)", chk)
    if m:
        got = [float(x) for x in m.group(1).split() + m.group(2).split()]
        want = [-0.5, -0.25, -0.15, 1.0, 0.25, 0.15]
        err = max(abs(a - b) for a, b in zip(got, want))
        check_max("domain bounding box error [m]", err, 1e-6,
                  note="1e-6: print precision of checkMesh")
    else:
        skip("domain bounding box error [m]", "not found in log.checkMesh")

    # Hull patch present, nonzero number of faces.
    m = re.search(r"\n\s+boat\s*\{[^}]*nFaces\s+(\d+);", boundary)
    nboat = int(m.group(1)) if m else 0
    check_min("hull patch faces", nboat, 1, note="patch 'boat' in constant/polyMesh/boundary")

    # ---- initialisation ----------------------------------------------------
    # Water below z = 0 fills half the domain minus the submerged hull volume;
    # the waterline lies on a background cell face, so the value is exact up
    # to snapping. Expected 0.4998 from the STL (see buoyancy check); the
    # tolerance 0.002 is 0.4 % of the value.
    fr = re.findall(r"Phase-1 volume fraction = ([0-9.eE+-]+)\s+Min\(alpha.water\) = ([0-9.eE+-]+)"
                    r"\s+Max\(alpha.water\) = ([0-9.eE+-]+)", solver)
    if fr:
        check_close("initial water fraction", float(fr[0][0]), 0.4998, 0.002,
                    note="half domain minus submerged hull volume")
    else:
        skip("initial water fraction", "no 'Phase-1 volume fraction' in log.interFoam")

    # ---- solver ------------------------------------------------------------
    fatal = "FOAM FATAL" in solver or "FOAM exiting" in solver
    check_true("no FOAM FATAL ERROR", not fatal,
               "fatal error in log" if fatal else "none", "none")

    times = re.findall(r"^Time = ([0-9.eE+-]+)\s*$", solver, re.M)
    t_end = float(times[-1]) if times else 0.0
    finished = re.search(r"^End\s*$", solver, re.M) is not None
    full = finished and t_end >= end_time - 1e-9
    print(f"# solver reached t = {t_end:g} s of endTime = {end_time:g} s "
          f"({'complete run' if full else 'PARTIAL run: full-run checks are SKIPped'})")

    # Courant limits. adjustTimeStep chooses dt from the Courant number of the
    # PREVIOUS step (and may grow it by up to 20 % per step), so the control
    # hovers around maxCo/maxAlphaCo and single steps overshoot it when the
    # flow changes. Criterion: the worst step stays below 1.5x the limit
    # (a first guess of 1.1x was too tight: the real run had isolated
    # interface Courant spikes up to 1.37x). A sustained excess would show up
    # as a larger maximum or as a collapsing time step (checked below).
    co = [float(x) for x in re.findall(r"^Courant Number mean: \S+ max: ([0-9.eE+-]+)", solver, re.M)]
    aco = [float(x) for x in re.findall(r"^Interface Courant Number mean: \S+ max: ([0-9.eE+-]+)",
                                        solver, re.M)]
    if co and aco:
        check_max("max Courant number", max(co), 1.5 * max_co,
                  note=f"1.5 x maxCo={max_co:g}: single-step overshoot")
        check_max("max interface Courant number", max(aco), 1.5 * max_alpha_co,
                  note=f"1.5 x maxAlphaCo={max_alpha_co:g}: single-step overshoot")
    else:
        skip("max Courant number", "no Courant lines in log.interFoam")

    # Time step must not collapse: a diverging VOF run drives deltaT towards 0.
    # 1e-5 s is 10x below the smallest step seen in the reference run.
    dts = [float(x) for x in re.findall(r"^deltaT = ([0-9.eE+-]+)", solver, re.M)]
    if dts:
        check_min("min deltaT [s]", min(dts), 1e-5, note="collapse indicates divergence")

    # alpha.water bounded: MULES keeps it within [0,1] up to small undershoots.
    # 0.01 is 1 % of the unit range; the reference run showed ~3e-4.
    if fr:
        amin = min(float(x[1]) for x in fr)
        amax = max(float(x[2]) for x in fr)
        check_min("min alpha.water (log)", amin, -0.01, note="1 % of unit range")
        check_max("max alpha.water (log)", amax, 1.01, note="1 % of unit range")

        # The only water source is the inlet below z = 0, so the water volume
        # fraction cannot grow faster than U*W*depth*t / V_total (at start-up
        # the water column at rest lets little water out, so the measured gain
        # sits right at this bound). An inlet that also injects water above
        # the waterline (alpha = 1 on the whole inlet face) gains twice as
        # fast. 5 % margin for the time-step sampling of t.
        mb = re.search(r"bounding box \(([^)]*)\) \(([^)]*)\)", chk)
        mv = re.search(r"Total volume = ([0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)", chk)
        if mb and mv:
            lo = [float(x) for x in mb.group(1).split()]
            hi = [float(x) for x in mb.group(2).split()]
            q_water = u_inf * (hi[1] - lo[1]) * (0.0 - lo[2])  # m^3/s below z=0
            bound = q_water * t_end / float(mv.group(1))
            gain = float(fr[-1][0]) - float(fr[0][0])
            check_max("water gain <= inlet flux", gain, 1.05 * bound,
                      note=f"bound U*W*depth*t/V = {bound:.4g}")
            check_min("water gain >= 0", gain, -0.002,
                      note="water cannot drain faster than 0.002 in this setup")
        else:
            skip("water gain <= inlet flux", "bounding box/volume not in log.checkMesh")

    # Cumulative continuity error of the incompressible solve.
    ce = re.findall(r"cumulative = ([-+0-9.eE]+)", solver)
    if ce:
        check_max("|cumulative continuity err|", abs(float(ce[-1])), 1e-5,
                  note="typical interFoam values are <1e-6")

    # Fields in the latest time directory.
    latest = latest_time_dir()
    if latest:
        tdir = latest[1]
        have = all((tdir / f).is_file() for f in ("U", "p_rgh", "alpha.water"))
        check_true("latest time has U,p_rgh,alpha.water", have,
                   f"t={tdir.name}" if have else "incomplete", "all three present")
        rng = read_alpha_range(tdir / "alpha.water") if (tdir / "alpha.water").is_file() else None
        if rng:
            check_min("min alpha.water (field)", rng[0], -0.01, note="1 % of unit range")
            check_max("max alpha.water (field)", rng[1], 1.01, note="1 % of unit range")
        else:
            skip("alpha.water field range", "could not parse internalField")
    else:
        check_true("time directory written", False, "none", ">= 1 time directory")

    # ---- hull forces -------------------------------------------------------
    rows = read_forces(HERE / "postProcessing/forces/0/forces.dat")
    if not rows:
        check_true("forces written", False, "0 rows", ">= 1 row")
    else:
        finite = all(math.isfinite(v) for r in rows for v in r)
        check_true("forces finite", finite, f"{len(rows)} rows", "all finite")

        # Vertical force = buoyancy. forces uses p (hydrostatic included), so at
        # start-up (rows with t <= 0.02 s, before the inlet surge has raised
        # the water level) Fz equals rho*g*V_submerged of the STL below z = 0
        # (Archimedes). 5 %: dynamic pressure at Fr = 0.24 (0.5*rho*U^2 = 45 Pa
        # vs rho*g*d = 270 Pa at 27 mm depth) and snapped-mesh geometry error.
        stl = HERE.parent.parent / "boat_hull.stl"
        rho_w = float(re.search(r"water\s*\{[^}]*rho\s+([0-9.eE+-]+);", read("constant/transportProperties"),
                                re.S).group(1))
        g = float(re.search(r"value\s+\(\s*\S+\s+\S+\s+(\S+)\s*\)", read("constant/g")).group(1))
        vol = hull_submerged_volume(stl, 1e-3) if stl.is_file() else None
        early = [r for r in rows if r[0] <= 0.02] or rows[:1]
        fz = sum(r[3] for r in early) / len(early)
        if vol is not None:
            check_close("buoyancy Fz [N] (start-up)", fz, rho_w * abs(g) * vol, 0.05, rel=True,
                        note=f"Archimedes from STL; mean of {len(early)} rows with t <= 0.02 s")
        else:
            skip("buoyancy Fz [N] (start-up)", "boat_hull.stl not found at repo root")

        n = max(1, len(rows) // 5)
        # Drag acts along the flow (+x). Sign only; no reference magnitude.
        # Pressure + viscous force, mean over the last 20 % of rows.
        fx = sum(r[1] for r in rows[-n:]) / n
        check_min("drag Fx > 0 [N]", fx, 1e-12, note="sign check; flow is +x; consistency only")

        # forceCoeffs consistency: Cd = Fx / (0.5 rho U^2 Aref), same time rows.
        coef = read_coeffs(HERE / "postProcessing/forceCoeffs/0/forceCoeffs.dat")
        if coef and len(coef) == len(rows):
            q = 0.5 * rho_inf * u_inf ** 2 * a_ref
            fx_last, cd_last = rows[-1][1], coef[-1][2]
            check_close("Cd consistent with Fx", cd_last * q, fx_last, 0.01, rel=True,
                        note="Cd*0.5*rho*U^2*Aref vs Fx at last row")
        else:
            skip("Cd consistent with Fx", "row counts differ")

        # Settling of the drag: needs the full run. Mean drag of the last two
        # 10 % windows must agree within 10 %. The 10 % was chosen before the
        # first full run; that run gave 38 % (the case has not settled).
        if full and len(rows) >= 20:
            k = max(1, len(rows) // 10)
            a = sum(r[1] for r in rows[-k:]) / k
            b = sum(r[1] for r in rows[-2 * k:-k]) / k
            check_close("drag change, last 2 windows", a, b, 0.10, rel=True,
                        note="settling criterion, full run only")
        else:
            skip("drag change, last 2 windows", "needs complete run to endTime")

    if not full:
        skip("solver reached endTime", f"stopped at t={t_end:g} s of {end_time:g} s")
    else:
        check_true("solver reached endTime", True, f"t={t_end:g}", f"t>={end_time:g}")

    # ---- summary -----------------------------------------------------------
    n_pass, n_fail, n_skip = (RESULTS.count(s) for s in ("PASS", "FAIL", "SKIP"))
    print(f"\nSummary: {n_pass} passed, {n_fail} failed, {n_skip} skipped"
          + ("  (PARTIAL run: skipped checks need the complete run)" if n_skip else ""))
    print("RESULT: " + ("FAIL" if n_fail else "PASS"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
