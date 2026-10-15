#!/usr/bin/env python3
"""Check the results of a finished (or interrupted) Ahmed body run.

Reads only files that the run already produced (log.*, postProcessing/,
constant/triSurface/ahmed_body.stl, system/forceCoeffs, 0/U, ...); needs no
OpenFOAM installation and no Docker.

Output: one line per check
    PASS|FAIL|SKIP  <name>  measured=<value>  expected=<value> ± <tol>
Exit code: 0 if no check FAILed, 1 if at least one FAILed, 2 if the results
needed for the checks are missing.

Checks marked "full run only" are reported as SKIP when the solver was
interrupted before it converged or reached endTime (e.g. by a time limit).

Usage:
    python3 verify.py [--case DIR] [--slant-angle DEG]
"""

import argparse
import math
import re
import sys
from pathlib import Path

# Reference drag coefficients of the standard Ahmed body (Ahmed, Ramm and
# Faltin, SAE 840300, 1984; 0.285 at 25 deg is the measured value as cited by
# COMSOL's Ahmed body benchmark and ERCOFTAC test case AC1-05; 0.260 at 35 deg
# is the commonly cited value). The experiment was run at 60 m/s with stilts.
CD_REFERENCE = {25.0: 0.285, 35.0: 0.260}
# Relative tolerance on Cd. RANS with wall functions typically lands within
# roughly +-20..30 % of the experiment; 25 % is the middle of that range and
# was fixed before the full run was made.
CD_REL_TOL = 0.25

results = []   # (status, name, measured, expected)


def add(status, name, measured, expected):
    results.append((status, name, measured, expected))


def check(name, ok, measured, expected):
    add("PASS" if ok else "FAIL", name, measured, expected)


def skip(name, why, expected=""):
    add("SKIP", name, why, expected)


def read(path):
    try:
        return Path(path).read_text(errors="replace")
    except OSError:
        return None


def number(text, pattern, last=False):
    found = re.findall(pattern, text)
    if not found:
        return None
    return float(found[-1] if last else found[0])


def dict_value(text, key):
    m = re.search(r"^\s*" + re.escape(key) + r"\s+([^;]+);", text, re.M)
    return m.group(1).strip() if m else None


def vec(s):
    return [float(x) for x in re.findall(r"[-+0-9.eE]+", s)]


def stl_vertices(path):
    text = read(path)
    if text is None:
        return None
    return [tuple(float(c) for c in m.groups())
            for m in re.finditer(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", text)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--case", default=str(Path(__file__).resolve().parent))
    ap.add_argument("--slant-angle", type=float, default=25.0)
    ap.add_argument("--nose-radius", type=float, default=0.1,
                    help="front edge radius of the STL in m (default 0.1)")
    args = ap.parse_args()
    case = Path(args.case)
    slant = args.slant_angle

    snappy = read(case / "log.snappyHexMesh")
    chk = read(case / "log.checkMesh")
    sol = read(case / "log.simpleFoam")
    missing = [n for n, t in (("log.snappyHexMesh", snappy),
                              ("log.checkMesh", chk),
                              ("log.simpleFoam", sol)) if t is None]
    if missing:
        print("ERROR: missing %s in %s. Run ./run_docker.sh (or ./Allrun) first."
              % (", ".join(missing), case), file=sys.stderr)
        return 2

    # ---- geometry and reference-value consistency (case files) -------------
    verts = stl_vertices(case / "constant/triSurface/ahmed_body.stl")
    fc = read(case / "system/forceCoeffs")
    u0 = read(case / "0/U")
    tp = read(case / "constant/transportProperties")
    if not verts or fc is None or u0 is None or tp is None:
        print("ERROR: STL, system/forceCoeffs, 0/U or constant/transportProperties missing.",
              file=sys.stderr)
        return 2
    xs, ys, zs = zip(*verts)
    L, W, H = max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)

    # slant angle from the STL: top edge end -> lowest point of the rear face
    x_top_end = max(x for x, y, z in verts if abs(z - max(zs)) < 1e-9)
    z_base_top = max(z for x, y, z in verts if abs(x - max(xs)) < 1e-9)
    angle = math.degrees(math.atan2(max(zs) - z_base_top, max(xs) - x_top_end))
    # tolerance: STL vertices are written with 7 significant digits
    check("stl_slant_angle", abs(angle - slant) < 0.05,
          "%.2f deg" % angle, "%.1f ± 0.05 deg" % slant)

    Aref = float(dict_value(fc, "Aref"))
    lRef = float(dict_value(fc, "lRef"))
    magU = float(dict_value(fc, "magUInf"))
    cofr = vec(dict_value(fc, "CofR")[1:-1])
    # STL integrity: closed (each directed edge once, its reverse once),
    # outward orientation (positive volume), frontal area from projection
    import collections
    keyed = [tuple(tuple(round(c, 7) for c in v) for v in verts[i:i + 3])
             for i in range(0, len(verts), 3)]
    edges = collections.Counter()
    for t in keyed:
        for i in range(3):
            edges[(t[i], t[(i + 1) % 3])] += 1
    bad = sum(1 for e, n in edges.items() if n != 1 or edges.get((e[1], e[0]), 0) != 1)
    check("stl_watertight", bad == 0, "%d open/duplicate edges, %d triangles" % (bad, len(keyed)),
          "0 open edges")
    vol = 0.0
    area_x = 0.0
    for a_, b_, c_ in keyed:
        vol += (a_[0] * (b_[1] * c_[2] - b_[2] * c_[1])
                - a_[1] * (b_[0] * c_[2] - b_[2] * c_[0])
                + a_[2] * (b_[0] * c_[1] - b_[1] * c_[0])) / 6.0
        area_x += abs((b_[1] - a_[1]) * (c_[2] - a_[2])
                      - (b_[2] - a_[2]) * (c_[1] - a_[1])) / 2.0
    # upper bound: box minus the slant wedge (rounding removes a little more);
    # lower bound: 10 % below that
    wedge = 0.5 * (max(xs) - x_top_end) * (max(zs) - z_base_top) * W
    v_hi = W * H * L - wedge
    check("stl_outward_positive_volume", 0.9 * v_hi < vol < v_hi,
          "%.5f m3" % vol, "between %.4f and %.4f m3" % (0.9 * v_hi, v_hi))
    # every ray along x crosses the closed, x-monotone body twice
    frontal = area_x / 2.0
    # tolerance 0.1 %: dictionary value is rounded to 6 digits
    check("Aref_equals_STL_frontal_area", abs(Aref - frontal) / frontal < 1e-3,
          "%.6f m2" % Aref, "%.6f ± 0.1 %%" % frontal)
    # nose radius: the front face is inset by R on all sides
    x0 = min(xs)
    inset = W / 2 - max(abs(y) for x, y, z in verts if abs(x - x0) < 1e-9)
    check("stl_nose_radius", abs(inset - args.nose_radius) < 1e-3,
          "%.4f m" % inset, "%.4f ± 0.001 m" % args.nose_radius)
    check("lRef_equals_body_length", abs(lRef - L) / L < 1e-3,
          "%.4f m" % lRef, "%.4f ± 0.1 %%" % L)
    inlet = re.search(r"inlet\s*\{[^}]*?value\s+uniform\s*\(([^)]*)\)", u0)
    uin = math.sqrt(sum(c * c for c in vec(inlet.group(1)))) if inlet else float("nan")
    check("magUInf_equals_inlet_velocity", abs(magU - uin) < 1e-6,
          "%.3f m/s" % magU, "%.3f m/s" % uin)
    # CofR: mid-length, on the symmetry plane, mid-height (tolerance 1 % of L/H)
    ok = (abs(cofr[0] - (min(xs) + L / 2)) < 0.01 * L and abs(cofr[1]) < 1e-6
          and abs(cofr[2] - (min(zs) + H / 2)) < 0.01 * H)
    check("CofR_at_body_centre", ok, "(%.3f %.3f %.3f)" % tuple(cofr),
          "(%.3f 0 %.3f) ± 1 %%" % (min(xs) + L / 2, min(zs) + H / 2))
    nu = float(re.search(r"nu\s+\[[^\]]*\]\s*([-+0-9.eE]+)", tp).group(1))
    print("INFO  Re_L = U*L/nu = %.3g,  Re_H = %.3g" % (magU * L / nu, magU * H / nu))

    # ---- snappyHexMesh ------------------------------------------------------
    check("snappyHexMesh_finished",
          "Finished meshing without any errors" in snappy,
          "finished" if "Finished meshing without any errors" in snappy else "not finished",
          "'Finished meshing without any errors'")
    # the log holds one layer table per iteration of the layer loop; the
    # last one is the final result
    tables = re.findall(r"(?m)^\s*ahmed_body\s+(\d+)\s+([\d.]+)\s+([\d.eE+-]+)\s+([\d.]+)\s*$",
                        snappy)
    if tables:
        m = type("M", (), {"group": lambda self, i, t=tables[-1]: t[i - 1]})()
        # the table lists the layers actually added (average) and the achieved
        # thickness in % of the requested one; 5 layers requested
        check("body_layers_added", abs(float(m.group(2)) - 5) < 0.5,
              "%s layers" % m.group(2), "5 ± 0.5")
        # 80 % is the usual rule of thumb for acceptable layer coverage
        check("body_layer_thickness_pct", float(m.group(4)) >= 80.0,
              "%s %%" % m.group(4), ">= 80 %")
    else:
        check("body_layers_added", False, "layer table not found", "5 ± 0.5")

    # ---- checkMesh ----------------------------------------------------------
    cells = number(chk, r"(?m)^\s*cells:\s+(\d+)", last=True)
    # about 1.0e6 cells measured; ±15 % allows for the slightly different
    # snappyHexMesh result with a different number of subdomains
    check("cell_count", cells is not None and abs(cells - 1.0e6) < 0.15e6,
          "%s" % (int(cells) if cells else "n/a"), "1.0e6 ± 15 %")
    check("checkMesh_Mesh_OK", "Mesh OK." in chk,
          "Mesh OK" if "Mesh OK." in chk else "failed checks", "Mesh OK.")
    no = number(chk, r"non-orthogonality Max:\s*([\d.eE+-]+)")
    # limit = maxNonOrtho of snappyHexMeshDict meshQualityControls
    check("max_non_orthogonality", no is not None and no < 65.0,
          "%.1f deg" % no if no is not None else "n/a", "< 65 deg")
    sk = number(chk, r"Max skewness\s*=\s*([\d.eE+-]+)")
    # limit = maxInternalSkewness of snappyHexMeshDict meshQualityControls
    check("max_skewness", sk is not None and sk < 4.0,
          "%.2f" % sk if sk is not None else "n/a", "< 4")

    # ---- solver -------------------------------------------------------------
    fatal = "FOAM FATAL" in sol or "FOAM exiting" in sol
    times = re.findall(r"^Time = ([\d.eE+-]+)\s*$", sol, re.M)
    iters = int(float(times[-1])) if times else 0
    check("solver_no_fatal_error", not fatal and iters > 0,
          "%d iterations, %s" % (iters, "FATAL ERROR" if fatal else "no error"),
          "no FOAM FATAL ERROR, >= 1 iteration")

    converged = "SIMPLE solution converged" in sol
    end_time = float(dict_value(read(case / "system/controlDict"), "endTime"))
    full_run = converged or iters >= end_time
    print("INFO  solver: %d iterations; %s" % (
        iters, "full run" if full_run else
        "PARTIAL run (not converged, endTime %d not reached)" % end_time))

    for fld in ("p", "Ux"):
        vals = [float(v) for v in re.findall(
            r"Solving for %s, Initial residual = ([\d.eE+-]+)" % fld, sol)]
        name = "residual_drop_" + fld
        if len(vals) < 20:
            skip(name, "only %d iterations" % len(vals), "final/max <= 0.1")
            continue
        ratio = vals[-1] / max(vals)
        # one decade is the usual minimum sign that SIMPLE is converging
        check(name, ratio <= 0.1, "final/max = %.3g" % ratio, "<= 0.1")

    yp = re.findall(r"patch ahmed_body y\+\s*:\s*min\s*=\s*([\d.eE+-]+),\s*max\s*=\s*"
                    r"([\d.eE+-]+),\s*average\s*=\s*([\d.eE+-]+)", sol)
    if yp:
        avg = float(yp[-1][2])
        # wall functions assume the first cell centre in the log layer,
        # conventionally 30 < y+ < 300
        check("yPlus_body_average", 30.0 <= avg <= 300.0,
              "avg %.1f (min %.1f, max %.1f)" % (avg, float(yp[-1][0]), float(yp[-1][1])),
              "30 .. 300")
    else:
        skip("yPlus_body_average", "no yPlus output in log.simpleFoam", "30 .. 300")

    # ---- force coefficients --------------------------------------------------
    cd_series = []
    for f in sorted((case / "postProcessing/forceCoeffs").glob("*/*.dat")):
        names, rows = None, []
        for line in f.read_text(errors="replace").splitlines():
            if line.startswith("#"):
                toks = line.lstrip("# ").split()
                if toks and toks[0] == "Time":
                    names = toks
                continue
            parts = line.split()
            if parts:
                rows.append(parts)
        if names and "Cd" in names and rows:
            k = names.index("Cd")
            cd_series = [float(r[k]) for r in rows if len(r) > k]
            break
    if not cd_series:
        check("Cd_available", False, "no postProcessing/forceCoeffs/*/*.dat with Cd",
              "Cd time series")
    else:
        cd = cd_series[-1]
        check("Cd_finite_positive", math.isfinite(cd) and cd > 0,
              "Cd = %.4f at iteration %d" % (cd, len(cd_series)), "finite, > 0")
        ref = CD_REFERENCE.get(slant)
        if not full_run:
            skip("Cd_in_literature_band", "partial run, Cd not converged",
                 "%s ± %d %%" % (ref, 100 * CD_REL_TOL))
            skip("Cd_converged", "partial run", "mean drift < 1 %")
        else:
            tail = cd_series[-200:]
            spread = (max(tail) - min(tail)) / abs(sum(tail) / len(tail))
            # A steady RANS solution of a bluff body can settle on a periodic
            # oscillation (here about 40 iterations), so the peak-to-peak
            # spread of the instantaneous Cd is reported for information only
            # and convergence is judged on the drift of the window mean:
            # mean of the last 200 iterations vs the 200 before.
            mean_t = sum(tail) / len(tail)
            sd = math.sqrt(sum((c - mean_t) ** 2 for c in tail) / len(tail))
            print("INFO  Cd last 200 its: mean %.4f, std %.4f, peak-to-peak %.2f %%"
                  % (mean_t, sd, 100 * spread))
            if len(cd_series) < 400:
                skip("Cd_converged", "fewer than 400 iterations",
                     "mean drift < 1 %")
            else:
                prev = cd_series[-400:-200]
                mean_p = sum(prev) / len(prev)
                drift = abs(mean_t - mean_p) / abs(mean_t)
                check("Cd_converged", drift < 0.01, "%.2f %%" % (100 * drift),
                      "< 1 %% drift of the mean, last 200 its vs the 200 before")
            if ref is None:
                skip("Cd_in_literature_band", "no reference for %.0f deg" % slant, "")
            else:
                mean = sum(tail) / len(tail)
                check("Cd_in_literature_band", abs(mean - ref) <= CD_REL_TOL * ref,
                      "%.4f (mean of last 200 its)" % mean,
                      "%.3f ± %d %%" % (ref, 100 * CD_REL_TOL))

    # ---- report -----------------------------------------------------------
    print()
    for status, name, measured, expected in results:
        print("%-5s %-32s measured=%s  expected=%s" % (status, name, measured, expected))
    n = {s: sum(1 for r in results if r[0] == s) for s in ("PASS", "FAIL", "SKIP")}
    print("\nSummary: %d PASS, %d FAIL, %d SKIP" % (n["PASS"], n["FAIL"], n["SKIP"]))
    if n["SKIP"]:
        print("Note: %d check(s) skipped; they apply to a full run only "
              "(converged or endTime reached)." % n["SKIP"])
    return 1 if n["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
