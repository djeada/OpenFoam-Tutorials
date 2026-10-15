#!/usr/bin/env python3
"""Check the results of a finished cylinder run (Re = 200, pimpleFoam).

Only reads files produced by ./Allrun (log.*, constant/polyMesh/boundary,
postProcessing/forceCoeffs, time directories). Needs no OpenFOAM and no
Docker; standard library only.

Usage:
    python3 verify.py [CASE_DIR]        (default: directory of this script)

Output: one line per check
    PASS|FAIL  <name>  measured=<v>  expected=<v> +- <tol>
Exit status: 0 if every check passes, 1 if any check fails, 2 if the
results needed for the checks are missing.
"""

import math
import os
import re
import sys

# Reference values for 2D laminar flow past a circular cylinder at Re = 200,
# as listed in secondary comparison tables (original articles not consulted):
#   Linnick & Fasel 2005       St 0.197  Cd 1.34 +- 0.044  Cl +-0.69
#   Taira & Colonius 2007      St 0.195  Cd 1.34 +- 0.047  Cl +-0.68
#   Liu, Zheng & Sung 1998     St 0.192  Cd 1.31 +- 0.049  Cl +-0.69
#   Williamson (experiment)    St 0.197
# (the Liu et al. row appears in two independent tables; only the first three
# rows and Williamson were seen in one table). Hence: St ~ 0.197, mean
# Cd ~ 1.34, Cl amplitude ~ 0.69; all approximate, not checked against the
# original articles.
# Tolerances are about +-10 %: published 2D codes scatter by ~3-5 %, and this
# domain (blockage 1/16, 8 D upstream, 25 D downstream, symmetry sides) and
# mesh differ from the unbounded-domain references.
REF_ST, TOL_ST = 0.197, 0.02
REF_CD, TOL_CD = 1.34, 0.13
REF_CL_AMP, TOL_CL_AMP = 0.69, 0.10
TOL_CL_MEAN = 0.05          # mean Cl is 0 by symmetry; finite window noise
EXPECTED_CELLS, TOL_CELLS = 22605, 0.10   # Gmsh version may change counts
MAX_NON_ORTHO = 70.0        # checkMesh error threshold
MAX_SKEWNESS = 4.0          # checkMesh error threshold (internal faces)
COURANT_SLACK = 1.05        # adjustable dt reacts one step late
MAX_CONTINUITY = 1e-5       # |cumulative continuity error|, incompressible
MIN_CL_AMP_SHEDDING = 0.1   # below this the wake did not start shedding

EXPECTED_PATCH_TYPES = {
    "inlet": "patch",
    "outlet": "patch",
    "top": "symmetryPlane",
    "bottom": "symmetryPlane",
    "cylinder": "wall",
    "frontAndBack": "empty",
}

results = []


class Missing(Exception):
    pass


def check(name, ok, measured, expected):
    results.append((bool(ok), name, measured, expected))
    print(f"{'PASS' if ok else 'FAIL'}  {name}  measured={measured}  "
          f"expected={expected}")


def read(case, rel):
    path = os.path.join(case, rel)
    if not os.path.isfile(path):
        raise Missing(f"missing {rel} (run ./run_docker.sh or ./Allrun first)")
    with open(path, errors="replace") as fh:
        return fh.read()


def foam_value(text, key):
    m = re.search(r"^\s*" + key + r"\s+([^;\s]+)\s*;", text, re.M)
    if not m:
        raise Missing(f"entry '{key}' not found")
    return m.group(1)


def search(pattern, text, what):
    m = re.search(pattern, text, re.M)
    if not m:
        raise Missing(f"could not find {what} in log")
    return m.group(1)


def mean(v):
    return sum(v) / len(v)


def upward_crossings(t, y, level):
    out = []
    for i in range(1, len(t)):
        if y[i - 1] < level <= y[i]:
            f = (level - y[i - 1]) / (y[i] - y[i - 1])
            out.append(t[i - 1] + f * (t[i] - t[i - 1]))
    return out


def latest_time_dir(case):
    times = []
    for d in os.listdir(case):
        if re.fullmatch(r"[0-9]+(\.[0-9]+)?", d) and os.path.isdir(
                os.path.join(case, d)):
            times.append((float(d), d))
    return max(times)[1] if times else None


def main():
    case = os.path.abspath(
        sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__))

    # ---- Mesh ----------------------------------------------------------
    cm = read(case, "log.checkMesh")
    check("checkMesh_ok", "Mesh OK" in cm, "Mesh OK" if "Mesh OK" in cm
          else "failed", "Mesh OK")

    cells = int(search(r"^\s*cells:\s+(\d+)", cm, "cell count in log.checkMesh"))
    check("cell_count",
          abs(cells - EXPECTED_CELLS) <= TOL_CELLS * EXPECTED_CELLS,
          cells, f"{EXPECTED_CELLS} +- {int(TOL_CELLS * 100)}%")

    no = float(search(r"non-orthogonality Max:\s*([0-9.eE+-]+)", cm,
                     "non-orthogonality in log.checkMesh"))
    check("max_non_orthogonality", no < MAX_NON_ORTHO, no, f"< {MAX_NON_ORTHO}")

    sk = float(search(r"Max skewness =\s*([0-9.eE+-]+)", cm,
                 "skewness in log.checkMesh"))
    check("max_skewness", sk < MAX_SKEWNESS, sk, f"< {MAX_SKEWNESS}")

    boundary = read(case, "constant/polyMesh/boundary")
    got = dict(re.findall(
        r"^\s{4}(\w+)\s*\{[^}]*?\btype\s+(\w+);", boundary, re.M | re.S))
    bad = {k: got.get(k) for k, v in EXPECTED_PATCH_TYPES.items()
           if got.get(k) != v}
    check("patch_types", not bad,
          "all match" if not bad else bad, "all match")

    # ---- Solver --------------------------------------------------------
    log = read(case, "log.pimpleFoam")
    check("no_fatal_error", "FOAM FATAL" not in log and "FOAM exiting" not in log,
          "none" if "FOAM FATAL" not in log else "FATAL ERROR found", "none")
    finished = re.search(r"^End\s*$", log, re.M) is not None
    check("solver_finished", finished, "End" if finished else "no 'End'", "End")

    ctrl = read(case, "system/controlDict")
    end_time = float(foam_value(ctrl, "endTime"))
    max_co = float(foam_value(ctrl, "maxCo"))
    times = [float(x) for x in re.findall(r"^Time = ([0-9.eE+-]+)\s*$", log, re.M)]
    last = times[-1] if times else float("nan")
    check("end_time_reached", abs(last - end_time) < 1e-3,
          last, f"{end_time} +- 0.001")

    cos = [float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
    cmax = max(cos) if cos else float("nan")
    check("max_courant", cmax <= max_co * COURANT_SLACK, cmax,
          f"<= {max_co * COURANT_SLACK:.3f} (maxCo {max_co} + 5%)")

    cont = [abs(float(x)) for x in re.findall(
        r"cumulative = ([0-9.eE+-]+)", log)]
    cc = max(cont) if cont else float("nan")
    check("continuity_error", cc < MAX_CONTINUITY, f"{cc:.2e}",
          f"< {MAX_CONTINUITY:.0e}")

    last_dir = latest_time_dir(case)
    have = []
    if last_dir:
        have = [f for f in ("U", "p", "UMean", "pMean", "UPrime2Mean",
                            "pPrime2Mean")
                if not os.path.isfile(os.path.join(case, last_dir, f))]
    check("last_time_fields", last_dir is not None and not have,
          f"{last_dir}: " + ("complete" if not have else f"missing {have}"),
          "U p UMean pMean UPrime2Mean pPrime2Mean present")

    # ---- Force coefficients -------------------------------------------
    fdir = os.path.join(case, "postProcessing", "forceCoeffs")
    if not os.path.isdir(fdir):
        raise Missing("missing postProcessing/forceCoeffs/")
    sub = sorted(d for d in os.listdir(fdir))
    fpath = None
    for d in sub:
        p = os.path.join(fdir, d, "forceCoeffs.dat")
        if os.path.isfile(p):
            fpath = p  # newest start-time directory wins
    if fpath is None:
        raise Missing("missing postProcessing/forceCoeffs/*/forceCoeffs.dat")
    txt = open(fpath).read()

    d_ref = float(re.search(r"#\s*lRef\s*:\s*(\S+)", txt).group(1))
    u_ref = float(re.search(r"#\s*magUInf\s*:\s*(\S+)", txt).group(1))
    t, cd, cl = [], [], []
    for line in txt.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c = line.split()
        t.append(float(c[0]))
        cd.append(float(c[2]))
        cl.append(float(c[3]))

    t_start = float(foam_value(
        re.search(r"fieldAverage\s*\{(.*?)\n    \}", ctrl, re.S).group(1),
        "timeStart"))
    idx = [i for i, x in enumerate(t) if x >= t_start]
    if len(idx) < 50:
        raise Missing(f"too few force samples after t = {t_start}")
    tw = [t[i] for i in idx]
    cdw = [cd[i] for i in idx]
    clw = [cl[i] for i in idx]

    cd_m = mean(cdw)
    check("mean_Cd", abs(cd_m - REF_CD) <= TOL_CD, f"{cd_m:.3f}",
          f"{REF_CD} +- {TOL_CD} (t >= {t_start:g})")

    cl_m = mean(clw)
    check("mean_Cl", abs(cl_m) <= TOL_CL_MEAN, f"{cl_m:.3f}",
          f"0 +- {TOL_CL_MEAN}")

    amp = (max(clw) - min(clw)) / 2.0
    shed = amp > MIN_CL_AMP_SHEDDING
    check("shedding_developed", shed, f"Cl amplitude {amp:.3f}",
          f"> {MIN_CL_AMP_SHEDDING}")
    check("Cl_amplitude", abs(amp - REF_CL_AMP) <= TOL_CL_AMP, f"{amp:.3f}",
          f"{REF_CL_AMP} +- {TOL_CL_AMP}")

    if shed:
        cross = upward_crossings(tw, clw, cl_m)
        if len(cross) >= 4:
            period = (cross[-1] - cross[0]) / (len(cross) - 1)
            st = d_ref / (u_ref * period)
            check("Strouhal", abs(st - REF_ST) <= TOL_ST,
                  f"{st:.3f} ({len(cross) - 1} periods)",
                  f"{REF_ST} +- {TOL_ST}")
        else:
            check("Strouhal", False, f"only {len(cross)} crossings",
                  ">= 4 crossings")
    else:
        check("Strouhal", False, "no shedding in window", "periodic Cl")

    # ---- Summary -------------------------------------------------------
    n_fail = sum(1 for r in results if not r[0])
    print(f"\n{len(results) - n_fail}/{len(results)} checks passed"
          + ("" if n_fail == 0 else f", {n_fail} FAILED"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Missing as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
