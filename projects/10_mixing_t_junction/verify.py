#!/usr/bin/env python3
"""Check the results of a finished run of the mixing T-junction case.

Usage:  python3 verify.py [case_dir]        (default: directory of this file)

Reads only files already produced by ./Allrun (log.*, the latest time
directory, constant/polyMesh) and the case dictionaries (0/, constant/).
Needs no OpenFOAM installation and only the Python standard library.
Exit code 0 only if every check passes.

Tolerances (justified next to each check) are acceptance thresholds, not
reference values: no experimental or benchmark data exist for this case, so
the checks are consistency checks (conservation, boundedness, convergence).
"""
import math
import os
import re
import sys

# Geometry of the inlets (system/blockMeshDict): channel height 0.1 m,
# branch width 0.05 m, depth (z) 0.01 m.
MAIN_AREA = 0.1 * 0.01
BRANCH_AREA = 0.05 * 0.01
EXPECTED_CELLS = 8800

results = []   # (ok, name, measured, expected, tol)


def record(ok, name, measured, expected, tol):
    results.append((bool(ok), name, measured, expected, tol))


def fail_missing(msg):
    print("FAIL  missing-input  %s" % msg)
    print("SUMMARY: cannot verify; run ./Allrun (or ./run_docker.sh) first")
    sys.exit(2)


def read(path):
    if not os.path.isfile(path):
        fail_missing("file not found: %s" % path)
    with open(path, errors="replace") as fh:
        return fh.read()


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def patch_block(text, patch):
    """Body of `patch { ... }` inside boundaryField / boundary."""
    m = re.search(r"\b%s\s*\{(.*?)\}" % re.escape(patch), text, flags=re.S)
    if not m:
        fail_missing("patch %s not found" % patch)
    return m.group(1)


def number_list(text, start):
    """Parse the first `N ( ... )` list at or after index `start`."""
    m = re.compile(r"(\d+)\s*\(", re.S).search(text, start)
    if not m:
        return None
    n = int(m.group(1))
    end = text.index("\n)", m.end())
    body = text[m.end():end]
    return n, body


def scalar_list(text):
    """All values of a `List<scalar> N (...)` or `List<label>` listing."""
    text = strip_comments(text)
    start = text.index("FoamFile")
    start = text.index("}", start)
    r = number_list(text, start)
    if r is None:
        fail_missing("cannot parse list")
    n, body = r
    vals = [float(v) for v in body.split()]
    if len(vals) != n:
        fail_missing("list length mismatch (%d != %d)" % (len(vals), n))
    return vals


def vector_list(text):
    text = strip_comments(text)
    start = text.index("}", text.index("FoamFile"))
    n, body = number_list(text, start)
    vals = [tuple(float(x) for x in ln.strip().strip("()").split())
            for ln in body.strip().split("\n")]
    if len(vals) != n:
        fail_missing("vector list length mismatch")
    return vals


def boundary_values(text, patch, nfaces):
    """Per-face values of a boundary patch (uniform or nonuniform scalar)."""
    block = patch_block(strip_comments(text[text.index("boundaryField"):]), patch)
    m = re.search(r"value\s+uniform\s+([-+0-9.eE]+)\s*;", block)
    if m:
        return [float(m.group(1))] * nfaces
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", block)
    if m:
        n = int(m.group(1))
        end = block.index(")", m.end())
        vals = [float(v) for v in block[m.end():end].split()]
        if len(vals) == n == nfaces:
            return vals
    fail_missing("cannot read boundary values of %s" % patch)


def latest_time(case):
    times = []
    for d in os.listdir(case):
        if os.path.isdir(os.path.join(case, d)) and re.fullmatch(r"[0-9.]+(e[-+]?\d+)?", d):
            try:
                t = float(d)
            except ValueError:
                continue
            if t > 0:
                times.append((t, d))
    if not times:
        fail_missing("no result time directory in %s" % case)
    return max(times)[1]


def main():
    case = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else \
        os.path.dirname(os.path.abspath(__file__))

    # ---- logs -----------------------------------------------------------
    chk = read(os.path.join(case, "log.checkMesh"))
    sol = read(os.path.join(case, "log.buoyantSimpleFoam"))

    # 1. Mesh check: exact structural facts, no tolerance.
    m = re.search(r"^\s*cells:\s+(\d+)", chk, flags=re.M)
    cells = int(m.group(1)) if m else -1
    record(cells == EXPECTED_CELLS, "mesh-cell-count", cells, EXPECTED_CELLS, "exact")
    record("Mesh OK." in chk, "checkMesh-ok", "Mesh OK" if "Mesh OK." in chk else "failed checks",
           "Mesh OK", "exact")

    # 2. Solver finished normally. (The "Radiation model not active ...
    # FOAM FATAL IO ERROR" text is a harmless message of the radiation
    # library, so look for the abort markers instead.)
    bad = re.search(r"FOAM exiting|FOAM aborting|Floating point exception|Segmentation fault", sol)
    ended = re.search(r"^End\s*$", sol, flags=re.M) is not None
    record(ended and not bad, "solver-finished",
           "End" if ended and not bad else (bad.group(0) if bad else "no End"),
           "End, no abort", "exact")

    # 3. Convergence by residualControl before endTime (2000 iterations).
    m = re.search(r"SIMPLE solution converged in (\d+) iterations", sol)
    iters = int(m.group(1)) if m else None
    end_time = 2000  # system/controlDict endTime
    cd = read(os.path.join(case, "system", "controlDict"))
    me = re.search(r"^\s*endTime\s+([0-9.eE+-]+);", cd, flags=re.M)
    if me:
        end_time = float(me.group(1))
    record(iters is not None and iters < end_time, "converged-before-endTime",
           iters if iters is not None else "not converged", "< %d" % end_time, "exact")

    # 4. All initial residuals of the last iteration below residualControl.
    # 1e-4 is the value set in system/fvSolution (residualControl).
    last = sol.rsplit("\nTime = ", 1)[-1]
    res = {}
    for fld, val in re.findall(r"Solving for (\w+), Initial residual = ([0-9.eE+-]+)", last):
        res[fld] = float(val)
    worst = max(res.values()) if res else float("nan")
    record(res and worst <= 1e-4, "final-initial-residuals", "%.2e" % worst, "<= 1e-4",
           "residualControl")

    # ---- fields ---------------------------------------------------------
    t = latest_time(case)
    tdir = os.path.join(case, t)
    T = scalar_list(read(os.path.join(tdir, "T")))
    owner = [int(v) for v in scalar_list(read(os.path.join(case, "constant", "polyMesh", "owner")))]
    bnd = strip_comments(read(os.path.join(case, "constant", "polyMesh", "boundary")))

    def patch_info(name):
        b = patch_block(bnd, name)
        return (int(re.search(r"nFaces\s+(\d+)", b).group(1)),
                int(re.search(r"startFace\s+(\d+)", b).group(1)))

    phi_txt = read(os.path.join(tdir, "phi"))
    flux = {}
    for p in ("mainInlet", "branchInlet", "outlet"):
        n, start = patch_info(p)
        flux[p] = (boundary_values(phi_txt, p, n), [owner[start + i] for i in range(n)])
    m_main = -sum(flux["mainInlet"][0])
    m_branch = -sum(flux["branchInlet"][0])
    m_out = sum(flux["outlet"][0])

    # density model and BC values from the case dictionaries
    thermo = strip_comments(read(os.path.join(case, "constant", "thermophysicalProperties")))
    rho0 = float(re.search(r"\brho0\s+([0-9.eE+-]+);", thermo).group(1))
    T0 = float(re.search(r"\bT0\s+([0-9.eE+-]+);", thermo).group(1))
    beta = float(re.search(r"\bbeta\s+([0-9.eE+-]+);", thermo).group(1))
    T_txt = strip_comments(read(os.path.join(case, "0", "T")))
    U_txt = strip_comments(read(os.path.join(case, "0", "U")))
    T_main = float(re.search(r"value\s+uniform\s+([0-9.eE+-]+);", patch_block(T_txt[T_txt.index("boundaryField"):], "mainInlet")).group(1))
    T_br = float(re.search(r"value\s+uniform\s+([0-9.eE+-]+);", patch_block(T_txt[T_txt.index("boundaryField"):], "branchInlet")).group(1))
    um = re.search(r"value\s+uniform\s+\(([^)]*)\)", patch_block(U_txt[U_txt.index("boundaryField"):], "mainInlet")).group(1).split()
    ub = re.search(r"value\s+uniform\s+\(([^)]*)\)", patch_block(U_txt[U_txt.index("boundaryField"):], "branchInlet")).group(1).split()
    U_main = math.hypot(float(um[0]), float(um[1]))
    U_br = math.hypot(float(ub[0]), float(ub[1]))

    def rho(temp):
        return rho0 * (1.0 - beta * (temp - T0))

    # 5. Inlet mass flow equals rho(T_in) * U * A. Both U and T are fixed
    # at the inlets, so only write precision (8 digits) can differ: 1e-4.
    for name, got, exp in (("mass-flow-mainInlet", m_main, rho(T_main) * U_main * MAIN_AREA),
                           ("mass-flow-branchInlet", m_branch,
                            rho(T_br) * U_br * BRANCH_AREA)):
        rel = abs(got - exp) / exp
        record(rel <= 1e-4, name, "%.6f kg/s (rel.err %.1e)" % (got, rel),
               "%.6f kg/s" % exp, "1e-4 rel")

    # 6. Global mass balance of the boundary fluxes. Tolerance 1e-4 of the
    # outflow: the p_rgh equation is converged to 1e-4 (residualControl), so
    # an imbalance above that level would indicate an unconverged or
    # inconsistent solution.
    imb = abs(m_main + m_branch - m_out) / m_out
    record(imb <= 1e-4, "mass-balance", "%.2e rel" % imb, "0", "1e-4 rel")

    # 7. Energy balance (adiabatic walls, constant Cp): the mass-flux
    # weighted outlet temperature must equal the inlet-flux weighted mean.
    # Tolerance 0.5 K = 1 % of the inlet difference: residualControl does
    # not bound the global enthalpy imbalance of the discretisation
    # (linearUpwind, steady SIMPLE), so this is an acceptance level.
    vals, cells_o = flux["outlet"]
    T_out_cells = [T[c] for c in cells_o]   # zeroGradient: face T = cell T
    T_out = sum(f * tc for f, tc in zip(vals, T_out_cells)) / sum(vals)
    T_mix = (m_main * T_main + m_branch * T_br) / (m_main + m_branch)
    record(abs(T_out - T_mix) <= 0.5, "outlet-temperature-energy-balance",
           "%.2f K" % T_out, "%.2f K" % T_mix, "0.5 K")

    # 8. Temperature stays within the inlet values. Central/linearUpwind
    # schemes are not bounded, so a small overshoot is allowed: 1.0 K =
    # 2 % of the 50 K inlet difference.
    lo, hi = min(T_br, T_main), max(T_br, T_main)
    tmin, tmax = min(T), max(T)
    record(tmin >= lo - 1.0 and tmax <= hi + 1.0, "temperature-bounded",
           "%.2f .. %.2f K" % (tmin, tmax), "%.0f .. %.0f K" % (lo, hi), "1.0 K")

    # 9. Streams are not fully mixed at the outlet and the cold stream is on
    # top (branch side): the outlet temperature spans T_mix, with the
    # hottest cell at the bottom and the coldest at the top. Face order of
    # the outlet patch follows y (blockMesh), so first = bottom.
    ordered = T_out_cells[0] > T_mix > T_out_cells[-1]
    record(ordered, "outlet-stratification",
           "%.1f K (bottom) / %.1f K (top)" % (T_out_cells[0], T_out_cells[-1]),
           "bottom > %.1f K > top" % T_mix, "exact")

    # ---- report ---------------------------------------------------------
    print("Case: %s   latest time: %s" % (case, t))
    nfail = 0
    for ok, name, measured, expected, tol in results:
        if not ok:
            nfail += 1
        print("%s  %s  measured=%s  expected=%s ± %s"
              % ("PASS" if ok else "FAIL", name, measured, expected, tol))
    print("SUMMARY: %d/%d checks passed" % (len(results) - nfail, len(results)))
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
