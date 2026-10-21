#!/usr/bin/env python3
"""Check completed runs; no OpenFOAM or third-party Python packages required.

Acceptance checks are numerical consistency checks, not experimental validation.
Missing/incomplete output is a failure, including an otherwise clean early stop.
"""
import math
from pathlib import Path
import re
import sys

from analyze import ROOT, SCENARIOS, history


def scalar_field(path):
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(), flags=re.S)
    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        values = [float(uniform.group(1))]
    else:
        match = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
        if not match:
            raise ValueError(f"cannot parse scalar internalField in {path}")
        values = [float(value) for value in match.group(2).split()]
        if len(values) != int(match.group(1)):
            raise ValueError(f"wrong scalar count in {path}")
    if not values or not all(map(math.isfinite, values)):
        raise ValueError(f"non-finite or empty internalField in {path}")
    return values


def verify_case(case):
    info, volume, rows = history(case)
    checks = []

    def check(ok, name):
        checks.append(bool(ok))
        print(f"{'PASS' if ok else 'FAIL'}  {case.name}: {name}")

    mesh = (case / "log.checkMesh").read_text()
    log = (case / "log.buoyantBoussinesqPimpleFoam").read_text()
    check("Mesh OK." in mesh and not re.search(r"Failed\s+[1-9]", mesh), "mesh quality")
    fatal = re.search(r"FOAM FATAL|FOAM exiting|FOAM aborting|Floating point exception|Segmentation fault", log)
    check(re.search(r"^End\s*$", log, re.M) is not None and not fatal, "solver completed normally")
    times = [float(time) for time in re.findall(r"^Time = (\S+)", log, re.M)]
    end = info["end_time"]
    check(bool(times) and math.isclose(times[-1], end, abs_tol=1e-7), "solver reached requested end time")
    check(math.isclose(rows[-1]["time_s"], end, abs_tol=1e-7), "statistics reach requested end time")
    expected_times = [i * info["write_interval"] for i in range(1, round(end / info["write_interval"]) + 1)]
    written_times = [row["time_s"] for row in rows]
    check(all(any(math.isclose(t, w, abs_tol=1e-7) for w in written_times) for t in expected_times),
          "all requested statistics were written")
    # Volume is an exact geometric property. Allow 0.1% for output rounding.
    check(math.isclose(volume, info["room_volume"], rel_tol=1e-3), "room zone excludes the exterior and pane")
    # Window area is the frame-plane area in both cases, not the top-gap area.
    flux_files = list((case / "postProcessing" / "windowFlux").glob("*/surfaceFieldValue.dat"))
    area = re.search(r"^#\s*Area\s*:\s*(\S+)", flux_files[0].read_text(), re.M | re.I)
    check(bool(area) and math.isclose(float(area.group(1)), 1.2 * 0.1, rel_tol=1e-3),
          "flux zone spans the window plane")
    result_dirs = []
    for path in case.iterdir():
        if path.is_dir():
            try:
                if math.isclose(float(path.name), end, abs_tol=1e-7):
                    result_dirs.append(path)
            except ValueError:
                pass
    if len(result_dirs) != 1:
        raise ValueError(f"{case.name}: missing final field directory at {end:g} seconds")
    final = result_dirs[0]
    c, T = scalar_field(final / "c"), scalar_field(final / "T")
    scalar_field(final / "p_rgh")
    check(min(c) >= -1e-5 and max(c) <= 1 + 1e-5, "local tracer bounded in [0,1] (1e-5 tolerance)")
    check(min(T) >= info["outdoor_temperature"] - 0.05 and
          max(T) <= info["radiator_temperature"] + 0.05, "local temperature bounded (0.05 K tolerance)")
    check(all(-0.001 <= row["fresh_air_percent"] <= 100.001 for row in rows), "room tracer averages bounded")
    courant = [float(value) for value in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)]
    check(bool(courant) and all(math.isfinite(value) and value <= 0.5 for value in courant),
          "maximum Courant number <= 0.5 (target 0.4)")
    check(all(math.isfinite(row["exchange_ACH"]) and row["exchange_ACH"] >= 0 for row in rows),
          "finite non-negative exchange rate")
    return all(checks)


def main():
    passed = True
    for scenario in SCENARIOS:
        try:
            passed = verify_case(ROOT / "cases" / scenario) and passed
        except (OSError, ValueError, KeyError) as error:
            print(f"FAIL  {scenario}: {error}")
            passed = False
    print("Consistency checks passed; no experimental validation is implied." if passed else
          "Verification failed. Run the cases to completion and inspect their logs.")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
