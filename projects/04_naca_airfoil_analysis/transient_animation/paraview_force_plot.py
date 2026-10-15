#!/usr/bin/env python3
"""Plot Cl and Cd against time for the NACA 0012 transient case (offscreen).

Usage (after ../run_docker.sh transient):
    pvbatch paraview_force_plot.py [--dat postProcessing/forceCoeffs/0/forceCoeffs.dat]
                                   [--out ../images/transient_forces.png]

Needs only ParaView (no matplotlib). The OpenFOAM .dat file is converted to a
temporary CSV with the standard library, then drawn in a line chart view.
"""

import argparse
import csv
import re
import tempfile
from pathlib import Path

from paraview.simple import (
    CSVReader,
    CreateView,
    SaveScreenshot,
    Show,
)

HERE = Path(__file__).resolve().parent

ap = argparse.ArgumentParser()
ap.add_argument("--dat", default=str(HERE / "postProcessing/forceCoeffs/0/forceCoeffs.dat"))
ap.add_argument("--out", default=str(HERE.parent / "images/transient_forces.png"))
ap.add_argument("--size", nargs=2, type=int, default=[1200, 600])
args, _ = ap.parse_known_args()

columns = None
rows = []
for line in Path(args.dat).read_text().splitlines():
    if line.startswith("#"):
        if re.match(r"#\s*Time\b", line):
            columns = line.lstrip("#").split()
        continue
    if line.strip():
        rows.append(line.split())
if columns is None or not rows:
    raise SystemExit("could not parse %s" % args.dat)

want = ("Time", "Cl", "Cd")
keep = [columns.index(c) for c in want]

tmp = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, newline="")
writer = csv.writer(tmp)
writer.writerow(want)
for r in rows:
    writer.writerow([r[i] for i in keep])
tmp.close()

reader = CSVReader(FileName=[tmp.name])
view = CreateView("XYChartView")
view.ViewSize = args.size
view.ChartTitle = "NACA 0012, 15 degrees, Re = 10 000: force coefficients"
view.BottomAxisTitle = "time [s]"
view.LeftAxisTitle = "Cl, Cd"
view.ShowLegend = 1

disp = Show(reader, view)
disp.UseIndexForXAxis = 0
disp.XArrayName = "Time"
disp.SeriesVisibility = ["Cl", "Cd"]
disp.SeriesColor = ["Cl", "0.15", "0.35", "0.85", "Cd", "0.85", "0.25", "0.10"]
disp.SeriesLineThickness = ["Cl", "2", "Cd", "2"]

SaveScreenshot(args.out, view, ImageResolution=args.size)
Path(tmp.name).unlink()
print("Wrote", args.out)
