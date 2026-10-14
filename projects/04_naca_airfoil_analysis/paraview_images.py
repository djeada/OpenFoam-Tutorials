#!/usr/bin/env python3
"""Render the README figures for the steady NACA 0012 case (offscreen).

Usage (after ./run_docker.sh):
    pvbatch paraview_images.py [--case airfoil.foam] [--out images]

The mesh images show the real cells (internalMesh seen from +z).
Writes mesh_domain.png, mesh_airfoil.png, pressure.png and velocity.png.
"""

import argparse
from pathlib import Path

from paraview.simple import (
    CellDatatoPointData,
    ColorBy,
    CreateView,
    GetColorTransferFunction,
    GetScalarBar,
    Hide,
    OpenFOAMReader,
    SaveScreenshot,
    Show,
    Slice,
    StreamTracer,
)

HERE = Path(__file__).resolve().parent

ap = argparse.ArgumentParser()
ap.add_argument("--case", default=str(HERE / "airfoil.foam"))
ap.add_argument("--out", default=str(HERE / "images"))
ap.add_argument("--size", nargs=2, type=int, default=[1200, 700])
ap.add_argument("--mesh-size", nargs=2, type=int, default=[960, 560],
                help="Size of the mesh images (keeps the PNGs small).")
args, _ = ap.parse_known_args()

out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)
size = args.size
mesh_size = args.mesh_size

# --- pipeline ---------------------------------------------------------------
flow = OpenFOAMReader(FileName=args.case)
flow.MeshRegions = ["internalMesh"]
flow.CellArrays = ["U", "p"]
flow.UpdatePipelineInformation()
t_last = flow.TimestepValues[-1] if flow.TimestepValues else 0.0
flow.UpdatePipeline(t_last)

pts = CellDatatoPointData(Input=flow)
mid = Slice(Input=pts)
mid.SliceType = "Plane"
mid.SliceType.Origin = [0.0, 0.0, 0.01]
mid.SliceType.Normal = [0.0, 0.0, 1.0]

wall = OpenFOAMReader(FileName=args.case)
wall.MeshRegions = ["patch/airfoil"]
wall.CellArrays = ["p"]
wall_line = Slice(Input=wall)
wall_line.SliceType = "Plane"
wall_line.SliceType.Origin = [0.0, 0.0, 0.01]
wall_line.SliceType.Normal = [0.0, 0.0, 1.0]

streams = StreamTracer(Input=pts, SeedType="Line")
streams.Vectors = ["POINTS", "U"]
streams.SeedType.Point1 = [-1.0, -0.5, 0.01]
streams.SeedType.Point2 = [-1.0, 0.5, 0.01]
streams.SeedType.Resolution = 25
streams.MaximumStreamlineLength = 6.0


def make_view(vsize=None):
    vsize = vsize or size
    view = CreateView("RenderView")
    view.ViewSize = vsize
    view.InteractionMode = "2D"
    view.CameraParallelProjection = 1
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.Background = [1.0, 1.0, 1.0]
    view.ViewTime = t_last
    return view


def camera(view, cx, cy, scale):
    view.CameraPosition = [cx, cy, 10.0]
    view.CameraFocalPoint = [cx, cy, 0.0]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelScale = scale


def add_wall(view, width=2.0):
    d = Show(wall_line, view)
    ColorBy(d, None)
    d.AmbientColor = [0.0, 0.0, 0.0]
    d.DiffuseColor = [0.0, 0.0, 0.0]
    d.LineWidth = width


def scalar_bar(view, display, array, title):
    lut = GetColorTransferFunction(array)
    display.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = title
    bar.ComponentTitle = ""
    bar.WindowLocation = "Any Location"
    bar.Position = [0.80, 0.15]
    bar.ScalarBarLength = 0.7
    bar.TitleColor = [0, 0, 0]
    bar.LabelColor = [0, 0, 0]
    return lut


# 1. mesh, whole domain -------------------------------------------------------
v = make_view(mesh_size)
d = Show(flow, v)
ColorBy(d, None)
d.Representation = "Surface With Edges"
d.AmbientColor = [0.2, 0.2, 0.2]
d.DiffuseColor = [0.85, 0.9, 1.0]
d.EdgeColor = [0.1, 0.1, 0.4]
add_wall(v, 1.5)
camera(v, 2.5, 0.0, 17.0)
SaveScreenshot(str(out / "mesh_domain.png"), v, ImageResolution=mesh_size)
Hide(flow, v)

# 2. mesh near the airfoil ------------------------------------------------------
v = make_view(mesh_size)
d = Show(flow, v)
ColorBy(d, None)
d.Representation = "Surface With Edges"
d.DiffuseColor = [0.85, 0.9, 1.0]
d.EdgeColor = [0.1, 0.1, 0.4]
add_wall(v, 2.0)
camera(v, 0.5, 0.0, 0.3)
SaveScreenshot(str(out / "mesh_airfoil.png"), v, ImageResolution=mesh_size)

# 3. pressure -------------------------------------------------------------------
v = make_view()
d = Show(mid, v)
ColorBy(d, ("POINTS", "p"))
lut = scalar_bar(v, d, "p", "p / rho  [m^2/s^2]")
lut.ApplyPreset("Cool to Warm", True)
lut.RescaleTransferFunction(-1000.0, 1000.0)
add_wall(v, 2.0)
camera(v, 0.5, 0.0, 0.5)
SaveScreenshot(str(out / "pressure.png"), v, ImageResolution=size)

# 4. velocity magnitude + streamlines -------------------------------------------
v = make_view()
d = Show(mid, v)
ColorBy(d, ("POINTS", "U", "Magnitude"))
lut = scalar_bar(v, d, "U", "|U|  [m/s]")
lut.ApplyPreset("Cool to Warm", True)
lut.RescaleTransferFunction(30.0, 70.0)
s = Show(streams, v)
ColorBy(s, None)
s.AmbientColor = [0.0, 0.0, 0.0]
s.DiffuseColor = [0.0, 0.0, 0.0]
s.LineWidth = 1.0
add_wall(v, 2.0)
camera(v, 1.0, 0.0, 1.2)
SaveScreenshot(str(out / "velocity.png"), v, ImageResolution=size)

mid.UpdatePipeline(t_last)
for name in ("p", "U"):
    print("%s range at t=%s: %s" % (name, t_last, mid.PointData[name].GetRange(-1 if name == "U" else 0)))
print("Wrote images to %s (time %s)" % (out, t_last))
