#!/usr/bin/env pvbatch
"""Create the README screenshots from a finished run (offscreen ParaView).

    pvbatch make_images.py            # run inside the project directory

Reads the OpenFOAM case in this directory (needs time directories and
log.interFoam) and writes PNG files to images/.
"""
import csv
import os
import re

from paraview.simple import *  # noqa: F401,F403

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "images")
os.makedirs(OUT, exist_ok=True)
Z_EQ = 0.01596   # equilibrium z of the centre of mass (see verify.py)

# ------------------------------------------------------------- motion data
NUM = r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?"
rows = []
t = None
z = None
v = None
with open(os.path.join(HERE, "log.interFoam"), errors="replace") as f:
    for line in f:
        m = re.match(r"^Time = (" + NUM + ")", line)
        if m:
            t = float(m.group(1))
            continue
        m = re.search(r"Centre of mass: \(\S+ \S+ (" + NUM + r")\)", line)
        if m:
            z = float(m.group(1))
            continue
        m = re.search(r"Linear velocity: \(\S+ \S+ (" + NUM + r")\)", line)
        if m:
            v = float(m.group(1))
            rows.append((t, z, v))          # last outer corrector wins below
last = {}
for r in rows:
    last[r[0]] = r
data = [last[k] for k in sorted(last)]
csv_path = os.path.join(os.environ.get("TMPDIR", "/tmp"), "boat_hull_heave.csv")
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["time", "z_CoM", "v_z", "z_equilibrium"])
    for r in data:
        w.writerow([r[0], r[1], r[2], Z_EQ])

# ------------------------------------------------------------ heave chart
rd = CSVReader(FileName=[csv_path])
lay = CreateLayout("heave")
v1 = CreateView("XYChartView")
v2 = CreateView("XYChartView")
AssignViewToLayout(view=v1, layout=lay, hint=0)
lay.SplitHorizontal(0, 0.5)
AssignViewToLayout(view=v2, layout=lay, hint=2)
for vw, cols, ylab in ((v1, ["z_CoM", "z_equilibrium"], "z of centre of mass [m]"),
                       (v2, ["v_z"], "vertical velocity [m/s]")):
    d = Show(rd, vw)
    d.UseIndexForXAxis = 0
    d.XArrayName = "time"
    d.SeriesVisibility = cols
    for c in cols:
        d.SeriesPlotCorner = [c, "0"]
    vw.BottomAxisTitle = "time [s]"
    vw.LeftAxisTitle = ylab
    vw.ShowLegend = 1 if len(cols) > 1 else 0
lay.SetSize(1500, 520)
SaveScreenshot(os.path.join(OUT, "heave.png"), lay, ImageResolution=[1500, 520])
for s in (v1, v2):
    Delete(s)
Delete(lay)

# ------------------------------------------------------------------- case
foam = os.path.join(HERE, "boat_hull_floating.foam")
open(foam, "w").close()
case = OpenFOAMReader(FileName=foam)
case.MeshRegions = ["internalMesh"]
case.CellArrays = ["alpha.water", "U", "p_rgh"]
case.UpdatePipelineInformation()
times = list(case.TimestepValues) if case.TimestepValues else [0.0]
if not isinstance(times, list):
    times = [times]
pick = [t for t in times if t <= 0.2 + 1e-9] or times   # 0.05 ... 0.2 s (t = 0 has no water field)
print("time steps:", times)

# Hull patch (separate reader instance so that the patch has its own pipeline)
hull = OpenFOAMReader(FileName=foam)
avail = list(hull.GetProperty("MeshRegions").Available)
patch = [a for a in avail if a.endswith("boat")]
hull.MeshRegions = patch
hull.UpdatePipelineInformation()
print("hull regions:", patch)

RES = [900, 430]
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()


def side_view(time, name):
    """y = 0 mid-plane coloured by alpha.water with the hull outline."""
    scene.AnimationTime = time
    vw = CreateView("RenderView")
    vw.ViewSize = RES
    vw.OrientationAxesVisibility = 0
    vw.ViewTime = time
    vw.UseColorPaletteForBackground = 0
    vw.Background = [1.0, 1.0, 1.0]
    sl = Slice(Input=case)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, 0.0, 0.0]
    sl.SliceType.Normal = [0.0, 1.0, 0.0]
    sd = Show(sl, vw)
    ColorBy(sd, ("CELLS", "alpha.water"))
    lut = GetColorTransferFunction("alpha.water")
    # air (alpha = 0) white, water (alpha = 1) blue
    lut.RGBPoints = [0.0, 1.0, 1.0, 1.0, 1.0, 0.15, 0.4, 0.85]
    lut.ColorSpace = "RGB"
    lut.RescaleTransferFunction(0.0, 1.0)
    sd.SetScalarBarVisibility(vw, False)
    # hull outline: the boat patch cut by the same plane
    hs = Slice(Input=hull)
    hs.SliceType = "Plane"
    hs.SliceType.Origin = [0.0, 0.0, 0.0]
    hs.SliceType.Normal = [0.0, 1.0, 0.0]
    hd = Show(hs, vw)
    hd.ColorArrayName = ["POINTS", ""]
    hd.AmbientColor = [0.9, 0.2, 0.0]
    hd.DiffuseColor = [0.9, 0.2, 0.0]
    hd.LineWidth = 4.0
    txt = Text(Text="alpha.water on y = 0,  t = %.3g s" % time)
    td = Show(txt, vw)
    td.FontSize = 16
    td.Color = [0.0, 0.0, 0.0]
    td.WindowLocation = "Upper Center"
    Render(vw)                  # first render resets the camera, so set it afterwards
    vw.InteractionMode = "2D"
    vw.CameraPosition = [0.0, -2.0, 0.015]
    vw.CameraFocalPoint = [0.0, 0.0, 0.015]
    vw.CameraViewUp = [0.0, 0.0, 1.0]
    vw.CameraParallelProjection = 1
    vw.CameraParallelScale = 0.065
    Render(vw)
    SaveScreenshot(os.path.join(OUT, name), vw, ImageResolution=RES)
    for src in (txt, hs, sl):
        Delete(src)
    Delete(vw)


for k, tm in enumerate(pick):
    side_view(tm, "free_surface_%d.png" % k)

# ------------------------------------------------------- hull and mesh (t=0)
lay2 = CreateLayout("mesh")
vm = CreateView("RenderView")
AssignViewToLayout(view=vm, layout=lay2, hint=0)
vm.OrientationAxesVisibility = 0
vm.ViewTime = times[0]   # first written time (0.05 s): the hull has hardly moved yet
sl = Slice(Input=case)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.0, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]
vm.UseColorPaletteForBackground = 0
vm.Background = [1.0, 1.0, 1.0]
sd = Show(sl, vm)
sd.ColorArrayName = ["POINTS", ""]
sd.Representation = "Surface With Edges"
sd.AmbientColor = [0.3, 0.3, 0.3]
sd.DiffuseColor = [0.85, 0.9, 1.0]
hd = Show(hull, vm)
hd.ColorArrayName = ["POINTS", ""]
hd.DiffuseColor = [0.9, 0.35, 0.1]
vm.InteractionMode = "3D"
vm.CameraPosition = [0.45, -0.55, 0.30]
vm.CameraFocalPoint = [0.02, 0.0, -0.01]
vm.CameraViewUp = [0.0, 0.0, 1.0]
vm.CameraParallelProjection = 0
vm.CameraViewAngle = 22
Render(vm)
lay2.SetSize(900, 490)
SaveScreenshot(os.path.join(OUT, "hull_mesh.png"), lay2, ImageResolution=[900, 490])
print("done")
