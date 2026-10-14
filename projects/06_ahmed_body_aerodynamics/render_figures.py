#!/usr/bin/env pvbatch
"""Render the figures used in README.md from a finished (or partial) run.

Usage (offscreen, from the case directory):
    pvbatch render_figures.py

Reads the case through the OpenFOAM reader (reconstructed case if a result
time exists in the case directory, otherwise the decomposed processor*
directories) and writes PNG files to images/. Needs only ParaView.
"""
import glob
import os
import re

from paraview.simple import *  # noqa: F401,F403

CASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE, "images")
os.makedirs(OUT, exist_ok=True)
U_INF = 40.0
CP_DEN = 0.5 * U_INF * U_INF   # kinematic p -> Cp

foam = os.path.join(CASE, "ahmed_body.foam")
open(foam, "a").close()


def numeric_times(pattern):
    out = []
    for d in glob.glob(pattern):
        name = os.path.basename(d)
        if re.fullmatch(r"[0-9.]+(e[-+]?[0-9]+)?", name):
            out.append(float(name))
    return sorted(out)


recon = [t for t in numeric_times(os.path.join(CASE, "*")) if t > 0]
decomp = [t for t in numeric_times(os.path.join(CASE, "processor0", "*")) if t > 0]
case_type = "Reconstructed Case" if recon else "Decomposed Case"
times = recon or decomp
print("case type: %s, result times: %s" % (case_type, times))

reader = OpenFOAMReader(FileName=foam)
reader.CaseType = case_type
reader.MeshRegions = ["internalMesh", "patch/ahmed_body"]
reader.CellArrays = ["U", "p", "k", "omega", "nut"]
reader.UpdatePipelineInformation()
if times:
    reader.UpdatePipeline(max(times))
    last = max(times)
else:
    last = 0.0
scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
scene.AnimationTime = last

view = CreateView("RenderView")
view.ViewSize = [1200, 700]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except Exception:
    pass


def save(name, size=(1200, 700)):
    Render(view)
    SaveScreenshot(os.path.join(OUT, name), view, ImageResolution=list(size))
    print("wrote", name)


def hide_all():
    for src in list(GetSources().values()):
        try:
            Hide(src, view)
        except Exception:
            pass


def style_text(disp, color=(0, 0, 0)):
    pass


# split the multiblock into internal mesh and body patch
internal = ExtractBlock(Input=reader, Selectors=["/Root/internalMesh"])
body = ExtractBlock(Input=reader, Selectors=["/Root/boundary/ahmed_body"])
body_surf = ExtractSurface(Input=MergeBlocks(Input=body))

# ---- 1. body surface mesh ----------------------------------------------
hide_all()
d = Show(body_surf, view)
d.Representation = "Surface With Edges"
ColorBy(d, None)
d.AmbientColor = d.DiffuseColor = [0.75, 0.78, 0.82]
d.EdgeColor = [0.1, 0.1, 0.1]
Render(view)   # the first render resets the camera
view.CameraPosition = [1.05, -0.42, 0.42]
view.CameraFocalPoint = [0.90, -0.19, 0.30]
view.CameraViewUp = [0, 0, 1]
view.CameraParallelProjection = 0
save("mesh_body_surface.png")

# ---- 2. mesh in the centre plane --------------------------------------
cut = Slice(Input=MergeBlocks(Input=internal))
cut.SliceType = "Plane"
cut.SliceType.Origin = [0.5, 0.0, 0.0]
cut.SliceType.Normal = [0.0, 1.0, 0.0]
hide_all()
d = Show(cut, view)
d.Representation = "Surface With Edges"
ColorBy(d, None)
d.AmbientColor = d.DiffuseColor = [0.85, 0.9, 1.0]
d.EdgeColor = [0.15, 0.15, 0.35]
b = Show(body_surf, view)
b.Representation = "Surface"
ColorBy(b, None)
b.AmbientColor = b.DiffuseColor = [0.8, 0.2, 0.2]
view.CameraParallelProjection = 1
view.CameraPosition = [0.5, -10.0, 0.5]
view.CameraFocalPoint = [0.5, 0.0, 0.5]
view.CameraViewUp = [0, 0, 1]
view.CameraParallelScale = 0.45
view.CameraFocalPoint = [0.75, 0.0, 0.3]
view.CameraPosition = [0.75, -10.0, 0.3]
save("mesh_centreplane.png", size=(900, 520))
Hide(b, view)

if not times:
    print("no result times found: only mesh figures written")
    raise SystemExit(0)

# ---- 3. velocity magnitude in the centre plane -------------------------
hide_all()
d = Show(cut, view)
d.Representation = "Surface"
ColorBy(d, ("CELLS", "U", "Magnitude"))
lut = GetColorTransferFunction("U")
lut.ApplyPreset("Turbo", True)
lut.RescaleTransferFunction(0.0, 55.0)
d.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = "|U| [m/s]"
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]
b = Show(body_surf, view)
b.Representation = "Surface"
ColorBy(b, None)
b.AmbientColor = b.DiffuseColor = [0.25, 0.25, 0.25]
view.CameraParallelProjection = 1
view.CameraPosition = [1.2, -10.0, 0.45]
view.CameraFocalPoint = [1.2, 0.0, 0.45]
view.CameraParallelScale = 1.1
save("velocity_centreplane.png")

# ---- 4. surface pressure coefficient on the body -----------------------
calc = Calculator(Input=body_surf)
calc.AttributeType = "Cell Data"
calc.Function = "p/%g" % CP_DEN
calc.ResultArrayName = "Cp"
hide_all()
d = Show(calc, view)
d.Representation = "Surface"
ColorBy(d, ("CELLS", "Cp"))
lut = GetColorTransferFunction("Cp")
lut.ApplyPreset("Cool to Warm", True)
lut.RescaleTransferFunction(-1.0, 1.0)
d.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = "Cp [-]"
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]
view.CameraParallelProjection = 0
view.CameraFocalPoint = [0.5, 0.0, 0.2]
view.CameraViewUp = [0, 0, 1]
view.CameraPosition = [-1.6, -1.4, 1.0]
save("cp_body_front.png")
view.CameraPosition = [2.2, -1.6, 1.1]
save("cp_body_rear.png")

# ---- 5. Cd history ------------------------------------------------------
dat = sorted(glob.glob(os.path.join(CASE, "postProcessing", "forceCoeffs", "*", "*.dat")))
if dat:
    names, rows = None, []
    for line in open(dat[0]):
        if line.startswith("#"):
            toks = line.lstrip("# ").split()
            if toks and toks[0] == "Time":
                names = toks
            continue
        if line.strip():
            rows.append(line.split())
    if names and "Cd" in names:
        k = names.index("Cd")
        csv = os.path.join(OUT, "_cd.csv")
        with open(csv, "w") as f:
            f.write("iteration,Cd\n")
            for r in rows:
                if float(r[0]) >= 25:   # skip the start-up transient
                    f.write("%s,%s\n" % (r[0], r[k]))
        tbl = CSVReader(FileName=csv)
        cv = CreateView("XYChartView")
        cv.ViewSize = [1200, 600]
        cv.BottomAxisTitle = "Iteration"
        cv.LeftAxisTitle = "Cd"
        cv.ChartTitle = "Drag coefficient history"
        try:
            cv.ChartTitleColor = [0, 0, 0]
        except Exception:
            pass
        td = Show(tbl, cv)
        td.UseIndexForXAxis = 0
        td.XArrayName = "iteration"
        td.SeriesVisibility = ["Cd"]
        td.SeriesColor = ["Cd", "0.1", "0.2", "0.7"]
        Render(cv)
        SaveScreenshot(os.path.join(OUT, "cd_history.png"), cv, ImageResolution=[1200, 600])
        os.remove(csv)
        print("wrote cd_history.png")
