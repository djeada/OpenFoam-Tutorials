# Render the figures in images/ from a finished run (offscreen).
#   pvbatch make_images.py            (run from the case directory)
# Needs ParaView; the case must contain at least one result time directory.
import os
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()
reader = OpenFOAMReader(FileName="t_junction.foam")
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["T", "U", "k", "p_rgh"]
reader.UpdatePipelineInformation()
times = reader.TimestepValues
last = times[-1] if hasattr(times, "__len__") else times
reader.UpdatePipeline(last)
c2p = CellDatatoPointData(Input=reader)
view = CreateView("RenderView")
view.ViewSize = [1200, 520]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except Exception:
    pass
view.CameraPosition = [0.25, 0.1, 1.2]
view.CameraFocalPoint = [0.25, 0.1, 0.0]
view.CameraViewUp = [0, 1, 0]
view.CameraParallelProjection = 1
view.CameraParallelScale = 0.118
view.ViewTime = last


def show_field(src, name, rng, cmap, label, fname, comp="Magnitude", edges=False):
    for s in list(GetSources().values()):
        Hide(s, view)
    d = Show(src, view)
    d.Representation = "Surface With Edges" if edges else "Surface"
    ColorBy(d, ("POINTS", name, comp) if comp else ("POINTS", name))
    if True:
        lut = GetColorTransferFunction(name)
        lut.ApplyPreset(cmap, True)
        lut.RescaleTransferFunction(rng[0], rng[1])
        d.SetScalarBarVisibility(view, True)
        bar = GetScalarBar(lut, view)
        bar.Title = label
        bar.ComponentTitle = ""
        bar.TitleColor = [0, 0, 0]
        bar.LabelColor = [0, 0, 0]
        bar.WindowLocation = "Any Location"
        bar.Orientation = "Horizontal"
        bar.Position = [0.52, 0.72]
        bar.ScalarBarLength = 0.40
    Render(view)
    SaveScreenshot(fname, view, ImageResolution=[1200, 520])
    return d


# 1. mesh
d = Show(reader, view)
d.Representation = "Surface With Edges"
d.ColorArrayName = ["POINTS", ""]
d.AmbientColor = [0.1, 0.1, 0.1]
d.DiffuseColor = [0.85, 0.92, 1.0]
d.EdgeColor = [0.35, 0.35, 0.35]
Render(view)
SaveScreenshot("images/mesh.png", view, ImageResolution=[1200, 520])
Hide(reader, view)

# 2. temperature
d = show_field(c2p, "T", (300, 350), "Cool to Warm", "T [K]", "images/temperature.png", comp="")

# 3. velocity magnitude with streamlines
d = show_field(c2p, "U", (0, 1.6), "Viridis (matplotlib)", "|U| [m/s]", "images/velocity.png")
st = StreamTracer(Input=c2p, SeedType="Line")
st.Vectors = ["POINTS", "U"]
st.SeedType.Point1 = [0.001, 0.002, 0.005]
st.SeedType.Point2 = [0.001, 0.098, 0.005]
st.SeedType.Resolution = 14
st.MaximumStreamlineLength = 1.0
st2 = StreamTracer(Input=c2p, SeedType="Line")
st2.Vectors = ["POINTS", "U"]
st2.SeedType.Point1 = [0.127, 0.198, 0.005]
st2.SeedType.Point2 = [0.173, 0.198, 0.005]
st2.SeedType.Resolution = 8
st2.MaximumStreamlineLength = 1.0
# lift the lines in front of the (z = 0.01) surface so the camera sees them
tr = Transform(Input=MergeBlocks(Input=AppendDatasets(Input=[st, st2])))
tr.Transform.Translate = [0, 0, 0.01]
for s in (tr,):
    sd = Show(s, view)
    sd.ColorArrayName = ["POINTS", ""]
    sd.AmbientColor = [0, 0, 0]
    sd.DiffuseColor = [0, 0, 0]
    sd.LineWidth = 1.5
Render(view)
SaveScreenshot("images/velocity_streamlines.png", view, ImageResolution=[1200, 520])
Hide(tr, view)

# 4. turbulent kinetic energy
show_field(c2p, "k", (0, 0.18), "Inferno (matplotlib)", "k [m2/s2]", "images/k.png", comp="")

# 5. temperature profile across the outlet (cell-centre data, line at x = 0.499)
try:
    for s in list(GetSources().values()):
        Hide(s, view)
    pol = PlotOverLine(Input=c2p)
    pol.Point1 = [0.499, 0.0, 0.005]
    pol.Point2 = [0.499, 0.1, 0.005]
    pol.Resolution = 100
    chart = CreateView("XYChartView")
    chart.ViewSize = [1200, 520]
    chart.ViewTime = last
    cd = Show(pol, chart)
    cd.UseIndexForXAxis = 0
    cd.XArrayName = "Points_Y"
    cd.SeriesVisibility = ["T"]
    chart.LeftAxisTitle = "T [K]"
    chart.BottomAxisTitle = "y [m] at x = 0.499 m (outlet)"
    chart.LeftAxisRangeMinimum = 300.0
    chart.LeftAxisRangeMaximum = 355.0
    chart.LeftAxisUseCustomRange = 1
    Render(chart)
    SaveScreenshot("images/outlet_temperature_profile.png", chart, ImageResolution=[1200, 520])
except Exception as exc:
    print("outlet profile chart skipped:", exc)

# 6. residual history from log.buoyantSimpleFoam (initial residuals per iteration)
try:
    import re
    import tempfile
    log = open("log.buoyantSimpleFoam").read()
    fields = ["Ux", "h", "p_rgh", "k", "epsilon"]
    hist = {f: [] for f in fields}
    for block in log.split("\nTime = ")[1:]:
        for f, v in re.findall(r"Solving for (\w+), Initial residual = ([0-9.eE+-]+)", block):
            if f in hist:
                hist[f].append(float(v))
    n = min(len(v) for v in hist.values())
    csv = os.path.join(tempfile.mkdtemp(), "residuals.csv")
    with open(csv, "w") as fh:
        fh.write("iteration," + ",".join(fields) + "\n")
        for i in range(n):
            fh.write("%d," % (i + 1) + ",".join("%.6e" % hist[f][i] for f in fields) + "\n")
    for s in list(GetSources().values()):
        Hide(s, view)
    rd = CSVReader(FileName=[csv])
    rview = CreateView("XYChartView")
    rview.ViewSize = [1200, 520]
    rdisp = Show(rd, rview)
    rdisp.UseIndexForXAxis = 0
    rdisp.XArrayName = "iteration"
    rdisp.SeriesVisibility = fields
    rview.LeftAxisLogScale = 1
    rview.LeftAxisTitle = "initial residual"
    rview.BottomAxisTitle = "iteration"
    rview.LeftAxisUseCustomRange = 1
    rview.LeftAxisRangeMinimum = 1e-5
    rview.LeftAxisRangeMaximum = 2.0
    rview.ChartTitle = "Initial residuals (residualControl: 1e-4)"
    Render(rview)
    SaveScreenshot("images/residuals.png", rview, ImageResolution=[1200, 520])
except Exception as exc:
    print("residual chart skipped:", exc)
