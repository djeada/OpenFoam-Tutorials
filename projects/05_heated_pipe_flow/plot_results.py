# Generate the README figures in images/ from a finished run.
# Usage (case directory, after ./Allrun):
#   buoyantSimpleFoam -postProcess -func wallHeatFlux -latestTime
#   pvbatch plot_results.py
# Field plots are stretched 5x in y so the 0.1 m high channel is visible.
import os
import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from vtk.numpy_interface import dataset_adapter as dsa

os.chdir(os.path.dirname(os.path.abspath(__file__)))
open("heated_pipe.foam", "a").close()
os.makedirs("images", exist_ok=True)

rdr = OpenFOAMReader(FileName="heated_pipe.foam")
rdr.MeshRegions = ["internalMesh"]
rdr.CellArrays = ["T", "U", "wallHeatFlux"]
last = max(rdr.TimestepValues)
rdr.UpdatePipeline(last)
GetAnimationScene().AnimationTime = last
print("time", last)

STRETCH = 5.0
SIZE = [1600, 560]
stretched = Transform(Input=rdr)
stretched.Transform.Scale = [1.0, STRETCH, 1.0]


def setup_view():
    v = CreateView("RenderView")
    v.ViewSize = SIZE
    v.OrientationAxesVisibility = 0
    v.UseColorPaletteForBackground = 0
    v.Background = [1, 1, 1]
    v.InteractionMode = "2D"
    return v


def set_camera(v):
    # set after the first render, which would otherwise reset the camera
    v.CameraPosition = [1.0, 0.05 * STRETCH, 6.0]
    v.CameraFocalPoint = [1.0, 0.05 * STRETCH, 0.0]
    v.CameraViewUp = [0, 1, 0]
    v.CameraParallelProjection = 1
    v.CameraParallelScale = 0.40


def field_view(name, fname, comp, title, preset, rng=None):
    v = setup_view()
    d = Show(stretched, v)
    d.Representation = "Surface"
    ColorBy(d, ("CELLS", name, comp) if comp else ("CELLS", name))
    lut = GetColorTransferFunction(name)
    lut.ApplyPreset(preset, True)
    if rng:
        lut.RescaleTransferFunction(*rng)
    else:
        d.RescaleTransferFunctionToDataRange(True, False)
    d.SetScalarBarVisibility(v, True)
    sb = GetScalarBar(lut, v)
    sb.Title = title
    sb.ComponentTitle = ""
    sb.Orientation = "Horizontal"
    sb.WindowLocation = "Any Location"
    sb.Position = [0.30, 0.03]
    sb.ScalarBarLength = 0.4
    sb.TitleColor = [0, 0, 0]
    sb.LabelColor = [0, 0, 0]
    sb.TitleFontSize = 20
    sb.LabelFontSize = 18
    sb.RangeLabelFormat = "%-#6.3g"
    Render(v)
    set_camera(v)
    Render(v)
    SaveScreenshot(fname, v, ImageResolution=SIZE)
    Delete(v)


field_view("T", "images/temperature.png", None, "T [K]", "Cool to Warm", [300, 350])
field_view("U", "images/velocity.png", "Magnitude", "|U| [m/s]", "Viridis (matplotlib)")

# mesh
v = setup_view()
d = Show(stretched, v)
d.Representation = "Surface With Edges"
d.ColorArrayName = [None, ""]
d.AmbientColor = [0.85, 0.92, 1.0]
d.DiffuseColor = [0.85, 0.92, 1.0]
d.EdgeColor = [0.1, 0.1, 0.1]
Render(v)
set_camera(v)
Render(v)
SaveScreenshot("images/mesh.png", v, ImageResolution=SIZE)
Delete(v)


# ---- numerical data for line charts -------------------------------------
def fetch_cells(src):
    c = MergeBlocks(Input=CellCenters(Input=src))
    c.UpdatePipeline(last)
    o = dsa.WrapDataObject(sm.Fetch(c))
    return np.array(o.Points), {k: np.array(o.PointData[k]) for k in o.PointData.keys()}


pts, f = fetch_cells(rdr)
Tf = f["T"]
Um = np.linalg.norm(f["U"], axis=1)
cols = [0, 25, 50, 99]
xs = sorted(set(np.round(pts[:, 0], 6)))
sel = [xs[i] for i in cols]


def column(x):
    m = np.isclose(pts[:, 0], x)
    o = np.argsort(pts[m, 1])
    return m, o


def table(columns):
    names = list(columns)
    n = len(columns[names[0]])
    script = (
        "import numpy as np\n"
        "from vtkmodules.vtkCommonCore import vtkFloatArray\n"
        "t = self.GetOutput()\n"
        "data = %r\n"
        "for k, v in data.items():\n"
        "    a = vtkFloatArray(); a.SetName(k)\n"
        "    for x in v: a.InsertNextValue(x)\n"
        "    t.AddColumn(a)\n" % {k: [float(x) for x in columns[k]] for k in names}
    )
    p = ProgrammableSource()
    p.OutputDataSetType = "vtkTable"
    p.Script = script
    p.UpdatePipeline()
    return p


def chart(fname, src, xname, series, ytitle, xtitle, size=(1000, 600)):
    cv = CreateView("XYChartView")
    cv.ViewSize = list(size)
    cv.ChartTitle = ""
    cv.LeftAxisTitle = ytitle
    cv.BottomAxisTitle = xtitle
    cv.LegendLocation = "TopRight"
    d = Show(src, cv)
    d.UseIndexForXAxis = 0
    d.XArrayName = xname
    d.SeriesVisibility = series
    Render(cv)
    SaveScreenshot(fname, cv, ImageResolution=list(size))
    Delete(cv)


m0, o0 = column(sel[0])
ycoord = pts[m0, 1][o0]
tcols = {"y": ycoord}
ucols = {"y": ycoord}
tnames = []
unames = []
for x in sel:
    m, o = column(x)
    lbl = "x = %.2f m" % x
    tcols[lbl] = Tf[m][o]
    ucols[lbl] = Um[m][o]
    tnames.append(lbl)
chart("images/profile_T.png", table(tcols), "y", tnames, "T [K]", "y [m]")
chart("images/profile_U.png", table(ucols), "y", tnames, "|U| [m/s]", "y [m]")

# wall heat flux (positive: heat into the fluid) and wall-adjacent T
w = OpenFOAMReader(FileName="heated_pipe.foam")
w.MeshRegions = ["patch/bottomWall", "patch/topWall"]
w.CellArrays = ["wallHeatFlux"]
w.UpdatePipeline(last)
wp, wf = fetch_cells(w)
q = wf["wallHeatFlux"]
bot = np.isclose(wp[:, 1], 0.0)
top = np.isclose(wp[:, 1], 0.1)
ob = np.argsort(wp[bot, 0])
ot = np.argsort(wp[top, 0])
chart(
    "images/wall_heat_flux.png",
    table({"x": wp[bot, 0][ob], "bottomWall": q[bot][ob], "topWall": q[top][ot]}),
    "x",
    ["bottomWall", "topWall"],
    "wall heat flux [W/m^2]",
    "x [m]",
)
print("q integral (W, 0.01 m slab):",
      q[bot].sum() * 0.02 * 0.01, q[top].sum() * 0.02 * 0.01)
