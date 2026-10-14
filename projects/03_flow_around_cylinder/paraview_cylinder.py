"""ParaView Python script — 2D cylinder vortex shedding visualisation.

Interactive (default): shows the Z-component of vorticity on the last time
step with a diverging colour map, the cylinder as a dark disk.

    paraview --script paraview_cylinder.py

Batch, single screenshot of the vorticity field:

    pvbatch paraview_cylinder.py --screenshot vorticity.png

Batch, the full image set used in the README (mesh, velocity, pressure,
vorticity, Cd/Cl history) written to a directory:

    pvbatch paraview_cylinder.py --images images

(add `--only mesh velocity pressure vorticity cd_cl` to select images).

Requires cylinder.foam (created by ./Allrun) and a finished or running case.
"""

import csv
import os
import sys
import tempfile

from paraview.simple import (
    AnnotateTimeFilter,
    CellDatatoPointData,
    ColorBy,
    CreateView,
    Cylinder,
    CSVReader,
    Delete,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetScalarBar,
    Gradient,
    Hide,
    OpenFOAMReader,
    Render,
    SaveScreenshot,
    Show,
    Slice,
    Transform,
)

# ── Configuration ─────────────────────────────────────────────────────────
CASE_DIR = os.path.dirname(os.path.abspath(__file__))
FOAM_FILE = os.path.join(CASE_DIR, "cylinder.foam")
FORCE_FILE = os.path.join(CASE_DIR, "postProcessing", "forceCoeffs", "0",
                          "forceCoeffs.dat")
SPAN = 0.1                     # must match --span of generate_cylinder_geo.py
RADIUS = 0.5                   # cylinder radius (D = 1)
VORTICITY_RANGE = [-2.0, 2.0]  # 1/s; the shear layers saturate, the street is visible
SPEED_RANGE = [0.0, 1.6]       # m/s
PRESSURE_RANGE = [-0.9, 0.6]   # m^2/s^2 (kinematic)
IMAGE_SIZE = [1200, 560]
DARK = [0.02, 0.02, 0.07]
LIGHT_TEXT = [0.9, 0.9, 0.9]

# ── Arguments ─────────────────────────────────────────────────────────────
screenshot_path = None
images_dir = None
if "--screenshot" in sys.argv:
    i = sys.argv.index("--screenshot")
    screenshot_path = sys.argv[i + 1] if i + 1 < len(sys.argv) else "vorticity.png"
if "--images" in sys.argv:
    i = sys.argv.index("--images")
    images_dir = sys.argv[i + 1] if i + 1 < len(sys.argv) else "images"
only = set()
if "--only" in sys.argv:   # e.g. --only mesh cd_cl  (with --images)
    only = set(a for a in sys.argv[sys.argv.index("--only") + 1:]
               if not a.startswith("--"))
batch = screenshot_path is not None or images_dir is not None

# ── Pipeline ──────────────────────────────────────────────────────────────
reader = OpenFOAMReader(FileName=FOAM_FILE)
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["U", "p"]
reader.UpdatePipelineInformation()

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
times = list(reader.TimestepValues) if reader.TimestepValues else []
if times:
    scene.AnimationTime = times[-1]
t_last = times[-1] if times else 0.0

c2p = CellDatatoPointData(Input=reader)
c2p.CellDataArraytoprocess = ["U", "p"]

# Mid-plane slice (the case is one cell thick in z); keep the quad cells
zslice = Slice(Input=c2p)
zslice.SliceType = "Plane"
zslice.SliceType.Origin = [0.0, 0.0, SPAN / 2]
zslice.SliceType.Normal = [0.0, 0.0, 1.0]
zslice.Triangulatetheslice = 0

# Vorticity from the point-data velocity (computed on the 3D mesh, then sliced)
grad3d = Gradient(Input=c2p)
grad3d.ScalarArray = ["POINTS", "U"]
grad3d.ComputeVorticity = 1
grad3d.VorticityArrayName = "Vorticity"
grad3d.ComputeGradient = 0

grad = Slice(Input=grad3d)
grad.SliceType = "Plane"
grad.SliceType.Origin = [0.0, 0.0, SPAN / 2]
grad.SliceType.Normal = [0.0, 0.0, 1.0]
grad.Triangulatetheslice = 0

# Cylinder drawn as a dark disk (a solid body inside the hole of the mesh)
# (ParaView's Cylinder source has its axis along y: rotate it to z. It is
# centred at z = 0 and 1.2 * SPAN tall, so it covers the slice at SPAN / 2.)
body = Cylinder(Radius=RADIUS, Height=SPAN * 1.2, Resolution=96)
body.Capping = 1
body_t = Transform(Input=body)
body_t.Transform = "Transform"
body_t.Transform.Rotate = [90.0, 0.0, 0.0]


def make_view(size):
    v = CreateView("RenderView") if batch else GetActiveViewOrCreate("RenderView")
    v.ViewSize = size
    v.UseColorPaletteForBackground = 0
    v.Background = DARK
    v.OrientationAxesVisibility = 0
    v.InteractionMode = "2D"
    v.CameraParallelProjection = 1
    return v


def camera(v, cx, cy, half_height):
    v.CameraPosition = [cx, cy, 20.0]
    v.CameraFocalPoint = [cx, cy, SPAN / 2]
    v.CameraViewUp = [0.0, 1.0, 0.0]
    v.CameraParallelScale = half_height


def show_body(v):
    d = Show(body_t, v)
    d.Representation = "Surface"
    d.ColorArrayName = ["POINTS", ""]
    d.AmbientColor = d.DiffuseColor = [0.35, 0.35, 0.35]
    return d


def scalar_bar(lut, v, title):
    lut.ApplyPreset("Cool to Warm", True)
    sb = GetScalarBar(lut, v)
    sb.Title = title
    sb.ComponentTitle = ""
    sb.TitleColor = LIGHT_TEXT
    sb.LabelColor = LIGHT_TEXT
    sb.TitleFontSize = 16
    sb.LabelFontSize = 14
    return sb


def time_label(v, source):
    ann = AnnotateTimeFilter(Input=source)
    ann.Format = "t = {time:.1f} s"
    d = Show(ann, v)
    d.FontSize = 14
    d.Color = LIGHT_TEXT
    d.WindowLocation = "Upper Right Corner"
    return ann


def field_view(source, array, assoc, comp, rng, title, cx=6.0, cy=0.0,
               half_height=4.0):
    v = make_view(IMAGE_SIZE)
    camera(v, cx, cy, half_height)
    d = Show(source, v)
    d.Representation = "Surface"
    if comp is None:
        ColorBy(d, (assoc, array))
    else:
        ColorBy(d, (assoc, array, comp))
    lut = GetColorTransferFunction(array)
    lut.VectorMode = "Magnitude" if comp is None else "Component"
    if comp is not None:
        lut.VectorComponent = {"X": 0, "Y": 1, "Z": 2}[comp]
    lut.RescaleTransferFunction(*rng)
    scalar_bar(lut, v, title)
    d.SetScalarBarVisibility(v, True)
    show_body(v)
    time_label(v, source)
    return v


def force_plot(path_out):
    """Cd and Cl versus time from forceCoeffs.dat (needs the case results)."""
    rows = []
    with open(FORCE_FILE) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            rows.append((float(c[0]), float(c[2]), float(c[3])))
    tmp = tempfile.mkdtemp()
    csv_path = os.path.join(tmp, "forces.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Time", "Cd", "Cl"])
        w.writerows(rows)
    src = CSVReader(FileName=[csv_path])
    v = CreateView("XYChartView")
    v.ViewSize = IMAGE_SIZE
    d = Show(src, v, "XYChartRepresentation")
    d.UseIndexForXAxis = 0
    d.XArrayName = "Time"
    d.SeriesVisibility = ["Cd", "Cl"]
    d.SeriesColor = ["Cd", "0.85", "0.2", "0.2", "Cl", "0.1", "0.35", "0.8"]
    d.SeriesLineThickness = ["Cd", "2", "Cl", "2"]
    v.BottomAxisTitle = "Time [s]"
    v.LeftAxisTitle = "Cd, Cl [-]"
    v.ChartTitle = "Drag and lift coefficients"
    v.ShowLegend = 1
    Render(v)
    SaveScreenshot(path_out, v, ImageResolution=IMAGE_SIZE)
    Delete(src)


# ── Batch: image set ──────────────────────────────────────────────────────
if images_dir:
    os.makedirs(images_dir, exist_ok=True)

    def out(name):
        return os.path.join(images_dir, name)

    # Mesh: whole domain and zoom on the cylinder (cells + boundary layer)
    def want(key):
        return not only or key in only

    for name, cx, cy, hh, size in (
            ("mesh_domain.png", 8.5, 0.0, 8.2, [800, 392]),
            ("mesh_cylinder.png", 0.0, 0.0, 0.9, [1000, 500])):
        if not want("mesh"):
            break
        v = make_view(size)
        v.Background = [1, 1, 1]
        d = Show(zslice, v)
        Render(v)          # the first render of a new view resets its camera
        camera(v, cx, cy, hh)
        d.Representation = "Surface With Edges"
        d.ColorArrayName = ["POINTS", ""]
        d.AmbientColor = d.DiffuseColor = [0.85, 0.9, 0.95]
        d.EdgeColor = [0.1, 0.1, 0.1]
        Render(v)
        SaveScreenshot(out(name), v, ImageResolution=size)
        Hide(zslice, v)

    if want("velocity"):
        v = field_view(zslice, "U", "POINTS", None, SPEED_RANGE, "|U| [m/s]")
        Render(v)
        SaveScreenshot(out("velocity_magnitude.png"), v,
                       ImageResolution=IMAGE_SIZE)

    if want("pressure"):
        v = field_view(zslice, "p", "POINTS", None, PRESSURE_RANGE,
                       "p/rho [m^2/s^2]")
        Render(v)
        SaveScreenshot(out("pressure.png"), v, ImageResolution=IMAGE_SIZE)

    if want("vorticity"):
        v = field_view(grad, "Vorticity", "POINTS", "Z", VORTICITY_RANGE,
                       "Vorticity z [1/s]")
        Render(v)
        SaveScreenshot(out("vorticity.png"), v, ImageResolution=IMAGE_SIZE)

    if want("cd_cl") and os.path.isfile(FORCE_FILE):
        force_plot(out("cd_cl_history.png"))
    print(f"Wrote images to {images_dir} (t = {t_last})")

# ── Batch: single vorticity screenshot / interactive ──────────────────────
else:
    v = field_view(grad, "Vorticity", "POINTS", "Z", VORTICITY_RANGE,
                   "Vorticity z [1/s]")
    Render(v)
    if screenshot_path:
        SaveScreenshot(screenshot_path, v, ImageResolution=IMAGE_SIZE)
        print(f"Saved screenshot to {screenshot_path}")
