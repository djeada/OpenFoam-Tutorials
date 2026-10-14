"""ParaView Python script for the OpenFOAM elbow tutorial.

Run from the case directory, after the case has been solved and elbow.foam
exists (open_paraview.sh creates it).

Default: show velocity magnitude at the last time step on a z-midplane slice.

    paraview --script paraview_elbow.py
    pvbatch paraview_elbow.py [--screenshot velocity.png]

Write the figures used in the README (mesh, velocity, pressure, streamlines):

    pvbatch paraview_elbow.py --images images
"""

import os
import sys

from paraview.simple import (
    CellDatatoPointData,
    ColorBy,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetScalarBar,
    Hide,
    HideAll,
    OpenFOAMReader,
    Render,
    SaveScreenshot,
    Show,
    Slice,
    StreamTracer,
)


def get_arg(name):
    """Return the value following a command-line option, or None."""
    if name in sys.argv:
        idx = sys.argv.index(name)
        if idx + 1 < len(sys.argv):
            return sys.argv[idx + 1]
    return None


screenshot_file = get_arg("--screenshot")
images_dir = get_arg("--images")
SIZE = [1200, 1100]

# ── Open the case and go to the last time step ────────────────────────
reader = OpenFOAMReader(FileName="elbow.foam")
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["U", "p"]

animation = GetAnimationScene()
animation.UpdateAnimationUsingDataTimeSteps()
animation.GoToLast()

# ── Cell data -> point data for smooth rendering ──────────────────────
c2p = CellDatatoPointData(Input=reader)
c2p.CellDataArraytoprocess = ["U", "p"]

# ── Z-normal slice through the midplane (2-D case) ────────────────────
slice1 = Slice(Input=c2p)
slice1.SliceType = "Plane"
slice1.SliceType.Origin = [0.0, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.SliceOffsetValues = [0.0]

# ── Render view: white background, 2-D parallel projection ────────────
view = GetActiveViewOrCreate("RenderView")
view.ViewSize = SIZE
if hasattr(view, "UseColorPaletteForBackground"):  # ParaView >= 5.10
    view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.InteractionMode = "2D"
view.CameraPosition = [0.0, 0.0, 1.0]
view.CameraFocalPoint = [0.0, 0.0, 0.0]
view.CameraViewUp = [0.0, 1.0, 0.0]
view.CameraParallelProjection = 1
view.OrientationAxesVisibility = 0 if images_dir else 1


def show_scalar(source, array, title, preset="Cool to Warm", component="Magnitude"):
    """Show `source` coloured by a point array and return its display."""
    display = Show(source, view)
    if component:
        ColorBy(display, ("POINTS", array, component))
    else:
        ColorBy(display, ("POINTS", array))
    lut = GetColorTransferFunction(array)
    lut.ApplyPreset(preset, True)
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = title
    bar.ComponentTitle = ""
    bar.TitleColor = [0.0, 0.0, 0.0]
    bar.LabelColor = [0.0, 0.0, 0.0]
    return display


def save(name):
    view.ResetCamera()
    Render()
    SaveScreenshot(os.path.join(images_dir, name), view, ImageResolution=SIZE)
    print("Saved", os.path.join(images_dir, name))


if images_dir:
    os.makedirs(images_dir, exist_ok=True)

    # 1. Mesh
    HideAll(view)
    mesh = Show(reader, view)
    mesh.SetRepresentationType("Surface With Edges")
    ColorBy(mesh, None)
    mesh.AmbientColor = [0.0, 0.0, 0.0]
    mesh.DiffuseColor = [0.85, 0.9, 1.0]
    mesh.EdgeColor = [0.1, 0.1, 0.1]
    save("mesh.png")

    # 2. Velocity magnitude
    HideAll(view)
    show_scalar(slice1, "U", "|U| (m/s)")
    save("velocity.png")

    # 3. Kinematic pressure
    HideAll(view)
    show_scalar(slice1, "p", "p (m^2/s^2)", component=None)
    save("pressure.png")

    # 4. Streamlines seeded across both inlets, over the |U| field
    HideAll(view)
    show_scalar(slice1, "U", "|U| (m/s)")
    seeds = [
        ([1.0, 0.5, 0.0], [1.0, 15.5, 0.0], 16),         # velocity-inlet-5
        ([52.8, -4.0, 0.0], [56.3, -4.0, 0.0], 6),       # velocity-inlet-6
    ]
    for p1, p2, res in seeds:
        tracer = StreamTracer(Input=c2p, SeedType="Line")
        tracer.Vectors = ["POINTS", "U"]
        tracer.MaximumStreamlineLength = 200.0
        tracer.SeedType.Point1 = p1
        tracer.SeedType.Point2 = p2
        tracer.SeedType.Resolution = res
        lines = Show(tracer, view)
        ColorBy(lines, None)
        lines.AmbientColor = [0.0, 0.0, 0.0]
        lines.DiffuseColor = [0.0, 0.0, 0.0]
        lines.LineWidth = 1.5
    save("streamlines.png")

elif screenshot_file:
    HideAll(view)
    show_scalar(slice1, "U", "U")
    view.ResetCamera()
    Render()
    SaveScreenshot(screenshot_file, view, ImageResolution=[1920, 1080])
    print(f"Screenshot saved to {screenshot_file}")

else:
    HideAll(view)
    show_scalar(slice1, "U", "U")
    view.ResetCamera()
    Render()
