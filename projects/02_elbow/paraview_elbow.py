"""ParaView Python script for the OpenFOAM elbow tutorial.

Visualises velocity magnitude on a z-midplane slice of the L-shaped
elbow pipe (two inlets merging, one outlet).

Usage:
    paraview --script paraview_elbow.py
    pvbatch paraview_elbow.py [--screenshot output.png]
"""

import sys
from paraview.simple import (
    CellDatatoPointData,
    ColorBy,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    GetScalarBar,
    Hide,
    OpenFOAMReader,
    Render,
    SaveScreenshot,
    Show,
    Slice,
)

# ── Parse optional --screenshot argument ──────────────────────────────
screenshot_file = None
if "--screenshot" in sys.argv:
    idx = sys.argv.index("--screenshot")
    if idx + 1 < len(sys.argv):
        screenshot_file = sys.argv[idx + 1]

# ── Open the case ─────────────────────────────────────────────────────
reader = OpenFOAMReader(FileName="elbow.foam")
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["U", "p"]

# Advance to the last time-step
animation = GetAnimationScene()
animation.UpdateAnimationUsingDataTimeSteps()
animation.GoToLast()

# ── Convert cell data to point data for smooth rendering ──────────────
c2p = CellDatatoPointData(Input=reader)
c2p.CellDataArraytoprocess = ["U", "p"]

# ── Z-normal slice through the midplane (2-D case) ───────────────────
slice1 = Slice(Input=c2p)
slice1.SliceType = "Plane"
slice1.SliceType.Origin = [0.0, 0.0, 0.0]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
slice1.SliceOffsetValues = [0.0]

# ── Set up the render view ───────────────────────────────────────────
view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1920, 1080]
view.Background = [1.0, 1.0, 1.0]

# Hide the raw reader; show only the slice
Hide(reader, view)
sliceDisplay = Show(slice1, view)

# ── Colour by velocity magnitude ─────────────────────────────────────
ColorBy(sliceDisplay, ("POINTS", "U", "Magnitude"))
sliceDisplay.RescaleTransferFunctionToDataRange(True, False)
sliceDisplay.SetScalarBarVisibility(view, True)

uLUT = GetColorTransferFunction("U")
uLUT.ApplyPreset("Cool to Warm", True)

uPWF = GetOpacityTransferFunction("U")

# ── Scalar bar ────────────────────────────────────────────────────────
scalarBar = GetScalarBar(uLUT, view)
scalarBar.Title = "U"
scalarBar.ComponentTitle = "Magnitude"
scalarBar.TitleColor = [0.0, 0.0, 0.0]
scalarBar.LabelColor = [0.0, 0.0, 0.0]

# ── Camera – 2-D parallel projection showing the full elbow ──────────
view.InteractionMode = "2D"
view.CameraPosition = [0.0, 0.0, 1.0]
view.CameraFocalPoint = [0.0, 0.0, 0.0]
view.CameraViewUp = [0.0, 1.0, 0.0]
view.CameraParallelScale = 0.12
view.CameraParallelProjection = 1
view.ResetCamera()

Render()

# ── Optional screenshot ──────────────────────────────────────────────
if screenshot_file:
    SaveScreenshot(screenshot_file, view, ImageResolution=[1920, 1080])
    print(f"Screenshot saved to {screenshot_file}")
