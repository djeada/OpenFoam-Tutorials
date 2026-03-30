"""ParaView Python script — 2D cylinder vortex shedding visualisation.

Displays the Z-component of vorticity with a diverging colourmap to
highlight the von Kármán vortex street.  The cylinder wall is shown as
a dark solid surface.

Usage (interactive):
    paraview --script paraview_cylinder.py

Usage (batch screenshot):
    pvbatch paraview_cylinder.py --screenshot vorticity.png
"""

import sys
import os

from paraview.simple import (
    AnnotateTimeFilter,
    CellDatatoPointData,
    ColorBy,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    GetScalarBar,
    GetTimeKeeper,
    Gradient,
    Hide,
    OpenFOAMReader,
    Render,
    SaveScreenshot,
    SetActiveSource,
    Show,
    Slice,
)

# ── Configuration ─────────────────────────────────────────────────────────
FOAM_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "cylinder.foam")
VORTICITY_RANGE = [-5.0, 5.0]
BACKGROUND_COLOR = [0.02, 0.02, 0.07]
SCREENSHOT_SIZE = [1920, 1080]

# ── Parse --screenshot flag ──────────────────────────────────────────────
screenshot_path = None
if "--screenshot" in sys.argv:
    idx = sys.argv.index("--screenshot")
    if idx + 1 < len(sys.argv):
        screenshot_path = sys.argv[idx + 1]
    else:
        screenshot_path = "vorticity.png"

# ── Open case ─────────────────────────────────────────────────────────────
reader = OpenFOAMReader(FileName=FOAM_FILE)
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["U", "p"]

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()

# Jump to last time step
tk = GetTimeKeeper()
times = reader.TimestepValues
if times:
    scene.AnimationTime = times[-1]

# ── Cell-to-Point interpolation ──────────────────────────────────────────
c2p = CellDatatoPointData(Input=reader)
c2p.CellDataArraytoprocess = ["U", "p"]

# ── Z-slice at midplane (z = 0.05) ──────────────────────────────────────
zslice = Slice(Input=c2p)
zslice.SliceType = "Plane"
zslice.SliceType.Origin = [0.0, 0.0, 0.05]
zslice.SliceType.Normal = [0.0, 0.0, 1.0]
zslice.Triangulatetheslice = False

# ── Compute vorticity via Gradient filter ────────────────────────────────
grad = Gradient(Input=zslice)
grad.ScalarArray = ["POINTS", "U"]
grad.ComputeVorticity = 1
grad.VorticityArrayName = "Vorticity"
grad.ComputeGradient = 0

# ── View setup ────────────────────────────────────────────────────────────
view = GetActiveViewOrCreate("RenderView")
view.ViewSize = SCREENSHOT_SIZE
view.Background = BACKGROUND_COLOR
view.InteractionMode = "2D"

# Camera: parallel projection centred on the wake
view.CameraParallelProjection = 1
view.CameraPosition = [5.0, 0.0, 20.0]
view.CameraFocalPoint = [5.0, 0.0, 0.05]
view.CameraViewUp = [0.0, 1.0, 0.0]
view.CameraParallelScale = 8.0

# ── Display vorticity field ──────────────────────────────────────────────
SetActiveSource(grad)
display = Show(grad, view)
ColorBy(display, ("POINTS", "Vorticity", "Z"))

vorticityLUT = GetColorTransferFunction("Vorticity")
vorticityLUT.ApplyPreset("Cool to Warm", True)
vorticityLUT.RescaleTransferFunction(VORTICITY_RANGE[0], VORTICITY_RANGE[1])

vorticityPWF = GetOpacityTransferFunction("Vorticity")
vorticityPWF.RescaleTransferFunction(VORTICITY_RANGE[0], VORTICITY_RANGE[1])

display.SetScalarBarVisibility(view, True)

# ── Scalar bar styling ───────────────────────────────────────────────────
scalarBar = GetScalarBar(vorticityLUT, view)
scalarBar.Title = "Vorticity Z"
scalarBar.ComponentTitle = ""
scalarBar.TitleColor = [0.9, 0.9, 0.9]
scalarBar.LabelColor = [0.9, 0.9, 0.9]
scalarBar.TitleFontSize = 16
scalarBar.LabelFontSize = 14

# ── Cylinder wall ────────────────────────────────────────────────────────
Hide(grad, view)

SetActiveSource(reader)
wallDisplay = Show(reader, view)
wallDisplay.Representation = "Surface"
wallDisplay.AmbientColor = [0.15, 0.15, 0.15]
wallDisplay.DiffuseColor = [0.15, 0.15, 0.15]

# Re-show the vorticity slice on top
SetActiveSource(grad)
Show(grad, view)

# ── Time annotation ──────────────────────────────────────────────────────
SetActiveSource(grad)
timeAnnotation = AnnotateTimeFilter(Input=grad)
timeAnnotation.Format = "Time: {time:.2f} s"

timeDisplay = Show(timeAnnotation, view)
timeDisplay.FontSize = 14
timeDisplay.Color = [0.9, 0.9, 0.9]
timeDisplay.WindowLocation = "Upper Right Corner"

# ── Render / Screenshot ──────────────────────────────────────────────────
Render()

if screenshot_path:
    SaveScreenshot(screenshot_path, view,
                   ImageResolution=SCREENSHOT_SIZE,
                   TransparentBackground=0)
    print(f"Saved screenshot to {screenshot_path}")
