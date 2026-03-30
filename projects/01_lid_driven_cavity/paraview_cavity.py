"""ParaView Python script for the lid-driven cavity case.

Visualizes velocity field with streamlines on the XY midplane.
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
    GetTimeKeeper,
    Hide,
    OpenFOAMReader,
    Render,
    SaveScreenshot,
    Show,
    Slice,
    StreamTracer,
)

# --- Read the OpenFOAM case ---
cavity = OpenFOAMReader(registrationName="cavity.foam", FileName="cavity.foam")
cavity.MeshRegions = ["internalMesh"]
cavity.CellArrays = ["U", "p"]

# Advance to the last timestep
animationScene = GetAnimationScene()
animationScene.UpdateAnimationUsingDataTimeSteps()
timeKeeper = GetTimeKeeper()
timesteps = timeKeeper.TimestepValues
if timesteps:
    animationScene.AnimationTime = timesteps[-1]

# --- Set up the render view ---
renderView = GetActiveViewOrCreate("RenderView")
renderView.ViewSize = [1280, 960]
renderView.Background = [1.0, 1.0, 1.0]

# 2D parallel camera looking down the Z axis at the XY midplane
renderView.InteractionMode = "2D"
renderView.CameraPosition = [0.05, 0.05, 0.5]
renderView.CameraFocalPoint = [0.05, 0.05, 0.005]
renderView.CameraViewUp = [0.0, 1.0, 0.0]
renderView.CameraParallelScale = 0.055
renderView.CameraParallelProjection = 1

# --- Convert cell data to point data ---
cellToPoint = CellDatatoPointData(registrationName="CellDatatoPointData", Input=cavity)
cellToPoint.CellDataArraytoprocess = ["U", "p"]

# --- Create a z-slice at the midplane (z = 0.005) ---
zSlice = Slice(registrationName="ZSlice", Input=cellToPoint)
zSlice.SliceType = "Plane"
zSlice.SliceType.Origin = [0.05, 0.05, 0.005]
zSlice.SliceType.Normal = [0.0, 0.0, 1.0]
zSlice.SliceOffsetValues = [0.0]

# Show slice colored by velocity magnitude
zSliceDisplay = Show(zSlice, renderView)
ColorBy(zSliceDisplay, ("POINTS", "U", "Magnitude"))
zSliceDisplay.SetScalarBarVisibility(renderView, True)

# Apply "Cool to Warm" color preset
uLUT = GetColorTransferFunction("U")
uLUT.ApplyPreset("Cool to Warm", True)
uLUT.RescaleTransferFunction(0.0, 1.0)

uPWF = GetOpacityTransferFunction("U")
uPWF.RescaleTransferFunction(0.0, 1.0)

# Configure the scalar bar
scalarBar = GetScalarBar(uLUT, renderView)
scalarBar.Title = "U"
scalarBar.ComponentTitle = "Magnitude"
scalarBar.TitleColor = [0.0, 0.0, 0.0]
scalarBar.LabelColor = [0.0, 0.0, 0.0]
scalarBar.Visibility = 1

# --- Create streamlines for flow visualization ---
streamTracer = StreamTracer(registrationName="StreamTracer", Input=zSlice, SeedType="Line")
streamTracer.Vectors = ["POINTS", "U"]
streamTracer.MaximumStreamlineLength = 0.2

# Seed line across the cavity
streamTracer.SeedType.Point1 = [0.005, 0.005, 0.005]
streamTracer.SeedType.Point2 = [0.095, 0.095, 0.005]
streamTracer.SeedType.Resolution = 20

streamTracerDisplay = Show(streamTracer, renderView)
ColorBy(streamTracerDisplay, ("POINTS", "U", "Magnitude"))
streamTracerDisplay.LineWidth = 2.0

# Hide the raw reader
Hide(cavity, renderView)

Render()

# --- Save screenshot if requested ---
if "--screenshot" in sys.argv:
    idx = sys.argv.index("--screenshot")
    if idx + 1 < len(sys.argv):
        output_file = sys.argv[idx + 1]
    else:
        output_file = "cavity.png"
    SaveScreenshot(output_file, renderView, ImageResolution=[1280, 960])
    print(f"Screenshot saved to {output_file}")
