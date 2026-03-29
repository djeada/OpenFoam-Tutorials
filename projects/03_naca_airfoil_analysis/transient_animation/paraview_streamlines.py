#!/usr/bin/env python3
import argparse
from pathlib import Path

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
    ResetCamera,
    SaveScreenshot,
    SetActiveSource,
    Show,
    Slice,
    StreamTracer,
)


def parse_args():
    script_file = globals().get("__file__")
    script_dir = Path(script_file).resolve().parent if script_file else Path.cwd()
    default_case = script_dir / "airfoil.foam"

    parser = argparse.ArgumentParser(
        description="Build a ParaView view for the transient NACA 0012 case."
    )
    parser.add_argument(
        "--case",
        default=str(default_case),
        help="Path to the .foam file.",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=None,
        help="Simulation time to display. Defaults to the latest available time.",
    )
    parser.add_argument(
        "--screenshot",
        default=None,
        help="Optional output image path for batch rendering.",
    )
    return parser.parse_args()


def choose_time(requested_time, reader):
    timesteps = list(getattr(reader, "TimestepValues", []))
    if not timesteps:
        return 0.0
    if requested_time is None:
        return float(timesteps[-1])
    return float(min(timesteps, key=lambda value: abs(value - requested_time)))


def configure_view(view, time_value):
    view.ViewTime = time_value
    view.InteractionMode = "2D"
    view.OrientationAxesVisibility = 1
    view.Background = [1.0, 1.0, 1.0]
    view.UseColorPaletteForBackground = 0
    view.CameraPosition = [0.55, 0.0, 6.0]
    view.CameraFocalPoint = [0.55, 0.0, 0.01]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 0.52


def show_slice(flow_points, view):
    slice_filter = Slice(Input=flow_points)
    slice_filter.SliceType = "Plane"
    slice_filter.SliceType.Origin = [0.5, 0.0, 0.01]
    slice_filter.SliceType.Normal = [0.0, 0.0, 1.0]

    display = Show(slice_filter, view)
    display.Representation = "Surface"
    ColorBy(display, ("POINTS", "p"))

    return slice_filter, display


def show_airfoil(case_file, view):
    airfoil = OpenFOAMReader(FileName=str(case_file))
    airfoil.MeshRegions = ["patch/airfoil"]
    airfoil.CellArrays = ["U", "p"]

    display = Show(airfoil, view)
    display.Representation = "Surface"
    display.DiffuseColor = [0.35, 0.36, 0.40]
    display.AmbientColor = [0.35, 0.36, 0.40]
    display.Specular = 0.15

    return airfoil, display


def show_streamlines(flow_points, view):
    tracers = []
    seeds = [
        ([-1.30, -0.55, 0.01], [-1.30, 0.65, 0.01], 120),
        ([-0.45, -0.30, 0.01], [-0.45, 0.35, 0.01], 80),
    ]

    for index, (point1, point2, resolution) in enumerate(seeds, start=1):
        tracer = StreamTracer(Input=flow_points, SeedType="Line")
        tracer.Vectors = ["POINTS", "U"]
        tracer.IntegrationDirection = "FORWARD"
        tracer.MaximumStreamlineLength = 5.0
        tracer.InitialStepLength = 0.005
        tracer.MaximumSteps = 4000
        tracer.SeedType.Point1 = point1
        tracer.SeedType.Point2 = point2
        tracer.SeedType.Resolution = resolution

        display = Show(tracer, view)
        display.Representation = "Surface"
        display.ColorArrayName = [None, ""]
        display.DiffuseColor = [1.0, 1.0, 1.0]
        display.AmbientColor = [1.0, 1.0, 1.0]
        display.LineWidth = 1.2 if index == 1 else 1.0
        display.Opacity = 0.78

        tracers.append((tracer, display))

    return tracers


def configure_coloring(view, slice_display):
    p_lut = GetColorTransferFunction("p")
    p_pwf = GetOpacityTransferFunction("p")

    slice_display.RescaleTransferFunctionToDataRange(True, False)
    slice_display.SetScalarBarVisibility(view, True)

    p_lut.ApplyPreset("Cool to Warm", True)

    scalar_bar = GetScalarBar(p_lut, view)
    scalar_bar.Title = "p"
    scalar_bar.ComponentTitle = ""

    # Match the opacity function range to the visible data range.
    data_range = slice_display.LookupTable.RGBPoints
    if len(data_range) >= 8:
        p_min = data_range[0]
        p_max = data_range[-4]
        p_pwf.RescaleTransferFunction(p_min, p_max)


def main():
    args = parse_args()
    case_file = Path(args.case).resolve()
    if not case_file.exists():
        raise FileNotFoundError(f"Could not find OpenFOAM case file: {case_file}")

    view = GetActiveViewOrCreate("RenderView")
    animation_scene = GetAnimationScene()
    time_keeper = GetTimeKeeper()

    flow = OpenFOAMReader(FileName=str(case_file))
    flow.MeshRegions = ["internalMesh"]
    flow.CellArrays = ["U", "p"]

    flow_points = CellDatatoPointData(Input=flow)
    flow_points.PassCellData = 1

    animation_scene.UpdateAnimationUsingDataTimeSteps()
    time_value = choose_time(args.time, flow)
    animation_scene.AnimationTime = time_value
    time_keeper.Time = time_value
    configure_view(view, time_value)

    slice_filter, slice_display = show_slice(flow_points, view)
    _airfoil, _airfoil_display = show_airfoil(case_file, view)
    _tracers = show_streamlines(flow_points, view)

    Hide(flow, view)
    Hide(flow_points, view)

    configure_coloring(view, slice_display)

    SetActiveSource(slice_filter)
    Render()
    ResetCamera(view)
    configure_view(view, time_value)
    Render()

    if args.screenshot:
        SaveScreenshot(args.screenshot, view, ImageResolution=[1600, 900])


if __name__ == "__main__":
    main()
