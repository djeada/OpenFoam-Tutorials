#!/usr/bin/env python3
import argparse
from pathlib import Path

from paraview.simple import (
    CellDatatoPointData,
    ColorBy,
    GetAnimationScene,
    GetColorTransferFunction,
    GetOpacityTransferFunction,
    GetScalarBar,
    GetTimeKeeper,
    GetTransferFunction2D,
    GetActiveViewOrCreate,
    Hide,
    OpenFOAMReader,
    Render,
    ResetCamera,
    SaveScreenshot,
    SetActiveSource,
    Show,
    Slice,
    StreamTracer,
    Tube,
)


def parse_args():
    script_file = globals().get("__file__")
    if script_file:
        script_dir = Path(script_file).resolve().parent
    else:
        script_dir = Path.cwd()
    default_case = script_dir / "boat_hull_fixed.foam"

    parser = argparse.ArgumentParser(
        description="Build a ParaView view for project 10 with a center slice, hull, and streamlines."
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
    view.Background = [0.22, 0.24, 0.35]
    view.UseColorPaletteForBackground = 0
    view.CameraPosition = [0.25, -1.5, 0.0]
    view.CameraFocalPoint = [0.25, 0.0, 0.0]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = 0.22


def show_slice(flow_points, view):
    slice_filter = Slice(Input=flow_points)
    slice_filter.SliceType = "Plane"
    slice_filter.SliceType.Origin = [0.0, 0.0, 0.0]
    slice_filter.SliceType.Normal = [0.0, 1.0, 0.0]

    display = Show(slice_filter, view)
    display.Representation = "Surface"
    display.Opacity = 0.35
    ColorBy(display, ("POINTS", "U", "Magnitude"))

    return slice_filter, display


def show_boat(case_file, view):
    boat = OpenFOAMReader(FileName=str(case_file))
    boat.MeshRegions = ["patch/boat"]
    boat.CellArrays = ["U", "alpha.water", "p_rgh"]

    display = Show(boat, view)
    display.Representation = "Surface"
    display.DiffuseColor = [0.18, 0.18, 0.22]
    display.AmbientColor = [0.18, 0.18, 0.22]
    display.Specular = 0.2

    return boat, display


def show_streamlines(flow_points, view):
    tracers = []
    seeds = [
        ([-0.22, 0.0, -0.05], [-0.22, 0.0, 0.01], 18),
        ([-0.22, 0.0, -0.08], [-0.22, 0.0, -0.03], 14),
    ]

    for point1, point2, resolution in seeds:
        tracer = StreamTracer(Input=flow_points, SeedType="Line")
        tracer.Vectors = ["POINTS", "U"]
        tracer.IntegrationDirection = "FORWARD"
        tracer.MaximumStreamlineLength = 2.0
        tracer.SeedType.Point1 = point1
        tracer.SeedType.Point2 = point2
        tracer.SeedType.Resolution = resolution

        tube = Tube(Input=tracer)
        tube.Radius = 0.0015
        tube.NumberofSides = 16

        display = Show(tube, view)
        display.Representation = "Surface"
        ColorBy(display, ("POINTS", "U", "Magnitude"))

        tracers.append((tracer, tube, display))

    return tracers


def configure_coloring(view, displays):
    u_lut = GetColorTransferFunction("U")
    u_pwf = GetOpacityTransferFunction("U")
    u_tf2d = GetTransferFunction2D("U")

    for display in displays:
        display.RescaleTransferFunctionToDataRange(True, False)
        display.SetScalarBarVisibility(view, True)

    u_lut.ApplyPreset("Cool to Warm", True)
    u_lut.RescaleTransferFunction(0.0, 2.2)
    u_pwf.RescaleTransferFunction(0.0, 2.2)
    u_tf2d.RescaleTransferFunction(0.0, 2.2, 0.0, 1.0)

    scalar_bar = GetScalarBar(u_lut, view)
    scalar_bar.Title = "U Magnitude"
    scalar_bar.ComponentTitle = ""


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
    flow.CellArrays = ["U", "alpha.water", "p_rgh"]

    flow_points = CellDatatoPointData(Input=flow)
    flow_points.PassCellData = 1

    animation_scene.UpdateAnimationUsingDataTimeSteps()
    time_value = choose_time(args.time, flow)
    animation_scene.AnimationTime = time_value
    time_keeper.Time = time_value
    configure_view(view, time_value)

    slice_filter, slice_display = show_slice(flow_points, view)
    boat, _boat_display = show_boat(case_file, view)
    tracers = show_streamlines(flow_points, view)

    Hide(flow, view)
    Hide(flow_points, view)

    displays = [slice_display] + [display for _tracer, _tube, display in tracers]
    configure_coloring(view, displays)

    SetActiveSource(slice_filter)
    Render()
    ResetCamera(view)
    configure_view(view, time_value)
    Render()

    if args.screenshot:
        SaveScreenshot(args.screenshot, view, ImageResolution=[1600, 900])


if __name__ == "__main__":
    main()
