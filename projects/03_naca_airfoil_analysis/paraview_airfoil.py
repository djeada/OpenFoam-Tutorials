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
    GetLayout,
    GetScalarBar,
    Hide,
    OpenFOAMReader,
    RenameSource,
    Render,
    ResetCamera,
    SetActiveSource,
    Show,
    Slice,
    StreamTracer,
)


SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
CASE_FILE = SCRIPT_DIR / "airfoil.foam"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a ParaView view for the main NACA 0012 airfoil case."
    )
    parser.add_argument(
        "--case",
        default=str(CASE_FILE),
        help="Path to the .foam file.",
    )
    parser.add_argument(
        "--screenshot",
        default=None,
        help="Optional screenshot output path.",
    )
    parser.add_argument(
        "--pressure-limit",
        type=float,
        default=300.0,
        help="Symmetric pressure color limit; uses [-limit, +limit].",
    )
    parser.add_argument(
        "--scalar",
        choices=("U", "p"),
        default="U",
        help="Primary scalar shown on the slice.",
    )
    parser.add_argument(
        "--u-min",
        type=float,
        default=20.0,
        help="Lower limit for U magnitude coloring.",
    )
    parser.add_argument(
        "--u-max",
        type=float,
        default=70.0,
        help="Upper limit for U magnitude coloring.",
    )
    return parser.parse_args()


def latest_time_value(case_dir: Path) -> float:
    values = []
    for child in case_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            values.append(float(child.name))
        except ValueError:
            continue
    if not values:
        return 0.0
    return max(values)


def configure_view(view):
    view.InteractionMode = "2D"
    view.CameraParallelProjection = 1
    view.UseColorPaletteForBackground = 0
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 0
    view.CameraPosition = [0.52, 0.01, 6.0]
    view.CameraFocalPoint = [0.52, 0.01, 0.01]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelScale = 0.17


def configure_pressure_coloring(view, slice_display, pressure_limit):
    p_lut = GetColorTransferFunction("p")
    p_pwf = GetOpacityTransferFunction("p")

    slice_display.SetScalarBarVisibility(view, True)
    p_lut.ApplyPreset("Cool to Warm", True)
    p_lut.RescaleTransferFunction(-pressure_limit, pressure_limit)
    p_pwf.RescaleTransferFunction(-pressure_limit, pressure_limit)

    scalar_bar = GetScalarBar(p_lut, view)
    scalar_bar.Title = "p"
    scalar_bar.ComponentTitle = ""
    scalar_bar.WindowLocation = "Lower Center"
    scalar_bar.HorizontalTitle = 1


def configure_velocity_coloring(view, slice_display, u_min, u_max):
    u_lut = GetColorTransferFunction("U")
    u_pwf = GetOpacityTransferFunction("U")

    slice_display.SetScalarBarVisibility(view, True)
    u_lut.ApplyPreset("Cool to Warm", True)
    u_lut.RescaleTransferFunction(u_min, u_max)
    u_pwf.RescaleTransferFunction(u_min, u_max)

    scalar_bar = GetScalarBar(u_lut, view)
    scalar_bar.Title = "U Magnitude"
    scalar_bar.ComponentTitle = ""
    scalar_bar.WindowLocation = "Lower Center"
    scalar_bar.HorizontalTitle = 1


def build_scene(case_file, screenshot=None, pressure_limit=300.0, scalar="U", u_min=20.0, u_max=70.0):
    view = GetActiveViewOrCreate("RenderView")
    layout = GetLayout()
    if layout is not None:
        layout.SetSize(1600, 900)

    animation_scene = GetAnimationScene()
    time_value = latest_time_value(case_file.parent)

    flow = OpenFOAMReader(FileName=str(case_file))
    RenameSource("flow", flow)
    flow.CaseType = "Reconstructed Case"
    flow.MeshRegions = ["internalMesh"]
    flow.CellArrays = ["U", "p", "nuTilda"]

    airfoil = OpenFOAMReader(FileName=str(case_file))
    RenameSource("airfoil_patch", airfoil)
    airfoil.CaseType = "Reconstructed Case"
    airfoil.MeshRegions = ["patch/airfoil"]
    airfoil.CellArrays = ["U", "p"]

    animation_scene.UpdateAnimationUsingDataTimeSteps()
    animation_scene.AnimationTime = time_value
    view.ViewTime = time_value

    configure_view(view)

    flow_points = CellDatatoPointData(Input=flow)
    RenameSource("flow_points", flow_points)
    flow_points.PassCellData = 1

    slice_filter = Slice(Input=flow_points)
    RenameSource("pressure_slice", slice_filter)
    slice_filter.SliceType = "Plane"
    slice_filter.SliceType.Origin = [0.5, 0.0, 0.01]
    slice_filter.SliceType.Normal = [0.0, 0.0, 1.0]

    slice_display = Show(slice_filter, view)
    slice_display.Representation = "Surface LIC"
    slice_display.Opacity = 0.92
    if scalar == "p":
        ColorBy(slice_display, ("POINTS", "p"))
    else:
        ColorBy(slice_display, ("POINTS", "U", "Magnitude"))
    slice_display.SelectInputVectors = ["POINTS", "U"]
    slice_display.InterpolateScalarsBeforeMapping = 0
    slice_display.NumberOfSteps = 32
    slice_display.StepSize = 0.2
    slice_display.LICIntensity = 0.9
    slice_display.EnhanceContrast = "LIC and Color"
    slice_display.LowLICContrastEnhancementFactor = 0.2
    slice_display.HighLICContrastEnhancementFactor = 0.1
    slice_display.NoiseTextureSize = 64
    slice_display.AntiAlias = 1

    airfoil_display = Show(airfoil, view)
    airfoil_display.Representation = "Surface"
    airfoil_display.DiffuseColor = [1.0, 1.0, 1.0]
    airfoil_display.AmbientColor = [1.0, 1.0, 1.0]
    airfoil_display.Ambient = 0.05
    airfoil_display.Diffuse = 0.95
    airfoil_display.Specular = 0.15

    Hide(flow, view)
    Hide(flow_points, view)

    if scalar == "p":
        configure_pressure_coloring(view, slice_display, pressure_limit)
    else:
        configure_velocity_coloring(view, slice_display, u_min, u_max)

    SetActiveSource(slice_filter)
    Render()
    ResetCamera(view)
    configure_view(view)
    Render()

    if screenshot:
        from paraview.simple import SaveScreenshot

        SaveScreenshot(screenshot, view, ImageResolution=[1600, 900])


args = parse_args()
case_file = Path(args.case).resolve()
if not case_file.exists():
    raise FileNotFoundError(f"Could not find case file: {case_file}")

build_scene(
    case_file,
    args.screenshot,
    args.pressure_limit,
    args.scalar,
    args.u_min,
    args.u_max,
)
