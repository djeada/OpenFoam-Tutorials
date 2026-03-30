#!/usr/bin/env python3
"""
ParaView vorticity animation for the NACA 0012 transient case.

Renders the z-component of vorticity on a dark background with the airfoil
shown as a solid shape.  Zero-vorticity regions blend seamlessly into the
background via a custom colormap.  In batch mode every written timestep is
exported as a numbered PNG so that ffmpeg can compose them into video.

Interactive (GUI):  paraview --script paraview_streamlines.py
Single snapshot:    pvbatch  paraview_streamlines.py --screenshot snap.png
Full animation:     pvbatch  paraview_streamlines.py --batch --output-dir frames
"""

import argparse
from pathlib import Path

from paraview.simple import (
    AnnotateTimeFilter,
    CellDatatoPointData,
    Clip,
    ColorBy,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
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
    StreamTracer,
    UpdatePipeline,
)


# Background colour -- used by both the view and the zero-vorticity colormap.
BG = [0.05, 0.05, 0.12]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    script_file = globals().get("__file__")
    script_dir = Path(script_file).resolve().parent if script_file else Path.cwd()
    default_case = script_dir / "airfoil.foam"

    ap = argparse.ArgumentParser(
        description="NACA 0012 transient vorticity animation.",
    )
    ap.add_argument("--case", default=str(default_case),
                    help="Path to .foam file.")
    ap.add_argument("--batch", action="store_true",
                    help="Export every timestep as a numbered PNG.")
    ap.add_argument("--output-dir", default="frames",
                    help="Directory for exported frames (default: frames).")
    ap.add_argument("--screenshot", default=None,
                    help="Save a single screenshot at --time.")
    ap.add_argument("--time", type=float, default=None,
                    help="Simulation time for single-frame mode (default: latest).")
    ap.add_argument("--skip-before", type=float, default=2.0,
                    help="Skip timesteps earlier than this (default: 2.0).")
    ap.add_argument("--resolution", nargs=2, type=int, default=[1920, 1080],
                    metavar=("W", "H"),
                    help="Frame resolution (default: 1920 1080).")
    ap.add_argument("--vort-range", nargs=2, type=float, default=[-30, 30],
                    metavar=("MIN", "MAX"),
                    help="Symmetric vorticity colormap range (default: -30 30).")
    ap.add_argument("--zoom", type=float, default=1.5,
                    help="Camera parallel scale -- smaller means more zoom.")
    ap.add_argument("--camera-x", type=float, default=2.0,
                    help="Camera focal-point x (default: 2.0).")
    ap.add_argument("--camera-y", type=float, default=-0.1,
                    help="Camera focal-point y (default: -0.1).")
    ap.add_argument("--streamlines", action="store_true",
                    help="Overlay semi-transparent streamlines.")
    return ap.parse_known_args()[0]


# ---------------------------------------------------------------------------
# Pipeline helpers
# ---------------------------------------------------------------------------

def open_case(case_file, time):
    """Read the OpenFOAM case and interpolate to point data."""
    reader = OpenFOAMReader(FileName=str(case_file))
    reader.MeshRegions = ["internalMesh"]
    reader.CellArrays = ["U", "p"]
    UpdatePipeline(time=time)

    pts = CellDatatoPointData(Input=reader)
    pts.PassCellData = 1
    UpdatePipeline(time=time)
    return reader, pts


def make_z_slice(source, time, z=0.01):
    """Slice through the pseudo-2D domain at constant z."""
    sl = Slice(Input=source)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.5, 0.0, z]
    sl.SliceType.Normal = [0.0, 0.0, 1.0]
    UpdatePipeline(time=time)
    return sl


def clip_region(source, time,
                xmin=-2.0, xmax=8.0, ymin=-2.0, ymax=2.0):
    """Clip the slice to show only the near-airfoil region."""
    planes = [
        ([xmin, 0, 0], [ 1, 0, 0]),
        ([xmax, 0, 0], [-1, 0, 0]),
        ([0, ymin, 0], [0,  1, 0]),
        ([0, ymax, 0], [0, -1, 0]),
    ]
    current = source
    clips = []
    for origin, normal in planes:
        c = Clip(Input=current)
        c.ClipType = "Plane"
        c.ClipType.Origin = origin
        c.ClipType.Normal = normal
        c.Invert = 0
        UpdatePipeline(time=time)
        clips.append(c)
        current = c
    return current, clips


def compute_vorticity(source, time):
    """Return the gradient filter with vorticity enabled."""
    grad = Gradient(Input=source)
    grad.ScalarArray = ["POINTS", "U"]
    grad.ComputeVorticity = 1
    grad.VorticityArrayName = "Vorticity"
    grad.ComputeGradient = 0
    UpdatePipeline(time=time)
    return grad


def open_airfoil_patch(case_file, view):
    """Show the airfoil wall as a light solid shape on the dark background."""
    af = OpenFOAMReader(FileName=str(case_file))
    af.MeshRegions = ["patch/airfoil"]
    af.CellArrays = []

    disp = Show(af, view)
    disp.Representation = "Surface"
    disp.DiffuseColor = [0.85, 0.85, 0.88]
    disp.AmbientColor = [0.85, 0.85, 0.88]
    disp.Specular = 0.2
    return af, disp


def add_streamlines(flow_points, view):
    """Optional thin streamlines for directional context."""
    seeds = [
        ([-1.5, -0.65, 0.01], [-1.5, 0.75, 0.01], 100),
        ([-0.5, -0.35, 0.01], [-0.5, 0.40, 0.01], 60),
    ]
    tracers = []
    for pt1, pt2, res in seeds:
        tr = StreamTracer(Input=flow_points, SeedType="Line")
        tr.Vectors = ["POINTS", "U"]
        tr.IntegrationDirection = "FORWARD"
        tr.MaximumStreamlineLength = 6.0
        tr.InitialStepLength = 0.005
        tr.MaximumSteps = 5000
        tr.SeedType.Point1 = pt1
        tr.SeedType.Point2 = pt2
        tr.SeedType.Resolution = res

        d = Show(tr, view)
        d.Representation = "Surface"
        d.ColorArrayName = [None, ""]
        d.DiffuseColor = [0.55, 0.65, 0.85]
        d.AmbientColor = [0.55, 0.65, 0.85]
        d.LineWidth = 0.8
        d.Opacity = 0.30
        tracers.append((tr, d))
    return tracers


# ---------------------------------------------------------------------------
# View / coloring
# ---------------------------------------------------------------------------

def configure_view(view, resolution, cam_x=2.0, cam_y=-0.1, zoom=1.5):
    """Dark background, 2-D parallel camera."""
    view.ViewSize = resolution
    view.InteractionMode = "2D"
    view.OrientationAxesVisibility = 0
    view.Background = BG
    view.UseColorPaletteForBackground = 0
    view.CameraPosition = [cam_x, cam_y, 6.0]
    view.CameraFocalPoint = [cam_x, cam_y, 0.01]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.CameraParallelProjection = 1
    view.CameraParallelScale = zoom


def apply_vorticity_coloring(view, display, vort_range):
    """Custom colormap: zero vorticity blends into the dark background."""
    ColorBy(display, ("POINTS", "Vorticity", "Z"))

    vmin, vmax = float(vort_range[0]), float(vort_range[1])

    lut = GetColorTransferFunction("Vorticity")
    lut.VectorMode = "Component"
    lut.VectorComponent = 2
    # Custom diverging colormap -- zero maps to the background colour so
    # the freestream disappears and only vortical structures are visible.
    lut.RGBPoints = [
        vmin,          0.15, 0.30, 0.90,
        vmin * 0.25,   0.10, 0.15, 0.45,
        vmin * 0.04,   BG[0], BG[1], BG[2],
        0.0,           BG[0], BG[1], BG[2],
        vmax * 0.04,   BG[0], BG[1], BG[2],
        vmax * 0.25,   0.55, 0.18, 0.10,
        vmax,          0.95, 0.30, 0.10,
    ]

    display.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = "\u03c9z  [1/s]"
    bar.ComponentTitle = ""
    bar.TitleColor = [0.90, 0.90, 0.90]
    bar.LabelColor = [0.90, 0.90, 0.90]
    bar.TitleFontSize = 18
    bar.LabelFontSize = 14


def add_time_annotation(source, view):
    """Overlay current simulation time in the upper-left corner."""
    ann = AnnotateTimeFilter(Input=source)
    ann.Format = "t = %.3f s"
    ann.Scale = 1.0
    ann.Shift = 0.0

    disp = Show(ann, view)
    disp.FontFamily = "Courier"
    disp.FontSize = 20
    disp.Color = [0.92, 0.92, 0.92]
    disp.WindowLocation = "Upper Left Corner"
    return ann, disp


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def pick_time(requested, timesteps):
    """Return the closest available timestep to *requested*, or the last one."""
    if not timesteps:
        return 0.0
    if requested is None:
        return float(timesteps[-1])
    return float(min(timesteps, key=lambda v: abs(v - requested)))


def export_frames(view, scene, timesteps, output_dir, resolution, skip_before):
    """Render each qualifying timestep to a numbered PNG."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    frames = [t for t in timesteps if t >= skip_before]
    n = len(frames)
    if n == 0:
        print("WARNING: no timesteps after skip-before threshold; nothing to export.")
        return
    print("Exporting {} frames to {}/".format(n, out))

    for i, t in enumerate(frames):
        scene.AnimationTime = t
        view.ViewTime = t
        Render()
        path = str(out / "frame_{:05d}.png".format(i))
        SaveScreenshot(path, view, ImageResolution=resolution)
        if i == 0 or (i + 1) % 20 == 0 or i == n - 1:
            pct = 100.0 * (i + 1) / n
            print("  [{:5.1f}%]  frame {:>5d}/{:<5d}  t={:.4f}".format(pct, i + 1, n, t))

    print("Done -- {} frames written.".format(n))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    case_file = Path(args.case).resolve()
    if not case_file.exists():
        raise FileNotFoundError("Case file not found: {}".format(case_file))

    view = GetActiveViewOrCreate("RenderView")
    scene = GetAnimationScene()
    tk = GetTimeKeeper()

    # --- build pipeline (UpdatePipeline at each stage for PV 5.11) ---
    reader, pts = open_case(case_file, time=0)
    scene.UpdateAnimationUsingDataTimeSteps()
    timesteps = sorted(getattr(reader, "TimestepValues", []))
    t0 = pick_time(args.time, timesteps)

    sl = make_z_slice(pts, time=t0)
    clipped, clip_list = clip_region(sl, time=t0)
    grad = compute_vorticity(clipped, time=t0)

    configure_view(view, args.resolution,
                   cam_x=args.camera_x, cam_y=args.camera_y, zoom=args.zoom)

    # Vorticity field
    vort_disp = Show(grad, view)
    vort_disp.Representation = "Surface"
    apply_vorticity_coloring(view, vort_disp, args.vort_range)

    # Airfoil solid surface
    af, _af_disp = open_airfoil_patch(case_file, view)

    # Optional streamlines
    if args.streamlines:
        add_streamlines(pts, view)

    # Hide intermediate pipeline objects
    for src in [reader, pts, sl] + clip_list:
        Hide(src, view)

    # Time annotation
    SetActiveSource(grad)
    add_time_annotation(grad, view)

    # Set display time
    scene.AnimationTime = t0
    view.ViewTime = t0
    tk.Time = t0
    Render()

    # --- output mode ---
    if args.batch:
        export_frames(view, scene, timesteps, args.output_dir,
                      args.resolution, args.skip_before)
    elif args.screenshot:
        SaveScreenshot(args.screenshot, view, ImageResolution=args.resolution)
        print("Screenshot saved: {}".format(args.screenshot))
    # else: interactive -- the ParaView GUI stays open


if __name__ == "__main__":
    main()
