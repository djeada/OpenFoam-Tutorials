#!/usr/bin/env python3
"""Offscreen screenshots of the fixed boat hull case (README figures).

Run from the project directory after the case has been solved:

    touch boat_hull_fixed.foam
    pvbatch paraview_screenshots.py [--time T] [--out images]

Needs constant/triSurface/boat_hull.stl (created by Allrun) and at least one
written time directory. Without --time the latest written time is used.
"""
import argparse
import contextlib
import sys
import traceback
from pathlib import Path

from paraview import simple as pv

HERE = Path(__file__).resolve().parent
FOAM = HERE / "boat_hull_fixed.foam"
STL = HERE / "constant" / "triSurface" / "boat_hull.stl"
SIZE = [1200, 675]
BACKGROUND = [0.22, 0.24, 0.35]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time", type=float, default=None,
                        help="time to render (default: latest written time)")
    parser.add_argument("--out", default=str(HERE / "images"),
                        help="output directory for the PNG files")
    return parser.parse_known_args()[0]


@contextlib.contextmanager
def step(name):
    """Run one figure; report a failure but keep going with the others."""
    try:
        yield
    except Exception:
        print("FAILED:", name)
        traceback.print_exc()


def new_view(time_value, size=None):
    view = pv.CreateView("RenderView")
    view.ViewSize = size or SIZE
    view.UseColorPaletteForBackground = 0
    view.Background = BACKGROUND
    view.OrientationAxesVisibility = 1
    view.ViewTime = time_value
    return view


def set_camera(view, position, focal, scale=None):
    view.CameraPosition = list(position)
    view.CameraFocalPoint = list(focal)
    view.CameraViewUp = [0.0, 0.0, 1.0]
    if scale is not None:
        view.CameraParallelProjection = 1
        view.CameraParallelScale = scale


def save(view, path, size=None):
    pv.Render(view)
    pv.SaveScreenshot(str(path), view, ImageResolution=size or SIZE)
    print("wrote", path)
    pv.Delete(view)


def color_by_point(source, view, array, component, lut_range, title):
    display = pv.Show(source, view)
    display.Representation = "Surface"
    pv.ColorBy(display, ("POINTS", array, component))
    lut = pv.GetColorTransferFunction(array)
    lut.ApplyPreset("Cool to Warm", True)
    if lut_range:
        lut.RescaleTransferFunction(*lut_range)
    display.SetScalarBarVisibility(view, True)
    pv.GetScalarBar(lut, view).Title = title
    pv.GetScalarBar(lut, view).ComponentTitle = ""
    return display, lut


def open_case(regions):
    reader = pv.OpenFOAMReader(FileName=str(FOAM))
    available = list(reader.GetProperty("MeshRegions").GetAvailable()) \
        if hasattr(reader.GetProperty("MeshRegions"), "GetAvailable") else []
    print("mesh regions:", available)
    reader.MeshRegions = regions
    reader.CellArrays = ["U", "p_rgh", "alpha.water"]
    return reader


def main():
    args = parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    flow = open_case(["internalMesh"])
    times = list(flow.TimestepValues) if flow.TimestepValues else [0.0]
    t = min(times, key=lambda v: abs(v - args.time)) if args.time is not None \
        else times[-1]
    print("rendering time", t)
    flow.UpdatePipeline(t)

    # Hull patch only (region name 'patch/boat' in ParaView >= 5.10)
    boat = open_case(["patch/boat"])
    boat.UpdatePipeline(t)

    # Ranges printed for the README
    print("p_rgh on hull [Pa]:", boat.CellData["p_rgh"].GetRange())
    print("|U| in domain [m/s]:", flow.CellData["U"].GetRange(-1))

    points = pv.CellDatatoPointData(Input=flow)
    points.PassCellData = 1

    # 1. Hull geometry (scaled STL) with the z = 0 waterline plane ---------
    with step('hull_stl'):
        view = new_view(t)
        stl = pv.STLReader(FileNames=[str(STL)])
        d = pv.Show(stl, view)
        d.Representation = "Surface With Edges"
        d.ColorArrayName = ["POINTS", ""]
        d.DiffuseColor = [0.85, 0.85, 0.88]
        d.EdgeColor = [0.25, 0.25, 0.3]
        plane = pv.Plane(Origin=[-0.12, -0.09, 0.0], Point1=[0.12, -0.09, 0.0],
                         Point2=[-0.12, 0.09, 0.0])
        pd = pv.Show(plane, view)
        pd.DiffuseColor = [0.2, 0.5, 0.9]
        pd.Opacity = 0.35
        set_camera(view, [0.22, -0.26, 0.16], [0.0, 0.0, 0.0])
        pv.ResetCamera(view)
        save(view, out / "hull_stl.png")

    # 2. Mesh on the centre plane: whole domain and close-up ----------------
    with step('mesh'):
        for name, size, focal, scale in (
            ("mesh_domain.png", [1400, 560], [0.25, 0.0, 0.0], 0.31),
            ("mesh_hull.png", SIZE, [0.1, 0.0, 0.0], 0.12),
        ):
            view = new_view(t, size)
            sl = pv.Slice(Input=flow)
            sl.SliceType = "Plane"
            sl.SliceType.Origin = [0.0, 0.0, 0.0]
            sl.SliceType.Normal = [0.0, 1.0, 0.0]
            d = pv.Show(sl, view)
            d.Representation = "Surface With Edges"
            d.ColorArrayName = ["POINTS", ""]
            d.DiffuseColor = [0.8, 0.85, 0.95]
            d.EdgeColor = [0.1, 0.1, 0.2]
            set_camera(view, [focal[0], -1.5, 0.0], focal, scale)
            save(view, out / name, size)

    # 3. Free surface: alpha.water = 0.5 isosurface coloured by height ------
    with step('free_surface'):
        view = new_view(t)
        iso = pv.Contour(Input=points)
        iso.ContourBy = ["POINTS", "alpha.water"]
        iso.Isosurfaces = [0.5]
        iso.ComputeNormals = 1
        height = pv.Calculator(Input=iso)
        height.Function = "coordsZ*1000"
        height.ResultArrayName = "height_mm"
        hd = pv.Show(height, view)
        pv.ColorBy(hd, ("POINTS", "height_mm"))
        lut = pv.GetColorTransferFunction("height_mm")
        lut.ApplyPreset("Cool to Warm", True)
        hd.RescaleTransferFunctionToDataRange(True, False)
        hd.SetScalarBarVisibility(view, True)
        pv.GetScalarBar(lut, view).Title = "free surface height z [mm]"
        bd = pv.Show(boat, view)
        bd.Representation = "Surface"
        bd.DiffuseColor = [0.15, 0.15, 0.18]
        set_camera(view, [0.75, -0.95, 0.5], [0.2, 0.0, 0.0])
        save(view, out / "free_surface.png")

    # 4. Velocity magnitude on the centre plane with the interface line -----
    with step('velocity_slice'):
        view = new_view(t)
        sl = pv.Slice(Input=points)
        sl.SliceType = "Plane"
        sl.SliceType.Origin = [0.0, 0.0, 0.0]
        sl.SliceType.Normal = [0.0, 1.0, 0.0]
        color_by_point(sl, view, "U", "Magnitude", (0.0, 0.5), "|U| [m/s]")
        line = pv.Contour(Input=sl)
        line.ContourBy = ["POINTS", "alpha.water"]
        line.Isosurfaces = [0.5]
        ld = pv.Show(line, view)
        ld.DiffuseColor = [0.0, 0.0, 0.0]
        ld.AmbientColor = [0.0, 0.0, 0.0]
        ld.LineWidth = 3.0
        ld.ColorArrayName = ["POINTS", ""]
        set_camera(view, [0.1, -1.5, 0.0], [0.1, 0.0, 0.0], 0.15)
        save(view, out / "velocity_slice.png")

    # 5. p_rgh on the hull surface ------------------------------------------
    with step('hull_pressure'):
        view = new_view(t)
        bd = pv.Show(boat, view)
        bd.Representation = "Surface"
        pv.ColorBy(bd, ("CELLS", "p_rgh"))
        lut = pv.GetColorTransferFunction("p_rgh")
        lut.ApplyPreset("Cool to Warm", True)
        bd.RescaleTransferFunctionToDataRange(True, False)
        bd.SetScalarBarVisibility(view, True)
        pv.GetScalarBar(lut, view).Title = "p_rgh [Pa]"
        set_camera(view, [-0.2, -0.22, 0.14], [0.0, 0.0, 0.0])
        pv.ResetCamera(view)
        save(view, out / "hull_pressure.png")

    # 6. Streamlines: reuse the interactive view builder --------------------
    with step('streamlines'):
        pv.ResetSession()
        sys.path.insert(0, str(HERE))
        sys.argv = ["paraview_streamlines.py", "--case", str(FOAM),
                    "--time", str(t), "--screenshot", str(out / "streamlines.png")]
        import paraview_streamlines
        paraview_streamlines.main()
        print("wrote", out / "streamlines.png")

main()
