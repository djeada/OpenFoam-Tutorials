# Render the result images used in README.md (offscreen).
# Run after a finished case:  pvpython make_images.py     (or pvbatch)
# Output: images/alpha_snapshots.png, images/velocity_pressure.png
# Needs ParaView (pvpython/pvbatch) and, for tiling, ImageMagick `montage`.
import os
import shutil
import subprocess
import tempfile

from paraview.simple import *

case = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(case, "images")
os.makedirs(out, exist_ok=True)
tmp = tempfile.mkdtemp()

open(os.path.join(case, "dam_break.foam"), "a").close()
reader = OpenFOAMReader(FileName=os.path.join(case, "dam_break.foam"))
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["alpha.water", "U", "p_rgh"]
reader.UpdatePipelineInformation()

view = CreateRenderView()
view.ViewSize = [400, 420]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.CameraPosition = [0.292, 0.292, 2.0]
view.CameraFocalPoint = [0.292, 0.292, 0.0]
view.CameraViewUp = [0, 1, 0]
view.CameraParallelScale = 0.31
view.Background = [1, 1, 1]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass

display = Show(reader, view)
display.Ambient = 1.0   # flat colours, no shading
display.Diffuse = 0.0
label = Text(Text="")
label_display = Show(label, view)
label_display.FontSize = 20
label_display.Color = [0, 0, 0]
label_display.WindowLocation = "Upper Center"


def frame(time, array, comp, rng, preset, name, title, bar=False):
    view.ViewTime = time
    ColorBy(display, ("CELLS", array, comp) if comp else ("CELLS", array))
    lut = GetColorTransferFunction(array)
    if preset is None:  # air = light grey, water = blue
        lut.RGBPoints = [0.0, 0.93, 0.93, 0.93, 1.0, 0.10, 0.30, 0.75]
        lut.ColorSpace = "RGB"
    else:
        lut.ApplyPreset(preset, True)
    lut.RescaleTransferFunction(*rng)
    display.SetScalarBarVisibility(view, bar)
    if bar:
        sb = GetScalarBar(lut, view)
        sb.Title = title
        sb.ComponentTitle = ""
        sb.TitleColor = [0, 0, 0]
        sb.LabelColor = [0, 0, 0]
        sb.WindowLocation = "Any Location"
        sb.Orientation = "Horizontal"
        sb.Position = [0.2, 0.02]
        sb.ScalarBarLength = 0.6
        sb.AddRangeLabels = 0
    label.Text = title if not bar else "t = %.2f s" % time
    Render(view)
    path = os.path.join(tmp, name + ".png")
    SaveScreenshot(path, view, ImageResolution=view.ViewSize)
    if bar:
        GetScalarBar(lut, view).Visibility = 0
    return path


def tile(paths, cols, target):
    if shutil.which("montage") is None:
        raise SystemExit("ImageMagick 'montage' not found; frames are in " + tmp)
    subprocess.check_call(
        ["montage"] + paths + ["-tile", "%dx" % cols, "-geometry", "+2+2",
                               "-background", "white", target])


times = [0.0, 0.2, 0.3, 0.4, 0.6, 1.0]
paths = [frame(t, "alpha.water", None, (0, 1), None, "alpha_%.2f" % t,
               "t = %.2f s" % t) for t in times]
tile(paths, 3, os.path.join(out, "alpha_snapshots.png"))

# velocity magnitude and p_rgh at the wall impact (t = 0.3 s)
view.ViewSize = [600, 640]
t = 0.3
u = frame(t, "U", "Magnitude", (0, 6), "Viridis (matplotlib)", "U", "|U| (m/s)", bar=True)
p = frame(t, "p_rgh", None, (0, 4000), "Inferno (matplotlib)", "p", "p_rgh (Pa)", bar=True)
tile([u, p], 2, os.path.join(out, "velocity_pressure.png"))
shutil.rmtree(tmp)
print("wrote", out)
