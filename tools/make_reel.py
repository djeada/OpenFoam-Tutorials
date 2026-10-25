"""Render a vertical (9:16) reel/short from any finished OpenFOAM case in this repo.

Run with ParaView's batch interpreter; ffmpeg (libx264 + AAC) encodes the video:

    pvbatch tools/make_reel.py projects/08_dam_break
    pvbatch tools/make_reel.py projects/13_window_ventilation          # uses its reel.json
    pvbatch tools/make_reel.py projects/02_elbow --field vorticity --streamlines 12
    pvbatch tools/make_reel.py --all                                   # every project with results

Videos land in output/artifacts/<project>/<reel>.mp4 (git-ignored).

Each reel is an intro card, the simulation played back smoothly (fields are
interpolated in time, so sparse output still gives 30 fps motion), and a closing
card, over a procedurally generated ambient soundtrack (tools/ambient_music.py) or
your own track (--music). Frames are rendered at twice the output size (one
render, no tiling) and downscaled for clean edges and text.

Without configuration the script finds the case(s), picks a field (alpha.water,
T, U, p_rgh, p), detects 2D meshes from their `empty` patches, slices 3D meshes
through their centre, fixes one colour range for the whole clip, and takes the
intro text from the project README. A project may add a `reel.json` next to its
README to choose fields, ranges, labels, several cases stacked for comparison,
live numbers read from postProcessing files, captions timed to the simulation,
and intro/outro text; see projects/13_window_ventilation/reel.json.
Command-line options override reel.json for every reel it defines.
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from paraview.simple import (AppendDatasets, Calculator, ColorBy, CreateRenderView,
                             GetColorTransferFunction, Gradient, GetScalarBar, Line, OpenFOAMReader,
                             Plane, ResetSession, SaveScreenshot, Show, Slice,
                             StreamTracerWithCustomSource, TemporalInterpolator, Text,
                             Transform)
from vtkmodules.vtkCommonCore import vtkObject
from vtkmodules.vtkRenderingCore import vtkTextProperty, vtkTextRenderer
import vtkmodules.vtkRenderingFreeType  # noqa: F401  (registers the FreeType text backend)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ambient_music  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
OUTPUT = REPO / "output" / "artifacts"
WIDTH, HEIGHT = 1080, 1920  # output size; rendering happens at SCALE times this
SCALE = 2  # set from --supersample; pixel sizes below are in output pixels
FPS = 30
BACKGROUND = [0.04, 0.05, 0.09]
TEXT = [0.95, 0.96, 0.98]
MUTED = [0.62, 0.66, 0.74]
CARD = [0.09, 0.11, 0.18]
FIELD_PRIORITY = ("alpha.water", "T", "U", "p_rgh", "p")
# A preset name, or RGB points (value r g b ...) on a 0..1 scale stretched to the range.
COLORMAPS = {"alpha.water": [0, 0.09, 0.11, 0.17, 0.5, 0.05, 0.42, 0.72, 1, 0.38, 0.86, 1.0],
             "T": "Turbo", "U": "Turbo", "c": "Viridis (matplotlib)",
             "k": "Plasma (matplotlib)",
             # Dark where the fluid does not spin, so vortices stand out on the background.
             "vorticity": [0, 0.55, 0.85, 1.0, 0.25, 0.13, 0.38, 0.85, 0.5, 0.06, 0.07, 0.12,
                           0.75, 0.85, 0.25, 0.12, 1, 1.0, 0.80, 0.45]}
# Plain-language colour-bar titles; units are appended from the field's dimensions.
FIELD_NAMES = {"alpha.water": "Water fraction", "T": "Temperature", "U": "Flow speed",
               "p": "Pressure", "p_rgh": "Pressure (hydrostatic part removed)",
               "k": "Turbulent kinetic energy", "nut": "Turbulent viscosity",
               "vorticity": "Spin (vorticity)"}
AXES = {"x": 0, "y": 1, "z": 2}
# OpenFOAM dimension sets [kg m s K mol A cd] -> unit label for the colour bar.
UNITS = {(0, 1, -1): "m/s", (1, -1, -2): "Pa", (0, 2, -2): "m²/s²", (0, 0, -1): "1/s",
         (0, 2, -1): "m²/s", (0, 2, -3): "m²/s³"}
# Normalised screen rows (0 = bottom, 1 = top) of the portrait layout.
TITLE_Y, SUBTITLE_Y = 0.905, 0.872
CONTENT_TOP, CONTENT_BOTTOM, CONTENT_WIDTH = 0.80, 0.37, 0.90
CAPTION_TOP, CLOCK_Y, SPEED_Y, BAR_Y, FOOTER_Y = 0.322, 0.178, 0.152, 0.05, 0.016
FADE = 0.6  # seconds for captions, cards and panel dimming


# --------------------------------------------------------------------------- case discovery

def time_dirs(case):
    """Written result times of a case (reconstructed or decomposed), excluding 0."""
    def numeric(path):
        try:
            return float(path.name) > 0 and path.is_dir()
        except ValueError:
            return False
    times = [p for p in case.iterdir() if numeric(p)]
    if not times and (case / "processor0").is_dir():
        times = [p for p in (case / "processor0").iterdir() if numeric(p)]
    return times


def find_cases(project):
    """The project itself if it is a case, else case folders below it (depth <= 2)."""
    if (project / "system" / "controlDict").is_file():
        return [project]
    found = sorted(p.parent.parent for p in project.glob("*/system/controlDict"))
    found += sorted(p.parent.parent for p in project.glob("*/*/system/controlDict"))
    return found


def read_boundary(case):
    """Patch name -> type from constant/polyMesh/boundary (an ASCII dictionary)."""
    path = case / "constant" / "polyMesh" / "boundary"
    if not path.is_file():
        path = case / "processor0" / "constant" / "polyMesh" / "boundary"
    if not path.is_file():
        return {}
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(errors="replace"), flags=re.S)
    return {name: kind for name, kind in re.findall(r"(\S+)\s*\{[^{}]*?\btype\s+(\w+)\s*;", text)}


def gravity_axis(case):
    """Index and sign of the axis that points 'up' (opposite to g), or None."""
    path = case / "constant" / "g"
    if not path.is_file():
        return None
    match = re.search(r"value\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)", path.read_text())
    if not match:
        return None
    g = [float(v) for v in match.groups()]
    axis = max(range(3), key=lambda i: abs(g[i]))
    return (axis, -1.0 if g[axis] > 0 else 1.0) if g[axis] else None


def field_units(case, field):
    """Unit label from the field's `dimensions` entry in the initial conditions."""
    if field == "vorticity":
        return "1/s"
    for folder in ("0", "0.orig"):
        path = case / folder / field
        if path.is_file():
            match = re.search(r"dimensions\s*\[([^\]]*)\]", path.read_text(errors="replace"))
            if match:
                dims = [float(v) for v in match.group(1).split()] + [0] * 7
                if dims[3] == 1 and not any(dims[:3]):
                    return "K"
                return UNITS.get(tuple(int(v) for v in dims[:3]), "")
    return ""


def project_title(project):
    readme = project / "README.md"
    if readme.is_file():
        for line in readme.read_text().splitlines():
            if line.startswith("# "):
                return re.sub(r"^(Project\s+)?\d+\s*[—–:-]\s*", "", line[2:]).strip()
    return re.sub(r"^\d+_", "", project.name).replace("_", " ").capitalize()


def project_blurb(project):
    """First sentence of the README's first prose paragraph, markdown stripped."""
    readme = project / "README.md"
    if not readme.is_file():
        return ""
    for block in re.split(r"\n\s*\n", readme.read_text()):
        block = block.strip()
        if not block or block[0] in "#!|<>`[-*" or block.startswith("$$"):
            continue
        text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", " ".join(block.split()))
        text = re.sub(r"[*_`]|\\\(|\\\)|\$", "", text)
        sentence = re.split(r"(?<=[.!?])\s", text)[0]
        return sentence if len(sentence) <= 180 else sentence[:177].rsplit(" ", 1)[0] + "…"
    return ""


# --------------------------------------------------------------------------- configuration

def load_reels(project, args):
    """Reel definitions from reel.json (if any) with CLI overrides applied."""
    path = project / "reel.json"
    config = json.loads(path.read_text()) if path.is_file() else {}
    # Top-level keys are defaults for every entry of "reels" (which may override them).
    base = {key: value for key, value in config.items() if key != "reels"}
    reels = config.get("reels") or [{}]
    overrides = {key: value for key, value in dict(
        field=args.field, range=args.range, colormap=args.colormap,
        lic=True if args.lic else (False if args.no_lic else None),
        streamlines=args.streamlines,
        normal=args.normal, up=args.up, duration=args.duration, skip_before=args.skip_before,
        skip_after=args.skip_after,
        zoom=args.zoom, center=args.center, music=args.music,
    ).items() if value is not None}
    if args.title:
        base["title"] = args.title
    if args.reel:
        reels = [r for r in reels if r.get("name") in args.reel]
    result = []
    for reel in reels:
        merged = {**base, **reel, **overrides}
        if args.name:
            merged["name"] = args.name
        result.append(merged)
    if args.field and len(result) > 1:  # one explicit field -> one reel
        result = result[:1]
    return result


def resolve_cases(project, reel):
    if reel.get("cases"):
        return [dict(entry, path=project / entry["path"]) for entry in reel["cases"]]
    cases = [c for c in find_cases(project) if time_dirs(c)]
    if not cases:
        return []
    return [dict(path=cases[0], label="")]


# --------------------------------------------------------------------------- pipeline

def open_case(case, include_zero=False):
    foam = case / f"{case.name}.foam"
    if not foam.exists():
        foam.touch()
    reader = OpenFOAMReader(FileName=str(foam))
    if include_zero:
        reader.SkipZeroTime = 0
        reader.Refresh()  # the time list is only rebuilt on a refresh
    if (case / "processor0").is_dir() and not any(
            p.parent == case for p in time_dirs(case)):
        reader.CaseType = "Decomposed Case"
    reader.MeshRegions = ["internalMesh"]
    reader.Createcelltopointfiltereddata = 1  # smooth point fields for colouring and LIC
    reader.UpdatePipelineInformation()
    return reader


def pick_field(reader, wanted):
    available = list(reader.CellArrays.Available)
    if wanted:
        if wanted not in available:
            raise SystemExit(f"field {wanted!r} not found; available: {', '.join(available)}")
        return wanted
    for name in FIELD_PRIORITY:
        if name in available:
            return name
    if not available:
        raise SystemExit("no cell fields found")
    return available[0]


def bounds_of(source, time=None):
    # A temporal interpolator produces nothing unless a time is requested.
    source.UpdatePipeline(time) if time is not None else source.UpdatePipeline()
    return source.GetDataInformation().GetBounds()


def orientation(case, bounds, reel):
    """(normal axis, up axis, is_2d): view direction and screen-up direction."""
    extent = [bounds[2 * i + 1] - bounds[2 * i] for i in range(3)]
    empty = "empty" in read_boundary(case).values()
    gravity = gravity_axis(case)
    if reel.get("normal"):
        normal = AXES[reel["normal"]]
    elif empty:
        normal = min(range(3), key=lambda i: extent[i])
    else:
        normal = None
    if reel.get("up"):
        up = AXES[reel["up"]]
    elif gravity and gravity[0] != normal:
        up = gravity[0]
    else:
        up = 1 if normal != 1 else 2
    if normal is None:  # 3D: side view through the centre, looking along the shorter axis
        normal = min((i for i in range(3) if i != up), key=lambda i: extent[i])
    if normal == up:
        raise SystemExit("--normal and --up must be different axes")
    return normal, up, empty


def build_panel(case_entry, reel, field_hint, interpolate):
    """Reader -> (time interpolation) -> (slice) -> (scaled field) for one case."""
    case = case_entry["path"]
    reader = open_case(case, reel.get("include_zero", True))
    vorticity = field_hint == "vorticity"  # derived from U: spin about the view axis
    field = pick_field(reader, "U" if vorticity else field_hint)
    needs_u = reel.get("lic") or reel.get("streamlines")
    arrays = [field] + (["U"] if needs_u and field != "U" else [])
    reader.CellArrays = [a for a in arrays if a in reader.CellArrays.Available]
    times = list(reader.TimestepValues or [])
    if not ((case / "0" / field).is_file() and reel.get("include_zero", True)):
        times = [t for t in times if t > 0]  # no initial field to show
    bounds = bounds_of(reader)
    normal, up, is_2d = orientation(case, bounds, reel)
    source = reader
    if interpolate and len(times) >= 2:
        source = TemporalInterpolator(Input=source)  # smooth motion between written times
    if not is_2d:
        origin = reel.get("origin") or [(bounds[2 * i] + bounds[2 * i + 1]) / 2 for i in range(3)]
        source = Slice(Input=source)
        source.SliceType = "Plane"
        source.SliceType.Origin = origin
        source.SliceType.Normal = [1.0 if i == normal else 0.0 for i in range(3)]
    array, component = field, None
    if vorticity:
        # Only in-plane derivatives enter the normal component, so a slice is enough.
        source = Gradient(Input=source)
        source.ScalarArray = ["POINTS", "U"]
        source.ComputeGradient, source.ComputeVorticity = 0, 1
        field, array, component = "vorticity", "Vorticity", "XYZ"[normal]
    elif reel.get("scale", 1) != 1 or reel.get("offset", 0) != 0:
        source = Calculator(Input=source)
        source.ResultArrayName = "display"
        source.Function = f'{reel.get("scale", 1)}*"{field}"+{reel.get("offset", 0)}'
        array = "display"
    elif field.startswith("U"):
        component = "Magnitude"
    return dict(entry=case_entry, case=case, reader=reader, source=source, field=field,
                array=array, component=component, times=times, normal=normal, up=up,
                is_2d=is_2d, bounds=bounds_of(source, times[0] if times else None))


def add_streamlines(view, panel, spec, toward_camera):
    """Thin white streamlines from fixed seeds, drawn over the coloured field.

    `spec` is a seed-grid size (N seeds across the panel) or a dict with "grid": N, or
    "lines": [[[x, y, z], [x, y, z], seeds], ...] for seeds along lines (e.g. across an
    inlet, or from a vortex centre to a wall so each closed loop is drawn once), plus
    optional "opacity" and "direction" ("BOTH", "FORWARD", "BACKWARD"). Fixed seeds keep
    the lines continuous from frame to frame. The lines are lifted just in front of the
    coloured surface so the mesh does not hide them.
    """
    spec = spec if isinstance(spec, dict) else {"grid": 14 if spec is True else int(spec)}
    b, normal, up = list(panel["bounds"]), panel["normal"], panel["up"]
    stack_shift = panel["top"] - b[2 * up + 1]  # panels below the first are moved down
    b[2 * up], b[2 * up + 1] = panel["bottom"], panel["top"]
    right = 3 - normal - up
    mid = [(b[2 * i] + b[2 * i + 1]) / 2 for i in range(3)]
    width, height = b[2 * right + 1] - b[2 * right], b[2 * up + 1] - b[2 * up]
    if spec.get("lines"):
        seeds = []
        for start, end, count in spec["lines"]:
            line = Line(Resolution=max(1, int(count) - 1))
            line.Point1, line.Point2 = [
                [mid[i] if i == normal else point[i] + (stack_shift if i == up else 0.0)
                 for i in range(3)] for point in (start, end)]
            seeds.append(line)
        seed = seeds[0] if len(seeds) == 1 else AppendDatasets(Input=seeds)
    else:
        spacing = max(width, height) / spec.get("grid", 14)
        seed = Plane()
        origin, p_right, p_up = list(mid), list(mid), list(mid)
        origin[right], origin[up] = b[2 * right] + spacing / 2, b[2 * up] + spacing / 2
        p_right[right], p_right[up] = b[2 * right + 1] - spacing / 2, origin[up]
        p_up[right], p_up[up] = origin[right], b[2 * up + 1] - spacing / 2
        seed.Origin, seed.Point1, seed.Point2 = origin, p_right, p_up
        seed.XResolution = max(1, int(round(width / spacing)) - 1)
        seed.YResolution = max(1, int(round(height / spacing)) - 1)
    lines = StreamTracerWithCustomSource(Input=panel["source"], SeedSource=seed)
    lines.Vectors = ["POINTS", "U"]
    lines.IntegrationDirection = spec.get("direction", "BOTH")
    # Short steps keep the curves smooth on coarse meshes.
    lines.InitialStepLength, lines.MaximumStepLength = 0.05, 0.1
    lines.MaximumSteps = 20000
    lines.MaximumStreamlineLength = spec.get("length", 2 * max(width, height))
    lifted = Transform(Input=lines)
    depth = b[2 * normal + 1] - b[2 * normal]
    shift = depth / 2 + 1e-3 * max(width, height)
    lifted.Transform.Translate = [shift * toward_camera if i == normal else 0.0 for i in range(3)]
    display = Show(lifted, view)
    display.ColorArrayName = ["POINTS", ""]  # solid colour
    display.AmbientColor = display.DiffuseColor = [1.0, 1.0, 1.0]
    display.Ambient, display.Diffuse = 1.0, 0.0
    display.Opacity = spec.get("opacity", 0.45)
    display.LineWidth = 1.5 * SCALE
    return lifted


def screen_axes(normal, up):
    right = 3 - normal - up
    # Put the camera on the side of the plane where `right` increases to the right.
    r, u = [0.0] * 3, [0.0] * 3
    r[right], u[up] = 1.0, 1.0
    cross = [r[1] * u[2] - r[2] * u[1], r[2] * u[0] - r[0] * u[2], r[0] * u[1] - r[1] * u[0]]
    return right, cross[normal]


def sample_range(panels, times):
    """Colour range from up to six time steps across all panels."""
    lo, hi = math.inf, -math.inf
    picks = sorted({times[int(round(i * (len(times) - 1) / 5))] for i in range(6)}) if times else [None]
    for panel in panels:
        for t in picks:
            panel["source"].UpdatePipeline(t) if t is not None else panel["source"].UpdatePipeline()
            info = panel["source"].GetPointDataInformation().GetArray(panel["array"])
            if info is None:
                continue
            a, b = info.GetRange({"Magnitude": -1, None: 0, "X": 0, "Y": 1, "Z": 2}[panel["component"]])
            lo, hi = min(lo, a), max(hi, b)
    if not math.isfinite(lo) or lo == hi:
        return [0.0, 1.0] if not math.isfinite(lo) else [lo - 0.5, hi + 0.5]
    return nice_range(lo, hi)


def nice_range(lo, hi):
    """Round outwards to readable bar labels; snap volume fractions to [0, 1]."""
    span = hi - lo
    if -1e-3 <= lo and hi <= 1 + 1e-3 and span > 0.5:
        return [0.0, 1.0]
    step = 10 ** math.floor(math.log10(span)) / 4
    lo = 0.0 if 0 <= lo < 0.02 * span or (lo < 0 and -lo < 1e-3 * span) else math.floor(lo / step) * step
    return [lo, math.ceil(hi / step) * step]


def bar_ticks(lo, hi):
    """4-6 round tick values spanning [lo, hi], readable on a phone screen."""
    span = hi - lo
    for factor in (1, 2, 2.5, 5, 10):
        step = factor * 10 ** math.floor(math.log10(span / 5))
        if span / step <= 6:
            break
    first = math.ceil(lo / step - 1e-9) * step
    ticks = [first + i * step for i in range(int((hi - first) / step + 1e-9) + 1)]
    for end in (lo, hi):  # label the bar's ends unless a tick is already close
        if min(abs(end - t) for t in ticks) >= 0.4 * step:
            ticks.append(end)
    return sorted(ticks)


# --------------------------------------------------------------------------- overlays

def read_series(path, column, scale, offset):
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            values = line.replace("(", " ").replace(")", " ").split()
            rows.append((float(values[0]), offset + scale * float(values[column])))
    return rows


def series_at(rows, t):
    if not rows:
        return None
    if t <= rows[0][0]:
        return rows[0][1]
    for (t0, v0), (t1, v1) in zip(rows, rows[1:]):
        if t <= t1:
            return v0 + (v1 - v0) * (t - t0) / (t1 - t0) if t1 > t0 else v1
    return rows[-1][1]


def overlay_series(panel):
    spec = panel["entry"].get("overlay")
    if not spec:
        return None
    files = sorted(panel["case"].glob(spec["file"]))
    if not files:
        print(f"  overlay file {spec['file']} not found in {panel['case']}; skipping it")
        return None
    rows = read_series(files[0], spec.get("column", 1), spec.get("scale", 1), spec.get("offset", 0))
    if "initial" in spec:  # value before the first written time, e.g. 0% fresh air at t=0
        rows.insert(0, (0.0, spec["initial"]))
    return rows


class Label:
    """Screen text anchored at its left edge, or centred on x when centre=True.

    ParaView positions text by the lower-left corner of its box, so centring and
    vertical anchoring (anchor="bottom" | "top" | "middle" of y) need the rendered
    size, measured with VTK's own text renderer (72 dpi = pixels at the view size).
    `box` is the background opacity of an optional rounded box; set() also takes
    an opacity for fades.
    """

    def __init__(self, view, string, size, color, x, y, bold=False, centre=False,
                 box=0.0, anchor="bottom", wrap=None):
        self.source = Text(Text="")
        self.display = Show(self.source, view)
        self.display.FontFamily = "Arial"
        self.display.FontSize = size * SCALE
        self.display.Bold = int(bold)
        self.display.Color = color
        self.display.WindowLocation = "Any Location"
        self.display.Justification = "Center" if centre else "Left"
        self.display.LineSpacing = 1.15
        self.pad = 0
        if box:
            self.pad = 22 * SCALE
            self.display.ShowBorder = "Always"
            self.display.BorderThickness = 1.5 * SCALE
            self.display.CornerRadius = 0.2
            self.display.Padding = self.pad
        self.box = box
        self.prop = vtkTextProperty()
        self.prop.SetFontFamilyToArial()
        self.prop.SetFontSize(size * SCALE)
        self.prop.SetBold(int(bold))
        self.prop.SetLineSpacing(1.15)
        self.x, self.y, self.centre, self.anchor, self.wrap = x, y, centre, anchor, wrap
        self.set(string)

    def measure(self, string):
        bbox = [0, 0, 0, 0]
        vtkTextRenderer.GetInstance().GetBoundingBox(self.prop, string, bbox, 72)
        return bbox[1] - bbox[0], bbox[3] - bbox[2]

    def wrapped(self, string):
        if not self.wrap:
            return string
        lines = []
        for paragraph in string.split("\n"):
            line = ""
            for word in paragraph.split(" "):
                trial = f"{line} {word}".strip()
                if line and self.measure(trial)[0] > self.wrap * WIDTH * SCALE - 2 * self.pad:
                    lines.append(line)
                    line = word
                else:
                    line = trial
            lines.append(line)
        return "\n".join(lines)

    def set(self, string, opacity=1.0):
        string = self.wrapped(string)
        width, height = self.measure(string) if string else (0, 0)
        x, y = self.x, self.y
        if self.centre:
            x -= (width / 2 + self.pad) / (WIDTH * SCALE)
        if self.anchor == "top":
            y -= (height + 2 * self.pad) / (HEIGHT * SCALE)
        elif self.anchor == "middle":
            y -= (height / 2 + self.pad) / (HEIGHT * SCALE)
        self.source.Text = string
        self.display.Position = [x, y]
        self.display.Opacity = opacity
        self.display.Visibility = int(bool(string) and opacity > 0.01)
        if self.box:
            self.display.BackgroundColor = CARD + [self.box * opacity]
            self.display.BorderColor = [c * opacity + b * (1 - opacity)
                                        for c, b in zip([0.25, 0.30, 0.42], BACKGROUND)]


def is_steady(case):
    """Steady solvers (simpleFoam, ...) count iterations, not seconds."""
    schemes = case / "system" / "fvSchemes"
    return schemes.is_file() and "steadyState" in schemes.read_text(errors="replace")


def format_time(t, end, steady=False):
    if steady:
        return f"iteration {t:.0f}"
    if end >= 120:
        minutes, seconds = divmod(int(round(t)), 60)
        return f"{minutes} min {seconds:02d} s"
    if end >= 10:
        return f"t = {t:.1f} s"
    return f"t = {t:.2f} s"


def speed_note(sim_span, video_span):
    if sim_span <= 0 or video_span <= 0:
        return ""
    ratio = sim_span / video_span
    if ratio >= 1.05:
        return f"played {ratio:.0f}× faster than real time" if ratio >= 10 else \
            f"played {ratio:.1f}× faster than real time"
    if ratio <= 0.95:
        return f"slowed down {1 / ratio:.0f}×"
    return "real time"


def fade(now, start, end):
    """0..1 visibility for something shown during [start, end) with FADE in/out."""
    if now < start or now >= end:
        return 0.0
    return max(0.0, min(1.0, (now - start) / FADE, (end - now) / FADE))


# --------------------------------------------------------------------------- audio/video

class Encoder:
    """Streams full-size PNG renders into ffmpeg, which downscales them."""

    def __init__(self, path):
        self.tmp = Path(tempfile.mkdtemp(prefix="reel_"))
        command = ["ffmpeg", "-y", "-loglevel", "error", "-f", "image2pipe",
                   "-framerate", str(FPS), "-c:v", "png", "-i", "-",
                   "-vf", f"scale={WIDTH}:{HEIGHT}:flags=lanczos,format=yuv420p",
                   "-c:v", "libx264", "-preset", "slow", "-crf", "16",
                   "-profile:v", "high", "-tune", "animation", "-r", str(FPS),
                   "-movflags", "+faststart", str(path)]
        self.process = subprocess.Popen(command, stdin=subprocess.PIPE)

    def add(self, view):
        frame = self.tmp / "frame.png"
        # ImageResolution equal to the view size: one render, no tile seams.
        SaveScreenshot(str(frame), view, CompressionLevel="1",
                       ImageResolution=[WIDTH * SCALE, HEIGHT * SCALE])
        self.process.stdin.write(frame.read_bytes())

    def close(self):
        self.process.stdin.close()
        if self.process.wait() != 0:
            raise SystemExit("ffmpeg failed while encoding frames")
        shutil.rmtree(self.tmp, ignore_errors=True)


def add_music(silent, output, seconds, music):
    """Mux a soundtrack (generated ambient, or a user file looped to length)."""
    tmp = Path(tempfile.mkdtemp(prefix="reel_audio_"))
    if music in (None, True, "ambient"):
        track = tmp / "ambient.wav"
        ambient_music.write_wav(track, ambient_music.compose(seconds + 0.5))
        source = ["-i", str(track)]
    else:
        source = ["-stream_loop", "-1", "-i", str(Path(music).expanduser())]
    fade_out = max(0.0, seconds - 3.0)
    command = ["ffmpeg", "-y", "-loglevel", "error", "-i", str(silent), *source,
               "-map", "0:v", "-map", "1:a", "-c:v", "copy",
               "-af", f"atrim=0:{seconds:.3f},loudnorm=I=-16:TP=-1.5:LRA=11,"
                      f"afade=t=in:d=1.5,afade=t=out:st={fade_out:.3f}:d=3,aresample=48000",
               "-c:a", "aac", "-b:a", "192k", "-t", f"{seconds:.3f}",
               "-movflags", "+faststart", str(output)]
    subprocess.run(command, check=True)
    shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------- render

def render_reel(project, reel, args):
    entries = resolve_cases(project, reel)
    if not entries:
        print(f"skip {project.name}: no finished case with written time steps")
        return None
    missing = [e["path"] for e in entries if not e["path"].is_dir() or not time_dirs(e["path"])]
    if missing:
        print(f"skip {project.name}: no results in {', '.join(map(str, missing))}")
        return None

    ResetSession()
    view = CreateRenderView()
    view.ViewSize = [WIDTH * SCALE, HEIGHT * SCALE]
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    view.UseColorPaletteForBackground = 0
    view.Background = BACKGROUND

    panels = [build_panel(entry, reel, reel.get("field"), not args.no_interpolate)
              for entry in entries]
    first = panels[0]
    normal, up = first["normal"], first["up"]
    right, sign = screen_axes(normal, up)

    # Stack the panels top to bottom along the screen-up axis.
    gap = (0.28 if len(panels) > 1 else 0.0) * max(
        p["bounds"][2 * up + 1] - p["bounds"][2 * up] for p in panels)
    cursor = first["bounds"][2 * up + 1]
    for panel in panels:
        height = panel["bounds"][2 * up + 1] - panel["bounds"][2 * up]
        shift = cursor - panel["bounds"][2 * up + 1]
        if shift:
            moved = Transform(Input=panel["source"])
            moved.Transform.Translate = [shift if i == up else 0.0 for i in range(3)]
            panel["source"] = moved
        panel["top"], panel["bottom"] = cursor, cursor - height
        cursor -= height + gap
    lo_up, hi_up = panels[-1]["bottom"], panels[0]["top"]
    lo_r = min(p["bounds"][2 * right] for p in panels)
    hi_r = max(p["bounds"][2 * right + 1] for p in panels)

    # Fit the stack into the content band of the portrait frame.
    aspect = WIDTH / HEIGHT
    band = CONTENT_TOP - CONTENT_BOTTOM
    scale = max((hi_up - lo_up) / (2 * band), (hi_r - lo_r) / (2 * CONTENT_WIDTH * aspect))
    # Wide, flat domains would be a thin strip in portrait: zoom in and crop the sides.
    zoom = reel.get("zoom") or min(3.0, max(1.0, 0.4 * band / ((hi_up - lo_up) / (2 * scale))))
    scale /= zoom
    centre_r, middle_up = (hi_r + lo_r) / 2, (hi_up + lo_up) / 2
    if reel.get("center"):
        centre_r, middle_up = reel["center"][right], reel["center"][up]
    centre_up = middle_up - ((CONTENT_TOP + CONTENT_BOTTOM) / 2 - 0.5) * 2 * scale
    depth_mid = (first["bounds"][2 * normal] + first["bounds"][2 * normal + 1]) / 2

    def camera(zoom=1.0):
        focal = [0.0] * 3
        focal[up], focal[right], focal[normal] = centre_up, centre_r, depth_mid
        position = list(focal)
        position[normal] += 10 * max(hi_r - lo_r, hi_up - lo_up, 1e-6) * (1 if sign > 0 else -1)
        view.CameraFocalPoint = focal
        view.CameraPosition = position
        view.CameraViewUp = [1.0 if i == up else 0.0 for i in range(3)]
        view.CameraParallelScale = scale / zoom

    def to_screen(value_up, value_r):
        return (0.5 + (value_r - centre_r) / (2 * scale * aspect),
                0.5 + (value_up - centre_up) / (2 * scale))

    times = first["times"]
    if reel.get("skip_before"):
        times = [t for t in times if t >= reel["skip_before"]]
    if reel.get("skip_after"):
        times = [t for t in times if t <= reel["skip_after"]]
    value_range = reel.get("range") or sample_range(panels, times)
    lut = GetColorTransferFunction(first["array"])
    lut.AutomaticRescaleRangeMode = "Never"
    colormap = reel.get("colormap") or COLORMAPS.get(first["field"], "Cool to Warm")
    if isinstance(colormap, str):
        lut.ApplyPreset(colormap, True)
    else:
        lut.RGBPoints = colormap

    for panel in panels:
        display = Show(panel["source"], view)
        display.Ambient, display.Diffuse = 1.0, 0.0  # flat colours
        if reel.get("lic"):
            display.SetRepresentationType("Surface LIC")
            display.SelectInputVectors = ["POINTS", "U"]
            # A light blend keeps the colour scale readable under the streak texture.
            # Contrast enhancement re-normalises every frame and makes video flicker.
            display.EnhanceContrast = "Off"
            display.ColorMode = "Blend"
            display.LICIntensity = 0.4
            display.NumberOfSteps = 60
            display.StepSize = 0.5
            display.NoiseGrainSize = 2 * SCALE
        if panel["component"]:
            ColorBy(display, ("POINTS", panel["array"], panel["component"]))
        else:
            ColorBy(display, ("POINTS", panel["array"]))
        display.LookupTable = lut
        display.SetScalarBarVisibility(view, False)
        if reel.get("streamlines") and not reel.get("lic"):
            add_streamlines(view, panel, reel["streamlines"], sign)
    lut.RescaleTransferFunction(*value_range)  # after ColorBy, which rescales to the data

    # Header, panel labels, captions, cards, clock, colour bar.
    title = reel.get("title") or project_title(project)
    Label(view, title, 52 if len(title) < 24 else 40, TEXT, 0.5, TITLE_Y, bold=True, centre=True)
    if reel.get("subtitle"):
        Label(view, reel["subtitle"], 26, MUTED, 0.5, SUBTITLE_Y, centre=True)
    labels = []
    for panel in panels:
        x, y = to_screen(panel["top"], lo_r)
        label = Label(view, panel["entry"].get("label", ""), 30, TEXT,
                      max(x, 0.04), min(y + 0.01, CONTENT_TOP + 0.03), bold=True)
        labels.append((label, panel["entry"].get("label", ""), overlay_series(panel),
                       panel["entry"].get("overlay", {}).get("format", "{:.1f}")))
    caption = Label(view, "", 34, TEXT, 0.5, CAPTION_TOP, centre=True, box=0.9, anchor="top",
                    wrap=0.88)
    # Cards sit over the panels; an opaque box rather than dimming the panels,
    # because Surface LIC stops rendering when its display is translucent.
    card = Label(view, "", 40, TEXT, 0.5, (CONTENT_TOP + CONTENT_BOTTOM) / 2, bold=True,
                 centre=True, box=0.96, anchor="middle", wrap=0.86)
    clock = Label(view, "", 46, TEXT, 0.5, CLOCK_Y, bold=True, centre=True)
    speed = Label(view, "", 22, MUTED, 0.5, SPEED_Y, centre=True)
    Label(view, "OpenFOAM · " + project.name, 20, MUTED, 0.5, FOOTER_Y, centre=True)

    bar = GetScalarBar(lut, view)
    bar.Visibility = 1
    bar.Orientation = "Horizontal"
    bar.WindowLocation = "Any Location"
    bar.Position = [0.1, BAR_Y]
    bar.ScalarBarLength = 0.8
    bar.ScalarBarThickness = 22 * SCALE
    if reel.get("bar_title"):
        bar.Title = reel["bar_title"]
    else:
        units = field_units(first["case"], first["field"])
        name = FIELD_NAMES.get(first["field"], first["field"])
        bar.Title = f"{name} ({units})" if units else name
    bar.ComponentTitle = ""
    bar.TitleFontSize, bar.LabelFontSize = 28 * SCALE, 22 * SCALE
    bar.TitleColor = bar.LabelColor = TEXT
    bar.TitleFontFamily = bar.LabelFontFamily = "Arial"
    bar.TitleBold = 1
    bar.RangeLabelFormat = bar.LabelFormat = reel.get("bar_format", "%-.3g")
    bar.UseCustomLabels = 1
    bar.CustomLabels = bar_ticks(*value_range)
    bar.AddRangeLabels = 0
    bar.DrawTickMarks = 0

    # Timeline: intro card -> simulation (or a slow zoom for steady results) -> outro card.
    intro = reel.get("intro", {"text": project_blurb(project), "seconds": 4})
    outro = reel.get("outro", {"text": "", "seconds": 3})
    intro_s = float(intro.get("seconds", 4)) if intro.get("text") else 0.0
    outro_s = float(outro.get("seconds", 3))
    transient = len(times) >= 2
    main_s = float(reel.get("duration") or (35 if transient else 10))
    t_start, t_end = (times[0], times[-1]) if times else (0.0, 0.0)
    main_n = max(2, int(round(main_s * FPS)))

    def sim_time(index):  # simulation time shown at main frame `index`
        return t_start + (t_end - t_start) * index / (main_n - 1) if transient else t_end

    def video_time(sim):  # inverse, for captions given in simulation time
        if not transient or t_end == t_start:
            return intro_s
        return intro_s + main_s * (sim - t_start) / (t_end - t_start)

    captions = [(video_time(c["from"]), video_time(c["to"]), c["text"])
                for c in reel.get("captions", [])]
    end_values = [series_at(series, t_end) for _, _, series, _ in labels if series]
    outro_text = outro.get("text", "")
    if outro_text:
        try:
            outro_text = outro_text.format(*end_values)
        except (IndexError, ValueError):
            print("  outro needs overlay values that are missing; showing it unformatted")
    steady = is_steady(first["case"])
    if steady:
        speed.set("solver converging towards a steady state" if transient else "")
    else:
        speed.set(speed_note(t_end - t_start, main_s) if transient else "")

    schedule = []  # (video seconds, simulation time, camera zoom)
    for i in range(int(round(intro_s * FPS))):
        schedule.append((i / FPS, t_start, 1.0))
    for i in range(main_n):
        zoom_now = 1.0 if transient else 1.0 + 0.12 * i / (main_n - 1)
        schedule.append((intro_s + i / FPS, sim_time(i), zoom_now))
    last_zoom = schedule[-1][2]
    for i in range(int(round(outro_s * FPS))):
        schedule.append((intro_s + main_s + i / FPS, t_end, last_zoom))
    total_s = len(schedule) / FPS

    name = reel.get("name") or first["field"].replace(".", "_")
    out_dir = OUTPUT / project.name
    out_dir.mkdir(parents=True, exist_ok=True)

    def draw(now, t, zoom_now):
        view.ViewTime = t
        camera(zoom_now)
        clock.set(format_time(t, t_end, steady) if transient else "")
        for label, string, series, fmt in labels:
            value = series_at(series, t) if series else None
            label.set(string + (f"  ·  {fmt.format(value)}" if value is not None else ""))
        shown = [(fade(now, a, b), text) for a, b, text in captions]
        alpha, text = max(shown, default=(0.0, ""))
        caption.set(text if alpha > 0 else "", alpha)
        intro_alpha = fade(now, -FADE, intro_s) if intro_s else 0.0
        outro_alpha = fade(now, intro_s + main_s, total_s + FADE) if outro_text else 0.0
        if intro_alpha:
            card.set(intro.get("text", ""), intro_alpha)
        else:
            card.set(outro_text if outro_alpha else "", outro_alpha)

    if args.still:  # a single PNG of the final state, without cards or captions
        still = out_dir / f"{name}.png"
        draw(intro_s + main_s - 1 / FPS - FADE, t_end, last_zoom)
        caption.set("")
        SaveScreenshot(str(still), view, ImageResolution=[WIDTH * SCALE, HEIGHT * SCALE])
        print(f"  wrote {still.relative_to(REPO)}")
        return still

    if args.segment:  # quick silent preview of part of the timeline
        schedule = [f for f in schedule if args.segment[0] <= f[0] < args.segment[1]]
        name += "_preview"
    print(f"{project.name}/{name}: {len(schedule)} frames ({total_s:.1f} s), field "
          f"{first['field']}, range {value_range[0]:.4g}..{value_range[1]:.4g}")
    video = out_dir / f"{name}.mp4"
    silent = out_dir / f".{name}.silent.mp4"
    encoder = Encoder(silent)
    # LIC complains about zero-velocity frames (e.g. the initial state); keep logs readable.
    vtkObject.GlobalWarningDisplayOff()
    for index, (now, t, zoom_now) in enumerate(schedule):
        draw(now, t, zoom_now)
        encoder.add(view)
        if index % 150 == 0:
            print(f"  frame {index + 1}/{len(schedule)}", flush=True)
    encoder.close()
    vtkObject.GlobalWarningDisplayOn()
    music = reel.get("music", "ambient")
    if music is False or args.no_music or args.segment:
        os.replace(silent, video)
    else:
        add_music(silent, video, total_s, music)
        silent.unlink()
    print(f"  wrote {video.relative_to(REPO)}")
    return video


# --------------------------------------------------------------------------- CLI

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("projects", nargs="*", type=Path, help="project directories")
    parser.add_argument("--all", action="store_true", help="every project under projects/")
    parser.add_argument("--reel", action="append", help="only the named reel(s) from reel.json")
    parser.add_argument("--field", help="cell field to colour by (overrides reel.json)")
    parser.add_argument("--range", nargs=2, type=float, metavar=("MIN", "MAX"))
    parser.add_argument("--colormap", help="ParaView preset name, e.g. 'Turbo'")
    parser.add_argument("--lic", action="store_true",
                        help="texture the field with U (Surface LIC); best for steady or slowly "
                             "changing flow, it shimmers and greys the colours in fast transients")
    parser.add_argument("--no-lic", action="store_true", help="plain colours, no flow texture")
    parser.add_argument("--streamlines", type=int, metavar="N",
                        help="overlay white streamlines seeded on an N-wide grid (0 = off)")
    parser.add_argument("--normal", choices=AXES, help="view/slice normal axis")
    parser.add_argument("--up", choices=AXES, help="axis pointing up on screen")
    parser.add_argument("--duration", type=float,
                        help="seconds of video for the simulation itself (default 35)")
    parser.add_argument("--zoom", type=float, help="magnify the view (crops the sides)")
    parser.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help="world point to centre the view on (useful with --zoom)")
    parser.add_argument("--skip-before", type=float, help="ignore times earlier than this")
    parser.add_argument("--skip-after", type=float, help="ignore times later than this")
    parser.add_argument("--title", help="headline text")
    parser.add_argument("--name", help="output file name (without .mp4)")
    parser.add_argument("--music", help="audio file to use instead of the generated ambient track")
    parser.add_argument("--no-music", action="store_true", help="silent video")
    parser.add_argument("--no-interpolate", action="store_true",
                        help="snap to written times (needed if the mesh topology changes)")
    parser.add_argument("--supersample", type=int, default=2,
                        help="render at N times the output size (default 2; 1 = fast draft)")
    parser.add_argument("--segment", nargs=2, type=float, metavar=("START", "END"),
                        help="render only this part of the video (seconds), silent, as *_preview.mp4")
    parser.add_argument("--still", action="store_true", help="write a PNG of the last time only")
    args = parser.parse_args(argv)
    global SCALE
    SCALE = max(1, args.supersample)
    if shutil.which("ffmpeg") is None:
        parser.error("ffmpeg is required on PATH")
    projects = sorted((REPO / "projects").glob("[0-9]*_*")) if args.all else args.projects
    if not projects:
        parser.error("give one or more project directories, or --all")
    videos = []
    for project in projects:
        project = project.resolve()
        for reel in load_reels(project, args):
            video = render_reel(project, reel, args)
            if video:
                videos.append(video)
    print(f"{len(videos)} reel(s) written to {OUTPUT.relative_to(REPO)}/")
    return 0 if videos else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
