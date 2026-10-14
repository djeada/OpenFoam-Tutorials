#!/usr/bin/env python3
"""Generate the Ahmed body STL (without stilts) for OpenFOAM snappyHexMesh.

Dimensions (standard Ahmed body, full scale, 1:1):
  Length:           1.044 m
  Width:            0.389 m
  Height:           0.288 m
  Ground clearance: 0.050 m
  Front radius:     0.100 m on all four front edges (top, bottom, both sides)
  Rear slant angle: 25 degrees (default), slant length 0.222 m

The body is centred in y, the nose is at x = 0 and the underside is at
z = 0.05 m. The four stilts of the experiment are omitted: the body is
free-standing above the ground plane.

Construction: the body is {z_bottom(x) <= z <= z_top(x), |y| <= half_width(x)},
i.e. a rectangular cross-section whose size depends on x only. Rounding the
front edges with radius R shrinks the cross-section by
R - sqrt(R^2 - (R - x)^2) on all four sides for 0 <= x <= R (quarter circle in
the side view and in the plan view); the rear slant lowers the top surface.
The longitudinal edges and the rear edges stay sharp, as on the real body.
The corner where two fillets meet is the intersection of two cylinders instead
of a sphere octant; the difference is below the mesh resolution.

Usage:
    python3 generate_ahmed_body.py
    python3 generate_ahmed_body.py --output constant/triSurface/ahmed_body.stl
    python3 generate_ahmed_body.py --slant-angle 35
    python3 generate_ahmed_body.py --nose-radius 0      # sharp front edges
"""

import argparse
import math
from pathlib import Path

L = 1.044       # total length
W = 0.389       # width
H = 0.288       # height
GC = 0.050      # ground clearance
SL = 0.222      # slant surface length
FILLET_SEGMENTS = 12   # 7.5 degrees per segment on the quarter circle


def stations(slant_angle_deg, radius):
    """Return a list of (x, half_width, z_bottom, z_top) cross-sections."""
    a = math.radians(slant_angle_deg)
    x_slant = L - SL * math.cos(a)           # start of the slant on the roof
    z_top = GC + H
    z_end = z_top - SL * math.sin(a)         # top of the rear base
    out = []

    def section(x, inset, z_t):
        return (x, W / 2 - inset, GC + inset, z_t - inset)

    if radius > 0:
        for k in range(FILLET_SEGMENTS + 1):
            th = (math.pi / 2) * k / FILLET_SEGMENTS
            x = radius * (1 - math.cos(th))
            inset = radius * (1 - math.sin(th))
            out.append(section(x, inset, z_top))
    else:
        out.append(section(0.0, 0.0, z_top))
    out.append(section(x_slant, 0.0, z_top))
    out.append((L, W / 2, GC, z_end))
    # the sections must not be collapsed or reversed
    assert all(b[0] > a_[0] for a_, b in zip(out, out[1:])), "radius too large"
    return out


def build_triangles(slant_angle_deg, radius):
    """Closed triangulated surface with outward-facing orientation."""
    secs = stations(slant_angle_deg, radius)

    def corners(s):
        x, h, zb, zt = s
        # loop order: bottom-left, bottom-right, top-right, top-left
        return [(x, -h, zb), (x, h, zb), (x, h, zt), (x, -h, zt)]

    rings = [corners(s) for s in secs]
    tris = []
    for r0, r1 in zip(rings, rings[1:]):
        for j in range(4):
            k = (j + 1) % 4
            a, b, c, d = r0[j], r0[k], r1[k], r1[j]
            tris.append((a, b, c))
            tris.append((a, c, d))
    # caps: fan, reversed loop at the front, forward loop at the rear
    f, r = rings[0], rings[-1]
    tris.append((f[0], f[3], f[2]))
    tris.append((f[0], f[2], f[1]))
    tris.append((r[0], r[1], r[2]))
    tris.append((r[0], r[2], r[3]))

    if signed_volume(tris) < 0:
        tris = [(a, c, b) for a, b, c in tris]
    return tris


def signed_volume(tris):
    v = 0.0
    for a, b, c in tris:
        v += (a[0] * (b[1] * c[2] - b[2] * c[1])
              - a[1] * (b[0] * c[2] - b[2] * c[0])
              + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0
    return v


def triangle_normal(v0, v1, v2):
    """Unit normal from the vertex order (right-hand rule)."""
    e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
    e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
    nx = e1[1] * e2[2] - e1[2] * e2[1]
    ny = e1[2] * e2[0] - e1[0] * e2[2]
    nz = e1[0] * e2[1] - e1[1] * e2[0]
    mag = math.sqrt(nx * nx + ny * ny + nz * nz)
    if mag < 1e-30:
        return (0.0, 0.0, 0.0)
    return (nx / mag, ny / mag, nz / mag)


def write_ascii_stl(filepath: Path, triangles: list, solid_name: str = "ahmed_body"):
    """Write triangles to an ASCII STL file."""
    with open(filepath, "w") as f:
        f.write(f"solid {solid_name}\n")
        for v0, v1, v2 in triangles:
            n = triangle_normal(v0, v1, v2)
            f.write(f"  facet normal {n[0]:.6e} {n[1]:.6e} {n[2]:.6e}\n")
            f.write("    outer loop\n")
            for vx in (v0, v1, v2):
                f.write(f"      vertex {vx[0]:.6e} {vx[1]:.6e} {vx[2]:.6e}\n")
            f.write("    endloop\n")
            f.write("  endfacet\n")
        f.write(f"endsolid {solid_name}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate the Ahmed body STL (no stilts) for OpenFOAM."
    )
    parser.add_argument(
        "--output",
        default="constant/triSurface/ahmed_body.stl",
        help="Output STL file path (default: constant/triSurface/ahmed_body.stl)",
    )
    parser.add_argument(
        "--slant-angle",
        type=float,
        default=25.0,
        help="Rear slant angle in degrees (default: 25)",
    )
    parser.add_argument(
        "--nose-radius",
        type=float,
        default=0.1,
        help="Radius of the four front edges in m (default: 0.1, 0 = sharp)",
    )
    args = parser.parse_args()

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    triangles = build_triangles(args.slant_angle, args.nose_radius)
    write_ascii_stl(outpath, triangles)

    print(f"Ahmed body STL written to {outpath}")
    print(f"  Slant angle: {args.slant_angle}°")
    print(f"  Nose radius: {args.nose_radius} m")
    print(f"  Triangles:   {len(triangles)}")
    print(f"  Volume:      {signed_volume(triangles):.6f} m3")


if __name__ == "__main__":
    main()
