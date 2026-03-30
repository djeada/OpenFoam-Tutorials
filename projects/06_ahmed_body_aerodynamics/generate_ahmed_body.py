#!/usr/bin/env python3
"""Generate a simplified Ahmed body STL for OpenFOAM snappyHexMesh.

Dimensions (standard Ahmed body at 25° slant):
  Length:           1.044 m
  Width:            0.389 m
  Height:           0.288 m
  Ground clearance: 0.050 m
  Rear slant angle: 25°
  Slant length:     0.222 m

The body is centred in y and sits at z = 0.05 m (ground clearance).
Front-edge rounding and stilts are omitted for simplicity.

Usage:
    python3 generate_ahmed_body.py
    python3 generate_ahmed_body.py --output constant/triSurface/ahmed_body.stl
    python3 generate_ahmed_body.py --slant-angle 35
"""

import argparse
import math
from pathlib import Path


def build_vertices(slant_angle_deg: float):
    """Return the 10 vertices of the simplified Ahmed body."""
    L = 1.044       # total length
    W = 0.389       # width
    H = 0.288       # height
    GC = 0.050      # ground clearance
    SL = 0.222      # slant surface length

    half_w = W / 2.0
    z_bot = GC
    z_top = GC + H

    # Rear slant geometry
    slant_rad = math.radians(slant_angle_deg)
    slant_dx = SL * math.cos(slant_rad)    # horizontal projection
    slant_dz = SL * math.sin(slant_rad)    # vertical drop

    x_slant = L - slant_dx                 # where slant begins on top
    z_slant = z_top - slant_dz             # z at bottom of slant

    return {
        # Bottom rectangle
        "BFL": (0.0,  -half_w, z_bot),
        "BFR": (0.0,   half_w, z_bot),
        "BRL": (L,    -half_w, z_bot),
        "BRR": (L,     half_w, z_bot),
        # Top rectangle (shortened by slant)
        "TFL": (0.0,  -half_w, z_top),
        "TFR": (0.0,   half_w, z_top),
        "TSL": (x_slant, -half_w, z_top),
        "TSR": (x_slant,  half_w, z_top),
        # Slant bottom edge
        "SBL": (L,    -half_w, z_slant),
        "SBR": (L,     half_w, z_slant),
    }


def triangle_normal(v0, v1, v2):
    """Compute outward normal (not necessarily unit) from three vertices."""
    e1 = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
    e2 = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
    nx = e1[1] * e2[2] - e1[2] * e2[1]
    ny = e1[2] * e2[0] - e1[0] * e2[2]
    nz = e1[0] * e2[1] - e1[1] * e2[0]
    mag = math.sqrt(nx * nx + ny * ny + nz * nz)
    if mag < 1e-30:
        return (0.0, 0.0, 0.0)
    return (nx / mag, ny / mag, nz / mag)


def build_triangles(v):
    """Return list of (v0, v1, v2) triangles with outward-facing normals."""
    tris = []

    # Bottom face (normal -z)
    tris.append((v["BFL"], v["BRR"], v["BRL"]))
    tris.append((v["BFL"], v["BFR"], v["BRR"]))

    # Top face (normal +z)
    tris.append((v["TFL"], v["TSL"], v["TSR"]))
    tris.append((v["TFL"], v["TSR"], v["TFR"]))

    # Slant face (normal up-and-rearward)
    tris.append((v["TSL"], v["SBL"], v["SBR"]))
    tris.append((v["TSL"], v["SBR"], v["TSR"]))

    # Front face (normal -x)
    tris.append((v["BFL"], v["TFL"], v["TFR"]))
    tris.append((v["BFL"], v["TFR"], v["BFR"]))

    # Rear base face (normal +x)
    tris.append((v["BRR"], v["SBL"], v["BRL"]))
    tris.append((v["BRR"], v["SBR"], v["SBL"]))

    # Left side (normal -y): pentagon BFL, BRL, SBL, TSL, TFL
    tris.append((v["BFL"], v["BRL"], v["SBL"]))
    tris.append((v["BFL"], v["SBL"], v["TSL"]))
    tris.append((v["BFL"], v["TSL"], v["TFL"]))

    # Right side (normal +y): pentagon BFR, TFR, TSR, SBR, BRR
    tris.append((v["BFR"], v["TFR"], v["TSR"]))
    tris.append((v["BFR"], v["TSR"], v["SBR"]))
    tris.append((v["BFR"], v["SBR"], v["BRR"]))

    return tris


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
        description="Generate a simplified Ahmed body STL for OpenFOAM."
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
    args = parser.parse_args()

    outpath = Path(args.output)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    vertices = build_vertices(args.slant_angle)
    triangles = build_triangles(vertices)
    write_ascii_stl(outpath, triangles)

    print(f"Ahmed body STL written to {outpath}")
    print(f"  Slant angle: {args.slant_angle}°")
    print(f"  Triangles:   {len(triangles)}")


if __name__ == "__main__":
    main()
