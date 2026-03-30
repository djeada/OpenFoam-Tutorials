import argparse
import math
import struct
from pathlib import Path
from typing import Optional


def cosine_spacing(chord_length: float, n_points: int):
    return [
        0.5 * chord_length * (1.0 - math.cos(beta))
        for beta in [math.pi * i / (n_points - 1) for i in range(n_points)]
    ]


def naca_0012_airfoil(t: float = 0.12, c: float = 1.0, n_points: int = 200):
    x_values = cosine_spacing(c, n_points)
    upper = []
    lower = []

    for x in x_values:
        xc = x / c
        yt = 5.0 * t * c * (
            0.2969 * math.sqrt(xc)
            - 0.1260 * xc
            - 0.3516 * xc**2
            + 0.2843 * xc**3
            - 0.1015 * xc**4
        )
        upper.append((x, yt))
        lower.append((x, -yt))

    return upper + list(reversed(lower))[1:-1]


def rotate_point(x, y, angle_deg, origin_x=0.25, origin_y=0.0):
    angle_rad = math.radians(angle_deg)
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    dx = x - origin_x
    dy = y - origin_y
    xr = origin_x + dx * cos_a - dy * sin_a
    yr = origin_y + dx * sin_a + dy * cos_a
    return xr, yr


def write_legacy_vtk(profile, vtk_path: Path) -> None:
    vtk_path.parent.mkdir(parents=True, exist_ok=True)
    with vtk_path.open("w", encoding="ascii") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("NACA 0012 profile\n")
        handle.write("ASCII\n")
        handle.write("DATASET POLYDATA\n")
        handle.write("POINTS {} float\n".format(len(profile)))
        for x, y in profile:
            handle.write("{:.9f} {:.9f} 0.0\n".format(x, y))
        handle.write("POLYGONS 1 {}\n".format(len(profile) + 2))
        handle.write(str(len(profile)))
        for idx in range(len(profile)):
            handle.write(" {}".format(idx))
        handle.write("\n")


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def normal(a, b, c):
    ab = sub(b, a)
    ac = sub(c, a)
    n = cross(ab, ac)
    mag = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
    if mag == 0.0:
        return (0.0, 0.0, 0.0)
    return (n[0] / mag, n[1] / mag, n[2] / mag)


def write_binary_stl(triangles, stl_path: Path) -> None:
    stl_path.parent.mkdir(parents=True, exist_ok=True)
    with stl_path.open("wb") as handle:
        handle.write(b"NACA 0012 airfoil".ljust(80, b" "))
        handle.write(struct.pack("<I", len(triangles)))
        for tri in triangles:
            n = normal(tri[0], tri[1], tri[2])
            payload = list(n)
            for point in tri:
                payload.extend(point)
            handle.write(struct.pack("<12fH", *(payload + [0])))


def triangulate_caps(front_ring, back_ring):
    triangles = []
    front_center = (
        sum(point[0] for point in front_ring) / len(front_ring),
        sum(point[1] for point in front_ring) / len(front_ring),
        front_ring[0][2],
    )
    back_center = (
        sum(point[0] for point in back_ring) / len(back_ring),
        sum(point[1] for point in back_ring) / len(back_ring),
        back_ring[0][2],
    )

    for idx in range(len(front_ring)):
        nxt = (idx + 1) % len(front_ring)
        triangles.append((front_center, front_ring[nxt], front_ring[idx]))
        triangles.append((back_center, back_ring[idx], back_ring[nxt]))

    return triangles


def triangulate_sides(front_ring, back_ring):
    triangles = []
    for idx in range(len(front_ring)):
        nxt = (idx + 1) % len(front_ring)
        a0 = front_ring[idx]
        a1 = front_ring[nxt]
        b0 = back_ring[idx]
        b1 = back_ring[nxt]
        triangles.append((a0, a1, b1))
        triangles.append((a0, b1, b0))
    return triangles


def create_naca_0012_airfoil(
    t=0.12,
    c=1.0,
    n_points=200,
    thickness=0.1,
    angle_deg=0.0,
    vtk_filename="naca0012.vtk",
    stl_filename="naca0012.stl",
):
    profile = naca_0012_airfoil(t, c, n_points)
    if angle_deg:
        profile = [rotate_point(x, y, angle_deg) for x, y in profile]
    half_span = 0.5 * thickness

    front_ring = [(x, y, -half_span) for x, y in profile]
    back_ring = [(x, y, half_span) for x, y in profile]

    triangles = triangulate_caps(front_ring, back_ring)
    triangles.extend(triangulate_sides(front_ring, back_ring))
    write_binary_stl(triangles, Path(stl_filename))

    if vtk_filename:
        write_legacy_vtk(profile, Path(vtk_filename))

    print(
        "3D NACA 0012 airfoil saved to {}".format(stl_filename)
        + (" and {}".format(vtk_filename) if vtk_filename else "")
    )


def main():
    parser = argparse.ArgumentParser(
        description="Generate a centered 3D NACA 0012 airfoil as STL and optional VTK."
    )
    parser.add_argument("--t", type=float, default=0.12, help="Airfoil thickness ratio.")
    parser.add_argument("--c", type=float, default=1.0, help="Chord length.")
    parser.add_argument(
        "--n_points",
        type=int,
        default=200,
        help="Number of cosine-spaced points on one surface.",
    )
    parser.add_argument(
        "--thickness",
        type=float,
        default=0.1,
        help="Extrusion thickness in z.",
    )
    parser.add_argument(
        "--angle_deg",
        type=float,
        default=0.0,
        help="Rotate the airfoil geometry by this angle in degrees.",
    )
    parser.add_argument(
        "--vtk_filename",
        type=str,
        default="naca0012.vtk",
        help="Optional legacy VTK filename. Pass an empty string to skip.",
    )
    parser.add_argument(
        "--stl_filename",
        type=str,
        default="naca0012.stl",
        help="Output STL filename.",
    )

    args = parser.parse_args()
    vtk_filename = args.vtk_filename or None

    create_naca_0012_airfoil(
        t=args.t,
        c=args.c,
        n_points=args.n_points,
        thickness=args.thickness,
        angle_deg=args.angle_deg,
        vtk_filename=vtk_filename,
        stl_filename=args.stl_filename,
    )


if __name__ == "__main__":
    main()
