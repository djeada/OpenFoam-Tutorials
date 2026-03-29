#!/usr/bin/env python3

import argparse
import math
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a Gmsh .geo file for a wall-resolved NACA 0012 case."
    )
    parser.add_argument(
        "--output",
        default="constant/airfoil_benchmark.geo",
        help="Output .geo file path.",
    )
    parser.add_argument(
        "--chord-points",
        type=int,
        default=81,
        help="Number of cosine-spaced chord points.",
    )
    parser.add_argument(
        "--angle-deg",
        type=float,
        default=10.0,
        help="Angle of attack applied by rotating the airfoil geometry about quarter chord.",
    )
    parser.add_argument(
        "--domain-xmin",
        type=float,
        default=-20.0,
        help="Farfield x-min extent in chord lengths.",
    )
    parser.add_argument(
        "--domain-xmax",
        type=float,
        default=40.0,
        help="Farfield x-max extent in chord lengths.",
    )
    parser.add_argument(
        "--domain-yhalf",
        type=float,
        default=20.0,
        help="Half-height of the farfield domain in chord lengths.",
    )
    parser.add_argument(
        "--span",
        type=float,
        default=0.02,
        help="Spanwise thickness for the pseudo-2D extrusion.",
    )
    parser.add_argument(
        "--mesh-size-airfoil",
        type=float,
        default=0.012,
        help="Base point size assigned to the airfoil spline points.",
    )
    parser.add_argument(
        "--mesh-size-farfield",
        type=float,
        default=1.5,
        help="Base point size assigned to the outer farfield rectangle.",
    )
    parser.add_argument(
        "--first-layer",
        type=float,
        default=2.0e-3,
        help="First boundary-layer thickness in chord lengths.",
    )
    parser.add_argument(
        "--bl-thickness",
        type=float,
        default=0.04,
        help="Total boundary-layer thickness in chord lengths.",
    )
    parser.add_argument(
        "--bl-ratio",
        type=float,
        default=1.15,
        help="Boundary-layer growth ratio.",
    )
    parser.add_argument(
        "--near-size",
        type=float,
        default=0.025,
        help="Target size near the airfoil outside the boundary layer.",
    )
    parser.add_argument(
        "--near-dist",
        type=float,
        default=2.0,
        help="Distance over which the near-body size transitions.",
    )
    parser.add_argument(
        "--wake-size",
        type=float,
        default=0.03,
        help="Target size in the wake refinement box.",
    )
    parser.add_argument(
        "--wake-xmax",
        type=float,
        default=12.0,
        help="Wake box downstream extent in chord lengths.",
    )
    parser.add_argument(
        "--wake-yhalf",
        type=float,
        default=2.5,
        help="Wake box half-height in chord lengths.",
    )
    parser.add_argument(
        "--split-farfield",
        action="store_true",
        help="Emit separate inlet/outlet/top/bottom patches instead of a single farfield patch.",
    )
    parser.add_argument(
        "--no-recombine",
        action="store_true",
        help="Leave the mesh triangular/tetra/prism-based instead of forcing recombination.",
    )
    return parser.parse_args()


def naca0012_thickness(x):
    return 5.0 * 0.12 * (
        0.2969 * math.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x * x
        + 0.2843 * x * x * x
        - 0.1036 * x * x * x * x
    )


def rotate_point(x, y, angle_rad, origin_x=0.25, origin_y=0.0):
    dx = x - origin_x
    dy = y - origin_y
    c = math.cos(angle_rad)
    s = math.sin(angle_rad)
    return (
        origin_x + c * dx - s * dy,
        origin_y + s * dx + c * dy,
    )


def cosine_spacing(count):
    return [0.5 * (1.0 - math.cos(math.pi * i / (count - 1))) for i in range(count)]


def format_point(point_id, x, y, z, lc):
    return f"Point({point_id}) = {{{x:.12g}, {y:.12g}, {z:.12g}, {lc:.12g}}};"


def main():
    args = parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    xs = cosine_spacing(args.chord_points)
    angle_rad = math.radians(args.angle_deg)

    upper_interior = []
    lower_interior = []
    for x in xs[1:-1]:
        yt = naca0012_thickness(x)
        upper_interior.append(rotate_point(x, yt, angle_rad))
        lower_interior.append(rotate_point(x, -yt, angle_rad))

    te = rotate_point(1.0, 0.0, angle_rad)
    le = rotate_point(0.0, 0.0, angle_rad)

    lines = [
        'SetFactory("OpenCASCADE");',
        "",
        "Mesh.Algorithm = 6;",
        "Mesh.RecombinationAlgorithm = 1;",
        f"Mesh.RecombineAll = {0 if args.no_recombine else 1};",
        "Mesh.CharacteristicLengthExtendFromBoundary = 0;",
        "Mesh.CharacteristicLengthFromCurvature = 0;",
        "Mesh.CharacteristicLengthFromPoints = 1;",
        "",
    ]

    point_id = 1
    te_id = point_id
    lines.append(format_point(point_id, te[0], te[1], 0.0, args.mesh_size_airfoil))
    point_id += 1

    upper_ids = [te_id]
    for x, y in reversed(upper_interior):
        upper_ids.append(point_id)
        lines.append(format_point(point_id, x, y, 0.0, args.mesh_size_airfoil))
        point_id += 1

    le_id = point_id
    upper_ids.append(le_id)
    lines.append(format_point(le_id, le[0], le[1], 0.0, args.mesh_size_airfoil))
    point_id += 1

    lower_ids = [le_id]
    for x, y in lower_interior:
        lower_ids.append(point_id)
        lines.append(format_point(point_id, x, y, 0.0, args.mesh_size_airfoil))
        point_id += 1
    lower_ids.append(te_id)

    p1 = point_id
    lines.append(
        format_point(
            p1,
            args.domain_xmin,
            -args.domain_yhalf,
            0.0,
            args.mesh_size_farfield,
        )
    )
    point_id += 1
    p2 = point_id
    lines.append(
        format_point(
            p2,
            args.domain_xmax,
            -args.domain_yhalf,
            0.0,
            args.mesh_size_farfield,
        )
    )
    point_id += 1
    p3 = point_id
    lines.append(
        format_point(
            p3,
            args.domain_xmax,
            args.domain_yhalf,
            0.0,
            args.mesh_size_farfield,
        )
    )
    point_id += 1
    p4 = point_id
    lines.append(
        format_point(
            p4,
            args.domain_xmin,
            args.domain_yhalf,
            0.0,
            args.mesh_size_farfield,
        )
    )
    lines.extend(
        [
            "",
            f"Spline(1) = {{{', '.join(str(i) for i in upper_ids)}}};",
            f"Spline(2) = {{{', '.join(str(i) for i in lower_ids)}}};",
            f"Line(11) = {{{p1}, {p2}}};",
            f"Line(12) = {{{p2}, {p3}}};",
            f"Line(13) = {{{p3}, {p4}}};",
            f"Line(14) = {{{p4}, {p1}}};",
            "",
            "Curve Loop(21) = {11, 12, 13, 14};",
            "Curve Loop(22) = {1, 2};",
            "Plane Surface(31) = {21, 22};",
            "",
            "Field[1] = BoundaryLayer;",
            "Field[1].CurvesList = {1, 2};",
            f"Field[1].Size = {args.first_layer:.12g};",
            f"Field[1].Thickness = {args.bl_thickness:.12g};",
            f"Field[1].Ratio = {args.bl_ratio:.12g};",
            f"Field[1].Quads = {0 if args.no_recombine else 1};",
            "",
            "Field[2] = Distance;",
            "Field[2].CurvesList = {1, 2};",
            "Field[2].Sampling = 200;",
            "",
            "Field[3] = Threshold;",
            "Field[3].InField = 2;",
            f"Field[3].SizeMin = {args.near_size:.12g};",
            f"Field[3].SizeMax = {args.mesh_size_farfield:.12g};",
            f"Field[3].DistMin = {args.bl_thickness:.12g};",
            f"Field[3].DistMax = {args.near_dist:.12g};",
            "",
            "Field[4] = Box;",
            f"Field[4].VIn = {args.wake_size:.12g};",
            f"Field[4].VOut = {args.mesh_size_farfield:.12g};",
            f"Field[4].XMin = 0.0;",
            f"Field[4].XMax = {args.wake_xmax:.12g};",
            f"Field[4].YMin = {-args.wake_yhalf:.12g};",
            f"Field[4].YMax = {args.wake_yhalf:.12g};",
            "Field[4].ZMin = -1;",
            "Field[4].ZMax = 1;",
            "",
            "Field[5] = Min;",
            "Field[5].FieldsList = {3, 4};",
            "Background Field = 5;",
            "BoundaryLayer Field = 1;",
            "",
            f"extruded[] = Extrude {{0, 0, {args.span:.12g}}} {{",
            "  Surface{31};",
            "  Layers{1};",
            "};",
            "",
            'Physical Surface("frontAndBack") = {31, extruded[0]};',
        ]
    )

    if not args.no_recombine:
        lines.insert(lines.index("};"), "  Recombine;")

    if args.split_farfield:
        lines.extend(
            [
                'Physical Surface("bottom") = {extruded[2]};',
                'Physical Surface("outlet") = {extruded[3]};',
                'Physical Surface("top") = {extruded[4]};',
                'Physical Surface("inlet") = {extruded[5]};',
            ]
        )
    else:
        lines.append(
            'Physical Surface("farfield") = {extruded[2], extruded[3], extruded[4], extruded[5]};'
        )

    lines.extend(
        [
            'Physical Surface("airfoil") = {extruded[6], extruded[7]};',
            'Physical Volume("fluid") = {extruded[1]};',
            "",
        ]
    )

    output.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
