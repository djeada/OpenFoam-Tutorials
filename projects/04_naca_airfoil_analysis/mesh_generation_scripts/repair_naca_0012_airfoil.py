import argparse
import sys

import bmesh
import bpy


def blender_args():
    return sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []


def import_stl(filepath: str) -> None:
    if hasattr(bpy.ops.wm, "stl_import"):
        bpy.ops.wm.stl_import(filepath=filepath)
    else:
        bpy.ops.import_mesh.stl(filepath=filepath)


def export_stl(filepath: str) -> None:
    if hasattr(bpy.ops.wm, "stl_export"):
        bpy.ops.wm.stl_export(filepath=filepath, export_selected_objects=True)
    else:
        bpy.ops.export_mesh.stl(filepath=filepath, use_selection=True)


def repair_stl(input_stl: str, repaired_stl: str) -> None:
    bpy.ops.wm.read_factory_settings(use_empty=True)
    import_stl(input_stl)

    obj = bpy.context.selected_objects[0]
    bpy.context.view_layer.objects.active = obj
    mesh = obj.data

    bpy.ops.object.mode_set(mode="EDIT")
    bm = bmesh.from_edit_mesh(mesh)

    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=1e-5)
    loose = [vert for vert in bm.verts if not vert.link_faces]
    if loose:
        bmesh.ops.delete(bm, geom=loose, context="VERTS")
    bmesh.ops.holes_fill(bm, edges=bm.edges[:], sides=0)
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])

    bmesh.update_edit_mesh(mesh)
    bpy.ops.object.mode_set(mode="OBJECT")
    obj.select_set(True)

    export_stl(repaired_stl)
    print(f"Repaired STL file saved to {repaired_stl}")


def main():
    parser = argparse.ArgumentParser(description="Repair an STL file using Blender.")
    parser.add_argument("--input_stl", default="naca0012.stl", help="Input STL path.")
    parser.add_argument(
        "--repaired_stl",
        default="naca0012_repaired.stl",
        help="Output repaired STL path.",
    )
    args = parser.parse_args(blender_args())
    repair_stl(args.input_stl, args.repaired_stl)


if __name__ == "__main__":
    main()
