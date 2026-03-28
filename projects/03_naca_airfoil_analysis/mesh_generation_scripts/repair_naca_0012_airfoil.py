import bpy
import bmesh
import argparse

# Function to repair STL file
def repair_stl(input_stl, repaired_stl):
    # Clear existing mesh data
    bpy.ops.wm.read_factory_settings(use_empty=True)

    # Import the STL file
    bpy.ops.import_mesh.stl(filepath=input_stl)

    # Get the imported mesh
    obj = bpy.context.selected_objects[0]
    mesh = obj.data

    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT')

    # Create a BMesh from the mesh
    bm = bmesh.from_edit_mesh(mesh)

    # Remove duplicate vertices
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.0001)

    # Remove loose geometry (vertices, edges, faces not part of any faces)
    to_remove = [v for v in bm.verts if not v.link_faces]
    bmesh.ops.delete(bm, geom=to_remove, context='VERTS')

    # Fill holes
    bmesh.ops.holes_fill(bm, edges=bm.edges)

    # Validate the mesh (removes non-manifold elements)
    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    to_remove = [v for v in bm.verts if not v.link_faces]
    bmesh.ops.delete(bm, geom=to_remove, context='VERTS')

    # Update the mesh and switch back to object mode
    bmesh.update_edit_mesh(mesh)
    bpy.ops.object.mode_set(mode='OBJECT')

    # Export the repaired STL file
    bpy.ops.export_mesh.stl(filepath=repaired_stl)

    print(f"Repaired STL file saved to {repaired_stl}")

def main():
    parser = argparse.ArgumentParser(description="Repair an STL file using Blender.")
    parser.add_argument("--input_stl", type=str, default="naca0012.stl", help="Input STL filename (default: naca0012.stl)")
    parser.add_argument("--repaired_stl", type=str, default="naca0012_repaired.stl", help="Output repaired STL filename (default: naca0012_repaired.stl)")

    args = parser.parse_args()

    repair_stl(args.input_stl, args.repaired_stl)

if __name__ == "__main__":
    main()
