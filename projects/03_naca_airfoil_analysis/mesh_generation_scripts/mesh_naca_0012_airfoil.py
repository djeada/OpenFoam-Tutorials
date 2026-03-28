import bpy
import bmesh
import argparse

# Function to mesh STL file
def mesh_stl(input_stl, output_stl):
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

    # Ensure the mesh is manifold
    bmesh.ops.triangulate(bm, faces=bm.faces)
    bmesh.ops.smooth_vert(bm, verts=bm.verts, factor=1.0, use_axis_x=True, use_axis_y=True, use_axis_z=True)

    # Update the mesh and switch back to object mode
    bmesh.update_edit_mesh(mesh)
    bpy.ops.object.mode_set(mode='OBJECT')

    # Export the meshed STL file
    bpy.ops.export_mesh.stl(filepath=output_stl)

    print(f"Meshed STL file saved to {output_stl}")

def main():
    parser = argparse.ArgumentParser(description="Mesh an STL file using Blender.")
    parser.add_argument("--input_stl", type=str, default="naca0012_repaired.stl", help="Input STL filename (default: naca0012_repaired.stl)")
    parser.add_argument("--output_stl", type=str, default="naca0012_meshed.stl", help="Output meshed STL filename (default: naca0012_meshed.stl)")

    args = parser.parse_args()

    mesh_stl(args.input_stl, args.output_stl)

if __name__ == "__main__":
    main()
