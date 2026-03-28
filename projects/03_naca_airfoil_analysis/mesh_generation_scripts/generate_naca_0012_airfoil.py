import numpy as np
import pyvista as pv
import argparse

# Function to generate NACA 0012 coordinates
def naca_0012_airfoil(t=0.12, c=1.0, n_points=100):
    x = np.linspace(0, c, n_points)
    yt = (t / 0.2) * c * (0.2969 * np.sqrt(x / c) - 0.1260 * (x / c) - 0.3516 * (x / c)**2 + 0.2843 * (x / c)**3 - 0.1036 * (x / c)**4)

    xu = x
    yu = yt
    xl = x
    yl = -yt

    x_coords = np.concatenate((xu, xl[::-1], [xu[0]]))
    y_coords = np.concatenate((yu, yl[::-1], [yu[0]]))
    return x_coords, y_coords

# Function to create 3D NACA 0012 airfoil and save it to files
def create_naca_0012_airfoil(t=0.12, c=1.0, n_points=100, thickness=0.1, vtk_filename="naca0012.vtk", stl_filename="naca0012.stl"):
    # Generate the airfoil coordinates
    x_coords, y_coords = naca_0012_airfoil(t, c, n_points)

    # Create the 2D profile as a continuous polyline in pyvista
    points = np.c_[x_coords, y_coords, np.zeros_like(x_coords)]
    n_points = points.shape[0]

    # Create the line array correctly
    lines = np.zeros((n_points, 3), dtype=int)
    lines[:, 0] = 2  # Number of points per line segment
    lines[:-1, 1] = np.arange(n_points - 1)
    lines[:-1, 2] = np.arange(1, n_points)
    lines[-1, 1] = n_points - 1
    lines[-1, 2] = 0  # Closing the loop

    profile_2d = pv.PolyData()
    profile_2d.points = points
    profile_2d.lines = lines.flatten()

    # Extrude the 2D profile to 3D
    profile_3d = profile_2d.extrude([0, 0, thickness], capping=True)

    # Save the 3D profile to VTK and STL files
    profile_3d.save(vtk_filename)
    profile_3d.save(stl_filename)

    print(f"3D NACA 0012 airfoil saved to {vtk_filename} and {stl_filename}")

def main():
    parser = argparse.ArgumentParser(description="Generate a 3D NACA 0012 airfoil profile and save to VTK and STL files.")
    parser.add_argument("--t", type=float, default=0.12, help="Thickness of the airfoil (default: 0.12)")
    parser.add_argument("--c", type=float, default=1.0, help="Chord length of the airfoil (default: 1.0)")
    parser.add_argument("--n_points", type=int, default=100, help="Number of points for the airfoil (default: 100)")
    parser.add_argument("--thickness", type=float, default=0.1, help="Extrusion thickness (default: 0.1)")
    parser.add_argument("--vtk_filename", type=str, default="naca0012.vtk", help="Output VTK filename (default: naca0012.vtk)")
    parser.add_argument("--stl_filename", type=str, default="naca0012.stl", help="Output STL filename (default: naca0012.stl)")

    args = parser.parse_args()

    create_naca_0012_airfoil(args.t, args.c, args.n_points, args.thickness, args.vtk_filename, args.stl_filename)

if __name__ == "__main__":
    main()
