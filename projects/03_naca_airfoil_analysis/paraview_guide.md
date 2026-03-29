# ParaView View For This Case

To get a view closer to the reference image:

1. Open `airfoil.foam` in ParaView.
2. In the mesh parts list, enable only `internalMesh` and `patch/airfoil`.
3. Click `Apply`.
4. Set the view direction to `+Z` so you look normal to the 2D plane.
5. Click `Reset Camera`, then zoom tightly around the airfoil.
6. Set the `internalMesh` representation to `Surface`.
7. Color by `U`, then switch the component selector to `Magnitude`.
8. Rescale to the current data range.
9. Apply `Filters -> Stream Tracer` with a line seed upstream of the airfoil.
10. Color the streamlines by `U Magnitude`.

Suggested line seed:
- Point 1: `(-1.5, -0.8, 0)`
- Point 2: `(-1.5, 0.8, 0)`
- Resolution: `120`

Suggested camera:
- Position: `(0.5, 0.0, 8.0)`
- Focal Point: `(0.5, 0.0, 0.0)`
- View Up: `(0, 1, 0)`

This case is now set to about `8°` nose-up angle of attack by rotating the airfoil geometry while keeping the freestream horizontal.
