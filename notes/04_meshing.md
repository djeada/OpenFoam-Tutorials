# Meshing in OpenFOAM

Meshing is arguably the most critical step in any CFD simulation. The mesh — the
discrete representation of your physical domain — directly controls solution
accuracy, convergence behaviour, and computational cost. A careless mesh can
produce results that look plausible but are quantitatively wrong; a well-crafted
mesh makes everything downstream easier.

```
    ╔══════════════════════════════════════════════════════════════════════╗
    ║                  WHY MESHING MATTERS                               ║
    ╠══════════════════════════════════════════════════════════════════════╣
    ║                                                                    ║
    ║   Physical Domain          Discretised Domain        Numerical     ║
    ║   (continuous)             (mesh / grid)             Solution      ║
    ║                                                                    ║
    ║   ┌─────────────┐         +-----+-----+-----+                     ║
    ║   │  fluid /    │         |     |     |     |       u, p, T …     ║
    ║   │  solid      │  ───►   +-----+-----+-----+  ───►  at each     ║
    ║   │  region     │         |     |     |     |       cell centre   ║
    ║   └─────────────┘         +-----+-----+-----+                     ║
    ║                                                                    ║
    ║   Mesh quality  ──►  Discretisation error  ──►  Solution accuracy  ║
    ╚══════════════════════════════════════════════════════════════════════╝
```

OpenFOAM provides two native meshing utilities — `blockMesh` for structured
hexahedral grids and `snappyHexMesh` for complex-geometry unstructured meshes —
as well as converters for meshes created in external tools.

---

## 1  Mesh Fundamentals

### 1.1  What Is a Mesh?

A mesh (or grid) partitions the computational domain into a finite number of
non-overlapping **cells** (control volumes). OpenFOAM stores values at **cell
centres** and fluxes through **cell faces** — this is the finite-volume method.

```
  Continuous domain               Discretised domain (mesh)

  ┌───────────────────┐           +---+---+---+---+---+
  │                   │           |   |   |   |   |   |
  │                   │    ───►   +---+---+---+---+---+
  │      Ω            │           |   |   |   |   |   |
  │                   │           +---+---+---+---+---+
  │                   │           |   |   |   |   |   |
  └───────────────────┘           +---+---+---+---+---+
                                  Each small box = one cell
```

### 1.2  Cell Types

OpenFOAM is polyhedral-capable, but most meshes use one or more of these types:

```
  Hexahedron              Tetrahedron           Prism               Polyhedron
  (6 faces, 8 vertices)   (4 faces, 4 verts)    (5 faces, 6 verts)  (N faces)

  +--------+                  /\                   /\                 +----+
  /|       /|                / | \                /  \               /  /   \
 / |      / |               /  |  \              /    \             +  |     |
+--------+  |              /   |   \            +------+            |  +-----+
|  +-----|--+             /____|____\           |      |             \  \   /
| /      | /                                    |      |              +----+
|/       |/                                     +------+
+--------+                                                          Arbitrary
                                                                    N-gon faces
```

| Property        | Hexahedra         | Tetrahedra        | Prisms           | Polyhedra         |
|-----------------|-------------------|-------------------|------------------|-------------------|
| Accuracy        | Highest           | Lowest            | Good             | Good              |
| Cells needed    | Fewest            | Most (~5-10×)     | Moderate         | Moderate          |
| Automation      | Hard (structured) | Easy (auto-mesh)  | Used for layers  | snappyHexMesh     |
| Typical use     | Simple geometry   | Complex geometry  | Boundary layers  | Adaptive meshes   |

### 1.3  Structured vs Unstructured vs Hybrid

| Feature            | Structured              | Unstructured             | Hybrid                 |
|--------------------|-------------------------|--------------------------|------------------------|
| Cell arrangement   | Regular i,j,k indexing  | Arbitrary connectivity   | Mix of both            |
| Cell type          | Hexahedra               | Tet / poly / mixed       | Hex core + tet/poly    |
| Geometry fit       | Simple shapes only      | Any complex shape        | Complex with quality   |
| OpenFOAM tool      | `blockMesh`             | `snappyHexMesh`          | `snappyHexMesh`        |
| Memory efficiency  | Best                    | Moderate                 | Good                   |

### 1.4  Cell Quality Metrics

Three metrics dominate mesh quality assessment in OpenFOAM:

**Non-orthogonality** — the angle between the face-normal vector and the vector
connecting the two cell centres that share that face.

```
  GOOD  (orthogonal, ~0°)            BAD  (non-orthogonal, > 70°)

  +-------+-------+                  +-------+
  |       |       |                  |      / \
  |   C1 -|-> C2  |                  |  C1/    \ C2
  |       |       |                  |   /  θ   \
  +-------+-------+                  +--+--------+

  face normal n  ∥  d = C2 − C1      face normal n  ∦  d = C2 − C1
  θ ≈ 0°   ✓ Excellent               θ > 70°  ✗ Poor accuracy
```

**Skewness** — how far the intersection of vector d with the face is from the
face centre. High skewness degrades gradient calculations.

**Aspect ratio** — ratio of longest to shortest cell dimension. High aspect
ratios slow convergence and can cause instability.

```
  Aspect ratio ≈ 1  (ideal)          Aspect ratio >> 1  (stretched)

  +------+                           +-------------------------------+
  |      |                           |                               |
  |      |                           +-------------------------------+
  +------+
  ~1:1 ratio                          ~10:1 ratio — avoid unless
                                      aligned with flow direction
```

> **Tip:** Run `checkMesh` after every mesh generation step. Fix quality issues
> before attempting a simulation — bad meshes waste far more time than good
> meshes take to create.

---

## 2  blockMesh — Structured Meshing

### 2.1  Overview

`blockMesh` reads a dictionary called `blockMeshDict` (in the `system/`
directory) and generates a **structured hexahedral mesh**. It is the go-to tool
when:

- The geometry is simple (boxes, channels, pipes, wedges).
- You need precise control over cell distribution and grading.
- You want the highest numerical accuracy for a given cell count.

### 2.2  blockMeshDict — Field by Field

```
  blockMeshDict Structure
  ════════════════════════

  ┌──────────────────────────┐
  │  FoamFile header         │   version, format, class, object
  ├──────────────────────────┤
  │  convertToMeters  0.1;   │   scale factor applied to all vertex coords
  ├──────────────────────────┤
  │  vertices ( ... );       │   list of (x y z) coordinates
  ├──────────────────────────┤
  │  blocks ( ... );         │   hex blocks referencing vertex labels
  ├──────────────────────────┤
  │  edges ( ... );          │   optional curved-edge definitions
  ├──────────────────────────┤
  │  boundary ( ... );       │   patch definitions with face lists
  ├──────────────────────────┤
  │  mergePatchPairs ( ... );│   optional: stitch blocks together
  └──────────────────────────┘
```

#### convertToMeters

A scalar multiplier applied to every vertex coordinate. If your coordinates are
in millimetres, set `convertToMeters 0.001;`.

#### vertices

An ordered list of 3D points. **The index in the list IS the vertex label** —
vertex 0 is the first entry, vertex 1 is the second, and so on.

```
  Vertex numbering for a single hex block:

        7 ─────────── 6             Vertices 0-3: bottom face (z = zMin)
       /|            /|             Vertices 4-7: top face    (z = zMax)
      / |           / |
     4 ─────────── 5  |             Right-hand rule:
     |  |          |  |             Bottom face: 0 → 1 → 2 → 3  (CCW from below)
     |  3 ─────────|── 2            Top face:    4 → 5 → 6 → 7  (CCW from below)
     | /           | /
     |/            |/               Each top vertex is directly above
     0 ─────────── 1                its corresponding bottom vertex:
                                    0↔4, 1↔5, 2↔6, 3↔7
```

#### blocks

Each block is defined as:

```
hex (v0 v1 v2 v3 v4 v5 v6 v7) (nx ny nz) grading
```

- `(v0 … v7)` — vertex labels following the numbering convention above
- `(nx ny nz)` — number of cells in each local direction
- `grading` — cell-size distribution (see § 2.4)

#### edges

Optional curved edge definitions between vertex pairs, e.g.
`arc 1 5 (0.5 0.1 0.5)` creates a circular arc through the midpoint.

#### boundary

Patch definitions. Each face is specified by four vertex labels ordered so the
face normal points **out of the block** (right-hand rule).

### 2.3  Real Example — Lid-Driven Cavity

From `projects/01_lid_driven_cavity/system/blockMeshDict`:

```c
convertToMeters 0.1;

vertices
(
    (0 0 0)       // 0
    (1 0 0)       // 1
    (1 1 0)       // 2
    (0 1 0)       // 3
    (0 0 0.1)     // 4
    (1 0 0.1)     // 5
    (1 1 0.1)     // 6
    (0 1 0.1)     // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (20 20 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    movingWall
    {
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }
    fixedWalls
    {
        type wall;
        faces
        (
            (0 4 7 3)
            (2 6 5 1)
            (1 5 4 0)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }
);
```

```
  Lid-Driven Cavity — mesh layout (2D view, z = const)

         movingWall (lid, slides →)
       3 ═══════════════════════ 2
       ║  +--+--+--+--+--+--+  ║
       ║  |  |  |  |  |  |  |  ║
  fixed║  +--+--+--+--+--+--+  ║fixed
  Walls║  |  |  |  |  |  |  |  ║Walls
       ║  +--+--+--+--+--+--+  ║
       ║  |  |  |  |  |  |  |  ║
       0 ═══════════════════════ 1
              fixedWalls

  20×20 cells, single cell deep (2D via "empty" patches)
  convertToMeters 0.1  →  physical domain = 0.1 m × 0.1 m
```

> **Note:** The `frontAndBack` patches are `type empty` — this tells OpenFOAM
> this is a 2D simulation. Only one cell is needed in the z-direction.

### 2.4  Real Example — Airfoil Background Mesh

From `projects/04_naca_airfoil_analysis/system/blockMeshDict` — used as the
background mesh for `snappyHexMesh`:

```c
convertToMeters 1;

vertices
(
    (-0.5 -0.5 -0.1)  // 0
    (2.0 -0.5 -0.1)   // 1
    (2.0  0.5 -0.1)   // 2
    (-0.5  0.5 -0.1)  // 3
    (-0.5 -0.5  0.1)  // 4
    (2.0 -0.5  0.1)   // 5
    (2.0  0.5  0.1)   // 6
    (-0.5  0.5  0.1)  // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (60 20 1) simpleGrading (1 1 1)
);

boundary
(
    inlet
    {
        type patch;
        faces ( (0 4 7 3) );
    }
    outlet
    {
        type patch;
        faces ( (1 5 6 2) );
    }
    walls
    {
        type wall;
        faces
        (
            (3 2 6 7)
            (0 1 5 4)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            (0 3 2 1)
            (4 7 6 5)
        );
    }
);
```

```
  Airfoil background mesh domain  (2.5 m × 1.0 m)

  y = 0.5   3 ─────────────────────────────── 2
             │  inlet          outlet          │
             │  ←              →               │
             │       ┌──airfoil──┐             │
  y = 0      │       └───────────┘             │
             │                                 │
  y = -0.5  0 ─────────────────────────────── 1
           x=-0.5                           x=2.0

  60 × 20 cells, uniform grading — snappyHexMesh refines near the airfoil
```

### 2.5  Grading (Cell Size Distribution)

Grading controls how cell sizes vary across a block. `simpleGrading (1 1 1)`
means uniform cells. Non-unity values create expansion ratios.

```
  simpleGrading (1 1 1)            simpleGrading (4 1 1)
  (uniform spacing)                (cells expand 4× in x-direction)

  +---+---+---+---+---+           +--+--+---+-----+--------+
  |   |   |   |   |   |           |  |  |   |     |        |
  +---+---+---+---+---+           +--+--+---+-----+--------+
  Equal Δx everywhere              Δx grows left → right

  Multi-grading example:
  simpleGrading ( ((0.5 0.5 4)(0.5 0.5 0.25)) 1 1 )

  +--+--+---+----+----+---+--+--+
  |  |  |   |    |    |   |  |  |
  +--+--+---+----+----+---+--+--+
  ◄─── expand ───►◄── compress ──►
  Fine at centre, coarse at edges
```

> **Tip:** Use grading to cluster cells near walls, inlets, or regions of high
> gradients. This gives better resolution where it matters without increasing
> the total cell count.

### 2.6  Multi-Block Meshes

Complex geometries can be built from **multiple blocks** that share vertices:

```
  Multi-block example: L-shaped domain

       +----------+----------+
       |          |          |
       | Block 1  | Block 2  |
       |          |          |
       +----------+----------+
       |          |
       | Block 3  |
       |          |
       +----------+

  Blocks share vertices along common edges.
  Use mergePatchPairs to stitch non-conforming blocks.
```

---

## 3  snappyHexMesh — Complex Geometry Meshing

### 3.1  Overview

`snappyHexMesh` creates unstructured hex-dominant meshes around complex
geometries defined by **STL surface files**. It is the primary tool for
real-world industrial cases where `blockMesh` alone cannot represent the
geometry.

### 3.2  The Three-Step Process

```
  ╔═══════════════════════════════════════════════════════════════════════════╗
  ║               snappyHexMesh  —  Three-Step Mesh Generation              ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║                                                                         ║
  ║  Step 0             Step 1               Step 2              Step 3     ║
  ║  BACKGROUND         CASTELLATED          SNAPPED             LAYER      ║
  ║  MESH               MESH                 MESH                ADDITION   ║
  ║                                                                         ║
  ║  +----------+       +--+--+--+--+        +--+--+-..          +--+--+-.. ║
  ║  |          |       |  |  |xx|  |        |  |  / ·\          |  | /===\ ║
  ║  |          | ───►  |--+--+--+--| ───►   |--+-·    · ───►   |--+·||| | ║
  ║  |          |       |  |  |xx|  |        |  |  \ ·/          |  | \===/ ║
  ║  +----------+       +--+--+--+--+        +--+--+-··          +--+--+-.. ║
  ║                                                                         ║
  ║  blockMesh           Cells inside         Surface mesh        Boundary  ║
  ║  creates this        geometry are         vertices snap       layer     ║
  ║                      removed (xx)         to STL surface      cells     ║
  ║                      + refinement                             extruded  ║
  ╚═══════════════════════════════════════════════════════════════════════════╝
```

**Step 1 — Castellated mesh:** Refines the background mesh near surfaces and
removes cells that fall inside the geometry (or outside, depending on
`locationInMesh`).

**Step 2 — Snapping:** Iteratively moves surface vertices onto the STL geometry
to create a body-fitted mesh.

**Step 3 — Layer addition:** Extrudes prismatic boundary-layer cells from
specified surfaces. (Optional — controlled by `addLayers` flag.)

### 3.3  snappyHexMeshDict — Structure

The dictionary has five major sections:

```
  snappyHexMeshDict Layout
  ═════════════════════════

  ┌──────────────────────────────────┐
  │  castellatedMesh  true/false;    │   Enable/disable step 1
  │  snap              true/false;   │   Enable/disable step 2
  │  addLayers         true/false;   │   Enable/disable step 3
  ├──────────────────────────────────┤
  │  geometry { ... }                │   STL files and other surfaces
  ├──────────────────────────────────┤
  │  castellatedMeshControls { ... } │   Refinement levels, cell limits
  ├──────────────────────────────────┤
  │  snapControls { ... }            │   Snapping iterations, tolerances
  ├──────────────────────────────────┤
  │  addLayersControls { ... }       │   Layer thickness, expansion ratio
  ├──────────────────────────────────┤
  │  meshQualityControls { ... }     │   Quality thresholds
  └──────────────────────────────────┘
```

### 3.4  Real Example — NACA Airfoil snappyHexMeshDict

From `projects/04_naca_airfoil_analysis/system/snappyHexMeshDict`:

```c
castellatedMesh true;
snap            true;
addLayers       false;

geometry
{
    airfoil.stl                     // STL file in constant/triSurface/
    {
        type triSurfaceMesh;
        name airfoil;
    }
};
```

> **Note:** The STL file lives at
> `projects/04_naca_airfoil_analysis/constant/triSurface/airfoil.stl`.
> snappyHexMesh expects surface files in the `constant/triSurface/` directory.

#### castellatedMeshControls

```c
castellatedMeshControls
{
    maxLocalCells       1000000;    // max cells per processor
    maxGlobalCells      2000000;    // max total cells
    minRefinementCells  0;          // min cells to trigger refinement
    maxLoadUnbalance    0.10;       // load balance tolerance
    nCellsBetweenLevels 3;         // buffer cells between refinement levels

    features ( );                   // feature edge refinement (none here)

    refinementSurfaces
    {
        airfoil
        {
            level (3 3);            // min and max refinement level
        }                           // level 3 = each cell split 2³ = 8×
    }

    resolveFeatureAngle 30;         // refine edges sharper than 30°

    refinementRegions { }           // volume refinement regions (none here)

    locationInMesh (0 0 0);         // point that must be INSIDE the mesh
                                    // (i.e., outside the airfoil body)
    allowFreeStandingZoneFaces true;
}
```

```
  Refinement level illustration:

  Level 0 (base)     Level 1 (2×)       Level 2 (4×)       Level 3 (8×)

  +--------+         +----+----+         +--+--+--+--+      +-+-+-+-+-+-+-+-+
  |        |         |    |    |         |  |  |  |  |      |||||||||||||||||
  |        |         |    |    |         +--+--+--+--+      +-+-+-+-+-+-+-+-+
  |        |         +----+----+         |  |  |  |  |      |||||||||||||||||
  +--------+         |    |    |         +--+--+--+--+      +-+-+-+-+-+-+-+-+
                     |    |    |         |  |  |  |  |      |||||||||||||||||
                     +----+----+         +--+--+--+--+      +-+-+-+-+-+-+-+-+
  1 cell              4 cells             16 cells           64 cells
```

#### snapControls

```c
snapControls
{
    nSmoothPatch 3;         // smoothing iterations on surface patch
    tolerance    2.0;       // distance tolerance (relative to cell size)
    nSolveIter   30;        // mesh motion solver iterations
    nRelaxIter   5;         // relaxation iterations
}
```

#### addLayersControls

```c
addLayersControls
{
    relativeSizes    true;
    layers           { }    // no layers added in this case
    expansionRatio   1.0;
    finalLayerThickness 0.3;
    minThickness     0.25;
    nGrow            0;
    featureAngle     30;
    nRelaxIter       5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals   3;
    nSmoothThickness 10;
    maxFaceThicknessRatio   0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 130;
    nBufferCellsNoExtrude 0;
    nLayerIter       50;
}
```

> **Warning:** Layer addition is set to `false` in this airfoil case. For
> accurate boundary-layer resolution (important for drag prediction), you
> would enable layers and configure thickness parameters carefully.

#### meshQualityControls

```c
meshQualityControls
{
    maxNonOrtho          65;
    maxBoundarySkewness  20;
    maxInternalSkewness  4;
    maxConcave           80;
    minFlatness          0.5;
    minVol               1e-13;
    minTetQuality        1e-9;
    minArea              -1;
    minTwist             0.05;
    minDeterminant       0.001;
    minFaceWeight        0.05;
    minVolRatio          0.01;
    minTriangleTwist     -1;
    nSmoothScale         4;
    errorReduction       0.75;
}
```

### 3.5  STL Geometry Handling

snappyHexMesh reads surface geometry from **STL files** (ASCII or binary)
placed in `constant/triSurface/`. The NACA airfoil project uses a binary STL
file (`airfoil.stl`, ~40 KB) exported from a CAD tool.

> **Tip:** Ensure your STL surface is **watertight** (no gaps or overlaps).
> Non-watertight surfaces cause snappyHexMesh to fail at the castellated mesh
> step because it cannot determine inside vs outside.

---

## 4  Mesh from External Sources

Not every mesh needs to be generated inside OpenFOAM. External tools often have
superior meshing capabilities for complex geometries.

### 4.1  Fluent Mesh Import

The `02_elbow` project imports a mesh created in ANSYS Fluent using the
`fluentMeshToFoam` converter.

From `projects/02_elbow/allrun`:

```bash
#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Get application directory
application=$(getApplication)

runApplication fluentMeshToFoam elbow.msh
runApplication "$application"
runApplication foamMeshToFluent
runApplication foamDataToFluent
```

The key line is `runApplication fluentMeshToFoam elbow.msh` — this reads the
Fluent `.msh` file and creates an OpenFOAM mesh in `constant/polyMesh/`.

```
  External mesh workflow:

  ┌──────────────┐     ┌────────────────────┐     ┌──────────────────┐
  │ External tool │     │  Converter utility  │     │  OpenFOAM case   │
  │              │     │                    │     │                  │
  │ ANSYS Fluent │────►│ fluentMeshToFoam   │────►│ constant/        │
  │ elbow.msh    │     │                    │     │   polyMesh/      │
  └──────────────┘     └────────────────────┘     │     points       │
                                                   │     faces        │
                                                   │     owner        │
                                                   │     neighbour    │
                                                   │     boundary     │
                                                   └──────────────────┘
```

### 4.2  Other Converters

| Converter              | Source Format           | Usage                        |
|------------------------|-------------------------|------------------------------|
| `fluentMeshToFoam`     | Fluent `.msh`           | `fluentMeshToFoam file.msh`  |
| `gmshToFoam`           | Gmsh `.msh`             | `gmshToFoam file.msh`        |
| `ideasUnvToFoam`       | I-DEAS `.unv`           | `ideasUnvToFoam file.unv`    |
| `starToFoam`           | STAR-CD `.ccm`          | `starToFoam file`            |
| `plot3dToFoam`         | Plot3D multi-block      | `plot3dToFoam file.xyz`      |
| `netgenNeutralToFoam`  | Netgen `.mesh`          | `netgenNeutralToFoam f.mesh` |

> **Tip:** After importing any external mesh, always run `checkMesh` to verify
> that the conversion produced a valid OpenFOAM mesh. Boundary names and types
> may need manual correction in `constant/polyMesh/boundary`.

---

## 5  Mesh Quality Assessment

### 5.1  The `checkMesh` Utility

Run `checkMesh` in your case directory to get a comprehensive quality report.
It checks geometry, topology, and quality metrics.

Typical output sections:
- **Mesh stats:** cell count, face count, point count
- **Overall domain bounding box**
- **Cell types:** hex, tet, prism, etc.
- **Topology checks:** boundary openness, face ordering
- **Geometry checks:** non-orthogonality, skewness, aspect ratio

### 5.2  Key Quality Thresholds

| Metric              | Ideal       | Acceptable    | Problematic   | What It Means                          |
|---------------------|-------------|---------------|---------------|----------------------------------------|
| Non-orthogonality   | < 40°       | < 70°         | > 80°         | Angle between face normal and d-vector |
| Max skewness        | < 2         | < 4           | > 6           | Face centre offset from d intersection |
| Max aspect ratio    | < 10        | < 100         | > 1000        | Longest-to-shortest cell dimension     |
| Min volume          | > 0         | > 0           | ≤ 0           | Negative volume = inverted cell        |
| Min determinant     | > 0.1       | > 0.001       | ≤ 0           | Cell shape regularity measure          |

### 5.3  Visualising Non-Orthogonality

```
  GOOD (orthogonal):                    BAD (non-orthogonal):

  +--------+--------+                   +--------+
  |        |        |                   |       / \
  |   C1 ──┼──► C2  |                   |  C1 /    \ C2
  |        |        |                   |    / θ    \
  +--------+--------+                   +--+/        \+--
                                             face
  d = C2 − C1                           d = C2 − C1
  n = face normal                       n = face normal
  n ∥ d  →  θ = 0°  ✓                  n ∦ d  →  θ >> 0°  ✗

  The finite-volume method assumes n ≈ d for gradient accuracy.
  Large non-orthogonality requires correctors (increasing cost).
```

### 5.4  Quality Controls in snappyHexMeshDict

The `meshQualityControls` section in the NACA airfoil case (§ 3.4) enforces
quality during mesh generation. Key settings from that file:

```
maxNonOrtho          65      // reject cells with non-orth > 65°
maxBoundarySkewness  20      // boundary face skewness limit
maxInternalSkewness  4       // internal face skewness limit
minVol               1e-13   // reject negative or zero-volume cells
minDeterminant       0.001   // reject degenerate cells
```

> **Warning:** If snappyHexMesh reports many quality failures, don't just raise
> the thresholds. Instead, improve the background mesh resolution, adjust the
> STL geometry, or tune `nSmoothPatch` / `nRelaxIter` in `snapControls`.

---

## 6  Mesh Refinement Strategies

### 6.1  Uniform Refinement

The simplest approach: increase cell count in all directions equally.

```
  20 × 20 mesh                    40 × 40 mesh

  +--+--+--+--+                   +-+-+-+-+-+-+-+-+
  |  |  |  |  |                   | | | | | | | | |
  +--+--+--+--+                   +-+-+-+-+-+-+-+-+
  |  |  |  |  |                   | | | | | | | | |
  +--+--+--+--+                   +-+-+-+-+-+-+-+-+
  |  |  |  |  |                   | | | | | | | | |
  +--+--+--+--+                   +-+-+-+-+-+-+-+-+
                                  | | | | | | | | |
  400 cells                       +-+-+-+-+-+-+-+-+
                                  1600 cells (4× cost)
```

Simple, but expensive. Doubling resolution in 3D increases cell count by **8×**.

### 6.2  Adaptive / Local Refinement

Refine only where needed — near walls, in wake regions, around features. This
is exactly what snappyHexMesh's `refinementSurfaces` and `refinementRegions`
achieve.

```
  Local refinement around an airfoil:

  +-----+-----+-----+-----+-----+-----+-----+
  |     |     |     |     |     |     |     |
  +-----+--+--+--+--+--+--+--+--+-----+-----+
  |     |  |  |  |  |  |  |  |  |     |     |
  +-----+--+-++-++-++-++-++--+--+-----+-----+
  |     |  |~~~~~~~~~airfoil|  |  |     |     |
  +-----+--+-++-++-++-++-++--+--+-----+-----+
  |     |  |  |  |  |  |  |  |  |     |     |
  +-----+--+--+--+--+--+--+--+--+-----+-----+
  |     |     |     |     |     |     |     |
  +-----+-----+-----+-----+-----+-----+-----+

  Coarse far-field  →  Fine near surface
```

### 6.3  Boundary Layer Meshing

Resolving the boundary layer requires thin, high-aspect-ratio cells near walls.
The first cell height is governed by the **y+** value:

```
  y+ ≈ (y · uτ) / ν

  Wall-resolved (y+ ≈ 1):      Wall-function (y+ ≈ 30-300):

  ════════ wall ════════        ════════ wall ════════
  |  | thin layers (y+≈1)      |       | coarser first cell
  |  |                          |       |
  |   |                         |        |
  |    |  expanding             |         | fewer layers needed
  |     |                       |          |
  |      | to freestream        |           |
```

The `addLayersControls` section in snappyHexMeshDict controls boundary layer
mesh generation. Key parameters:

- `nSurfaceLayers` — number of layers per surface patch
- `expansionRatio` — growth ratio between successive layers
- `finalLayerThickness` — thickness of outermost layer (relative or absolute)

### 6.4  Grid Independence Study

A simulation result is only trustworthy if it does not change significantly when
the mesh is refined. A grid independence study uses 3+ meshes:

```
  Grid Independence Study Workflow:

  ┌────────────────┐    ┌────────────────┐    ┌────────────────┐
  │  Coarse mesh   │    │  Medium mesh   │    │  Fine mesh     │
  │  (N cells)     │    │  (~2× N)       │    │  (~4× N)       │
  └───────┬────────┘    └───────┬────────┘    └───────┬────────┘
          │                     │                     │
          ▼                     ▼                     ▼
       Solve                 Solve                 Solve
          │                     │                     │
          ▼                     ▼                     ▼
       Result₁               Result₂              Result₃
          │                     │                     │
          └─────────┬───────────┘─────────────────────┘
                    ▼
          Compare key quantity (drag, pressure drop, etc.)
          If Result₂ ≈ Result₃ → medium mesh is sufficient
```

> **Tip:** Always report your grid independence study in publications. A single
> mesh result without convergence evidence is not credible.

---

## 7  Practical Meshing Workflow

```
  ╔═══════════════════════════════════════════════════════════════════════════╗
  ║                 COMPLETE MESHING WORKFLOW                                ║
  ╠═══════════════════════════════════════════════════════════════════════════╣
  ║                                                                         ║
  ║  1. Prepare Geometry                                                    ║
  ║     │  • Clean CAD model, export STL (if using snappyHexMesh)           ║
  ║     │  • Ensure watertight surfaces, no self-intersections              ║
  ║     ▼                                                                   ║
  ║  2. Create Background Mesh                                              ║
  ║     │  • Write blockMeshDict (define domain extent, base resolution)    ║
  ║     │  • Run: blockMesh                                                 ║
  ║     │  • Run: checkMesh                                                 ║
  ║     ▼                                                                   ║
  ║  3. Refine (if needed)                                                  ║
  ║     │  • Write snappyHexMeshDict (refinement levels, snapping, layers)  ║
  ║     │  • Run: snappyHexMesh                                             ║
  ║     │  • Run: checkMesh                                                 ║
  ║     ▼                                                                   ║
  ║  4. Assess Quality                                                      ║
  ║     │  • Check non-orthogonality, skewness, aspect ratio                ║
  ║     │  • Visualise in ParaView — look for distorted cells               ║
  ║     │  • If quality is poor → go back to step 2 or 3 and adjust        ║
  ║     ▼                                                                   ║
  ║  5. Run Simulation                                                      ║
  ║     │  • If solver diverges → suspect mesh quality first                ║
  ║     ▼                                                                   ║
  ║  6. Grid Independence                                                   ║
  ║        • Repeat with 2-3 mesh densities                                 ║
  ║        • Confirm results are mesh-independent                           ║
  ║                                                                         ║
  ╚═══════════════════════════════════════════════════════════════════════════╝
```

---

## 8  Comparison: blockMesh vs snappyHexMesh

| Feature                  | blockMesh                         | snappyHexMesh                        |
|--------------------------|-----------------------------------|--------------------------------------|
| Mesh type                | Structured hexahedral             | Unstructured hex-dominant            |
| Geometry complexity      | Simple (blocks, wedges)           | Arbitrary (STL surfaces)             |
| Input                    | `blockMeshDict`                   | `blockMeshDict` + `snappyHexMeshDict`|
| Surface geometry         | Not needed                        | STL files required                   |
| Boundary layers          | Via grading                       | Automatic layer addition             |
| Local refinement         | Limited (multi-block grading)     | Flexible (surface + volume regions)  |
| Cell quality             | Excellent (all hex)               | Good (mostly hex + some poly)        |
| Setup complexity         | Low                               | Moderate–High                        |
| Run time                 | Seconds                           | Minutes–Hours                        |
| Typical projects         | Cavity, channel, pipe, wedge      | Airfoil, vehicle, building, turbine  |
| This repo                | `01_lid_driven_cavity`            | `04_naca_airfoil_analysis`           |
| External mesh equivalent | —                                 | `02_elbow` (imported from Fluent)    |

---

## 9  Common Meshing Pitfalls & Tips

> **Pitfall:** Running a simulation on a mesh you haven't checked.
> **Fix:** Always run `checkMesh` before solving. Always.

> **Pitfall:** Non-watertight STL surfaces causing snappyHexMesh to crash.
> **Fix:** Inspect STL in ParaView or MeshLab. Look for open edges,
> self-intersections, and duplicate triangles.

> **Pitfall:** `locationInMesh` point falls inside the geometry body.
> **Fix:** The point must be in the **fluid region** (outside the solid body).
> Check coordinates carefully against your STL bounding box.

> **Pitfall:** Too-coarse background mesh for snappyHexMesh.
> **Fix:** The background mesh should have at least a few cells across the
> smallest geometric feature. Increase `blockMeshDict` resolution.

> **Pitfall:** Ignoring mesh independence — reporting results from a single mesh.
> **Fix:** Run at least 3 mesh densities. Compare a key output quantity.

> **Pitfall:** Extremely high aspect ratio cells in boundary layers.
> **Fix:** Use gradual expansion ratios (1.1–1.3). Sudden jumps in cell size
> degrade solution accuracy and can cause convergence issues.

> **Pitfall:** 2D simulations with more than one cell in the z-direction.
> **Fix:** For 2D cases, use exactly 1 cell in z and `type empty` patches on
> front and back faces (see the cavity example in § 2.3).

> **Pitfall:** Boundary patches not matching between mesh and field files.
> **Fix:** Patch names in `blockMeshDict` / `snappyHexMeshDict` must exactly
> match the patch names used in `0/U`, `0/p`, and other field files.

---

## 10  Quick Reference

```
  ┌──────────────────────────────────────────────────────────────┐
  │                   MESHING COMMANDS CHEATSHEET                │
  ├──────────────────────────────────────────────────────────────┤
  │                                                              │
  │  blockMesh                Generate structured mesh           │
  │  snappyHexMesh            Generate complex-geometry mesh     │
  │  checkMesh                Check mesh quality                 │
  │  checkMesh -allGeometry   Extended geometry checks           │
  │  checkMesh -allTopology   Extended topology checks           │
  │  fluentMeshToFoam f.msh   Import Fluent mesh                │
  │  gmshToFoam f.msh         Import Gmsh mesh                  │
  │  transformPoints -scale '(0.001 0.001 0.001)'  Scale mesh   │
  │  renumberMesh             Optimise cell ordering             │
  │  refineMesh               Refine existing mesh               │
  │                                                              │
  └──────────────────────────────────────────────────────────────┘
```
