# Flow Around a Cylinder

`Difficulty: Intermediate` · `Solver: simpleFoam / pimpleFoam` · `Turbulence: k-ε` · `Mesh: snappyHexMesh`

> **What this project covers:** Setting up a complete external aerodynamics simulation of
> incompressible flow past a circular cylinder using OpenFOAM. This guide walks through the
> physics, geometry creation, mesh generation, boundary conditions, solver configuration,
> post-processing, and validation against published data. Every configuration file is provided
> in full so you can build the case from scratch.

---

## Problem Description

Flow around a circular cylinder is one of the most extensively studied problems in fluid
mechanics. Despite the simple geometry, it produces a remarkably rich set of flow phenomena:
boundary layer separation, recirculation zones, periodic vortex shedding (the Von Kármán
vortex street), and turbulent wake dynamics. This makes it an ideal benchmark for validating
CFD codes and learning simulation techniques.

```
          Free stream U∞ →→→→→→→→→→→→→→→→→
                →→→→→→→→→→→→→→→→→→→→→→→→→
                →→→→→→→ ╭─────╮ →→→→→→→→→
                →→→→→→→ │     │    wake   →
                →→→→→→→ │  ●  │  ~~~~~~~ →
                →→→→→→→ │     │    region →
                →→→→→→→ ╰─────╯ →→→→→→→→→
                →→→→→→→→→→→→→→→→→→→→→→→→→
          →→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→

    Inlet                Cylinder          Outlet
```

**Key flow features:**

- **Stagnation point** — The front of the cylinder where the flow velocity drops to zero and
  pressure reaches its maximum value. This is the highest-pressure point on the surface.
- **Separation point** — Where the boundary layer detaches from the cylinder surface due to the
  adverse pressure gradient. For laminar flow at moderate Re, separation occurs at roughly
  80°–90° from the front stagnation point.
- **Wake region** — The low-pressure recirculation zone downstream of the cylinder. The wake
  structure depends strongly on Reynolds number and drives the drag force.
- **Von Kármán vortex street** — At moderate Reynolds numbers (roughly $Re = 47$–$200$), vortices
  shed alternately from the top and bottom of the cylinder, producing a periodic oscillating
  wake pattern. This generates oscillating lift forces and a dominant shedding frequency.
- **Drag and lift forces** — The cylinder experiences a mean drag force (pressure drag dominates
  at moderate-to-high Re) and, during vortex shedding, an oscillating lift force at the
  shedding frequency.

> **💡 Tip:** For background on CFD fundamentals and the governing Navier-Stokes equations,
> see [`notes/01_short_intro_to_cfd.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md).

---

## Physics and Theory

### Reynolds Number

The Reynolds number is the single most important parameter governing this flow:

$$Re = U_\infty \cdot D / \nu$$

where `U∞` is the free-stream velocity, `D` is the cylinder diameter, and `ν` is the
kinematic viscosity of the fluid.

### Reynolds Number Regimes

| Re Range | Flow Regime | Description |
|----------|------------|-------------|
| Re < 5 | Creeping flow | Viscous forces dominate. Flow remains attached and symmetric fore-and-aft (nearly). No separation. |
| 5 – 40 | Attached vortex pair | A symmetric pair of steady recirculation vortices forms behind the cylinder. Flow is still steady and laminar. |
| 40 – 200 | Laminar vortex shedding | The wake becomes unstable. Vortices shed alternately from each side — the Von Kármán vortex street. Flow is periodic and laminar. |
| 200 – 300,000 | Subcritical turbulent | The wake is turbulent but the boundary layer on the cylinder remains laminar. Separation occurs before the top of the cylinder (~80°). Cd $\approx$ 1.0–1.2. |
| > 300,000 | Supercritical | The boundary layer transitions to turbulent before separation, which delays separation to ~120°. Dramatic drop in drag (the "drag crisis"). Cd drops to ~0.3. |

```
  Re < 5             5 < Re < 40          40 < Re < 200
  (Creeping)         (Steady vortex pair)  (Vortex shedding)

   →→  ○  →→         →→  ○ ◜◝  →→         →→  ○  ∿∿∿  →→
   →→     →→         →→    ◟◞  →→         →→     ∿∿∿  →→
  No separation      Symmetric wake        Alternating vortices
```

### Strouhal Number

For vortex shedding flows, the Strouhal number relates the shedding frequency to the
flow velocity and cylinder diameter:

$$St = f \cdot D / U_\infty$$

where `f` is the vortex shedding frequency. For a circular cylinder in the subcritical
regime, $St \approx 0.2$ is remarkably constant over a wide range of Reynolds numbers
(roughly $300 < Re < 300{,}000$). This is one of the most well-established empirical results
in fluid mechanics.

### Drag Coefficient

The drag coefficient is defined as:

$$Cd = F_D / (0.5 \cdot \rho \cdot U_\infty^2 \cdot A)$$

where `F_D` is the drag force, `ρ` is the fluid density, and `A` is the projected area
(A = D · L for a cylinder of span L). Key reference values:

| Re | Cd (approximate) | Source |
|----|-------------------|--------|
| 1 | 10.0 | Analytical / experimental |
| 10 | 2.8 | Experimental |
| 100 | 1.3 | Experimental |
| 1,000 | 1.0 | Experimental |
| 10,000 | 1.2 | Experimental |
| 100,000 | 1.2 | Experimental |
| 500,000 | 0.3 | Drag crisis — supercritical |

> **💡 Tip:** For details on turbulence modeling and when to use RANS vs LES, see
> [`notes/06_turbulence_models.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md).

---

## Geometry — Creating the Cylinder STL

Since we are using `snappyHexMesh` for mesh generation, we need a surface geometry file
(STL format) representing the cylinder. There are several ways to obtain one.

### Option 1: Python Script

The following Python script generates a simple cylinder STL using `numpy` and `numpy-stl`.
The cylinder has diameter `D = 1.0 m` and span `L = 1.0 m`, centered at the origin.

```python
#!/usr/bin/env python3
"""Generate a cylinder STL for OpenFOAM snappyHexMesh."""

import numpy as np
from stl import mesh

def create_cylinder_stl(filename, radius=0.5, length=1.0, n_facets=64):
    """Create a closed cylinder STL with the axis along the z-direction."""
    angles = np.linspace(0, 2 * np.pi, n_facets, endpoint=False)
    triangles = []

    for i in range(n_facets):
        j = (i + 1) % n_facets
        c0, s0 = radius * np.cos(angles[i]), radius * np.sin(angles[i])
        c1, s1 = radius * np.cos(angles[j]), radius * np.sin(angles[j])

        # Side face — two triangles per segment
        triangles.append([[c0, s0, 0], [c1, s1, 0], [c1, s1, length]])
        triangles.append([[c0, s0, 0], [c1, s1, length], [c0, s0, length]])

        # Bottom cap
        triangles.append([[0, 0, 0], [c1, s1, 0], [c0, s0, 0]])

        # Top cap
        triangles.append([[0, 0, length], [c0, s0, length], [c1, s1, length]])

    data = np.zeros(len(triangles), dtype=mesh.Mesh.dtype)
    for i, tri in enumerate(triangles):
        data["vectors"][i] = tri

    cylinder = mesh.Mesh(data)
    cylinder.save(filename)
    print(f"Saved {filename} with {len(triangles)} triangles")

if __name__ == "__main__":
    create_cylinder_stl("cylinder.stl")
```

Install the dependency and run:

```bash
pip install numpy-stl
python generate_cylinder.py
```

### Option 2: OpenSCAD

Create a file `cylinder.scad`:

```
cylinder(h=1, r=0.5, center=true, $fn=64);
```

Export to STL from the OpenSCAD GUI or command line:

```bash
openscad -o cylinder.stl cylinder.scad
```

### Option 3: FreeCAD

Open FreeCAD, create a cylinder primitive (Part → Cylinder, radius 0.5, height 1.0),
and export as STL via File → Export.

> **💡 Tip:** Whichever method you use, ensure the STL is watertight (no gaps or
> self-intersections). You can verify with `surfaceCheck cylinder.stl` in OpenFOAM.

---

## Complete Case Setup Guide

### Directory Structure

Create the full case directory:

```bash
mkdir -p $FOAM_RUN/cylinder
cd $FOAM_RUN/cylinder
mkdir -p 0 constant/triSurface system
```

Copy your cylinder STL into the geometry directory:

```bash
cp /path/to/cylinder.stl constant/triSurface/
```

The complete directory tree when finished:

```
cylinder/
├── 0/
│   ├── U
│   ├── p
│   ├── k
│   ├── epsilon
│   └── nut
├── constant/
│   ├── transportProperties
│   ├── turbulenceProperties
│   └── triSurface/
│       └── cylinder.stl
├── system/
│   ├── controlDict
│   ├── fvSchemes
│   ├── fvSolution
│   ├── blockMeshDict
│   ├── snappyHexMeshDict
│   └── decomposeParDict
└── Allrun
```

> **💡 Tip:** For a detailed explanation of the OpenFOAM case directory structure, see
> [`notes/02_openfoam_cases.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md).

---

### blockMeshDict — Background Mesh

The background mesh is a simple rectangular box that encloses the entire computational domain.
`snappyHexMesh` will later refine this mesh around the cylinder. The domain extends from
-5D upstream to 20D downstream, and $\pm$5D in the cross-stream direction, which is standard
practice for external flow simulations.

Create `system/blockMeshDict`:

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

// Cylinder diameter D = 1.0 m, centered at origin
// Domain: x = [-5, 20], y = [-5, 5], z = [0, 1]

scale   1;

vertices
(
    (-5  -5  0)    // 0
    (20  -5  0)    // 1
    (20   5  0)    // 2
    (-5   5  0)    // 3
    (-5  -5  1)    // 4
    (20  -5  1)    // 5
    (20   5  1)    // 6
    (-5   5  1)    // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (100 40 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }
    top
    {
        type symmetryPlane;
        faces
        (
            (3 7 6 2)
        );
    }
    bottom
    {
        type symmetryPlane;
        faces
        (
            (0 1 5 4)
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

mergePatchPairs
(
);
```

> **💡 Tip:** The `(100 40 1)` cell count gives 100 cells in the streamwise direction and
> 40 in the cross-stream direction. The single cell in z makes this a 2D simulation (with
> `empty` boundary condition on the front and back faces). For details on `blockMesh`, see
> [`notes/04_meshing.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/04_meshing.md).

---

### snappyHexMeshDict — Mesh Refinement Around the Cylinder

`snappyHexMesh` refines the background mesh around the STL geometry and snaps the mesh to
the cylinder surface. This dictionary controls refinement levels, surface layers, and mesh
quality.

Create `system/snappyHexMeshDict`:

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}

castellatedMesh true;
snap            true;
addLayers       true;

geometry
{
    cylinder.stl
    {
        type triSurfaceMesh;
        name cylinder;
    }

    refinementBox
    {
        type searchableBox;
        min (-1 -1.5 -0.5);
        max (8  1.5  1.5);
    }
}

castellatedMeshControls
{
    maxLocalCells       100000;
    maxGlobalCells      2000000;
    minRefinementCells  10;
    maxLoadUnbalance    0.10;
    nCellsBetweenLevels 3;
    features            ();

    refinementSurfaces
    {
        cylinder
        {
            level (4 4);

            patchInfo
            {
                type wall;
            }
        }
    }

    resolveFeatureAngle 30;

    refinementRegions
    {
        refinementBox
        {
            mode inside;
            levels ((1E15 2));
        }
    }

    locationInMesh (5 0 0.5);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch            3;
    tolerance               2.0;
    nSolveIter              100;
    nRelaxIter              5;
    nFeatureSnapIter        10;
    implicitFeatureSnap     false;
    explicitFeatureSnap     true;
    multiRegionFeatureSnap  false;
}

addLayersControls
{
    relativeSizes       true;

    layers
    {
        "cylinder.*"
        {
            nSurfaceLayers 5;
        }
    }

    expansionRatio      1.2;
    finalLayerThickness 0.3;
    minThickness        0.1;
    nGrow               0;
    featureAngle        60;
    slipFeatureAngle    30;
    nRelaxIter          3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals      3;
    nSmoothThickness    10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle  90;
    nBufferCellsNoExtrude 0;
    nLayerIter          50;
}

meshQualityControls
{
    maxNonOrtho         65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave          80;
    minVol              1e-13;
    minTetQuality       1e-15;
    minArea             -1;
    minTwist            0.02;
    minDeterminant      0.001;
    minFaceWeight       0.05;
    minVolRatio         0.01;
    minTriangleTwist    -1;
    nSmoothScale        4;
    errorReduction      0.75;
}

writeFlags
(
    scalarLevels
    layerSets
    layerFields
);
```

> **💡 Tip:** The `level (4 4)` on the cylinder surface means 4 levels of refinement
> (each level halves the cell size). With a background cell size of ~0.25 m, level 4 gives
> cells of ~0.016 m near the cylinder. The 5 boundary layers with expansion ratio 1.2 help
> resolve the boundary layer. Adjust these for your target Reynolds number.

---

### Boundary Conditions

All field files go in the `0/` directory. Each file specifies initial and boundary conditions
for one variable.

> **💡 Tip:** For a comprehensive reference on boundary condition types, see
> [`notes/05_boundary_conditions.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md).

#### Velocity — `0/U`

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (1 0 0);

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (1 0 0);    // U∞ = 1 m/s in x-direction
    }

    outlet
    {
        type            zeroGradient;
    }

    cylinder
    {
        type            noSlip;             // Wall: zero velocity
    }

    top
    {
        type            symmetryPlane;
    }

    bottom
    {
        type            symmetryPlane;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

#### Pressure — `0/p`

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            zeroGradient;
    }

    outlet
    {
        type            fixedValue;
        value           uniform 0;          // Reference pressure at outlet
    }

    cylinder
    {
        type            zeroGradient;
    }

    top
    {
        type            symmetryPlane;
    }

    bottom
    {
        type            symmetryPlane;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

#### Turbulent Kinetic Energy — `0/k`

The inlet values for `k` and `epsilon` are estimated from a turbulence intensity of 1%
and a turbulent length scale of 0.07D:

$$k = 1.5 \cdot (U_\infty \cdot I)^2 = 1.5 \cdot (1 \cdot 0.01)^2 = 1.5 \times 10^{-4} \text{ m}^2/\text{s}^2$$

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 1.5e-04;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 1.5e-04;
    }

    outlet
    {
        type            zeroGradient;
    }

    cylinder
    {
        type            kqRWallFunction;
        value           uniform 1.5e-04;
    }

    top
    {
        type            symmetryPlane;
    }

    bottom
    {
        type            symmetryPlane;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

#### Turbulent Dissipation Rate — `0/epsilon`

$$\varepsilon = C_\mu^{0.75} \cdot k^{1.5} / l = 0.09^{0.75} \cdot (1.5 \times 10^{-4})^{1.5} / 0.07 \approx 4.85 \times 10^{-7} \text{ m}^2/\text{s}^3$$

For simplicity we use a rounded value:

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      epsilon;
}

dimensions      [0 2 -3 0 0 0 0];

internalField   uniform 5e-07;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 5e-07;
    }

    outlet
    {
        type            zeroGradient;
    }

    cylinder
    {
        type            epsilonWallFunction;
        value           uniform 5e-07;
    }

    top
    {
        type            symmetryPlane;
    }

    bottom
    {
        type            symmetryPlane;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

#### Turbulent Viscosity — `0/nut`

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 0;
    }

    outlet
    {
        type            calculated;
        value           uniform 0;
    }

    cylinder
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    top
    {
        type            symmetryPlane;
    }

    bottom
    {
        type            symmetryPlane;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

---

### Transport and Turbulence Properties

#### `constant/transportProperties`

Set the kinematic viscosity to achieve your target Reynolds number. For $Re = 1000$ with
$U_\infty = 1$ m/s and $D = 1$ m, we need $\nu = U_\infty \cdot D / Re = 10^{-3}$ m$^2$/s:

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 1e-03;
```

#### `constant/turbulenceProperties`

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  RAS;

RAS
{
    RASModel        kEpsilon;
    turbulence      on;
    printCoeffs     on;
}
```

> **💡 Tip:** For a comparison of available turbulence models (k-$\varepsilon$, k-$\omega$ SST,
> Spalart-Allmaras) and guidance on which to choose, see
> [`notes/06_turbulence_models.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md).

---

### Control and Solver Settings

#### `system/controlDict`

This configuration includes `forceCoeffs` function objects to automatically compute drag
and lift coefficients during the simulation:

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     simpleFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1000;

deltaT          1;

writeControl    timeStep;
writeInterval   100;

purgeWrite      3;

writeFormat     ascii;
writePrecision  8;
writeCompression off;

timeFormat      general;
timePrecision   6;

runTimeModifiable true;

functions
{
    forceCoeffs
    {
        type            forceCoeffs;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   1;

        patches         (cylinder);
        rho             rhoInf;
        rhoInf          1;          // Density [kg/m³]
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        CofR            (0 0 0);    // Center of rotation
        pitchAxis       (0 0 1);
        magUInf         1;          // Free-stream velocity magnitude
        lRef            1;          // Reference length (diameter)
        Aref            1;          // Reference area (D * span)
    }

    fieldAverage
    {
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        timeStart       500;        // Start averaging after initial transient

        fields
        (
            U
            {
                mean        on;
                prime2Mean  on;
                base        time;
            }
            p
            {
                mean        on;
                prime2Mean  on;
                base        time;
            }
        );
    }
}
```

#### `system/fvSchemes`

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
    grad(p)         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

wallDist
{
    method          meshWave;
}
```

> **💡 Tip:** `linearUpwind` for the momentum equation provides a good balance between
> accuracy and stability. For details on discretization schemes and their accuracy, see
> [`notes/03_openfoam_dictionaries.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md).

#### `system/fvSolution`

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.1;
        smoother        GaussSeidel;
    }

    "(U|k|epsilon)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-06;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 1;
    consistent      yes;

    residualControl
    {
        p               1e-5;
        U               1e-5;
        "(k|epsilon)"   1e-5;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        k               0.7;
        epsilon         0.7;
    }
}
```

> **💡 Tip:** For an in-depth explanation of linear solvers (GAMG, PCG, smoothSolver)
> and how to choose tolerances, see
> [`notes/09_linear_solvers.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/09_linear_solvers.md).

#### `system/decomposeParDict`

Required only for parallel runs:

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

numberOfSubdomains 4;

method          scotch;
```

> **💡 Tip:** For parallel decomposition strategies and MPI configuration, see
> [`notes/07_parallelization.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/07_parallelization.md).

---

## Running the Simulation

### Serial Execution

Run the full simulation sequence:

```bash
cd $FOAM_RUN/cylinder

# 1. Generate background mesh
blockMesh

# 2. Refine mesh around cylinder using STL geometry
snappyHexMesh -overwrite

# 3. Verify mesh quality
checkMesh

# 4a. Run steady-state solver (RANS, time-averaged solution)
simpleFoam

# 4b. OR run transient solver (captures vortex shedding)
#     Change ddtSchemes to "Euler" and application to "pimpleFoam" first
# pimpleFoam
```

### Parallel Execution

For larger meshes, decompose and run in parallel:

```bash
cd $FOAM_RUN/cylinder

blockMesh
snappyHexMesh -overwrite
checkMesh

# Decompose the domain
decomposePar

# Run on 4 processors
mpirun -np 4 simpleFoam -parallel

# Reconstruct the results
reconstructPar
```

### Transient Simulation (Vortex Shedding)

To capture vortex shedding, switch to a transient solver. Modify `system/controlDict`:

```
application     pimpleFoam;
endTime         100;            // Physical seconds
deltaT          0.005;          // Time step — check CFL < 1
writeInterval   200;            // Write every 200 time steps
```

And change `system/fvSchemes` time discretization:

```
ddtSchemes
{
    default         Euler;      // First-order implicit for transient
}
```

> **💡 Tip:** For vortex shedding at $Re = 100$ (laminar), you can use `pimpleFoam` with
> `simulationType laminar` in `turbulenceProperties`. The Strouhal number should converge
> to $St \approx 0.164$. For guidance on choosing time steps using the CFL condition, see
> [`notes/08_cfl_number.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md).

---

## Post-Processing

### Force Coefficients

If you included the `forceCoeffs` function object in `controlDict` (as shown above), the
drag and lift coefficients are written to `postProcessing/forceCoeffs/0/coefficient.dat`.

Plot the convergence history:

```bash
gnuplot -persist -e "
    set xlabel 'Iteration';
    set ylabel 'Coefficient';
    set title 'Force Coefficients';
    plot 'postProcessing/forceCoeffs/0/coefficient.dat' using 1:2 with lines title 'Cd',
         '' using 1:3 with lines title 'Cl'
"
```

### Residuals

Monitor convergence by plotting residuals from the solver log:

```bash
foamLog simpleFoam.log
gnuplot -persist -e "
    set logscale y;
    set xlabel 'Iteration';
    set ylabel 'Residual';
    plot 'logs/p_0' with lines title 'p',
         'logs/Ux_0' with lines title 'Ux',
         'logs/Uy_0' with lines title 'Uy',
         'logs/k_0' with lines title 'k',
         'logs/epsilon_0' with lines title 'epsilon'
"
```

### ParaView Visualization

Open the case directly in ParaView (no conversion needed with modern ParaView):

```bash
paraview cylinder.foam &
# Create an empty .foam file if it doesn't exist:
# touch cylinder.foam
```

**Recommended visualizations:**

1. **Pressure contours** — Apply a "Surface" representation colored by `p`. Observe the
   high-pressure stagnation region at the front and the low-pressure wake behind the cylinder.

2. **Velocity streamlines** — Use Filters → Stream Tracer, seeded from a line upstream of
   the cylinder. This reveals the flow separation and recirculation zone.

3. **Vorticity / Q-criterion** — For transient simulations, compute Q-criterion to visualize
   vortex structures:
   - Filters → Calculator: `0.5*(mag(Vorticity)^2 - mag(Strain)^2)`
   - Or use the built-in Q-Criterion filter if available.

4. **Wake profile** — Use Filters → Plot Over Line to extract velocity profiles at various
   downstream locations (x/D = 1, 2, 5, 10) for comparison with experimental data.

---

## Expected Results

### What to Look For

| Quantity | Expected Value (Re $\approx$ 1000) | Notes |
|----------|---------------------------|-------|
| Drag coefficient (Cd) | ~1.0 – 1.2 | Compare with Wieselsberger (1922) data |
| Separation angle | ~80° from stagnation point | Laminar boundary layer separation |
| Wake length (steady) | ~1–2 D | Recirculation zone behind cylinder |
| Strouhal number (if transient) | ~0.2 | Vortex shedding frequency |
| Lift coefficient (Cl) | ~0 (mean), oscillating | Mean $\approx$ 0 by symmetry; amplitude depends on Re |

### Convergence Criteria

For a steady-state `simpleFoam` run:

- All residuals (p, U, k, epsilon) should drop below 1e-5.
- The drag coefficient should stabilize to a constant value.
- The simulation typically converges within 500–1000 iterations for this mesh.

### Common Issues

- **Divergence at startup** — Reduce relaxation factors (e.g., p → 0.2, U → 0.5) and
  increase `nNonOrthogonalCorrectors` if mesh quality is poor.
- **Asymmetric steady solution** — At $Re > 47$, the flow naturally wants to shed vortices.
  A steady solver will find one of the two possible asymmetric solutions or oscillate.
  Switch to `pimpleFoam` if you need the physical transient behavior.
- **High y+ values** — Check with `yPlus` function object. For wall functions, aim for
  30 < y+ < 300. If y+ < 30, add more boundary layers or switch to a low-Re model.

---

## Exercises

1. **Reynolds Number Study** — Run the simulation at $Re = 20, 100, 200, 1000$, and $10000$
   by adjusting the kinematic viscosity (`nu` in `transportProperties`). Compare drag
   coefficients with published data and observe the different flow regimes.

2. **Mesh Refinement Study** — Run with 3 progressively finer meshes (e.g., refinement
   levels 3, 4, and 5 in `snappyHexMeshDict`). Plot Cd vs. cell count to assess mesh
   independence. The Richardson extrapolation can give an estimate of the grid-converged value.

3. **Turbulence Model Comparison** — Run the same case with k-$\varepsilon$, k-$\omega$ SST, and
   Spalart-Allmaras models. Compare predicted Cd, separation angle, and wake structure.

4. **3D Extension** — Remove the `empty` boundary condition on `frontAndBack`, extend the
   domain to 4D in the spanwise direction, and add cells in z. Observe 3D wake instabilities
   (mode A and mode B) that appear at $Re > 190$.

5. **Vortex Shedding Frequency** — Run a transient simulation at $Re = 100$ using `pimpleFoam`
   with laminar flow. Measure the vortex shedding frequency from the lift coefficient
   time history and compute the Strouhal number. Compare with $St \approx 0.164$ from literature.

6. **Domain Size Sensitivity** — Test the effect of domain boundaries by running with
   smaller (10D downstream) and larger (40D downstream) domains. How does the outlet
   distance affect the predicted Cd?

---

## Allrun Script

Save the following as `Allrun` in the case directory and make it executable:

```bash
#!/bin/bash
#------------------------------------------------------------------------------
# Allrun script for flow around a cylinder
# Usage: ./Allrun [serial|parallel]
#------------------------------------------------------------------------------

cd "${0%/*}" || exit 1
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"

MODE="${1:-serial}"
NPROCS=4

echo "============================================="
echo "  Flow Around a Cylinder — OpenFOAM Case"
echo "  Mode: $MODE"
echo "============================================="

# Clean previous results
echo "Cleaning case..."
foamCleanCase 2>/dev/null
rm -rf postProcessing logs

# Generate background mesh
echo "Running blockMesh..."
runApplication blockMesh

# Refine mesh around cylinder
echo "Running snappyHexMesh..."
runApplication snappyHexMesh -overwrite

# Check mesh quality
echo "Running checkMesh..."
runApplication checkMesh

if [ "$MODE" = "parallel" ]; then
    echo "Decomposing for $NPROCS processors..."
    runApplication decomposePar

    echo "Running simpleFoam in parallel on $NPROCS processors..."
    runParallel simpleFoam

    echo "Reconstructing..."
    runApplication reconstructPar
else
    echo "Running simpleFoam..."
    runApplication simpleFoam
fi

echo "============================================="
echo "  Simulation complete!"
echo "  Check postProcessing/forceCoeffs/ for Cd/Cl"
echo "  Run: paraview cylinder.foam"
echo "============================================="
```

Make it executable:

```bash
chmod +x Allrun
```

Run:

```bash
./Allrun              # Serial
./Allrun parallel     # Parallel on 4 processors
```

---

## References

1. Schlichting, H., & Gersten, K. (2017). *Boundary-Layer Theory* (9th ed.). Springer.
   The definitive reference for boundary layer separation and wake flow theory.

2. Zdravkovich, M. M. (1997). *Flow Around Circular Cylinders, Vol. 1: Fundamentals*.
   Oxford University Press. Comprehensive review of circular cylinder aerodynamics.

3. Wieselsberger, C. (1922). "New data on the laws of fluid resistance." *NACA TN-84*.
   Classic experimental data for drag coefficient vs. Reynolds number.

4. Roshko, A. (1954). "On the development of turbulent wakes from vortex streets."
   *NACA Report 1191*. Foundational work on vortex shedding and Strouhal number.

5. Williamson, C. H. K. (1996). "Vortex dynamics in the cylinder wake."
   *Annual Review of Fluid Mechanics*, 28, 477–539. Modern review of wake instabilities.

6. OpenFOAM Foundation. *OpenFOAM User Guide*.
   [https://www.openfoam.com/documentation/user-guide](https://www.openfoam.com/documentation/user-guide)

7. Versteeg, H. K., & Malalasekera, W. (2007). *An Introduction to Computational Fluid
   Dynamics: The Finite Volume Method* (2nd ed.). Pearson. Accessible CFD textbook covering
   turbulence modeling and discretization schemes used in this tutorial.
