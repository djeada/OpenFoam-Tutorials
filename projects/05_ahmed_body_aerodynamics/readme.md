# Ahmed Body Aerodynamics

![Difficulty: Advanced](https://img.shields.io/badge/difficulty-advanced-red)
![Solver: simpleFoam](https://img.shields.io/badge/solver-simpleFoam-blue)
![Turbulence: k--ω SST](https://img.shields.io/badge/turbulence-k--ω_SST-green)
![Mesh: snappyHexMesh](https://img.shields.io/badge/mesh-snappyHexMesh-orange)
![Parallel: MPI](https://img.shields.io/badge/parallel-MPI-purple)

Steady-state RANS simulation of turbulent flow around the Ahmed reference body
using **simpleFoam** with **k-ω SST** turbulence modelling. This tutorial walks
through every file you need to reproduce the classic automotive-aerodynamics
benchmark and validate against published experimental data.

---

## Table of Contents

1. [Problem Description](#problem-description)
2. [Physics and Theory](#physics-and-theory)
3. [Case Directory Structure](#case-directory-structure)
4. [Geometry — Obtaining the Ahmed Body STL](#geometry--obtaining-the-ahmed-body-stl)
5. [Mesh Strategy](#mesh-strategy)
6. [Boundary Conditions](#boundary-conditions)
7. [Turbulence Setup — k-ω SST](#turbulence-setup--k-ω-sst)
8. [Solver Configuration](#solver-configuration)
9. [Parallel Decomposition](#parallel-decomposition)
10. [Force Coefficient Monitoring](#force-coefficient-monitoring)
11. [Running the Case](#running-the-case)
12. [Expected Results](#expected-results)
13. [Post-Processing Guide](#post-processing-guide)
14. [Exercises](#exercises)
15. [References](#references)

---

## Problem Description

The **Ahmed body** is a simplified, generic car shape introduced by Ahmed,
Ramm, and Faltin in 1984. It was designed to capture the essential aerodynamic
features of a road vehicle — particularly the rear-end flow separation and wake
structure — while being simple enough to manufacture accurately and simulate
reproducibly.

```
                    25° slant angle
                         ╱│
    ┌───────────────────╱ │
    │                  ╱  │ 288mm
    │   Ahmed Body    ╱   │
    │                │    │
    │   1044mm       │    │
    └────────────────┴────┘
    ────────────────────────── Ground plane
         50mm ground clearance

    Front view:              Side view (3D perspective):
    ┌──────────┐             ┌──────────────────╲
    │          │ 389mm       │                   ╲
    │          │             │                    │
    └──────────┘             └────────────────────┘
      389mm
```

### Standard Dimensions

| Parameter | Value |
|---|---|
| Overall length | 1044 mm |
| Width | 389 mm |
| Height | 288 mm |
| Ground clearance | 50 mm |
| Front-body radius (pillars) | 100 mm |
| Rear slant angle (baseline) | 25° |
| Stilts (4 cylindrical legs) | Ø 30 mm |

### Why the Ahmed Body Matters

The Ahmed body is **the** standard benchmark for automotive computational fluid
dynamics (CFD) validation. Virtually every turbulence model, meshing strategy,
and solver setting used in vehicle aerodynamics has been tested against Ahmed
body data at some point. Its value comes from:

- **Geometric simplicity** — the shape is fully defined by a handful of
  dimensions, so every research group is solving the same problem.
- **Rich flow physics** — despite the simple geometry, the flow contains
  separation, reattachment, longitudinal vortices, and a complex three-
  dimensional wake.
- **Extensive experimental database** — Ahmed et al. (1984) provided drag
  measurements across a range of slant angles, and Lienhart, Stoots, and
  Becker (2003) later added detailed LDA velocity profiles and Reynolds-stress
  measurements that are the gold standard for RANS validation.
- **Critical slant-angle sensitivity** — at approximately 30°, the rear-slant
  flow undergoes a dramatic topological change from an attached, vortex-
  dominated regime to a fully separated wake, providing a stringent test of
  turbulence-model accuracy.

---

## Physics and Theory

### Flow Regime

The standard test condition uses a freestream velocity of **U∞ = 40 m/s** in
air (ν = 1.5 × 10⁻⁵ m²/s). Based on the body height H = 0.288 m:

```
Re_H = U∞ · H / ν = 40 × 0.288 / 1.5e-5 ≈ 768,000
```

This is a **fully turbulent external flow** in the subcritical Reynolds-number
range for bluff bodies.

### Wake Structure

The Ahmed body wake is dominated by three interacting vortex systems:

1. **A-pillar vortices** — generated along the rounded front edges of the body,
   these are relatively weak and dissipate quickly.
2. **C-pillar (longitudinal) vortices** — the most prominent feature. These
   counter-rotating vortex pairs roll up from the lateral edges of the rear
   slant and trail far downstream. They are the dominant source of induced drag
   at slant angles below 30°.
3. **Recirculation zone** — a near-wake region of reversed flow immediately
   behind the vertical base. Its size and shape depend strongly on the slant
   angle.

### Drag Breakdown

Total aerodynamic drag on the Ahmed body comes from two contributions:

- **Pressure drag (form drag)** — dominates at roughly 70–80 % of total drag.
  Comes from the low-pressure wake behind the base and the suction on the rear
  slant.
- **Skin-friction drag** — accounts for the remaining 20–30 %. Due to the bluff
  shape, this is much smaller than for a streamlined body.

### The Critical Slant Angle

At slant angles below roughly 30°, the flow over the rear slant stays partially
attached. Strong C-pillar vortices form and actually help maintain attachment
through downwash. As the slant angle increases toward 30°, vortex strength
grows and drag increases.

At approximately **30°**, the rear-slant flow **separates completely**. The
C-pillar vortices collapse, and the wake transitions from a three-dimensional,
vortex-dominated topology to a quasi-two-dimensional, fully separated wake.
Drag drops abruptly.

This sharp transition is extremely difficult for RANS turbulence models to
capture accurately, which is what makes the Ahmed body such a demanding
validation case.

### Ground Effect

The 50 mm ground clearance creates an accelerated flow passage underneath the
body. This under-body flow interacts with the wake, and the ground boundary
condition (stationary vs. moving wall) has a noticeable effect on the predicted
drag and wake symmetry.

### Turbulence Model Selection

**k-ω SST** (Menter 1994) is the recommended turbulence model for this case:

- Blends k-ω near walls (superior for adverse-pressure-gradient boundary layers)
  with k-ε in the freestream (avoids the freestream sensitivity of standard k-ω).
- Includes a shear-stress-transport limiter that prevents over-prediction of
  turbulent shear stress in adverse-pressure-gradient flows — critical for
  predicting separation on the rear slant.
- Extensively validated for external vehicle aerodynamics.
- **k-ε models** tend to over-predict the attachment on the slant (under-predict
  drag) because they over-estimate turbulent mixing.

---

## Case Directory Structure

```
ahmed_body/
├── 0/
│   ├── U                          # Velocity field
│   ├── p                          # Pressure field (kinematic)
│   ├── k                          # Turbulent kinetic energy
│   ├── omega                      # Specific dissipation rate
│   └── nut                        # Turbulent viscosity
├── constant/
│   ├── transportProperties        # Fluid properties (nu)
│   ├── turbulenceProperties       # RANS / k-omega SST selection
│   └── triSurface/
│       └── ahmed_body.stl         # Surface geometry for snappyHexMesh
└── system/
    ├── controlDict                # Solver control + forceCoeffs
    ├── fvSchemes                  # Discretisation schemes
    ├── fvSolution                 # Linear solvers & SIMPLE settings
    ├── blockMeshDict              # Background hex mesh
    ├── snappyHexMeshDict          # Body-fitted mesh with refinement
    ├── decomposeParDict           # Domain decomposition for MPI
    ├── meshQualityDict            # Mesh quality controls
    └── surfaceFeatureExtractDict  # Feature-edge extraction
```

---

## Geometry — Obtaining the Ahmed Body STL

You need an STL surface of the Ahmed body placed at the correct location in
your computational domain.

### Option 1 — Download

Several public repositories provide Ahmed body STL files:

- The original OpenFOAM tutorial collection on GitLab
- GrabCAD / Thingiverse (search "Ahmed body CFD")
- University course material repositories

Make sure the dimensions match the standard specification above and that the
body is oriented with the long axis along **x**, width along **y**, and height
along **z**.

### Option 2 — Create in FreeCAD

1. Open FreeCAD → Part workbench.
2. Create a box: 1044 × 389 × 288 mm.
3. Fillet the front face edges with radius = 100 mm.
4. Cut the rear-top edge at 25° to form the slant.
5. Add four cylindrical stilts (Ø 30 mm, length 50 mm).
6. Export as STL (binary, fine resolution).

### Option 3 — Parametric Script (OpenSCAD)

```openscad
// Simplified Ahmed body (no stilts, no front radius)
module ahmed_body() {
    difference() {
        cube([1044, 389, 288], center=false);
        // Rear slant cut at 25 degrees
        translate([1044, 0, 288])
            rotate([0, 25, 0])
                cube([400, 389, 400], center=false);
    }
}
translate([0, 0, 50])  // Ground clearance
    ahmed_body();
```

### Positioning the STL

Place the body so that:

- The nose is at approximately x = 1.0 m from the inlet (adequate approach
  length).
- The body is centred in y (symmetry plane at y = 0 if using half-model).
- The bottom of the stilts sits on the ground plane (z = 0).

---

## Mesh Strategy

Meshing is the single most important factor for accuracy in external automotive
CFD. The Ahmed body requires careful attention to:

- Adequate domain size
- Surface and volume refinement
- Boundary-layer resolution
- Wake refinement

### Domain Sizing

Use the body length L = 1.044 m as the reference:

| Boundary | Distance from body |
|---|---|
| Inlet | 5 L upstream of nose |
| Outlet | 10 L downstream of base |
| Lateral walls | 5 L from body centreline |
| Top | 5 L above ground |
| Ground | z = 0 (coincides with floor) |

This gives a domain of approximately:

- x: -5.0 m to +12.0 m (17 m total)
- y: -5.0 m to +5.0 m (10 m total)
- z: 0.0 m to 5.0 m (5 m total)

### blockMeshDict — Background Mesh

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale   1;

vertices
(
    ( -5.0  -5.0   0.0 )   // 0
    ( 12.0  -5.0   0.0 )   // 1
    ( 12.0   5.0   0.0 )   // 2
    ( -5.0   5.0   0.0 )   // 3
    ( -5.0  -5.0   5.0 )   // 4
    ( 12.0  -5.0   5.0 )   // 5
    ( 12.0   5.0   5.0 )   // 6
    ( -5.0   5.0   5.0 )   // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (85 50 25) simpleGrading (1 1 1)
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

    ground
    {
        type wall;
        faces
        (
            (0 3 2 1)
        );
    }

    top
    {
        type symmetryPlane;
        faces
        (
            (4 5 6 7)
        );
    }

    side_left
    {
        type symmetryPlane;
        faces
        (
            (0 1 5 4)
        );
    }

    side_right
    {
        type symmetryPlane;
        faces
        (
            (3 7 6 2)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
```

This gives a base mesh of roughly 106,000 cells (85 × 50 × 25). The cell size
is approximately 0.2 m, which snappyHexMesh will refine down to the required
resolution.

### surfaceFeatureExtractDict

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      surfaceFeatureExtractDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ahmed_body.stl
{
    extractionMethod    extractFromSurface;

    extractFromSurfaceCoeffs
    {
        includedAngle   150;
    }

    writeObj            yes;
}

// ************************************************************************* //
```

### snappyHexMeshDict

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

castellatedMesh true;
snap            true;
addLayers       true;

geometry
{
    ahmed_body.stl
    {
        type triSurfaceMesh;
        name ahmed_body;
    }

    // Refinement box around the body
    refinementBox_body
    {
        type searchableBox;
        min ( -0.5  -0.4   0.0 );
        max (  2.0   0.4   0.6 );
    }

    // Refinement box for the wake region
    refinementBox_wake
    {
        type searchableBox;
        min (  1.0  -0.6   0.0 );
        max (  4.0   0.6   0.6 );
    }

    // Fine refinement close to the body
    refinementBox_near
    {
        type searchableBox;
        min ( -0.2  -0.3   0.0 );
        max (  1.5   0.3   0.5 );
    }
}

castellatedMeshControls
{
    maxLocalCells       1000000;
    maxGlobalCells      5000000;
    minRefinementCells  10;
    maxLoadUnbalance    0.10;
    nCellsBetweenLevels 3;

    features
    (
        {
            file "ahmed_body.eMesh";
            level 6;
        }
    );

    refinementSurfaces
    {
        ahmed_body
        {
            level (5 6);

            patchInfo
            {
                type wall;
            }
        }
    }

    resolveFeatureAngle 30;

    refinementRegions
    {
        refinementBox_near
        {
            mode inside;
            levels ((1e15 4));
        }

        refinementBox_body
        {
            mode inside;
            levels ((1e15 3));
        }

        refinementBox_wake
        {
            mode inside;
            levels ((1e15 3));
        }
    }

    locationInMesh (3.0 0.0 2.5);  // Point outside the body, inside the domain

    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch    3;
    tolerance       2.0;
    nSolveIter      100;
    nRelaxIter      5;

    nFeatureSnapIter 10;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap false;
}

addLayersControls
{
    relativeSizes   true;

    layers
    {
        "ahmed_body.*"
        {
            nSurfaceLayers 5;
        }
        ground
        {
            nSurfaceLayers 3;
        }
    }

    expansionRatio          1.2;
    finalLayerThickness     0.5;
    minThickness            0.1;
    nGrow                   0;
    featureAngle            130;
    slipFeatureAngle        30;
    nRelaxIter              5;
    nSmoothSurfaceNormals   1;
    nSmoothNormals          3;
    nSmoothThickness        10;
    maxFaceThicknessRatio   0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle      90;
    nBufferCellsNoExtrude   0;
    nLayerIter              50;
}

meshQualityControls
{
    maxNonOrtho     65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave      80;
    minVol          1e-13;
    minTetQuality   -1e30;
    minArea         -1;
    minTwist        0.02;
    minDeterminant  0.001;
    minFaceWeight   0.02;
    minVolRatio     0.01;
    minTriangleTwist -1;

    nSmoothScale    4;
    errorReduction  0.75;
}

writeFormat ascii;

mergeTolerance 1e-6;

// ************************************************************************* //
```

### Mesh Quality Targets

| Metric | Target |
|---|---|
| y+ on body | ~30 (wall functions) or ~1 (resolved) |
| Surface refinement | Level 5–6 (~3–6 mm cells) |
| Wake refinement | Level 3 (~25 mm cells) |
| Boundary layers on body | 5 layers, expansion ratio 1.2 |
| Total cell count | ~2–5 million |

For wall-function meshes (y+ ≈ 30), the first cell height should be
approximately 0.8 mm at Re = 768,000. The 5-layer prism stack with expansion
ratio 1.2 captures the inner portion of the boundary layer while relying on
`omegaWallFunction` and `kqRWallFunction` for the near-wall treatment.

---

## Boundary Conditions

### Velocity — `0/U`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (40 0 0);

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (40 0 0);
    }

    outlet
    {
        type            zeroGradient;
    }

    ground
    {
        // Moving wall to simulate road — matches freestream velocity
        // to eliminate ground boundary-layer buildup.
        // Use "noSlip" for a stationary-ground simulation instead.
        type            fixedValue;
        value           uniform (40 0 0);
    }

    top
    {
        type            symmetryPlane;
    }

    side_left
    {
        type            symmetryPlane;
    }

    side_right
    {
        type            symmetryPlane;
    }

    ahmed_body
    {
        type            noSlip;
    }
}

// ************************************************************************* //
```

### Pressure — `0/p`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
        value           uniform 0;
    }

    ground
    {
        type            zeroGradient;
    }

    top
    {
        type            symmetryPlane;
    }

    side_left
    {
        type            symmetryPlane;
    }

    side_right
    {
        type            symmetryPlane;
    }

    ahmed_body
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //
```

### Turbulent Kinetic Energy — `0/k`

Inlet values are calculated from freestream turbulence intensity I = 1 % and
length scale l = 0.01 m:

```
k = 1.5 * (U∞ * I)² = 1.5 * (40 * 0.01)² = 0.24 m²/s²
```

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0.24;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 0.24;
    }

    outlet
    {
        type            zeroGradient;
    }

    ground
    {
        type            kqRWallFunction;
        value           uniform 0.24;
    }

    top
    {
        type            symmetryPlane;
    }

    side_left
    {
        type            symmetryPlane;
    }

    side_right
    {
        type            symmetryPlane;
    }

    ahmed_body
    {
        type            kqRWallFunction;
        value           uniform 0.24;
    }
}

// ************************************************************************* //
```

### Specific Dissipation Rate — `0/omega`

```
omega = k^0.5 / (C_mu^0.25 * l) = 0.24^0.5 / (0.09^0.25 * 0.01) ≈ 89.4 s⁻¹
```

Rounding to **90 s⁻¹** for the inlet value.

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 -1 0 0 0 0];

internalField   uniform 90;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 90;
    }

    outlet
    {
        type            zeroGradient;
    }

    ground
    {
        type            omegaWallFunction;
        value           uniform 90;
    }

    top
    {
        type            symmetryPlane;
    }

    side_left
    {
        type            symmetryPlane;
    }

    side_right
    {
        type            symmetryPlane;
    }

    ahmed_body
    {
        type            omegaWallFunction;
        value           uniform 90;
    }
}

// ************************************************************************* //
```

### Turbulent Viscosity — `0/nut`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    ground
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    top
    {
        type            symmetryPlane;
    }

    side_left
    {
        type            symmetryPlane;
    }

    side_right
    {
        type            symmetryPlane;
    }

    ahmed_body
    {
        type            nutkWallFunction;
        value           uniform 0;
    }
}

// ************************************************************************* //
```

### Boundary Condition Summary

| Patch | U | p | k | omega | nut |
|---|---|---|---|---|---|
| inlet | fixedValue (40 0 0) | zeroGradient | fixedValue 0.24 | fixedValue 90 | calculated |
| outlet | zeroGradient | fixedValue 0 | zeroGradient | zeroGradient | calculated |
| ground | fixedValue (40 0 0) | zeroGradient | kqRWallFunction | omegaWallFunction | nutkWallFunction |
| top | symmetryPlane | symmetryPlane | symmetryPlane | symmetryPlane | symmetryPlane |
| side_left | symmetryPlane | symmetryPlane | symmetryPlane | symmetryPlane | symmetryPlane |
| side_right | symmetryPlane | symmetryPlane | symmetryPlane | symmetryPlane | symmetryPlane |
| ahmed_body | noSlip | zeroGradient | kqRWallFunction | omegaWallFunction | nutkWallFunction |

---

## Turbulence Setup — k-ω SST

### Why k-ω SST for External Aerodynamics

The Menter k-ω SST model is preferred over the standard k-ε model for the
Ahmed body for several reasons:

1. **Adverse-pressure-gradient performance** — the shear-stress-transport
   limiter prevents the over-prediction of turbulent shear stress that causes
   k-ε to delay separation on the rear slant.
2. **Near-wall treatment** — k-ω formulation in the inner boundary layer is
   more accurate than k-ε wall functions for non-equilibrium boundary layers.
3. **Freestream behaviour** — the blending function switches to k-ε in the
   outer flow, avoiding the known freestream sensitivity of the Wilcox k-ω
   model.

### `constant/turbulenceProperties`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}

// ************************************************************************* //
```

### `constant/transportProperties`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 1.5e-05;

// ************************************************************************* //
```

### Inlet Turbulence Calculation

For freestream turbulence intensity I = 1 % and turbulent length scale
l = 0.01 m (approximately 3 % of body height):

| Quantity | Formula | Value |
|---|---|---|
| k | 1.5 × (U∞ × I)² | 0.24 m²/s² |
| ω | k^0.5 / (C_μ^0.25 × l) | ≈ 90 s⁻¹ |
| ν_t | k / ω | ≈ 2.7 × 10⁻³ m²/s |
| ν_t / ν | — | ≈ 178 |

These values give a turbulent viscosity ratio of about 178, which is reasonable
for a low-turbulence wind-tunnel freestream.

---

## Solver Configuration

### `system/controlDict`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     simpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         3000;

deltaT          1;

writeControl    timeStep;

writeInterval   500;

purgeWrite      3;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

functions
{
    #include "forceCoeffs"
}

// ************************************************************************* //
```

### `system/fvSchemes`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
    grad(k)         cellLimited Gauss linear 1;
    grad(omega)     cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;

    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
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

// ************************************************************************* //
```

### `system/fvSolution`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.01;
        smoother        GaussSeidel;
        nPreSweeps      0;
        nPostSweeps     2;
        cacheAgglomeration on;
        agglomerator    faceAreaPair;
        nCellsInCoarsestLevel 200;
        mergeLevels     1;
    }

    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.01;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 1;

    consistent      yes;    // SIMPLEC for faster convergence

    residualControl
    {
        p               1e-5;
        U               1e-5;
        "(k|omega)"     1e-5;
    }
}

relaxationFactors
{
    equations
    {
        U               0.7;
        k               0.7;
        omega           0.7;
    }
    fields
    {
        p               0.3;
    }
}

// ************************************************************************* //
```

### Solver Notes

- **SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations) is the
  standard algorithm for steady-state incompressible flow.
- **SIMPLEC** (`consistent yes`) allows higher relaxation factors and faster
  convergence than standard SIMPLE.
- **GAMG** (Geometric Agglomerated Algebraic Multigrid) for pressure is much
  faster than Krylov solvers for this type of elliptic equation.
- **Relaxation factors**: 0.7 for momentum/turbulence and 0.3 for pressure are
  good starting values. Reduce if convergence is unstable.
- **Residual targets**: 1 × 10⁻⁵ for all fields should give well-converged
  force coefficients (< 1 % variation in last 200 iterations).

---

## Parallel Decomposition

The Ahmed body mesh (2–5 million cells) is too large for efficient serial
execution on most workstations. MPI-parallel execution across 4–8 cores is
strongly recommended.

### `system/decomposeParDict`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains  6;

method          scotch;

// ************************************************************************* //
```

**Scotch** is the recommended decomposition method because it automatically
determines the optimal partition without requiring user-specified direction
weights. It minimises the inter-processor communication surface area, which is
especially important for complex snappyHexMesh geometries.

### `Allrun` Script

```bash
#!/bin/bash
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions

nProcs=$(foamDictionary -entry numberOfSubdomains -value system/decomposeParDict)

# 1. Generate background mesh
runApplication blockMesh

# 2. Extract surface features
runApplication surfaceFeatureExtract

# 3. Run snappyHexMesh in parallel
runApplication decomposePar
runParallel snappyHexMesh -overwrite

# 4. Check mesh quality
runParallel checkMesh -latestTime -allTopology -allGeometry

# 5. Run the solver
runParallel $(getApplication)

# 6. Reconstruct
runApplication reconstructParMesh -constant
runApplication reconstructPar -latestTime

echo ""
echo "Done. Force coefficients written to postProcessing/forceCoeffs/"
echo ""
```

### `Allclean` Script

```bash
#!/bin/bash
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions

cleanCase
rm -rf processor*
rm -rf constant/extendedFeatureEdgeMesh
rm -f constant/triSurface/*.eMesh
rm -rf postProcessing

echo "Case cleaned."
```

---

## Force Coefficient Monitoring

Create a file `system/forceCoeffs` that is `#include`d from `controlDict`:

### `system/forceCoeffs`

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

forceCoeffs
{
    type            forceCoeffs;
    libs            (forces);
    writeControl    timeStep;
    writeInterval   1;

    patches         (ahmed_body);

    rho             rhoInf;     // incompressible
    rhoInf          1.225;      // air density [kg/m³]

    // Freestream velocity magnitude
    magUInf         40;

    // Reference length (body length)
    lRef            1.044;

    // Reference area (frontal area = width x height)
    // A_ref = 0.389 m × 0.288 m = 0.112032 m²
    Aref            0.112032;

    // Drag direction (streamwise)
    dragDir         (1 0 0);

    // Lift direction (vertical)
    liftDir         (0 0 1);

    // Pitch axis
    pitchAxis       (0 1 0);

    // Centre of rotation
    CofR            (0.522 0 0.194);
}

// ************************************************************************* //
```

### Expected Force Coefficients

| Source | Slant Angle | Cd |
|---|---|---|
| Ahmed et al. (1984) experiment | 25° | 0.285 |
| Ahmed et al. (1984) experiment | 35° | 0.260 |
| k-ω SST RANS (typical) | 25° | 0.28 – 0.30 |
| k-ω SST RANS (typical) | 35° | 0.25 – 0.27 |

The k-ω SST model typically predicts Cd within 5 % of the experimental value
for the 25° slant case. Accuracy for the 35° case (fully separated) depends
heavily on mesh resolution in the wake.

---

## Running the Case

### Prerequisites

- OpenFOAM v2312 or later (ESI or Foundation version)
- An Ahmed body STL file in `constant/triSurface/`
- At least 8 GB of RAM (16 GB recommended for > 3 M cells)
- 4–8 CPU cores

### Step-by-Step

```bash
# 1. Source OpenFOAM environment
source /opt/openfoam/etc/bashrc    # adjust path as needed

# 2. Navigate to the case directory
cd ahmed_body

# 3. Create background mesh
blockMesh

# 4. Extract surface features for snapping
surfaceFeatureExtract

# 5. Decompose for parallel meshing
decomposePar

# 6. Run snappyHexMesh in parallel
mpirun -np 6 snappyHexMesh -overwrite -parallel

# 7. Check mesh quality
mpirun -np 6 checkMesh -latestTime -allTopology -allGeometry -parallel

# 8. Run the solver
mpirun -np 6 simpleFoam -parallel | tee log.simpleFoam

# 9. Reconstruct for post-processing
reconstructParMesh -constant
reconstructPar -latestTime

# 10. Visualise
paraFoam
```

Or simply:

```bash
chmod +x Allrun Allclean
./Allclean && ./Allrun
```

### Monitoring Convergence

During the run, you can monitor residuals and force coefficients:

```bash
# In a separate terminal — watch residuals
foamMonitor -l postProcessing/forceCoeffs/0/coefficient.dat

# Or plot with gnuplot
gnuplot -persist -e "
    set xlabel 'Iteration';
    set ylabel 'Cd';
    plot 'postProcessing/forceCoeffs/0/coefficient.dat' using 1:2 with lines title 'Cd'
"
```

Convergence criteria:

- All residuals below 1 × 10⁻⁵
- Force coefficients stable to within 0.1 % over the last 200 iterations
- Typically requires 1500–3000 iterations

---

## Expected Results

### Drag Coefficient

For a well-resolved mesh (~3 M cells) with k-ω SST:

- **25° slant**: Cd ≈ 0.285 ± 0.015 (experimental: 0.285)
- **35° slant**: Cd ≈ 0.255 ± 0.020 (experimental: 0.260)

### Wake Structure

The simulation should capture:

1. **C-pillar vortices** — a strong counter-rotating vortex pair trailing from
   the lateral edges of the rear slant. These are clearly visible in streamline
   plots and Q-criterion iso-surfaces.
2. **Recirculation zone** — a closed recirculation bubble on the rear vertical
   base, with a saddle point separating the upper (slant-influenced) and lower
   (under-body) recirculation regions.
3. **Horseshoe vortex** — at the junction of the body and the ground plane
   (visible primarily in fine-mesh simulations).

### Pressure Distribution

- **Slant surface**: strong suction (negative Cp) near the top edge of the
  slant, recovering toward ambient at the trailing edge.
- **Rear base**: large region of negative Cp (the primary source of pressure
  drag).
- **Front stagnation**: Cp ≈ 1.0 at the nose centreline.

### Wall Shear Stress

The surface streamlines on the slant should show:

- Attachment line along the centreline (for 25° case)
- Separation lines curving outward toward the C-pillar vortex cores
- Complex 3D pattern on the rear base

---

## Post-Processing Guide

### ParaView Visualisation

Open the case in ParaView:

```bash
paraFoam -builtin
# or
touch case.foam && paraview case.foam
```

#### Pressure Contours

1. Apply the `OpenFOAMReader` source.
2. Select the latest time step.
3. Color by `p` → select the "Cool to Warm" colour map.
4. Add a **Slice** filter at y = 0 (centreline) to see the symmetry-plane
   pressure field.

#### Streamlines

1. Add a **Stream Tracer** filter.
2. Set seed type to "Point Cloud" or "Line".
3. Place seeds upstream of the body and in the wake.
4. Color streamlines by velocity magnitude.

#### Vortex Identification (Q-Criterion)

1. Apply **Filters → Alphabetical → GradientOfUnstructuredDataSet** on `U`.
2. Use a **Calculator** to compute Q-criterion:

   ```
   Q = 0.5 * (||Ω||² - ||S||²)
   ```

   Or use the built-in `Q Criterion` filter if available.

3. Add a **Contour** filter on Q at a positive iso-value (e.g., Q = 5000).
4. Color the iso-surface by pressure to visualise the vortex cores.

#### Surface Cp Distribution

1. Apply an **ExtractBlock** filter to isolate the `ahmed_body` patch.
2. Use a **Calculator** to compute Cp:

   ```
   Cp = p / (0.5 * 40 * 40)
   ```

   (Note: OpenFOAM pressure `p` is kinematic, p/ρ, in m²/s².)

3. Color the body surface by Cp.

### Extracting Force Coefficients from Log

```bash
# Extract Cd vs iteration from the forceCoeffs output
grep "Cd" log.simpleFoam | tail -20

# Or plot directly from the postProcessing directory
cat postProcessing/forceCoeffs/0/coefficient.dat | tail -20
```

### Centreline Cp Extraction

Use the `sample` utility or ParaView's **Plot Over Line** filter:

1. In ParaView: **Filters → Data Analysis → Plot Over Line**
2. Set line from (−0.5, 0, 0.194) to (2.0, 0, 0.194) — centreline at
   mid-height of the body.
3. Plot Cp along this line and compare with Lienhart et al. (2003) data.

---

## Exercises

### Exercise 1 — Slant Angle Study

Compare the aerodynamics of the Ahmed body at different rear slant angles:

1. Run the baseline 25° case as described above.
2. Modify the STL geometry to use a **35°** slant angle.
3. Run the same simulation with identical mesh settings and boundary conditions.
4. Compare:
   - Drag coefficient (Cd)
   - Wake structure (streamlines)
   - C-pillar vortex strength
   - Pressure distribution on the slant
5. You should observe a **drop** in Cd at 35° due to full flow separation
   on the slant and collapse of the C-pillar vortices.

### Exercise 2 — Mesh Sensitivity Study

Assess the influence of mesh resolution on the predicted drag:

1. Create three meshes with increasing refinement:
   - Coarse: ~1 M cells (reduce refinement levels by 1)
   - Medium: ~3 M cells (baseline settings above)
   - Fine: ~6 M cells (increase refinement levels by 1)
2. Run all three and compare Cd.
3. Plot Cd vs. cell count and assess whether the solution is mesh-independent.

### Exercise 3 — Turbulence Model Comparison

Repeat the 25° simulation with different turbulence models:

1. **k-ε (standard)** — `kEpsilon` in turbulenceProperties
2. **k-ε (realizable)** — `realizableKE`
3. **k-ω SST** — baseline
4. **Spalart-Allmaras** — `SpalartAllmaras` (single-equation model)

Compare the predicted Cd, slant separation pattern, and wake structure. You
should find that k-ε models tend to under-predict drag (delayed separation)
while k-ω SST gives the best agreement with experiments.

> **Note**: When switching turbulence models, you must update the `0/` boundary
> condition files to match the new model's required fields (e.g., `epsilon`
> instead of `omega`, or `nuTilda` for Spalart-Allmaras).

### Exercise 4 — Ground Clearance Effect

Study how ground clearance affects the aerodynamics:

1. Run the baseline case (50 mm clearance, moving ground).
2. Modify to **25 mm** ground clearance.
3. Modify to **100 mm** ground clearance.
4. Compare the drag coefficient and under-body flow acceleration for each.
5. Also try a **stationary ground** (change ground BC in `0/U` to `noSlip`)
   and observe the ground boundary-layer effect.

### Exercise 5 — Half-Model Simulation

For faster turnaround, exploit the geometric symmetry:

1. Slice the domain and STL at y = 0.
2. Add a `symmetryPlane` boundary condition on the cut face.
3. Halve the number of processors.
4. Verify that the full-model and half-model Cd values agree.

---

## References

1. **Ahmed, S.R., Ramm, G., and Faltin, G.** (1984). "Some Salient Features
   of the Time-Averaged Ground Vehicle Wake." *SAE Technical Paper 840300*.
   DOI: [10.4271/840300](https://doi.org/10.4271/840300).
   — *The original paper defining the Ahmed body geometry and presenting drag
   measurements across slant angles from 0° to 40°.*

2. **Lienhart, H., Stoots, C., and Becker, S.** (2003). "Flow and Turbulence
   Structures in the Wake of a Simplified Car Model (Ahmed Model)." *New
   Results in Numerical and Experimental Fluid Mechanics III*, Notes on
   Numerical Fluid Mechanics, Vol. 77, Springer, pp. 323–330.
   DOI: [10.1007/978-3-540-45466-3_39](https://doi.org/10.1007/978-3-540-45466-3_39).
   — *Detailed LDA measurements of mean velocity and Reynolds stresses in the
   wake. The standard dataset for RANS and LES validation.*

3. **Menter, F.R.** (1994). "Two-Equation Eddy-Viscosity Turbulence Models for
   Engineering Applications." *AIAA Journal*, 32(8), pp. 1598–1605.
   DOI: [10.2514/3.12149](https://doi.org/10.2514/3.12149).
   — *The paper introducing the k-ω SST turbulence model used in this tutorial.*

4. **Hinterberger, C., Garcia-Villalba, M., and Rodi, W.** (2004). "Large Eddy
   Simulation of Flow around the Ahmed Body." *The Aerodynamics of Heavy
   Vehicles: Trucks, Buses, and Trains*, Lecture Notes in Applied and
   Computational Mechanics, Vol. 19, Springer, pp. 77–87.
   — *LES reference showing the unsteady vortex dynamics not captured by RANS.*

5. **Guilmineau, E.** (2008). "Computational Study of Flow around a Simplified
   Car Body." *Journal of Wind Engineering and Industrial Aerodynamics*, 96(6–7),
   pp. 1207–1217.
   DOI: [10.1016/j.jweia.2007.06.041](https://doi.org/10.1016/j.jweia.2007.06.041).
   — *Systematic comparison of k-ω SST and k-ε for the Ahmed body at 25° and
   35° slant angles.*

---

## Troubleshooting

### Mesh Issues

| Problem | Solution |
|---|---|
| snappyHexMesh crashes | Reduce `maxGlobalCells`, check STL is watertight |
| Poor layer addition | Increase `nRelaxIter`, reduce `featureAngle` |
| High non-orthogonality | Add `nNonOrthogonalCorrectors 2` in fvSolution |
| Cells inside the body | Check `locationInMesh` is outside the STL |

### Solver Divergence

| Problem | Solution |
|---|---|
| Pressure diverges early | Reduce relaxation factors (p → 0.2, U → 0.5) |
| k/omega blow up | Use `bounded Gauss upwind` for turbulence, add `bound` utility |
| Oscillating Cd | Increase iterations, check mesh quality near separation |
| Very slow convergence | Use SIMPLEC (`consistent yes`), check mesh quality |

### Unexpected Cd Values

| Problem | Likely Cause |
|---|---|
| Cd too low (< 0.25) | Insufficient wake refinement, delayed separation |
| Cd too high (> 0.35) | Mesh too coarse on body surface, poor layer quality |
| Cd keeps changing | Not converged — run more iterations or check mesh |

---

*This tutorial is part of the [OpenFOAM Tutorials](../../README.md) collection.*
