# NACA 0012 Airfoil Analysis

> **Level:** Intermediate–Advanced · **Solver:** `simpleFoam` · **Flow:** Steady-state, incompressible, turbulent

---

## Table of Contents

1. [Problem Description](#problem-description)
2. [Physics & Theory](#physics--theory)
3. [Case Structure](#case-structure)
4. [Mesh Generation Pipeline](#mesh-generation-pipeline)
5. [Boundary Conditions](#boundary-conditions)
6. [Turbulence Setup](#turbulence-setup)
7. [Solver Configuration](#solver-configuration)
8. [How to Run](#how-to-run)
9. [Expected Results](#expected-results)
10. [Post-Processing Guide](#post-processing-guide)
11. [Exercises](#exercises)
12. [Comparison: Laminar vs Turbulent](#comparison-laminar-vs-turbulent)
13. [References](#references)

---

## Problem Description

This case simulates steady-state turbulent airflow over a **NACA 0012 airfoil** — one of
the most widely studied airfoil profiles in aerodynamics and a standard validation case
for CFD codes worldwide.

```
  Free stream U∞ →→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→
  →→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→
  →→→→→→→→→→    ╱‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾╲   →→→→→→→→→
  →→→→→→→→→→  ╱    NACA 0012        ╲  →→→→→→→→→
  →→→→→→→→→→ ╱________________________╲ →→→→→→→→→
  →→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→
  →→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→

  Inlet ──→                                 ──→ Outlet
             Stagnation    Suction side
             point         (low pressure)
```

### What NACA 0012 Means

The four-digit NACA naming convention encodes the airfoil geometry:

| Digit(s) | Meaning                         | NACA 0012 Value             |
|-----------|----------------------------------|-----------------------------|
| **00**    | Maximum camber (% of chord)     | 0 % → **symmetric airfoil** |
| **12**    | Maximum thickness (% of chord)  | **12 %** at 30 % chord      |

Because the camber is zero the profile is perfectly symmetric about its chord line, which
means it produces **zero lift at zero angle of attack** — an ideal property for
validation studies.

### Key Aerodynamic Concepts

- **Lift** — the net force perpendicular to the free-stream direction. For a symmetric
  airfoil at 0° angle of attack, the theoretical lift coefficient Cl = 0.
- **Drag** — the net force parallel to the free-stream direction, arising from both
  pressure (form drag) and viscous (skin friction) contributions.
- **Angle of Attack (α)** — the angle between the chord line and the incoming flow. This
  case is set up at α = 0°; see [Exercises](#exercises) for varying α.
- **Pressure Coefficient (Cp)** — the non-dimensional pressure distribution over the
  airfoil surface, a primary quantity for validation.

### Why This Is a Standard Validation Case

The NACA 0012 has been exhaustively tested in wind tunnels since the 1930s (Abbott & Von
Doenhoff, 1959). Extensive experimental data for Cp distributions, lift curves, and drag
polars are publicly available from NASA, making it the de facto benchmark for verifying
CFD turbulence models, meshing strategies, and numerical schemes.

---

## Physics & Theory

### Governing Equations

`simpleFoam` solves the **steady-state Reynolds-Averaged Navier-Stokes (RANS)** equations
for incompressible flow:

```
∇ · U = 0                                          (continuity)
(U · ∇)U = −∇p + ∇ · [(ν + νt) ∇U]               (momentum)
```

where `νt` is the turbulent (eddy) viscosity provided by the turbulence model.

### Steady-State Assumption

For angles of attack below the stall angle (~12–15° for NACA 0012 at this Reynolds
number), the flow is predominantly steady and well-attached. The `steadyState` ddt
scheme drops all time-derivative terms, making the solver iterate toward a converged
steady solution using the SIMPLE pressure-velocity coupling algorithm.

### Reynolds Number

```
Re = U∞ × c / ν = 1.0 × 1.0 / 1×10⁻⁵ = 100,000
```

| Quantity           | Symbol | Value        |
|--------------------|--------|--------------|
| Free-stream speed  | U∞     | 1.0 m/s      |
| Chord length       | c      | 1.0 m        |
| Kinematic viscosity| ν      | 1 × 10⁻⁵ m²/s |
| Reynolds number    | Re     | 100,000      |

At Re = 100,000 the flow is in the **transitional-to-turbulent** regime. A fully
turbulent RANS model (k-epsilon) is used here as a first approximation; see the
[Exercises](#exercises) section for exploring laminar or transition-sensitive models.

### Turbulent Boundary Layer

On the airfoil surface the no-slip condition creates a thin boundary layer where viscous
effects dominate near the wall. The k-epsilon model combined with wall functions
(`nutkWallFunction`, `kqRWallFunction`, `epsilonWallFunction`) bridges the viscous
sublayer, avoiding the need for an extremely fine near-wall mesh.

### Expected Aerodynamic Coefficients

For NACA 0012 at α = 0° and Re ≈ 100,000:

| Coefficient              | Expected Range | Notes                          |
|--------------------------|----------------|--------------------------------|
| Lift coefficient (Cl)    | ≈ 0.0          | Symmetric airfoil at 0° AoA   |
| Drag coefficient (Cd)    | 0.010 – 0.020  | Depends on turbulence model    |
| Moment coefficient (Cm)  | ≈ 0.0          | Symmetric about quarter-chord  |

### Related Notes

- [Turbulence Models](../../notes/06_turbulence_models.md) — k-epsilon theory and wall
  functions
- [Boundary Conditions](../../notes/05_boundary_conditions.md) — inlet/outlet treatment
- [Meshing](../../notes/04_meshing.md) — blockMesh and snappyHexMesh details
- [Linear Solvers](../../notes/09_linear_solvers.md) — PCG, smoothSolver configuration
- [OpenFOAM Dictionaries](../../notes/03_openfoam_dictionaries.md) — controlDict, fvSchemes, fvSolution

---

## Case Structure

```
03_naca_airfoil_analysis/
├── README.md                          ← This file
├── 0/                                 ← Initial & boundary conditions
│   ├── U                              ← Velocity field (fixedValue inlet, noSlip walls)
│   ├── p                              ← Pressure field (fixedValue outlet)
│   ├── k                              ← Turbulent kinetic energy
│   ├── epsilon                        ← Turbulent dissipation rate
│   └── nut                            ← Turbulent viscosity (wall functions)
├── constant/                          ← Physical properties & mesh
│   ├── transportProperties            ← Newtonian fluid, ν = 1×10⁻⁵ m²/s
│   ├── turbulenceProperties           ← RAS → kEpsilon model
│   ├── triSurface/                    ← STL geometry files
│   │   └── airfoil.stl               ← NACA 0012 surface mesh
│   └── polyMesh/                      ← Mesh data (generated)
│       ├── boundary                   ← Patch definitions
│       ├── points                     ← Vertex coordinates
│       ├── faces                      ← Face connectivity
│       ├── owner                      ← Face-to-cell ownership
│       ├── neighbour                  ← Face-to-neighbour-cell
│       ├── cellLevel                  ← snappyHexMesh refinement levels
│       ├── pointLevel                 ← Point refinement levels
│       ├── surfaceIndex               ← Surface intersection data
│       ├── cellZones                  ← Cell zone definitions
│       ├── faceZones                  ← Face zone definitions
│       ├── pointZones                 ← Point zone definitions
│       └── level0Edge                 ← Base-level edge data
├── system/                            ← Solver & discretization settings
│   ├── controlDict                    ← simpleFoam, 1000 iters, Δt = 0.1
│   ├── fvSchemes                      ← steadyState, linearUpwind(U), upwind(k,ε)
│   ├── fvSolution                     ← SIMPLE, relaxation, residual targets
│   ├── blockMeshDict                  ← Background mesh: 60×20×1
│   └── snappyHexMeshDict             ← Surface refinement: level 3
└── mesh_generation_scripts/           ← Geometry & mesh utilities
    ├── generate_naca_0012_airfoil.py  ← Generate airfoil STL (NumPy + PyVista)
    ├── repair_naca_0012_airfoil.py    ← Repair STL (Blender/bmesh)
    ├── mesh_naca_0012_airfoil.py      ← Mesh/smooth STL (Blender/bmesh)
    └── process_naca0012.sh            ← End-to-end pipeline script
```

---

## Mesh Generation Pipeline

The meshing proceeds in two major phases: geometry preparation and OpenFOAM meshing.

### Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     GEOMETRY PREPARATION                                │
│                                                                         │
│  generate_naca_0012_airfoil.py                                         │
│      │  (NumPy + PyVista)                                              │
│      │  NACA 0012 parametric equation → 100-point profile              │
│      │  Extrude to 3D slab (thickness = 0.1 m)                        │
│      ▼                                                                  │
│  naca0012.stl  +  naca0012.vtk                                         │
│      │                                                                  │
│      ▼                                                                  │
│  repair_naca_0012_airfoil.py  (Blender)                                │
│      │  Remove duplicate vertices (tol = 0.0001)                       │
│      │  Fill holes, recalculate normals                                │
│      ▼                                                                  │
│  naca0012_repaired.stl                                                  │
│      │                                                                  │
│      ▼                                                                  │
│  mesh_naca_0012_airfoil.py  (Blender)                                  │
│      │  Triangulate, smooth                                            │
│      ▼                                                                  │
│  naca0012_meshed.stl  →  copy to  constant/triSurface/airfoil.stl      │
└─────────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                     OPENFOAM MESHING                                    │
│                                                                         │
│  blockMesh                                                              │
│      │  Background hex mesh: 60 × 20 × 1 cells                        │
│      │  Domain: (−0.5, −0.5, −0.1) → (2.0, 0.5, 0.1)                 │
│      │  Uniform grading (1 1 1)                                        │
│      ▼                                                                  │
│  snappyHexMesh -overwrite                                              │
│      │  Castellated mesh + snap (no boundary layers)                   │
│      │  Refinement level 3 on airfoil surface                          │
│      │  Max cells: 2,000,000 global / 1,000,000 local                 │
│      │  3 cells between refinement levels                              │
│      │  Location in mesh: (0 0 0) — outside the airfoil               │
│      ▼                                                                  │
│  checkMesh                                                              │
│      │  Verify mesh quality (non-ortho < 65°, skewness, etc.)         │
│      ▼                                                                  │
│  Ready for simpleFoam                                                   │
└─────────────────────────────────────────────────────────────────────────┘
```

### Background Mesh — blockMeshDict

The background mesh is a single hexahedral block:

```
Domain extents:
    x: −0.5  →  2.0   (2.5 m — 0.5c upstream, 1.0c for airfoil, 1.0c downstream)
    y: −0.5  →  0.5   (1.0 m — ±0.5c above/below)
    z: −0.1  →  0.1   (0.2 m — single cell depth for 2D)

Resolution: 60 × 20 × 1  (1,200 background cells)
Grading:    uniform (1 1 1)
```

| Boundary     | Type   | Faces              |
|--------------|--------|--------------------|
| `inlet`      | patch  | x = −0.5 plane     |
| `outlet`     | patch  | x = 2.0 plane      |
| `walls`      | wall   | y = −0.5, y = 0.5  |
| `frontAndBack` | empty | z = −0.1, z = 0.1 |

### Surface Refinement — snappyHexMeshDict

snappyHexMesh reads `constant/triSurface/airfoil.stl` and refines the background mesh:

| Parameter               | Value     | Purpose                                   |
|--------------------------|-----------|-------------------------------------------|
| `castellatedMesh`        | true      | Cut cells around STL surface              |
| `snap`                   | true      | Snap points onto the STL surface          |
| `addLayers`              | false     | No prismatic boundary layers (see note)   |
| Refinement level         | (3 3)     | Min/max refinement on airfoil surface     |
| `nCellsBetweenLevels`    | 3         | Gradual transition between levels         |
| `resolveFeatureAngle`    | 30°       | Refine sharp features                     |
| `maxGlobalCells`         | 2,000,000 | Upper bound on total cell count           |
| `locationInMesh`         | (0 0 0)   | Point outside airfoil (flow region)       |

> **Note:** Boundary layers (`addLayers`) are disabled. For production-quality results,
> enabling boundary layers with appropriate y+ values is recommended (see
> [Exercises](#exercises)).

### Mesh Generation Scripts

| Script                              | Tool     | Purpose                                      |
|--------------------------------------|----------|----------------------------------------------|
| `generate_naca_0012_airfoil.py`     | Python   | NACA 0012 parametric generation → STL/VTK    |
| `repair_naca_0012_airfoil.py`       | Blender  | Remove duplicates, fill holes, fix normals   |
| `mesh_naca_0012_airfoil.py`         | Blender  | Triangulate and smooth the STL surface       |
| `process_naca0012.sh`               | Bash     | Runs all three scripts in sequence           |

> 📖 See also: [Meshing Notes](../../notes/04_meshing.md)

---

## Boundary Conditions

All five field variables must have consistent boundary conditions across every patch.

### Complete Boundary Condition Table

| Boundary       | U                              | p                  | k                      | epsilon                    | nut                     |
|----------------|--------------------------------|--------------------|------------------------|----------------------------|-------------------------|
| **inlet**      | `fixedValue` (1 0 0) m/s      | `zeroGradient`     | `fixedValue` 0.01 m²/s² | `fixedValue` 0.01 m²/s³   | `calculated` 0          |
| **outlet**     | `zeroGradient`                 | `fixedValue` 0     | `zeroGradient`         | `zeroGradient`             | `calculated` 0          |
| **walls**      | `noSlip`                       | `zeroGradient`     | `kqRWallFunction`      | `epsilonWallFunction`      | `nutkWallFunction`      |
| **airfoil**    | `noSlip`                       | `zeroGradient`     | `kqRWallFunction`      | `epsilonWallFunction`      | `nutkWallFunction`      |
| **frontAndBack** | `empty`                      | `empty`            | `empty`                | `empty`                    | `empty`                 |

### Boundary Condition Rationale

- **Inlet:** Fixed uniform velocity drives the flow. Turbulence quantities are prescribed
  at low free-stream levels. Pressure is extrapolated from the interior.
- **Outlet:** Zero-gradient velocity allows the flow to leave freely. Fixed pressure
  (p = 0) provides the reference pressure level.
- **Walls / Airfoil:** No-slip enforces zero velocity at solid surfaces. Wall functions
  bridge the viscous sublayer without requiring extremely fine near-wall cells.
- **frontAndBack:** The `empty` type enforces 2D behaviour by removing the z-direction
  from the solution — OpenFOAM's approach to pseudo-2D simulations.

### Internal Field Initialization

| Field    | Initial Value         | Units   |
|----------|-----------------------|---------|
| U        | uniform (0 0 0)      | m/s     |
| p        | uniform 0            | m²/s²   |
| k        | uniform 0.01         | m²/s²   |
| epsilon  | uniform 0.01         | m²/s³   |
| nut      | uniform 0            | m²/s    |

> 📖 See also: [Boundary Conditions Notes](../../notes/05_boundary_conditions.md)

---

## Turbulence Setup

### Why k-epsilon?

The standard k-epsilon model is chosen for this introductory turbulent case because:

1. **Robustness** — k-epsilon is one of the most numerically stable two-equation RANS
   models and converges reliably for external aerodynamic flows.
2. **Simplicity** — Only two additional transport equations (k and epsilon) are solved,
   making it computationally inexpensive.
3. **Wall function compatibility** — The model pairs naturally with standard wall
   functions, avoiding the need for a very fine near-wall mesh.
4. **Baseline comparison** — It provides a well-understood baseline from which to compare
   more advanced models (k-omega SST, Spalart-Allmaras) in the exercises.

### Model Coefficients (from `turbulenceProperties`)

```
RASModel    kEpsilon

Cmu         0.09        Eddy viscosity coefficient
C1          1.44        Production term coefficient
C2          1.92        Dissipation term coefficient
sigmak      1.0         Turbulent Prandtl number for k
sigmaEps    1.3         Turbulent Prandtl number for epsilon
```

### Inlet Turbulence Values

The inlet values of k and epsilon are set to low free-stream turbulence levels:

```
k       = 0.01 m²/s²
epsilon = 0.01 m²/s³
```

These correspond to:
- **Turbulence intensity:**  I = √(2k/3) / U∞ = √(2×0.01/3) / 1.0 ≈ 8.2 %
- **Turbulent viscosity ratio:**  νt/ν = Cmu × k² / (ε × ν) = 0.09 × 0.0001 / (0.01 × 1e-5) = 90

> **Note:** For a more physically realistic setup, calculate k and epsilon from measured
> or estimated turbulence intensity and turbulent length scale:
> ```
> k = 1.5 × (U∞ × I)²
> epsilon = Cmu^0.75 × k^1.5 / l_t
> ```
> where I is the turbulence intensity (typically 0.01–0.05 for clean wind tunnels) and
> l_t is the turbulent length scale (often ~0.07 × c).

### Wall Function Approach

| Wall Function            | Applied To      | Purpose                                     |
|--------------------------|-----------------|---------------------------------------------|
| `nutkWallFunction`       | nut (airfoil, walls) | Computes νt at first cell using log-law |
| `kqRWallFunction`        | k (airfoil, walls)   | Applies appropriate k BC at wall        |
| `epsilonWallFunction`    | epsilon (airfoil, walls) | Computes ε at first cell from k and y+ |

Wall functions are valid when the first cell centre is in the log-law region
(30 < y+ < 300). For the mesh used here, you should verify y+ with:

```bash
simpleFoam -postProcess -func yPlus
```

> 📖 See also: [Turbulence Models Notes](../../notes/06_turbulence_models.md)

---

## Solver Configuration

### SIMPLE Algorithm

`simpleFoam` uses the **SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations)
algorithm for steady-state pressure-velocity coupling:

```
┌──────────────────────────────────────────────┐
│              SIMPLE Iteration Loop            │
│                                               │
│  1. Solve momentum (U) with current p        │
│  2. Solve pressure correction equation       │
│  3. Correct U with pressure correction       │
│  4. Solve turbulence (k, epsilon)            │
│  5. Check residuals → converged? → STOP      │
│              ↓ (no)                           │
│         Back to step 1                        │
└──────────────────────────────────────────────┘
```

The `consistent yes` option enables the **SIMPLEC** variant, which uses a more accurate
velocity correction and often allows higher relaxation factors.

### Relaxation Factors

Under-relaxation is essential for SIMPLE stability:

| Variable  | Factor | Effect                                                |
|-----------|--------|-------------------------------------------------------|
| p         | 0.3    | Conservative — pressure is the most sensitive field   |
| U         | 0.7    | Standard for steady-state with SIMPLEC                |
| k         | 0.7    | Matches velocity relaxation                           |
| epsilon   | 0.7    | Matches velocity relaxation                           |

> **Why relaxation matters:** Without under-relaxation, the SIMPLE algorithm is
> inherently unstable. The factors limit how much each field changes per iteration,
> preventing divergence at the cost of slower convergence. Values that are too low waste
> iterations; values that are too high cause oscillation or divergence.

### Linear Solvers

| Field    | Solver          | Preconditioner / Smoother | Tolerance | Relative Tol |
|----------|-----------------|---------------------------|-----------|--------------|
| p        | PCG             | DIC                       | 1×10⁻⁶   | 0.05         |
| U        | smoothSolver    | symGaussSeidel            | 1×10⁻⁵   | 0.1          |
| k        | smoothSolver    | symGaussSeidel            | 1×10⁻⁵   | 0.1          |
| epsilon  | smoothSolver    | symGaussSeidel            | 1×10⁻⁵   | 0.1          |

- **PCG** (Preconditioned Conjugate Gradient) is optimal for symmetric pressure systems.
- **smoothSolver** with symmetric Gauss-Seidel is efficient for the non-symmetric
  momentum and turbulence systems.

### Convergence Criteria

The solver stops when **all** residuals drop below their targets:

| Field            | Residual Target |
|------------------|-----------------|
| p                | 1 × 10⁻³       |
| U                | 1 × 10⁻⁴       |
| k                | 1 × 10⁻⁴       |
| epsilon          | 1 × 10⁻⁴       |

Or when the maximum number of iterations (endTime / deltaT = 1000 / 0.1 = **10,000
iterations**) is reached.

### Numerical Schemes (fvSchemes)

| Term               | Scheme                          | Notes                           |
|--------------------|----------------------------------|---------------------------------|
| ddt                | `steadyState`                   | No time derivative              |
| grad               | `Gauss linear`                  | Second-order, central           |
| div(phi,U)         | `Gauss linearUpwind grad(U)`    | Second-order, bounded           |
| div(phi,k)         | `Gauss upwind`                  | First-order, very stable        |
| div(phi,epsilon)   | `Gauss upwind`                  | First-order, very stable        |
| laplacian          | `Gauss linear corrected`        | Second-order with correction    |
| snGrad             | `corrected`                     | Non-orthogonal correction       |
| interpolation      | `linear`                        | Central differencing            |

> Upwind schemes for k and epsilon prioritize stability over accuracy, which is standard
> practice since turbulence quantities can be stiff and prone to unphysical oscillations.

---

## How to Run

### Prerequisites

- OpenFOAM (v2012 or compatible)
- Python 3 with NumPy and PyVista (for geometry generation)
- Blender with Python API (for STL repair/meshing — optional)
- ParaView (for post-processing)

### Step-by-Step

```bash
# Navigate to the case directory
cd projects/03_naca_airfoil_analysis

# ──────────────────────────────────────────────
# STEP 1: Generate airfoil geometry (if needed)
# ──────────────────────────────────────────────
# Option A: Run the full pipeline (requires Blender)
cd mesh_generation_scripts
bash process_naca0012.sh
cp naca0012_meshed.stl ../constant/triSurface/airfoil.stl
cd ..

# Option B: Generate STL only (requires PyVista, no Blender)
cd mesh_generation_scripts
python generate_naca_0012_airfoil.py \
    --stl_filename ../constant/triSurface/airfoil.stl
cd ..

# ──────────────────────────────────────────────
# STEP 2: Generate the computational mesh
# ──────────────────────────────────────────────
blockMesh                    # Create background hex mesh (60×20×1)
snappyHexMesh -overwrite     # Refine around airfoil surface
checkMesh                    # Verify mesh quality

# ──────────────────────────────────────────────
# STEP 3: Run the solver
# ──────────────────────────────────────────────
simpleFoam | tee log.simpleFoam   # Run and save log

# ──────────────────────────────────────────────
# STEP 4: Post-process
# ──────────────────────────────────────────────
simpleFoam -postProcess -func yPlus    # Check y+ values
foamToVTK                              # Convert to VTK format
paraview                               # Visualize results
```

### Monitoring Convergence

While the solver runs, you can monitor residuals in a separate terminal:

```bash
# Plot residuals (if gnuplot is available)
foamMonitor -l postProcessing/residuals/0/residuals.dat

# Or simply watch the log
tail -f log.simpleFoam
```

---

## Expected Results

### Pressure Field

The pressure distribution around the airfoil should show:
- **Stagnation point** at the leading edge (highest pressure, Cp ≈ 1.0)
- **Suction peaks** on both upper and lower surfaces near the leading edge
  (since the airfoil is symmetric at 0° AoA, both peaks are equal)
- **Pressure recovery** toward the trailing edge
- **Wake region** downstream of the trailing edge with slightly reduced pressure

### Velocity Field

- Flow accelerates over the airfoil surface (thinnest streamtubes at max thickness)
- Thin boundary layer on the airfoil surface
- Possible small separation bubble near the trailing edge (depending on mesh quality)
- Undisturbed free-stream velocity far from the airfoil

### Pressure Coefficient Distribution

For NACA 0012 at α = 0° and Re = 100,000:
- Cp should be symmetric about the chord line
- Maximum Cp ≈ 1.0 at the stagnation point
- Minimum Cp ≈ −0.4 to −0.6 near the leading edge suction peak

```
     Cp
  1.0 ┤ ●                             ← Stagnation point
      │  ╲
  0.5 ┤   ╲
      │    ╲
  0.0 ┤─────●───────────────────●──── ← Trailing edge
      │      ╲                 ╱
 -0.5 ┤       ●───────────────●      ← Suction peak (symmetric)
      │
 -1.0 ┤
      └────────────────────────────
      0     0.25    0.5    0.75    1.0
                    x/c
```

### Force Coefficients

At convergence you should observe:
- **Cl ≈ 0** (symmetric airfoil, zero angle of attack)
- **Cd ≈ 0.010 – 0.020** (primarily skin friction drag)
- Residuals below the targets defined in `fvSolution`

---

## Post-Processing Guide

### Opening Results in ParaView

```bash
foamToVTK
paraview
# Open the VTK directory: File → Open → select VTK/...vtk
```

### Pressure Contours

1. Load the case and select the last time step
2. Apply the **Surface** representation
3. Color by `p` (pressure)
4. Use a diverging colormap (blue-white-red) centered at 0
5. The stagnation point (red) and suction regions (blue) should be clearly visible

### Velocity Streamlines

1. Apply **Filters → Stream Tracer**
2. Set seed type to **Line** positioned upstream of the airfoil
3. Adjust the number of seed points (50–100 gives good coverage)
4. Color streamlines by velocity magnitude
5. Observe flow acceleration over the airfoil surface

### Cp Distribution Extraction

1. Apply **Filters → Plot Over Line** or **Slice** on the airfoil surface
2. Calculate Cp: `Cp = (p - p∞) / (0.5 × ρ × U∞²)`
   Since `simpleFoam` uses kinematic pressure (p/ρ): `Cp = p / (0.5 × U∞²) = p / 0.5`
3. Plot Cp vs x/c — compare with NASA experimental data
4. Convention: invert the y-axis so negative Cp (suction) is on top

### Force Coefficients

To compute lift and drag coefficients, add a `forceCoeffs` function object to
`controlDict`:

```c
functions
{
    forceCoeffs
    {
        type            forceCoeffs;
        libs            ("libforces.so");
        writeControl    timeStep;
        writeInterval   1;
        patches         (airfoil);
        rho             rhoInf;
        rhoInf          1.0;
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        CofR            (0.25 0 0);  // Quarter-chord
        pitchAxis       (0 0 1);
        magUInf         1.0;
        lRef            1.0;         // Chord length
        Aref            0.2;         // Chord × span (1.0 × 0.2)
    }
}
```

---

## Exercises

### Exercise 1: Angle of Attack Study

Rotate the inlet velocity to simulate different angles of attack:

```c
// In 0/U, change inlet fixedValue:
// For α = 5°:
inlet
{
    type    fixedValue;
    value   uniform (0.9962 0.0872 0);   // (cos5°, sin5°, 0)
}
```

Run cases at α = 0°, 2°, 5°, 8°, 10°, 12° and plot the lift curve (Cl vs α).
Compare with experimental data from Abbott & Von Doenhoff.

### Exercise 2: Mesh Refinement Study

Test mesh independence by varying snappyHexMesh refinement:

| Case   | Refinement Level | Approx. Cells | Cd      |
|--------|-----------------|---------------|---------|
| Coarse | 2               | ~5,000        | —       |
| Medium | 3 (current)     | ~15,000       | —       |
| Fine   | 4               | ~50,000       | —       |
| V.Fine | 5               | ~200,000      | —       |

Plot Cd vs cell count — the solution should become mesh-independent.

### Exercise 3: Different Turbulence Models

Compare the following models by changing `turbulenceProperties`:

| Model             | RASModel          | Expected Behaviour                     |
|-------------------|-------------------|----------------------------------------|
| k-epsilon         | `kEpsilon`        | Current setup (baseline)               |
| k-omega SST       | `kOmegaSST`       | Better near-wall, requires different BCs |
| Spalart-Allmaras  | `SpalartAllmaras` | Single-equation, good for aero flows   |
| Laminar           | (simulationType laminar) | No turbulence model              |

> **Tip:** When switching to k-omega SST, replace the `epsilon` boundary file with
> `omega` and update wall functions accordingly.

### Exercise 4: Enable Boundary Layers

Edit `snappyHexMeshDict` to add prismatic layers on the airfoil:

```c
addLayers true;

layers
{
    airfoil
    {
        nSurfaceLayers 5;
    }
}
```

Compare results with and without boundary layers — the force coefficients should
become more accurate with proper near-wall resolution.

### Exercise 5: Compare with Experimental Data

Download NACA 0012 experimental Cp data from the NASA Turbulence Modeling Resource:
https://turbmodels.larc.nasa.gov/naca0012_val.html

Overlay your computed Cp distribution with experimental data at matching Re and α.

---

## Comparison: Laminar vs Turbulent

Comparing this project with a laminar case (e.g., the lid-driven cavity):

| Aspect                  | Laminar (Cavity)                | Turbulent (NACA 0012)                   |
|--------------------------|----------------------------------|-----------------------------------------|
| **Solver**              | `icoFoam`                       | `simpleFoam`                            |
| **Time treatment**      | Transient                       | Steady-state                            |
| **Algorithm**           | PISO                            | SIMPLE (SIMPLEC)                        |
| **Turbulence model**    | None (laminar)                  | k-epsilon (RAS)                         |
| **Extra fields**        | U, p                            | U, p, k, epsilon, nut                   |
| **Wall treatment**      | Direct resolution               | Wall functions                          |
| **Mesh complexity**     | Simple blockMesh                | blockMesh + snappyHexMesh               |
| **Geometry**            | Box (no STL)                    | STL surface (external)                  |
| **Reynolds number**     | O(100) – O(1000)                | 100,000                                 |
| **Relaxation factors**  | Not needed (PISO)               | Required (p=0.3, U/k/ε=0.7)            |
| **Boundary conditions** | 2 fields                        | 5 fields with wall functions            |
| **Convergence**         | Time-accurate evolution         | Residual-based iteration                |
| **Difficulty**          | Beginner                        | Intermediate–Advanced                   |

### Key Takeaways

1. **Turbulence modelling adds complexity** — three extra fields (k, epsilon, nut)
   each requiring their own boundary conditions, solvers, and discretization schemes.
2. **Wall functions are a practical compromise** — they avoid the extreme mesh refinement
   needed to resolve the viscous sublayer directly (y+ ≈ 1) but introduce modelling
   assumptions.
3. **Steady-state solvers iterate differently** — SIMPLE uses under-relaxation and
   residual control rather than time-stepping, requiring careful tuning.
4. **External aerodynamics needs snappyHexMesh** — complex geometries cannot be meshed
   with blockMesh alone; the STL-based refinement workflow is essential.

---

## References

1. **Abbott, I.H. & Von Doenhoff, A.E.** (1959). *Theory of Wing Sections*.
   Dover Publications. — The classic reference for NACA airfoil data including extensive
   wind tunnel results for the 0012 profile.

2. **NASA Turbulence Modeling Resource** — NACA 0012 validation case:
   https://turbmodels.larc.nasa.gov/naca0012_val.html

3. **OpenFOAM User Guide** — simpleFoam, snappyHexMesh, and turbulence modelling:
   https://www.openfoam.com/documentation/user-guide

4. **Ladson, C.L.** (1988). *Effects of Independent Variation of Mach and Reynolds
   Numbers on the Low-Speed Aerodynamic Characteristics of the NACA 0012 Airfoil
   Section*. NASA TM-4074. — Comprehensive experimental Cp data.

5. **Gregory, N. & O'Reilly, C.L.** (1970). *Low-Speed Aerodynamic Characteristics of
   NACA 0012 Aerofoil Section*. ARC R&M 3726. — Classic UK wind tunnel data.

6. **Versteeg, H.K. & Malalasekera, W.** (2007). *An Introduction to Computational
   Fluid Dynamics: The Finite Volume Method*. Pearson. — Excellent introduction to
   SIMPLE, turbulence models, and discretization schemes.

### Course Notes

- [Short Intro to CFD](../../notes/01_short_intro_to_cfd.md)
- [OpenFOAM Cases](../../notes/02_openfoam_cases.md)
- [OpenFOAM Dictionaries](../../notes/03_openfoam_dictionaries.md)
- [Meshing](../../notes/04_meshing.md)
- [Boundary Conditions](../../notes/05_boundary_conditions.md)
- [Turbulence Models](../../notes/06_turbulence_models.md)
- [Parallelization](../../notes/07_parallelization.md)
- [CFL Number](../../notes/08_cfl_number.md)
- [Linear Solvers](../../notes/09_linear_solvers.md)
