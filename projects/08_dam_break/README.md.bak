# 08 — Dam Break Simulation 🌊

| Property       | Value                                            |
|----------------|--------------------------------------------------|
| **Difficulty** | ⭐⭐⭐ Intermediate                              |
| **Solver**     | `interFoam` (transient, incompressible, VOF)     |
| **Physics**    | Two-phase flow (water + air), free surface        |
| **Turbulence** | Laminar                                          |
| **Mesh**       | blockMesh (structured, 2D)                       |
| **Key Concept**| Volume of Fluid (VOF) method                     |
| **Run Time**   | ~5–15 minutes (depending on hardware)            |

---

## Table of Contents

1. [Problem Description](#1-problem-description)
2. [Physics and Theory](#2-physics-and-theory)
3. [The Volume of Fluid (VOF) Method](#3-the-volume-of-fluid-vof-method)
4. [Case Structure](#4-case-structure)
5. [Mesh Generation](#5-mesh-generation)
6. [Initial Conditions and setFields](#6-initial-conditions-and-setfields)
7. [Boundary Conditions](#7-boundary-conditions)
8. [Solver Configuration](#8-solver-configuration)
9. [Numerical Schemes](#9-numerical-schemes)
10. [How to Run](#10-how-to-run)
11. [Expected Results](#11-expected-results)
12. [Post-Processing](#12-post-processing)
13. [Exercises](#13-exercises)
14. [References](#14-references)

---

## 1. Problem Description

The **dam break** problem is one of the most classic benchmarks in computational fluid
dynamics for testing free-surface and multiphase flow solvers. A column of water, initially
held in place at one end of a rectangular tank, is instantaneously released and collapses
under gravity. The water surges across the tank floor, impacts the far wall, and generates
complex splashing, wave reflections, and air entrainment.

This tutorial teaches you:

- The **Volume of Fluid (VOF)** method for tracking fluid interfaces
- **Multiphase flow** simulation with two immiscible fluids (water and air)
- **Transient free-surface dynamics** with large interface deformations
- **Adjustable time stepping** controlled by the Courant number
- The `setFields` utility for initializing non-uniform fields

### Problem Setup

```
                    atmosphere (open top boundary)
    ┌─────────────────────────────────────────────────────────┐
    │                                                         │
    │              AIR  (alpha.water = 0)                      │
    │              ρ = 1 kg/m³                                 │
    │              ν = 1.48×10⁻⁵ m²/s                         │
    │                                                         │
    │   0.292 m                                               │ 0.584 m
    │   ┌──────┐                                              │
    │   │██████│                                              │
    │   │██████│  WATER  (alpha.water = 1)                     │
    │   │██████│  ρ = 1000 kg/m³                               │
    │   │██████│  ν = 1×10⁻⁶ m²/s                             │
    │   │██████│                                              │
    └───┴──────┴──────────────────────────────────────────────┘
    ↑   0.146 m                                               ↑
 leftWall              lowerWall (no-slip)               rightWall
                         0.584 m
```

### Physical Parameters

| Parameter          | Value                    | Units   |
|--------------------|--------------------------|---------|
| Domain width       | 0.584                    | m       |
| Domain height      | 0.584                    | m       |
| Domain depth       | 0.01                     | m       |
| Water column width | 0.146                    | m       |
| Water column height| 0.292                    | m       |
| Water density      | 1000                     | kg/m³   |
| Air density        | 1                        | kg/m³   |
| Water viscosity    | 1 × 10⁻⁶                | m²/s    |
| Air viscosity      | 1.48 × 10⁻⁵             | m²/s    |
| Surface tension    | 0.07                     | N/m     |
| Gravity            | 9.81 (downward)          | m/s²    |

### Time Evolution

The simulation captures the following stages:

```
  t = 0.0 s              t = 0.2 s              t = 0.4 s
  (Initial state)        (Column collapse)       (Surge reaches wall)

  ┌──┬──────────────┐    ┌─────────────────┐    ┌─────────────────┐
  │  │              │    │                 │    │           ░░░░  │
  │  │              │    │                 │    │          ░░░░░  │
  │██│              │    │ ~~              │    │        ░░░░░░░  │
  │██│              │    │ ████            │    │  ~~~ ░░░░░░░░░  │
  │██│              │    │ ██████          │    │  ████████░░░░░  │
  │██│              │    │ ████████        │    │  ███████████░░  │
  │██│              │    │ ██████████      │    │  █████████████  │
  └──┴──────────────┘    └─────────────────┘    └─────────────────┘


  t = 0.6 s              t = 0.8 s              t = 1.0 s
  (Wall impact)          (Splash-up)            (Wave reflection)

  ┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
  │              ░░ │    │            ░ ██ │    │  ░░         ░░  │
  │             ░░░ │    │            ░░██ │    │  ░░░       ░░░  │
  │            ░░░░ │    │           ░░░██ │    │  ░░░░     ░░░░  │
  │  ~~       ░░░░░ │    │  ~~      ░░░░██ │    │  ░░░░░  ░░░░░  │
  │  ████████░░░░░░ │    │  ██████░░░░░░██ │    │  ░░░░░░░░░░░░  │
  │  ██████████████ │    │  ██████████████ │    │  ████████████░  │
  │  ██████████████ │    │  ██████████████ │    │  █████████████  │
  └─────────────────┘    └─────────────────┘    └─────────────────┘

  Legend:  ██ = water (alpha.water ≈ 1)
           ░░ = interface region (0 < alpha.water < 1)
           ~~ = free surface
```

---

## 2. Physics and Theory

### Governing Equations

The dam break problem involves two immiscible, incompressible Newtonian fluids
(water and air) sharing a computational domain. The flow is governed by the
**single-fluid Navier-Stokes equations** where fluid properties vary according to
the local volume fraction:

**Continuity equation:**

```
∇ · U = 0
```

**Momentum equation:**

```
∂(ρU)/∂t + ∇ · (ρUU) = -∇p_rgh - g · x ∇ρ + ∇ · (μ ∇U) + σ κ ∇α
```

where:
- `ρ` — mixture density (varies with α)
- `U` — velocity field (shared by both phases)
- `p_rgh` — dynamic pressure (p - ρ g · x), removes the hydrostatic component
- `g` — gravitational acceleration vector (0, -9.81, 0) m/s²
- `μ` — mixture dynamic viscosity
- `σ` — surface tension coefficient
- `κ` — interface curvature
- `α` — volume fraction (alpha.water)

### Mixture Properties

The fluid properties at any point in the domain are computed as weighted averages
of the two phases using the volume fraction α:

```
ρ = α · ρ_water + (1 - α) · ρ_air
  = α · 1000 + (1 - α) · 1

μ = α · μ_water + (1 - α) · μ_air
```

### Dynamic Pressure (p_rgh)

OpenFOAM's multiphase solvers use the **modified pressure** `p_rgh` instead of
the absolute pressure `p`. This is defined as:

```
p_rgh = p - ρ g · x
```

This formulation removes the hydrostatic pressure gradient from the momentum
equation, improving numerical accuracy — especially for flows dominated by gravity
where the hydrostatic pressure is much larger than the dynamic pressure variations.

### Why Laminar?

This small-scale dam break (0.584 m domain) operates at relatively low Reynolds
numbers during most of the simulation. While some turbulent features may develop
during the splash phase, a laminar simulation captures the essential physics of
the collapse and wave propagation. This keeps the case simple and focused on the
VOF method itself. For larger-scale dam breaks, turbulence modelling (e.g., k-ε
or LES) would be appropriate.

See also: [Turbulence Models](../../notes/06_turbulence_models.md)

---

## 3. The Volume of Fluid (VOF) Method

### Overview

The **Volume of Fluid** method, introduced by Hirt and Nichols (1981), is the
most widely used technique for simulating flows with sharp interfaces between
immiscible fluids. The key idea is elegantly simple:

> Track the **volume fraction** of each fluid in every cell, rather than
> explicitly tracking the interface geometry.

### The Alpha Field

The volume fraction field `alpha.water` (denoted α) represents the fraction of
each computational cell occupied by water:

```
  α = 1          α ≈ 0.5         α = 0
  ┌──────────┐   ┌──────────┐   ┌──────────┐
  │██████████│   │██████    │   │          │
  │██████████│   │██████    │   │          │
  │██████████│   │████      │   │          │
  │██████████│   │████      │   │          │
  │██████████│   │██        │   │          │
  └──────────┘   └──────────┘   └──────────┘
   Pure water     Interface       Pure air
                  cell
```

| Value of α | Meaning                     |
|------------|-----------------------------|
| α = 1      | Cell is fully water         |
| α = 0      | Cell is fully air           |
| 0 < α < 1  | Cell contains the interface |

### Transport Equation for Alpha

The volume fraction is advected by the flow velocity using a transport equation:

```
∂α/∂t + ∇ · (U α) + ∇ · (U_r α(1-α)) = 0
```

The last term `∇ · (U_r α(1-α))` is the **interface compression** term, unique
to OpenFOAM's VOF implementation. Here `U_r` is an artificial compression velocity
that acts only at the interface (where `α(1-α) ≠ 0`) to maintain a sharp interface
and prevent numerical diffusion from smearing it out.

The `cAlpha` parameter in `fvSolution` controls the strength of this compression:

| cAlpha | Effect                                           |
|--------|--------------------------------------------------|
| 0      | No compression (interface will smear)            |
| 1      | Conservative compression (recommended)           |
| > 1    | Enhanced compression (sharper but may oscillate) |

### MULES: Maintaining Boundedness

A critical requirement for VOF is that alpha must remain bounded between 0 and 1.
OpenFOAM uses the **MULES** (Multidimensional Universal Limiter with Explicit
Solution) algorithm to ensure this boundedness. MULES applies flux limiting to the
alpha transport equation, preventing overshoots and undershoots that would create
non-physical negative volume fractions or fractions exceeding unity.

```
  Without MULES:                  With MULES:
  α can exceed bounds             α stays bounded

  1.2 ─ ─ ─ ─ ─ ─ ─ ─            1.0 ─────────────────
                  ╱ ╲                           ┌─────
  1.0 ──────────╱   ╲            0.8           ╱
               ╱     ╲           0.6          ╱
  0.5         ╱       ╲          0.4         ╱
             ╱         ╲         0.2        ╱
  0.0 ──────╱           ╲────   0.0 ───────┘
       ‐0.1                          (bounded 0 to 1)
```

### Surface Tension: The Continuum Surface Force (CSF) Model

Surface tension at the interface is modelled using the **Brackbill CSF method**
(Brackbill et al., 1992). Instead of applying a force at the discrete interface,
the surface tension is converted to a volume force distributed across interface
cells:

```
F_σ = σ κ ∇α
```

where:
- `σ` = surface tension coefficient (0.07 N/m for water-air)
- `κ` = interface curvature, computed as: `κ = -∇ · n̂` where `n̂ = ∇α / |∇α|`
- `∇α` = gradient of volume fraction (non-zero only at the interface)

```
  Surface tension force at the interface:

  Air  (α = 0)
  ─────────────────────────────────
                                        F_σ acts normal to interface,
    ← ← ← F_σ → → →                   pointing toward center of
                                        curvature
  ─────────────────────────────────
  Water (α = 1)
```

For this dam break case, surface tension is relatively minor compared to gravity
and inertia forces, but it is included for physical completeness and helps
maintain a sharp interface.

---

## 4. Case Structure

```
08_dam_break/
├── 0/                              # Initial conditions (time = 0)
│   ├── U                           #   Velocity field
│   ├── p_rgh                       #   Dynamic pressure
│   └── alpha.water                 #   Volume fraction (water phase)
├── constant/                       # Physical properties and mesh
│   ├── g                           #   Gravitational acceleration
│   ├── transportProperties         #   Phase properties (water, air)
│   ├── turbulenceProperties        #   Turbulence model (laminar)
│   └── polyMesh/                   #   Mesh (generated by blockMesh)
├── system/                         # Solver and discretization settings
│   ├── blockMeshDict               #   Mesh generation dictionary
│   ├── controlDict                 #   Run control parameters
│   ├── fvSchemes                   #   Discretization schemes
│   ├── fvSolution                  #   Linear solver settings
│   └── setFieldsDict              #   Field initialization (water column)
├── Allrun                          # Run script
├── Allclean                        # Cleanup script
└── README.md                       # This file
```

### Key Files Explained

| File                   | Purpose                                                      |
|------------------------|--------------------------------------------------------------|
| `0/alpha.water`        | Volume fraction — **the defining field** of VOF simulations |
| `0/p_rgh`              | Dynamic pressure (not absolute pressure)                     |
| `0/U`                  | Shared velocity field for both phases                        |
| `constant/g`           | Gravity vector — drives the entire flow                      |
| `constant/transportProperties` | Density, viscosity for each phase; surface tension  |
| `system/setFieldsDict` | Defines where water starts (the "dam")                      |
| `system/controlDict`   | Adjustable time stepping with CFL control                   |

See also: [OpenFOAM Case Structure](../../notes/02_openfoam_cases.md),
[OpenFOAM Dictionaries](../../notes/03_openfoam_dictionaries.md)

---

## 5. Mesh Generation

The mesh is generated with `blockMesh` as a simple rectangular domain subdivided
into a uniform structured grid.

### Domain and Resolution

```
                      0.584 m
    ┌─────────────────────────────────────────┐
    │                                         │
    │     46 cells across (Δx ≈ 0.0127 m)    │
    │                                         │
    │     46 cells high   (Δy ≈ 0.0127 m)    │  0.584 m
    │                                         │
    │     1 cell deep     (Δz = 0.01 m)       │
    │                                         │
    └─────────────────────────────────────────┘
    Total cells: 46 × 46 × 1 = 2,116

    Cell aspect ratio: ~1.0 (square cells — ideal)
```

### Vertex Numbering

The block is defined by 8 vertices forming a single hexahedral block:

```
           7 ─────────── 6          z
          ╱│            ╱│          │  y
         ╱ │           ╱ │          │ ╱
        4 ─────────── 5  │          │╱
        │  │          │  │          └───── x
        │  3 ─────────│─ 2
        │ ╱           │ ╱
        │╱            │╱
        0 ─────────── 1

    Back face  (z=0):    0, 1, 2, 3
    Front face (z=0.01): 4, 5, 6, 7
```

### Boundary Patches

| Patch         | Type   | Face Vertices    | Physical Meaning               |
|---------------|--------|------------------|--------------------------------|
| `leftWall`    | wall   | (0 4 7 3)        | Left solid wall                |
| `rightWall`   | wall   | (1 2 6 5)        | Right solid wall               |
| `lowerWall`   | wall   | (0 1 5 4)        | Bottom solid wall              |
| `atmosphere`  | patch  | (3 7 6 2)        | Open top (atmospheric)         |
| `frontAndBack`| empty  | front & back     | 2D constraint (no z-direction) |

The `empty` type on `frontAndBack` tells OpenFOAM this is a 2D simulation — no
equations are solved in the z-direction.

See also: [Meshing](../../notes/04_meshing.md)

---

## 6. Initial Conditions and setFields

### The Problem: Non-Uniform Initial Fields

Unlike most OpenFOAM simulations where initial fields are uniform, the dam break
requires **spatially varying initial conditions** — the water column occupies only
a portion of the domain. The `0/alpha.water` file sets the default (air everywhere):

```
internalField   uniform 0;   // Everything starts as air
```

### The setFields Utility

The `setFields` utility modifies field values within specified geometric regions
**after** the mesh is created. In our case, it sets `alpha.water = 1` inside a
box corresponding to the water column:

```
  setFieldsDict defines:

  ┌──────────────────────────────────────────────────────────┐
  │                                                          │
  │   alpha.water = 0 (air) everywhere by default            │
  │                                                          │
  │   ┌────────┐                                             │
  │   │████████│ ← boxToCell region:                         │
  │   │████████│   box (0, 0, -1) to (0.1461, 0.292, 1)     │
  │   │████████│   alpha.water = 1 (water)                   │
  │   │████████│                                             │
  │   └────────┘                                             │
  └──────────────────────────────────────────────────────────┘
```

The z-bounds of (-1, 1) extend well beyond the mesh depth (0 to 0.01 m) to
ensure all cells in the z-direction are captured.

### How setFields Works

1. Reads the mesh from `constant/polyMesh/`
2. Reads `system/setFieldsDict`
3. For each region, finds cells whose **centre** falls within the specified geometry
4. Overwrites the field values for those cells
5. Writes the modified field back to the `0/` directory

> **Important:** `setFields` modifies the field files in place. You must run
> `blockMesh` before `setFields` (the mesh must exist), and if you re-run
> `setFields`, it will overwrite the previous initialization.

---

## 7. Boundary Conditions

### Velocity (U)

| Patch          | BC Type                          | Value           | Meaning                                        |
|----------------|----------------------------------|-----------------|-------------------------------------------------|
| `leftWall`     | `noSlip`                         | —               | Zero velocity at wall                           |
| `rightWall`    | `noSlip`                         | —               | Zero velocity at wall                           |
| `lowerWall`    | `noSlip`                         | —               | Zero velocity at wall                           |
| `atmosphere`   | `pressureInletOutletVelocity`    | uniform (0 0 0) | Allows free outflow; applies zero gradient for outgoing flow |
| `defaultFaces` | `noSlip`                         | —               | Catch-all for unmeshed faces                    |
| `frontAndBack` | `empty`                          | —               | 2D constraint                                   |

### Pressure (p_rgh)

| Patch          | BC Type              | Value        | Meaning                                 |
|----------------|----------------------|--------------|-----------------------------------------|
| `leftWall`     | `fixedFluxPressure`  | uniform 0    | Adjusts pressure gradient to match velocity BC |
| `rightWall`    | `fixedFluxPressure`  | uniform 0    | Same as above                           |
| `lowerWall`    | `fixedFluxPressure`  | uniform 0    | Same as above                           |
| `atmosphere`   | `totalPressure`      | p0 = 0       | Atmospheric reference pressure          |
| `defaultFaces` | `fixedFluxPressure`  | uniform 0    | Consistent with wall velocity           |
| `frontAndBack` | `empty`              | —            | 2D constraint                           |

### Volume Fraction (alpha.water)

| Patch          | BC Type        | Value              | Meaning                                |
|----------------|----------------|--------------------|----------------------------------------|
| `leftWall`     | `zeroGradient` | —                  | No flux of alpha through wall          |
| `rightWall`    | `zeroGradient` | —                  | No flux of alpha through wall          |
| `lowerWall`    | `zeroGradient` | —                  | No flux of alpha through wall          |
| `atmosphere`   | `inletOutlet`  | inletValue = 0     | Air enters if flow reverses inward     |
| `defaultFaces` | `zeroGradient` | —                  | Default wall treatment                 |
| `frontAndBack` | `empty`        | —                  | 2D constraint                          |

### Why fixedFluxPressure?

The `fixedFluxPressure` BC is specifically designed for multiphase flows with
gravity. It adjusts the pressure gradient at walls to be consistent with the
velocity BC while correctly accounting for the body force (gravity) term. This
prevents spurious velocities ("parasitic currents") near walls in
gravity-dominated multiphase flows.

See also: [Boundary Conditions](../../notes/05_boundary_conditions.md)

---

## 8. Solver Configuration

### controlDict — Time Control

The dam break uses **adjustable time stepping**, which is critical for VOF
simulations:

```
deltaT          0.001;          // Initial time step
adjustTimeStep  yes;            // Let OpenFOAM adjust dt
maxCo           0.5;            // Maximum Courant number
maxAlphaCo      0.5;            // Maximum Courant for alpha transport
maxDeltaT       0.01;           // Upper limit for time step
endTime         1;              // Simulate 1 second
writeInterval   0.05;           // Write output every 0.05 s (20 frames)
```

### Why Adjustable Time Stepping?

The dam break involves dramatically different velocities at different times:

```
  Velocity timeline:

  Speed
  (m/s)
    3 │            ╱╲     Impact!
      │           ╱  ╲
    2 │          ╱    ╲ ╱╲  Reflection
      │         ╱      ╳  ╲
    1 │    ╱───╱      ╱ ╲  ╲───
      │   ╱                    ──
    0 │──╱
      └──┬───┬───┬───┬───┬───┬──→ Time (s)
         0  0.2 0.4 0.6 0.8 1.0
         ↑              ↑
       Collapse      Splash
       begins        impact
```

A fixed time step small enough for the splash phase would waste computation
during the slow initial collapse. Adjustable time stepping adapts automatically.

### The Alpha Courant Number

For VOF simulations, there are **two** Courant numbers that matter:

```
Standard Courant:     Co = |U| × Δt / Δx

Alpha Courant:        Co_α = |U_interface| × Δt / Δx
```

The alpha Courant number is based on the velocity at the interface between fluids.
Because the interface advection equation requires high accuracy, `maxAlphaCo` is
often set equal to or lower than `maxCo`. In our case, both are 0.5.

See also: [CFL Number](../../notes/08_cfl_number.md)

### fvSolution — PIMPLE Algorithm

```
PIMPLE
{
    momentumPredictor   no;
    nOuterCorrectors    1;      // 1 = PISO mode
    nCorrectors         3;      // Pressure correction loops
    nNonOrthogonalCorrectors 0; // Not needed for orthogonal mesh
}
```

With `nOuterCorrectors = 1`, PIMPLE degenerates to the **PISO** algorithm, which
is standard for transient VOF simulations. The 3 corrector steps help ensure
pressure-velocity coupling convergence.

### Linear Solvers

| Field         | Solver        | Preconditioner  | Tolerance | Purpose                  |
|---------------|---------------|-----------------|-----------|--------------------------|
| alpha.water   | smoothSolver  | symGaussSeidel  | 1e-8      | Volume fraction (tight!) |
| pcorr         | PCG           | DIC             | 1e-5      | Flux correction          |
| p_rgh         | PCG           | DIC             | 1e-7      | Pressure (most expensive)|
| U             | smoothSolver  | symGaussSeidel  | 1e-6      | Velocity                 |

The alpha solver uses a very tight tolerance (1e-8) because small errors in the
volume fraction field can lead to mass conservation issues and interface smearing.

See also: [Linear Solvers](../../notes/09_linear_solvers.md)

---

## 9. Numerical Schemes

### fvSchemes — Discretization

| Category               | Scheme                        | Reason                                           |
|------------------------|-------------------------------|--------------------------------------------------|
| `ddtSchemes`           | `Euler`                       | First-order implicit time; robust for VOF        |
| `gradSchemes`          | `Gauss linear`                | Second-order gradient (standard)                 |
| `div(rhoPhi,U)`        | `Gauss linearUpwind grad(U)`  | Upwind-biased for momentum; stable & accurate    |
| `div(phi,alpha)`       | `Gauss vanLeer`               | TVD limiter for alpha; prevents overshoots       |
| `div(phirb,alpha)`     | `Gauss linear`                | Interface compression term                       |
| `laplacianSchemes`     | `Gauss linear corrected`      | Second-order diffusion with non-orthogonal correction |
| `snGradSchemes`        | `corrected`                   | Consistent with laplacian scheme                 |

### Why vanLeer for Alpha?

The volume fraction transport requires special care. Standard schemes can cause:

- **Upwind:** Excessive numerical diffusion → interface smears
- **Central/linear:** Oscillations → alpha goes below 0 or above 1
- **vanLeer (TVD):** Bounded, sharp interface with minimal diffusion

```
  Interface profile comparison:

  α
  1.0 ─────╮            ╭─── Exact (step function)
            │            │
            │╲          ╱│
  0.5       │ ╲   ╱╲  ╱ │    ── Central (oscillates!)
            │  ╲─╱  ╲╱  │    ── Upwind  (too diffuse)
            │    ────    │    ── vanLeer (good compromise)
  0.0 ──────╯            ╰───
         interface location
```

The `div(phirb,alpha)` term handles the interface compression velocity `U_r`,
which is OpenFOAM's approach to maintaining a sharp interface without geometric
reconstruction.

---

## 10. How to Run

### Prerequisites

- OpenFOAM 6+ installed and sourced (`source /opt/openfoam6/etc/bashrc` or similar)
- `blockMesh`, `setFields`, and `interFoam` available in your `$PATH`

### Quick Start

```bash
cd projects/08_dam_break
./Allrun
```

### Step-by-Step

```bash
# 1. Generate the mesh
blockMesh

# 2. Initialize the water column
setFields

# 3. Run the simulation
interFoam

# 4. Visualize results
paraFoam
```

### What Each Step Does

```
  blockMesh:                  setFields:                  interFoam:
  Creates mesh                Sets alpha.water=1          Solves the
  in polyMesh/                in water column             multiphase flow

  ┌─────────────────┐        ┌───┬─────────────┐        ┌─────────────────┐
  │ │ │ │ │ │ │ │ │ │        │███│             │        │                 │
  │─┼─┼─┼─┼─┼─┼─┼─│        │███│             │        │ ~~              │
  │ │ │ │ │ │ │ │ │ │        │███│             │        │ ████~           │
  │─┼─┼─┼─┼─┼─┼─┼─│        │███│             │        │ ████████~       │
  │ │ │ │ │ │ │ │ │ │        │███│             │        │ ██████████████  │
  └─┴─┴─┴─┴─┴─┴─┴─┘        └───┴─────────────┘        └─────────────────┘
  Empty grid                  Water initialized           Flow solved!
```

### Monitoring Progress

While `interFoam` runs, you can monitor:

```bash
# Watch the Courant number and time step in real-time
tail -f log.interFoam | grep "Courant"

# Check which time step the solver is at
tail -f log.interFoam | grep "^Time ="

# Monitor interface sharpness (alpha should stay 0-1)
tail -f log.interFoam | grep "alpha.water"
```

### Cleaning Up

To remove all generated files and start fresh:

```bash
./Allclean
```

---

## 11. Expected Results

### Console Output

During the run, you should see output like:

```
Time = 0.001

Courant Number mean: 0.000345 max: 0.123
Interface Courant Number mean: 0.000123 max: 0.089
deltaT = 0.00108
MULES: Solving for alpha.water
Phase-1 volume fraction = 0.0625  Min(alpha.water) = 0  Max(alpha.water) = 1
...
```

### Key Things to Verify

| Check                    | Expected Value        | Concern if Wrong                    |
|--------------------------|-----------------------|-------------------------------------|
| Min(alpha.water)         | 0 (exactly)           | Negative values = boundedness issue |
| Max(alpha.water)         | 1 (exactly)           | > 1 = MULES not working            |
| Phase volume fraction    | ~0.0625 (constant)    | Changing = mass not conserved       |
| Courant number           | ≤ 0.5                 | Higher = may need smaller maxCo     |
| Time step                | 0.001–0.01            | Very small = possible mesh issue    |

### Physical Behavior Timeline

| Time (s) | Event                              | Velocity Scale |
|----------|-------------------------------------|----------------|
| 0.0      | Water column released               | 0 m/s          |
| 0.0–0.2  | Column collapses, base surges right | ~1 m/s         |
| 0.2–0.4  | Water front crosses the tank floor  | ~2 m/s         |
| 0.4–0.6  | Front reaches right wall, climbs up | ~2–3 m/s       |
| 0.6–0.8  | Splash impact, water jets upward    | ~3 m/s (peak)  |
| 0.8–1.0  | Wave reflects, sloshing begins      | ~1–2 m/s       |

---

## 12. Post-Processing

### ParaView Visualization

```bash
# Open the case in ParaView
paraFoam
# Or, create a .foam file and open in ParaView
touch case.foam
paraview case.foam
```

### Recommended Visualizations

#### 1. Alpha Contour Animation (Primary Result)

- Select `alpha.water` as the display variable
- Set colour range: 0 to 1
- Use a blue-white-red or blue-red colour map
- Press Play to animate the dam break

#### 2. Velocity Magnitude

- Display `mag(U)` to see how the water accelerates
- Overlay velocity vectors (Glyph filter) to visualize flow direction
- Notice the high velocities at the wave front and during wall impact

#### 3. Pressure Field

- Display `p_rgh` to see the dynamic pressure distribution
- Observe the pressure spike during wall impact
- Note the hydrostatic-like distribution in the water at rest

#### 4. Interface Sharpness

- Plot alpha along a horizontal line (Plot Over Line filter)
- A sharp interface shows a step from 0 to 1 over 1–3 cells
- A smeared interface indicates excessive numerical diffusion

### Quantitative Post-Processing

#### Wave Front Position vs Time

Extract the leading edge of the water front over time:

```bash
# Create a postProcessing directory with sample data
# Add to controlDict under functions:
#
# functions
# {
#     frontPosition
#     {
#         type            surfaces;
#         libs            ("libsampling.so");
#         writeControl    writeTime;
#         surfaceFormat   raw;
#         fields          ( alpha.water );
#         surfaces
#         (
#             yNormal
#             {
#                 type        cuttingPlane;
#                 planeType   pointAndNormal;
#                 point       (0 0.001 0.005);
#                 normal      (0 1 0);
#             }
#         );
#     }
# }
```

#### Validation Against Experiment

The dam break problem has been experimentally studied by Martin and Moyce (1952).
Their non-dimensional results for the wave front position are:

```
  Non-dimensional front position Z/Z₀ vs non-dimensional time T*:

  Z/Z₀
  5.0 ┤                                           ○
      │                                       ○
  4.0 ┤                                   ○
      │                              ○ ○
  3.0 ┤                          ○
      │                     ○ ○
  2.0 ┤               ○ ○
      │          ○ ○
  1.0 ┤─ ○ ○ ○
      │
  0.0 ┤
      └───┬───┬───┬───┬───┬───┬───┬───┬──
          0  0.5 1.0 1.5 2.0 2.5 3.0 3.5    T*

  ○ = Martin & Moyce (1952) experimental data

  Where:
    Z₀ = initial column width (0.146 m)
    T* = t × √(2g / Z₀)
```

---

## 13. Exercises

### Exercise 1: Add an Obstacle

Add a small rectangular obstacle on the tank floor to study wave-structure interaction.
Modify `blockMeshDict` to include an obstacle at x = 0.292 m with dimensions
0.024 × 0.048 m. This requires splitting the domain into multiple blocks.

```
  ┌──────────────────────────────────────────────────────┐
  │                                                      │
  │   WATER  │     AIR                                   │
  │   COLUMN │                                           │
  │          │              ┌────┐                        │
  │          │              │ OB │  ← 0.024 × 0.048 m    │
  └──────────┴──────────────┴────┴────────────────────────┘
```

### Exercise 2: Change Water Column Size

Experiment with different initial water column aspect ratios:

| Case | Width (m) | Height (m) | Aspect Ratio | Expected Behavior            |
|------|-----------|------------|--------------|------------------------------|
| A    | 0.146     | 0.292      | 1:2          | Standard (this tutorial)     |
| B    | 0.146     | 0.146      | 1:1          | Shorter column, slower flow  |
| C    | 0.292     | 0.292      | 1:1          | More water, stronger impact  |
| D    | 0.146     | 0.584      | 1:4          | Tall column, violent splash  |

Modify the `setFieldsDict` box coordinates for each case and compare results.

### Exercise 3: Mesh Refinement Study

Test the sensitivity of results to mesh resolution:

| Mesh   | Cells    | Δx (m)   | Expected Outcome                       |
|--------|----------|----------|----------------------------------------|
| Coarse | 23×23    | 0.0254   | Fast but diffused interface            |
| Medium | 46×46    | 0.0127   | Good balance (this tutorial)           |
| Fine   | 92×92    | 0.0063   | Sharper interface, slower              |
| VFine  | 184×184  | 0.0032   | Detailed splashing features visible    |

Modify the cell counts in `blockMeshDict` and compare interface sharpness.

### Exercise 4: 3D Extension

Convert this 2D case to a full 3D simulation:

1. Change `frontAndBack` from `empty` to `wall` in `blockMeshDict`
2. Increase z-depth and z-cells: e.g., 0.146 m with 12 cells
3. Change `frontAndBack` boundary conditions from `empty` to `noSlip`/`fixedFluxPressure`/`zeroGradient`
4. Consider using parallel decomposition for the larger mesh

See also: [Parallelization](../../notes/07_parallelization.md)

### Exercise 5: Turbulence Model

For higher Reynolds number dam breaks, add turbulence modelling:

1. Change `turbulenceProperties` to use `RAS` with `kOmegaSST`
2. Add `0/k`, `0/omega`, and `0/nut` initial condition files
3. Set appropriate wall functions on boundaries

See also: [Turbulence Models](../../notes/06_turbulence_models.md)

### Exercise 6: Surface Tension Effects

Study the effect of surface tension by varying the `sigma` value in
`transportProperties`:

| σ (N/m) | Physical Meaning                | Effect on Flow                       |
|---------|----------------------------------|--------------------------------------|
| 0.0     | No surface tension               | Maximum interface fragmentation      |
| 0.07    | Water-air (standard)             | Slight interface smoothing           |
| 0.7     | 10× enhanced                     | Noticeable rounding of interface     |
| 7.0     | 100× enhanced                    | Strong surface tension dominance     |

---

## 14. References

### Experimental Validation

1. **Martin, J.C. and Moyce, W.J.** (1952). "Part IV. An Experimental Study of
   the Collapse of Liquid Columns on a Rigid Horizontal Plane."
   *Philosophical Transactions of the Royal Society of London. Series A,
   Mathematical and Physical Sciences*, 244(882), pp. 312–324.
   - The classic experimental benchmark for dam break validation
   - Provides non-dimensional front position and column height data

### VOF Method

2. **Hirt, C.W. and Nichols, B.D.** (1981). "Volume of Fluid (VOF) Method for
   the Dynamics of Free Boundaries." *Journal of Computational Physics*,
   39(1), pp. 201–225.
   - The original VOF paper introducing the method
   - Describes the concept of volume fraction tracking

### Surface Tension Model

3. **Brackbill, J.U., Kothe, D.B. and Zemach, C.** (1992). "A Continuum Method
   for Modeling Surface Tension." *Journal of Computational Physics*, 100(2),
   pp. 335–354.
   - The Continuum Surface Force (CSF) model used in OpenFOAM
   - Converts surface tension to a volume force using ∇α

### OpenFOAM Documentation

4. **OpenFOAM User Guide** — Chapter on multiphase flows and interFoam solver:
   [https://www.openfoam.com/documentation/user-guide](https://www.openfoam.com/documentation/user-guide)

5. **OpenFOAM Wiki — Dam Break Tutorial:**
   [https://openfoamwiki.net/index.php/InterFoam](https://openfoamwiki.net/index.php/InterFoam)

### Related Notes in This Repository

- [Introduction to CFD](../../notes/01_short_intro_to_cfd.md) — Fundamental CFD concepts
- [OpenFOAM Cases](../../notes/02_openfoam_cases.md) — Case directory structure
- [OpenFOAM Dictionaries](../../notes/03_openfoam_dictionaries.md) — Dictionary syntax
- [Meshing](../../notes/04_meshing.md) — Mesh generation with blockMesh
- [Boundary Conditions](../../notes/05_boundary_conditions.md) — BC types and usage
- [Turbulence Models](../../notes/06_turbulence_models.md) — For Exercise 5
- [Parallelization](../../notes/07_parallelization.md) — For Exercise 4 (3D)
- [CFL Number](../../notes/08_cfl_number.md) — Time step control and stability
- [Linear Solvers](../../notes/09_linear_solvers.md) — Solver algorithms

---

*This tutorial is part of the [OpenFOAM Tutorials](../../README.md) repository.*
