# 🚢 Project 11 — Fixed Boat Hull in Calm Water

```
  ╔══════════════════════════════════════════════════════════════════════╗
  ║                                                                      ║
  ║      ~~~~~ ~~~~~~~~ ~~~~~~~~~~  FIXED HULL  ~~~~~~~~~ ~~~~~~~~~      ║
  ║                         ___                                          ║
  ║     ~~~~~  ~~~  _______/   \___   ~~~~~  ~~~~~  ~~~~~  ~~~~~         ║
  ║                /               \                                     ║
  ║     --------==/   BOAT  HULL    \==---------------------------       ║
  ║               \_______________  /                                    ║
  ║     ~~~~~  ~~~ ~~~  ~~~~  ~~  \/  ~~~~~  ~~~~  ~~~~~  ~~~~  ~~       ║
  ║                                                                      ║
  ║      Two-Phase VOF  ·  interFoam  ·  Water + Air  ·  Laminar         ║
  ║                                                                      ║
  ╚══════════════════════════════════════════════════════════════════════╝
```

> **Simulate steady flow around a fixed toy boat hull submerged at the
> waterline, capturing the free-surface deformation (bow wave, stern wake)
> using the Volume of Fluid (VOF) method in OpenFOAM's `interFoam` solver.**

---

## Table of Contents

1.  [Overview](#overview)
2.  [Physics Background](#physics-background)
3.  [Boat Geometry](#boat-geometry)
4.  [Domain Sizing](#domain-sizing)
5.  [File-by-File Walkthrough](#file-by-file-walkthrough)
6.  [STL Scaling: Millimeters to Meters](#stl-scaling-millimeters-to-meters)
7.  [Meshing Strategy](#meshing-strategy)
8.  [Boundary Conditions](#boundary-conditions)
9.  [Phase Initialization (setFields)](#phase-initialization-setfields)
10. [VOF Theory](#vof-theory)
11. [The interFoam Solver](#the-interfoam-solver)
12. [Running the Simulation](#running-the-simulation)
13. [Post-Processing](#post-processing)
14. [Expected Results](#expected-results)
15. [Troubleshooting](#troubleshooting)
16. [Exercises](#exercises)
17. [Cross-References](#cross-references)

---

## Overview

This project simulates **steady flow around a fixed toy boat hull** sitting
at the air-water interface. The hull is held in place (no rigid-body motion)
while an inlet velocity of **0.3 m/s** drives water past it, creating:

- A **bow wave** where water piles up at the front of the hull
- A **stern wake** and wave pattern trailing behind
- **Pressure forces** on the hull (drag and lift)
- A **free surface** that deforms around the waterline

We use OpenFOAM's **`interFoam`** solver, which implements the **Volume of
Fluid (VOF)** method for tracking two immiscible, incompressible fluids
(water and air) sharing the same computational domain.

### Why "Fixed Hull"?

In naval CFD, you typically start with a **fixed hull** simulation before
moving to a **free-floating** hull. Fixing the hull:

- Eliminates the complexity of 6-DOF rigid body motion
- Gives a stable baseline for mesh and solver settings
- Still captures the essential wave-making physics
- Is the natural precursor to Project 12 (floating hull)

### Key Parameters

| Parameter            | Value                                |
|----------------------|--------------------------------------|
| Solver               | `interFoam` (VOF two-phase)          |
| Phases               | Water ($\rho$=998.2 kg/m$^3$) + Air ($\rho$=1.225 kg/m$^3$) |
| Flow speed           | 0.3 m/s (towing speed)               |
| Boat length          | 0.155 m (~6 inches)                  |
| Reynolds number      | Re $\approx$ 46,500 (based on boat length)   |
| Froude number        | Fr $\approx$ 0.243 (based on boat length)    |
| Turbulence model     | Laminar (first pass)                 |
| Surface tension      | $\sigma$ = 0.072 N/m                        |
| Gravity              | g = -9.81 m/s$^2$ (z-direction)         |

The **Froude number** Fr = U / $\sqrt{gL}$ $\approx$ 0.3 / $\sqrt{9.81 × 0.155}$ $\approx$ 0.243,
which is in the wave-making regime — we should see visible surface waves.

---

## Physics Background

### Two-Phase Flow

When a hull moves through water, it creates disturbances at the air-water
interface. Modeling this requires tracking **two fluids** simultaneously:

```
  ATMOSPHERE (air)         ρ = 1.225 kg/m³,  ν = 1.48e-5 m²/s
  ─────────────────────────────────────────────────────────────
  ~~~~~~~~~~~~ FREE SURFACE (interface) ~~~~~~~~~~~~~~~~~~~~~~~~
  ─────────────────────────────────────────────────────────────
  WATER                    ρ = 998.2 kg/m³,  ν = 1.00e-6 m²/s
```

The interface between water and air is not prescribed — it **evolves** as
part of the solution. This is what makes free-surface problems interesting
and computationally challenging.

### Volume of Fluid (VOF)

The VOF method tracks the interface using a **phase fraction** field
`alpha.water` ($\alpha$):

- $\alpha$ = 1: cell is fully water
- $\alpha$ = 0: cell is fully air
- 0 < $\alpha$ < 1: cell contains the interface

A single set of Navier-Stokes equations is solved for the **mixture**, with
material properties that vary according to $\alpha$:

```
ρ_mixture = α · ρ_water + (1 - α) · ρ_air
μ_mixture = α · μ_water + (1 - α) · μ_air
```

### Surface Tension

Surface tension $\sigma$ acts at the interface and is modeled using the
**Continuum Surface Force (CSF)** method. The surface tension force is
converted to a volume force:

```
F_σ = σ κ ∇α
```

where $\kappa$ is the interface curvature estimated from the $\alpha$ field.

---

## Boat Geometry

### Source STL File

The hull geometry is defined in `boat_hull.stl` at the repository root.
This is a **binary STL** file containing **960 triangles**.

> ⚠️ **The STL is in millimeters!** It must be scaled to meters before use.
> The `Allrun` script handles this automatically.

### Bounding Box

```
                        Units: millimeters (as-stored in STL)
     ┌─────────────────────────────────────────────────────────┐
     │                                                         │
     │         ┌───────────────────────────────────┐           │
     │         │                                   │           │
     │         │          BOAT HULL STL            │  Height   │
     │  Width  │                                   │  54.50 mm │
     │ 81.70mm │         Length: 155.06 mm          │           │
     │         │                                   │           │
     │         │           (centered at origin)    │           │
     │         └───────────────────────────────────┘           │
     │                                                         │
     └─────────────────────────────────────────────────────────┘

     Bounding Box (mm):
       X: [-77.53, +77.53]   (length axis, bow to stern)
       Y: [-40.85, +40.85]   (beam/width axis, port to starboard)
       Z: [-27.25, +27.25]   (vertical axis, keel to deck)

     Bounding Box (m) — after scaling ×0.001:
       X: [-0.0775, +0.0775]
       Y: [-0.0409, +0.0409]
       Z: [-0.0273, +0.0273]
```

### Orientation

The STL is centered at the origin with:

- **X-axis**: longitudinal (bow → stern direction, flow direction)
- **Y-axis**: lateral (port ↔ starboard)
- **Z-axis**: vertical (keel → deck, gravity acts in -Z)

The hull is **symmetric about the XZ-plane** (Y=0) and roughly symmetric
about the XY-plane (Z=0), meaning it sits half-above, half-below the
waterline at z = 0.

---

## Domain Sizing

### Rationale

The computational domain must be large enough that:
1. The inlet flow is undisturbed by the hull
2. Wake and waves can develop without hitting boundaries
3. Lateral boundaries don't constrain the flow

Rule of thumb for external flow around a body of length L:

| Direction    | Distance         | Value       | Reason                     |
|-------------|------------------|-------------|----------------------------|
| Upstream    | ~3L from hull    | 0.5 m       | Undisturbed inflow         |
| Downstream  | ~6L from hull    | 1.0 m       | Wake/wave development      |
| Lateral     | ~3W each side    | $\pm$0.25 m     | No wall interference       |
| Below       | ~3H below hull   | 0.15 m      | No bottom effects          |
| Above       | ~3H above water  | 0.15 m      | Air space for waves        |

### Domain Diagram

```
                           1.5 m total (x-direction)
            ◄─────────────────────────────────────────────────────────►

            x = -0.5                                          x = 1.0
         ┌────────────────────────────────────────────────────────────┐ z = 0.15
         │  ATMOSPHERE  (air, α = 0)                                  │
         │                                                            │ ▲
         │               ___......___                                 │ │ 0.15m
  INLET  │          ____/           \____                     OUTLET  │ │
  U=0.3→ │~~~~~~~~/     BOAT HULL      \~~~~~~~~~~~~~~~~~~~~~~~~→    │ z = 0
  m/s    │       \______           _____/                     (open)  │
         │              ‾‾‾......‾‾‾                                  │ │
         │  WATER  (water, α = 1)                                     │ │ 0.15m
         │                                                            │ ▼
         └────────────────────────────────────────────────────────────┘ z = -0.15
                                                                      
         y = -0.25  ◄──── 0.5 m (y-direction) ────►  y = +0.25

         ┌────────────────────────────────────────────────────────────┐
         │               SIDE (symmetryPlane)                          │
         │                                                            │
         │                    TOP VIEW                                │
         │             ___......___                                   │
  INLET  │        ____/    HULL    \____                      OUTLET  │
  ────►  │       ‾‾‾‾\___........___/‾‾‾‾                     ────►  │
         │                                                            │
         │               SIDE (symmetryPlane)                          │
         └────────────────────────────────────────────────────────────┘
```

---

## File-by-File Walkthrough

### Directory Structure

```
11_boat_hull_fixed/
├── Allrun                              ← Master run script
├── Allclean                            ← Cleanup script
├── README.md                           ← This file
│
├── 0/                                  ← Initial & boundary conditions
│   ├── U                               ← Velocity field
│   ├── p_rgh                           ← Pressure (hydrostatic removed)
│   └── alpha.water                     ← Phase fraction (water/air)
│
├── constant/                           ← Physical properties
│   ├── transportProperties             ← Fluid properties (water + air)
│   ├── turbulenceProperties            ← Turbulence model (laminar)
│   ├── g                               ← Gravity vector
│   └── triSurface/                     ← STL geometry (populated by Allrun)
│       └── boat_hull.stl               ← (copied & scaled at runtime)
│
└── system/                             ← Solver & mesh configuration
    ├── controlDict                     ← Time control, I/O, function objects
    ├── blockMeshDict                   ← Background mesh definition
    ├── snappyHexMeshDict               ← Mesh refinement around hull
    ├── surfaceFeatureExtractDict        ← Edge feature extraction
    ├── setFieldsDict                   ← Phase initialization
    ├── fvSchemes                       ← Discretization schemes
    └── fvSolution                      ← Linear solver settings
```

### 0/U — Velocity Field

Defines the velocity boundary conditions and initial internal field.

- **inlet**: `fixedValue (0.3 0 0)` — uniform towing speed
- **outlet**: `pressureInletOutletVelocity` — allows outflow, prevents inflow
- **atmosphere**: `pressureInletOutletVelocity` — open top
- **bottom**: `noSlip` — stationary tank floor
- **port**: `symmetryPlane` — lateral symmetry boundary at y = -0.25
- **starboard**: `symmetryPlane` — lateral symmetry boundary at y = +0.25
- **boat**: `noSlip` — no-slip on hull surface
- **internalField**: `uniform (0 0 0)` — fluid initially at rest

### 0/p_rgh — Pressure Field

The pressure variable `p_rgh = p - ρgh` has the hydrostatic component
removed, which improves numerical conditioning for free-surface flows.

- **inlet/bottom/boat**: `fixedFluxPressure` — adjusts pressure to match
  the velocity boundary condition
- **port/starboard**: `symmetryPlane`
- **outlet/atmosphere**: `totalPressure` with p0 = 0 — reference pressure

### 0/alpha.water — Phase Fraction

The water volume fraction. Initialized to 0 everywhere (all air), then
`setFields` fills everything below z = 0 with water ($\alpha$ = 1).

- **inlet**: `fixedValue 1` — water enters from the left
- **outlet/atmosphere**: `inletOutlet` — prevents air from flowing in as water
- **boat/bottom**: `zeroGradient` — interface meets walls naturally
- **port/starboard**: `symmetryPlane`

### constant/transportProperties

Defines the two-phase material properties:

```
Water:  ν = 1.00 × 10⁻⁶ m²/s,  ρ = 998.2 kg/m³   (at 20°C)
Air:    ν = 1.48 × 10⁻⁵ m²/s,  ρ = 1.225 kg/m³    (at 20°C, 1 atm)
Surface tension: σ = 0.072 N/m
```

### constant/turbulenceProperties

Set to `laminar` for this initial study. The Reynolds number (~46,500) is
in the transitional regime, but laminar gives a good first approximation
and avoids the complexity of turbulence modeling at the free surface.

### constant/g — Gravity

Gravity points in the -z direction: `(0 0 -9.81)`. This is essential for:
- Hydrostatic pressure distribution
- Free surface position (water settles under gravity)
- Wave dynamics

### system/blockMeshDict — Background Mesh

Creates the base hexahedral mesh that `snappyHexMesh` will refine.

- **Domain**: [-0.5, 1.0] $\times$ [-0.25, 0.25] $\times$ [-0.15, 0.15]
- **Cells**: 60 $\times$ 20 $\times$ 12 = 14,400 cells
- **Cell size**: ~25 mm (coarse — snappy will refine near hull)

### system/controlDict — Run Control

| Setting          | Value           | Rationale                           |
|-----------------|-----------------|-------------------------------------|
| deltaT          | 0.0001 s        | Small initial step for stability    |
| endTime         | 2 s             | ~4 flow-throughs of the domain      |
| adjustTimeStep  | yes             | Adaptive time stepping              |
| maxCo           | 0.5             | Courant number limit                |
| maxAlphaCo      | 0.3             | Interface Courant limit (stricter)  |
| writeInterval   | 0.05 s          | 40 output frames                    |

The `forceCoeffs` function object computes drag and lift coefficients on
the hull patch at every 10th time step.

### system/fvSchemes — Discretization

Key entries for interFoam:

| Term                                          | Scheme              | Notes                        |
|----------------------------------------------|---------------------|------------------------------|
| `div(rhoPhi,U)`                              | linearUpwindV       | Momentum convection          |
| `div(phi,alpha)`                             | vanLeer             | Phase fraction (TVD)         |
| `div(phirb,alpha)`                           | linear              | Interface compression        |
| `div(((rho*nuEff)*dev2(T(grad(U)))))`        | linear              | Viscous stress tensor        |

### system/fvSolution — Linear Solvers

- **alpha.water**: `smoothSolver` with MULES correction for boundedness
- **p_rgh**: `PCG` with GAMG preconditioner (efficient for pressure)
- **U**: `smoothSolver` with Gauss-Seidel
- **PIMPLE**: 1 outer, 3 pressure corrections, 1 non-orthogonal correction

### system/setFieldsDict — Phase Initialization

Uses `boxToCell` to set $\alpha$ = 1 in a huge box from (-10,-10,-10) to (10,10,0).
Since z = 0 is the waterline, everything below gets filled with water.

### system/surfaceFeatureExtractDict

Extracts sharp edges (angle > 150°) from the STL to create a `.eMesh`
file. These edges guide `snappyHexMesh` for accurate feature snapping.

### system/snappyHexMeshDict — Mesh Refinement

Three-stage mesh generation:

1. **Castellated mesh**: Refines cells near the hull surface
   - Surface refinement: levels 2-3 (cells 4-8$\times$ smaller than background)
   - Region refinement: level 2 in a box around the waterline
2. **Snap**: Moves cell vertices onto the STL surface
3. **Add layers**: Disabled for first pass (set `addLayers false`)

---

## STL Scaling: Millimeters to Meters

The original `boat_hull.stl` is authored in millimeters. OpenFOAM works in
SI units (meters). The `Allrun` script handles conversion:

```bash
# Copy from repo root to case triSurface directory
cp ../../boat_hull.stl constant/triSurface/boat_hull.stl

# Scale by 0.001 (mm → m) IN-PLACE
surfaceTransformPoints -scale '(0.001 0.001 0.001)' \
    constant/triSurface/boat_hull.stl \
    constant/triSurface/boat_hull.stl
```

**Before scaling** (mm):
```
X: [-77.53, 77.53]  →  Length = 155.06 mm
Y: [-40.85, 40.85]  →  Width  =  81.70 mm
Z: [-27.25, 27.25]  →  Height =  54.50 mm
```

**After scaling** (m):
```
X: [-0.0775, 0.0775]  →  Length = 0.155 m
Y: [-0.0409, 0.0409]  →  Width  = 0.082 m
Z: [-0.0273, 0.0273]  →  Height = 0.055 m
```

---

## Meshing Strategy

### Stage 1: Background Mesh (blockMesh)

`blockMesh` creates a uniform hexahedral grid filling the entire domain:

```
  ┌──┬──┬──┬──┬──┬──┬──┬──┬──┬──┐
  │  │  │  │  │  │  │  │  │  │  │     60 × 20 × 12 = 14,400 cells
  ├──┼──┼──┼──┼──┼──┼──┼──┼──┼──┤     Cell size ≈ 25 mm
  │  │  │  │  │  │  │  │  │  │  │
  ├──┼──┼──┼──┼──┼──┼──┼──┼──┼──┤
  │  │  │  │  │  │  │  │  │  │  │
  └──┴──┴──┴──┴──┴──┴──┴──┴──┴──┘
```

### Stage 2: Surface Feature Extraction

`surfaceFeatureExtract` identifies sharp edges on the STL surface (where
the angle between adjacent triangles exceeds 150°). These edges are saved
as `.eMesh` files and used by `snappyHexMesh` to ensure the mesh conforms
accurately to hull features like the bow, keel line, and transom.

### Stage 3: snappyHexMesh

`snappyHexMesh` transforms the coarse background mesh into a body-fitted
mesh around the hull:

```
  ┌──────────────────────────────────────────────────────────────────┐
  │                                                                  │
  │  Coarse background cells (level 0, ~25mm)                        │
  │                                                                  │
  │           ┌────────────────────────────┐                         │
  │           │  Refinement box (level 2)  │                         │
  │           │  ~6mm cells                │                         │
  │           │      ┌──────────┐          │                         │
  │           │      │Hull surf.│          │                         │
  │           │      │(level 3) │          │                         │
  │           │      │ ~3mm     │          │                         │
  │           │      └──────────┘          │                         │
  │           └────────────────────────────┘                         │
  │                                                                  │
  └──────────────────────────────────────────────────────────────────┘
```

The `locationInMesh` point `(0.3, 0, 0.05)` tells snappy which side of the
hull surface is the "fluid" side (the point must be clearly outside the
boat and inside the domain).

**Approximate cell counts after snappyHexMesh**: ~100k-200k cells
(depends on refinement level details).

---

## Boundary Conditions

### Domain Patch Layout

```
                    atmosphere (top, z = 0.15)
         ┌────────────────────────────────────────────┐
         │                 ↑ open air                  │
         │                                             │
  inlet  │  U = 0.3 m/s →        → → →        outlet  │
  (left) │  α = 1 (water)    ___/‾‾‾\___    (right)   │
  x=-0.5 │  p: fixedFlux    / BOAT HULL \   p: total  │  x=1.0
         │                  \___________/    α: outflow│
         │                                             │
         │  ══════════ water (α = 1) ═══════════       │
         │                                             │
         └────────────────────────────────────────────┘
                     bottom (wall, z = -0.15)

         port/starboard: symmetryPlane (y = ±0.25)
```

### Boundary Condition Summary Table

| Patch      | U                             | p_rgh              | alpha.water       |
|-----------|-------------------------------|--------------------|--------------------|
| inlet      | fixedValue (0.3 0 0)          | fixedFluxPressure  | fixedValue 1       |
| outlet     | pressureInletOutletVelocity   | totalPressure 0    | inletOutlet 0      |
| atmosphere | pressureInletOutletVelocity   | totalPressure 0    | inletOutlet 0      |
| bottom     | noSlip                        | fixedFluxPressure  | zeroGradient       |
| port       | symmetryPlane                 | symmetryPlane      | symmetryPlane      |
| starboard  | symmetryPlane                 | symmetryPlane      | symmetryPlane      |
| boat       | noSlip                        | fixedFluxPressure  | zeroGradient       |

### Why These Choices?

- **fixedFluxPressure**: Adjusts the pressure gradient at the boundary to
  be consistent with the velocity BC. Essential for gravity-driven flows.
- **totalPressure**: Sets p$_0$ = 0 at the outlet/atmosphere, allowing
  pressure to vary with the velocity (Bernoulli-like condition).
- **pressureInletOutletVelocity**: A smart BC that applies zero-gradient
  for outflow and fixedValue for inflow (preventing reverse flow).
- **inletOutlet for $\alpha$**: Prevents air from being "sucked in" as water at
  outflow boundaries.

---

## Phase Initialization (setFields)

After meshing, the entire domain is initialized as air ($\alpha$ = 0). The
`setFields` utility reads `setFieldsDict` and fills the water region:

```
  BEFORE setFields:                    AFTER setFields:
  ┌──────────────────────┐             ┌──────────────────────┐
  │                      │             │                      │
  │   α = 0 (all air)   │             │   α = 0 (air)        │  z > 0
  │                      │             │                      │
  │                      │             ├══════════════════════╡  z = 0
  │                      │             │                      │
  │                      │             │   α = 1 (water)      │  z < 0
  │                      │             │                      │
  └──────────────────────┘             └──────────────────────┘
```

The `boxToCell` selector fills a giant box from (-10,-10,-10) to (10,10,0)
with $\alpha$ = 1. The box is intentionally oversized — only cells inside the
mesh domain are affected. The z = 0 cutoff creates the waterline.

**Inside the hull**: Cells inside the boat geometry were removed by
snappyHexMesh, so `setFields` naturally skips them. The hull acts as a
solid obstacle.

---

## VOF Theory

### The Volume Fraction Equation

The VOF method solves a transport equation for the phase fraction $\alpha$:

```
∂α/∂t + ∇·(U α) + ∇·(U_r α(1-α)) = 0
```

where:
- `α`: volume fraction of water (0 = air, 1 = water)
- `U`: velocity field (shared by both phases)
- `U_r`: compression velocity (artificial, sharpens the interface)

The third term `∇·(U_r α(1-α))` is the **interface compression** term,
unique to OpenFOAM's VOF implementation. The factor `α(1-α)` ensures
compression only acts at the interface (where 0 < $\alpha$ < 1), not in pure
fluid regions.

### MULES Algorithm

OpenFOAM uses the **Multidimensional Universal Limiter with Explicit
Solution (MULES)** algorithm to solve the $\alpha$ equation. MULES ensures:

1. **Boundedness**: $\alpha$ stays between 0 and 1 (no unphysical overshoot)
2. **Sharpness**: The interface remains thin (1-3 cells wide)
3. **Conservation**: Total water volume is conserved

The `cAlpha` parameter (set to 1 in `fvSolution`) controls the strength
of interface compression. Higher values give sharper interfaces but may
introduce artifacts.

### Material Properties

The mixture properties are computed cell-by-cell:

```
ρ = α · 998.2 + (1-α) · 1.225        [kg/m³]
μ = α · 998.2e-6 + (1-α) · 1.225·1.48e-5  [Pa·s]
```

This creates a smooth transition across the interface, avoiding sharp
discontinuities that would destabilize the numerics.

---

## The interFoam Solver

### What It Solves

`interFoam` solves the following coupled system at each time step:

1. **Volume fraction equation** ($\alpha$ transport):
   ```
   ∂α/∂t + ∇·(Uα) + ∇·(U_c α(1-α)) = 0
   ```

2. **Momentum equation** (Navier-Stokes for the mixture):
   ```
   ∂(ρU)/∂t + ∇·(ρUU) = -∇p_rgh - g·x∇ρ + ∇·(μ_eff ∇U) + σκ∇α
   ```

3. **Pressure equation** (from continuity):
   ```
   ∇·(∇p_rgh / ρ) = ∇·U*  (approximate)
   ```

The variables are:
- `U`: velocity vector field
- `p_rgh = p - ρgh`: modified pressure (hydrostatic part removed)
- `α`: water volume fraction
- `ρ, μ`: mixture density and viscosity (functions of $\alpha$)
- `σκ∇α`: surface tension force (CSF model)

### Solution Algorithm (PIMPLE)

The PIMPLE algorithm (merged PISO-SIMPLE) proceeds:

```
For each time step:
  1. Solve α equation (MULES)
  2. Update ρ and μ from α
  3. Outer correction loop (nOuterCorrectors = 1):
     a. Solve momentum predictor → U*
     b. Inner correction loop (nCorrectors = 3):
        i.   Assemble and solve pressure equation → p_rgh
        ii.  Non-orthogonal corrections (nNonOrthogonalCorrectors = 1)
        iii. Correct velocity U from pressure gradient
  4. Advance time
```

### Adaptive Time Stepping

With `adjustTimeStep yes`, interFoam automatically adjusts $\Delta t$ to satisfy:

- **maxCo = 0.5**: Maximum Courant number for velocity
- **maxAlphaCo = 0.3**: Maximum Courant number for the interface

The interface Courant limit is stricter because the $\alpha$ equation is solved
explicitly and is more sensitive to time step size. See
[Note 08 — The CFL Number](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md) for details.

---

## Running the Simulation

### Prerequisites

- OpenFOAM installed and environment sourced (`. /opt/openfoam*/etc/bashrc`)
- The `boat_hull.stl` file at the repository root

### Step-by-Step Execution

The `Allrun` script automates everything:

```bash
cd projects/11_boat_hull_fixed
./Allrun
```

To reproduce the exact containerized run used for the checked-in results:

```bash
docker run --rm -u "$(id -u):$(id -g)" \
  -v "$PWD":/work \
  -w /work/projects/11_boat_hull_fixed \
  cfdengine/openfoam \
  bash -lc 'source /opt/openfoam6/etc/bashrc && ./Allclean && ./Allrun'
```

Here's what happens at each step:

#### Step 1: Copy and Scale STL

```bash
cp ../../boat_hull.stl constant/triSurface/boat_hull.stl
surfaceTransformPoints -scale '(0.001 0.001 0.001)' \
    constant/triSurface/boat_hull.stl constant/triSurface/boat_hull.stl
```

Copies the STL from the repo root and scales from millimeters to meters.
Check the log: `log.surfaceTransformPoints`

#### Step 2: Background Mesh

```bash
blockMesh
```

Creates 14,400 hexahedral cells. Check: `log.blockMesh`

#### Step 3: Feature Extraction

```bash
surfaceFeatureExtract
```

Reads `boat_hull.stl`, identifies sharp edges, writes `boat_hull.eMesh`.
Check: `log.surfaceFeatureExtract`

#### Step 4: Mesh Refinement

```bash
snappyHexMesh -overwrite
```

Refines and snaps the mesh around the hull. The `-overwrite` flag writes
the final mesh to `constant/polyMesh` instead of creating time directories.
Check: `log.snappyHexMesh`

#### Step 5: Mesh Quality Check

```bash
checkMesh
```

Reports mesh statistics and quality metrics. Look for:
- "Mesh OK" at the end
- No "***" warnings about severely non-orthogonal faces
Check: `log.checkMesh`

#### Step 6: Phase Initialization

```bash
setFields
```

Sets $\alpha$ = 1 below z = 0 (water region). Check: `log.setFields`

#### Step 7: Run Solver

```bash
interFoam
```

Runs the VOF solver for 2 seconds of simulation time.
Check: `log.interFoam`

### Monitoring Progress

While `interFoam` runs, you can monitor:

```bash
# Watch the log in real time
tail -f log.interFoam

# Check Courant numbers
grep "Courant Number" log.interFoam | tail -20

# Check interface Courant numbers
grep "Interface Courant" log.interFoam | tail -20

# Check continuity errors
grep "continuity errors" log.interFoam | tail -10

# Monitor time step
grep "^deltaT" log.interFoam | tail -20
```

### Parallel Execution

For faster runs, decompose and run in parallel:

```bash
# Add to Allrun (before interFoam):
runApplication decomposePar
runParallel interFoam
runApplication reconstructPar

# Create system/decomposeParDict:
# method scotch; numberOfSubdomains 4;
```

---

## Post-Processing

### ParaView Visualization

Open the results in ParaView:

```bash
paraFoam
# or
paraview --data=boat_hull_fixed.foam
```

For a prepared view with a center slice, hull surface, and streamlines:

```bash
./open_paraview_streamlines.sh
```

#### Visualize the Free Surface

1. Load the case and click **Apply**
2. Add a **Clip** filter:
   - Scalars: `alpha.water`
   - Value: `0.5`
   - Inside Out: unchecked
3. This clips away the air, showing only the water volume
4. Color by `p_rgh` to see pressure distribution on the water surface

#### Water Surface Contour

1. Apply a **Contour** filter on `alpha.water` at value `0.5`
2. This creates an isosurface at the air-water interface
3. Color by `U` magnitude to see flow speed along the waterline

#### Pressure on Hull

1. Extract the hull surface: **Extract Block** → select `boat` patch
2. Color by `p_rgh` to see pressure distribution
3. High pressure at bow (stagnation), low pressure along sides

#### Velocity Field

1. Apply a **Slice** filter at y = 0 (centerplane)
2. Color by `U` magnitude
3. Add **Glyph** filter for velocity vectors

#### Animation

1. Set the time toolbar to the desired range
2. File → Save Animation → choose format and resolution
3. This creates a frame-by-frame animation of wave evolution

### Force Extraction

The `forceCoeffs` function object writes data to:

```
postProcessing/forceCoeffs/0/forceCoeffs.dat
```

Plot drag coefficient vs. time:

```bash
# Using gnuplot:
gnuplot -e "
set xlabel 'Time (s)';
set ylabel 'Cd';
plot 'postProcessing/forceCoeffs/0/forceCoeffs.dat' u 1:3 w l title 'Drag';
pause -1
"

# Or extract drag force:
cat postProcessing/forces/0/forces.dat | tail -20
```

---

## Expected Results

### Flow Features

At Fr $\approx$ 0.24 and Re $\approx$ 46,500, you should observe:

```
  SIDE VIEW (centerplane, y = 0):

                     bow wave
                       /\
                      /  \
  ─── ─── ─── ──────/    \___     ___/\___/\___  ─── ─── ─── ───
                    │ BOW  │  \__/                    undisturbed
  ══════════════════│      │══════════════════════════════════════
                    │HULL  │     wake region
  ══════════════════│      │══════════════════════════════════════
                    │STERN │  ___
  ─── ─── ─── ─────\______/\/   \___/\___ ─── ─── ─── ─── ─── 
                       stern wave
```

1. **Bow wave**: Water piles up at the front, creating a standing wave
2. **Wake pattern**: V-shaped wave pattern behind the hull (Kelvin wake)
3. **Stern depression**: Water level drops behind the hull
4. **Boundary layer**: Thin layer of slow-moving water on the hull surface

### Pressure Distribution

```
  HULL SURFACE PRESSURE (approximate):

  Bow (front):    HIGH pressure  (stagnation, ~½ρU² ≈ 45 Pa)
  Sides:          LOW pressure   (accelerated flow)
  Stern (back):   MODERATE       (pressure recovery, some separation)
```

### Drag Components

The total drag on the hull has two main components:

1. **Pressure drag (form drag)**: Due to pressure difference bow vs. stern
2. **Friction drag**: Due to viscous shear stress on the hull surface

For a bluff body like a toy boat at this Reynolds number, **pressure drag
dominates** (~70-80% of total drag).

### Typical Values

| Quantity           | Expected Range      |
|-------------------|---------------------|
| Total drag force   | 0.01–0.05 N         |
| Drag coefficient   | Cd $\approx$ 0.1–0.5        |
| Bow wave height    | 1–5 mm              |
| Max Courant number | $\leq$ 0.5 (controlled)  |
| Time step (steady) | ~0.001–0.01 s       |

---

## Troubleshooting

### Common Issues

#### 1. Interface Smearing

**Symptom**: The water surface becomes diffuse ($\alpha$ spreads over many cells).

**Causes and fixes**:
- Mesh too coarse at waterline → increase refinement level in `snappyHexMeshDict`
- cAlpha too low → increase `cAlpha` in `fvSolution` (try 1.5 or 2)
- Time step too large → reduce `maxAlphaCo` (try 0.2)
- Wrong div scheme → ensure `div(phi,alpha)` uses a TVD scheme (vanLeer)

#### 2. Divergence (Solver Crash)

**Symptom**: "FOAM FATAL ERROR" or extremely large velocity/pressure values.

**Causes and fixes**:
- Time step too large → reduce `maxCo` and `maxAlphaCo`
- Bad mesh quality → run `checkMesh` and fix issues
- Incorrect BCs → verify boundary conditions match the physics
- locationInMesh inside the boat → move it clearly outside

#### 3. Mesh Quality Warnings

**Symptom**: `checkMesh` reports non-orthogonality > 70° or negative volumes.

**Fixes**:
- Increase `nNonOrthogonalCorrectors` in PIMPLE (try 2 or 3)
- Relax `maxNonOrtho` in `meshQualityControls` of `snappyHexMeshDict`
- Reduce refinement level or increase `nCellsBetweenLevels`

#### 4. Slow Convergence

**Symptom**: Simulation runs but residuals don't drop.

**Fixes**:
- Increase `nOuterCorrectors` (try 2-3 for better accuracy)
- Check that `adjustTimeStep` is enabled
- Ensure initial conditions are reasonable

#### 5. STL Problems

**Symptom**: snappyHexMesh fails or produces strange mesh.

**Fixes**:
- Verify STL is watertight: `surfaceCheck constant/triSurface/boat_hull.stl`
- Check STL bounds after scaling: `surfaceMeshInfo constant/triSurface/boat_hull.stl`
- Ensure the STL is inside the blockMesh domain

#### 6. Boundary Layer Not Resolved

**Symptom**: Flow doesn't slow down near the hull surface.

**Fix**: Enable `addLayers true` in `snappyHexMeshDict` and configure
appropriate layer parameters. This adds thin cells near the hull surface
for boundary layer resolution.

---

## Exercises

### Exercise 1: Mesh Sensitivity Study

Run the simulation with different refinement levels:

1. Change `refinementSurfaces` level from (2 3) to (1 2) — coarser
2. Change to (3 4) — finer
3. Compare: free-surface shape, drag force, and computation time
4. At what refinement does the drag converge to within 5%?

### Exercise 2: Speed Variation

Modify the inlet velocity and observe how waves change:

1. Run at U = 0.15 m/s (Fr $\approx$ 0.12) — minimal wave making
2. Run at U = 0.3 m/s (Fr $\approx$ 0.24) — moderate waves (baseline)
3. Run at U = 0.6 m/s (Fr $\approx$ 0.49) — strong waves
4. Plot drag vs. speed — does it follow a U$^2$ relationship?

### Exercise 3: Enable Turbulence

Switch from laminar to a turbulence model:

1. Change `turbulenceProperties` to `simulationType RAS`
2. Add `RAS { RASModel kOmegaSST; ... }` settings
3. Create `0/k`, `0/omega`, and `0/nut` boundary condition files
4. Compare drag with the laminar result — how much does it change?

### Exercise 4: Add Boundary Layers

Enable the layer addition in `snappyHexMeshDict`:

1. Set `addLayers true`
2. Configure 3-5 layers on the boat patch with expansion ratio 1.2
3. Check y+ values on the hull after running
4. How does the velocity profile near the hull change?

### Exercise 5: Domain Size Study

Test whether the domain is large enough:

1. Double the downstream length (xMax = 2.0 m)
2. Halve the lateral width (yMax = 0.125 m)
3. Compare results — if they differ significantly, the original domain
   was too small in that direction

### Exercise 6: Transition to Floating Hull

This project is the stepping stone to Project 12 (floating hull):

1. Review the current fixed-hull results
2. Read about `sixDoFRigidBodyMotion` in OpenFOAM
3. Estimate the hull's buoyancy: does the waterline position make sense?
4. What additional dictionary entries would you need for a floating hull?
5. Proceed to [Project 12](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/12_boat_hull_floating) to see the answer

---

## Cross-References

### Related Notes

| Note | Topic | Relevance |
|------|-------|-----------|
| [01 — CFD Introduction](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md) | Fundamentals of CFD | Navier-Stokes equations, discretization basics |
| [04 — Meshing](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/04_meshing.md) | Meshing in OpenFOAM | blockMesh, snappyHexMesh theory and practice |
| [05 — Boundary Conditions](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md) | BCs in OpenFOAM | Detailed explanation of all BC types used here |
| [08 — CFL Number](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md) | Stability & time stepping | Why maxCo and maxAlphaCo matter for interFoam |
| [09 — Linear Solvers](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/09_linear_solvers.md) | Solver settings | PCG, GAMG, smoothSolver configuration |

### Related Projects

| Project | Description | Connection |
|---------|-------------|------------|
| [08 — Dam Break](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/08_dam_break) | Two-phase dam break simulation | Also uses interFoam/VOF — simpler geometry, same physics |
| [12 — Floating Hull](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/12_boat_hull_floating) | Free-floating boat hull | Next step — adds 6-DOF rigid body motion to this case |
| [04 — NACA Airfoil](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis) | External aerodynamics | Similar external flow concepts, single-phase |
| [06 — Ahmed Body](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/06_ahmed_body_aerodynamics) | Bluff body aerodynamics | Similar drag analysis concepts |

### Learning Path

```
  ┌─────────────────────┐     ┌──────────────────────────┐
  │  08 — Dam Break     │ ──► │  11 — Fixed Boat Hull    │
  │  (intro to VOF)     │     │  (THIS PROJECT)          │
  └─────────────────────┘     └──────────┬───────────────┘
                                         │
                                         ▼
                              ┌──────────────────────────┐
                              │  12 — Floating Boat Hull  │
                              │  (adds 6-DOF motion)      │
                              └──────────────────────────┘
```

---

## References

- OpenFOAM User Guide: [interFoam](https://www.openfoam.com/documentation/guides/latest/doc/guide-applications-solvers-multiphase-interFoam.html)
- OpenFOAM Tutorial: `$FOAM_TUTORIALS/multiphase/interFoam/laminar/damBreak`
- Berberović, E. et al. (2009). "Drop impact onto a liquid layer of finite
  thickness." *Physical Review E*, 79(3).
- Hirt, C.W. & Nichols, B.D. (1981). "Volume of Fluid (VOF) Method for
  the Dynamics of Free Boundaries." *Journal of Computational Physics*,
  39(1), 201-225.

---

*Part of the [OpenFOAM Tutorials](https://github.com/djeada/OpenFoam-Tutorials/blob/main/README.md) learning series.*
