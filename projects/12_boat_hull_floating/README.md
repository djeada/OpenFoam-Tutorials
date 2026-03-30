# Project 12 — Floating Boat Hull (6-DoF Rigid Body Motion)

```
  ╔══════════════════════════════════════════════════════════════════════╗
  ║                                                                      ║
  ║   ██████╗  ██████╗  █████╗ ████████╗    ██╗  ██╗██╗   ██╗██╗         ║
  ║   ██╔══██╗██╔═══██╗██╔══██╗╚══██╔══╝    ██║  ██║██║   ██║██║         ║
  ║   ██████╔╝██║   ██║███████║   ██║       ███████║██║   ██║██║         ║
  ║   ██╔══██╗██║   ██║██╔══██║   ██║       ██╔══██║██║   ██║██║         ║
  ║   ██████╔╝╚██████╔╝██║  ██║   ██║       ██║  ██║╚██████╔╝███████╗    ║
  ║   ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝       ╚═╝  ╚═╝ ╚═════╝╚══════╝     ║
  ║                                                                      ║
  ║   F L O A T I N G   —   6 - D o F   R i g i d   B o d y              ║
  ║                                                                      ║
  ║   Solver:   interDyMFoam  (dynamic mesh + VOF + 6-DoF)               ║
  ║   Physics:  Two-phase free surface + rigid body motion               ║
  ║   Motion:   Heave · Pitch · Roll  (translation + rotation)           ║
  ║                                                                      ║
  ╚══════════════════════════════════════════════════════════════════════╝
```

---

## Table of Contents

1.  [Overview](#1-overview)
2.  [Prerequisites — Run Project 11 First](#2-prerequisites--run-project-11-first)
3.  [What Makes This Different from Project 11](#3-what-makes-this-different-from-project-11)
4.  [Physics Background](#4-physics-background)
    - [4.1 Volume of Fluid (VOF) — Recap](#41-volume-of-fluid-vof--recap)
    - [4.2 6-DoF Rigid Body Motion — The New Ingredient](#42-6-dof-rigid-body-motion--the-new-ingredient)
    - [4.3 Coupling: Fluid Forces ↔ Body Motion](#43-coupling-fluid-forces--body-motion)
5.  [Dynamic Mesh — How the Mesh Moves](#5-dynamic-mesh--how-the-mesh-moves)
6.  [Boat Geometry and Physical Properties](#6-boat-geometry-and-physical-properties)
7.  [Domain and Mesh Design](#7-domain-and-mesh-design)
8.  [File-by-File Walkthrough](#8-file-by-file-walkthrough)
    - [8.1  0/U](#81-0u)
    - [8.2  0/p_rgh](#82-0p_rgh)
    - [8.3  0/alpha.water](#83-0alphawater)
    - [8.4  0/pointDisplacement](#84-0pointdisplacement)
    - [8.5  constant/transportProperties](#85-constanttransportproperties)
    - [8.6  constant/turbulenceProperties](#86-constantturbulenceproperties)
    - [8.7  constant/g](#87-constantg)
    - [8.8  constant/dynamicMeshDict — The Key File](#88-constantdynamicmeshdict--the-key-file)
    - [8.9  system/blockMeshDict](#89-systemblockmeshdict)
    - [8.10 system/controlDict](#810-systemcontroldict)
    - [8.11 system/fvSchemes](#811-systemfvschemes)
    - [8.12 system/fvSolution](#812-systemfvsolution)
    - [8.13 system/setFieldsDict](#813-systemsetfieldsdict)
    - [8.14 system/surfaceFeatureExtractDict](#814-systemsurfacefeatureextractdict)
    - [8.15 system/snappyHexMeshDict](#815-systemsnappyhexmeshdict)
9.  [Running the Simulation](#9-running-the-simulation)
10. [Post-Processing](#10-post-processing)
11. [Expected Results](#11-expected-results)
12. [Stability and Troubleshooting](#12-stability-and-troubleshooting)
13. [Exercises](#13-exercises)
14. [Cross-References](#14-cross-references)

---

## 1. Overview

This project simulates a **toy boat hull floating freely** on a water surface.
Unlike Project 11 — where the hull was *fixed* in the domain and the water
flowed past it — here the boat is an unconstrained (partially constrained)
rigid body that **responds to fluid forces in real time**.

The solver is **`interDyMFoam`**, which combines three capabilities:

| Capability              | What it does                                      |
|-------------------------|---------------------------------------------------|
| **VOF** (Volume of Fluid)    | Tracks the air / water interface               |
| **Dynamic mesh**             | Deforms the mesh as the body moves             |
| **6-DoF rigid body motion**  | Solves Newton / Euler equations for the hull   |

```
  interDyMFoam  =  interFoam  +  dynamic mesh  +  sixDoFRigidBodyMotion
                   (VOF solver)   (mesh morph)    (rigid body dynamics)
```

The simulation loop each timestep:

```
  ┌─────────────────────────────────────────────────────────────┐
  │  1. Solve Navier-Stokes + VOF on current mesh               │
  │  2. Compute pressure + viscous forces on boat surface       │
  │  3. Pass forces → 6-DoF solver → new position / rotation    │
  │  4. Move mesh to conform to new body position               │
  │  5. Go to next timestep                                     │
  └─────────────────────────────────────────────────────────────┘
```

The boat starts at the water surface and should settle under gravity to its
equilibrium draft — the depth at which the buoyancy force equals its weight
(Archimedes' principle). You will observe damped heave oscillations as the
boat "bobs" up and down before reaching a steady state.

---

## 2. Prerequisites — Run Project 11 First

**Project 11** (`11_boat_hull_fixed/`) is the *fixed-hull* version of this
case. It uses `interFoam` — the same VOF two-phase solver, but with a
**static mesh**. The hull sits motionless while water flows around it.

If you haven't run Project 11 yet, do that first. It teaches:

- How to import and scale an STL surface (mm → m)
- Background mesh creation with `blockMesh`
- Surface snapping with `snappyHexMesh`
- Two-phase VOF initialisation with `setFields`
- Interpreting free-surface results in ParaView

Project 12 builds on *all* of those skills and adds:

- Dynamic mesh motion (`dynamicMotionSolverFvMesh`)
- Rigid body dynamics (`sixDoFRigidBodyMotion`)
- New boundary conditions (`movingWallVelocity`, `pointDisplacement`)
- Much tighter timestep control for mesh-motion stability

```
  Progression:   Project 08  →  Project 11  →  Project 12
                 (dam break)    (fixed hull)   (floating hull)
                  basic VOF      STL + VOF      STL + VOF + 6-DoF
```

---

## 3. What Makes This Different from Project 11

| Aspect                     | Project 11 (Fixed)           | Project 12 (Floating)              |
|----------------------------|------------------------------|-------------------------------------|
| Solver                     | `interFoam`                  | `interDyMFoam`                      |
| Mesh                       | Static                       | Dynamic — deforms every timestep    |
| Hull motion                | None (fixed in space)        | Heave + Pitch + Roll (3-DoF)        |
| Boat wall BC for U         | `noSlip`                     | `movingWallVelocity`                |
| Extra BC file              | —                            | `0/pointDisplacement`               |
| Extra constant file        | —                            | `constant/dynamicMeshDict`          |
| fvSolution                 | PIMPLE                       | PIMPLE + `moveMeshOuterCorrectors`  |
| fvSchemes                  | Standard interFoam           | + `cellDisplacement` laplacian      |
| Timestep                   | ~1e-4                        | ~5e-5 (smaller for mesh stability)  |
| Runtime                    | ~2 s                         | ~3 s (needs settling time)          |
| Computational cost         | Moderate                     | High — mesh motion is expensive     |

The **dynamicMeshDict** is the most important new file. It tells OpenFOAM:
- Which patches belong to the rigid body (`boat`)
- The body's mass, centre of mass, and moments of inertia
- How far the mesh morphing extends from the body
- What constraints apply (e.g., suppress surge/sway/yaw)

---

## 4. Physics Background

### 4.1 Volume of Fluid (VOF) — Recap

The VOF method tracks the interface between two immiscible fluids (water and
air) using a scalar field **α** (alpha):

```
  α = 1  →  pure water
  α = 0  →  pure air
  0 < α < 1  →  interface region
```

The transport equation for α:

```
  ∂α/∂t + ∇·(U α) + ∇·(U_r α(1-α)) = 0
```

where `U_r` is a compression velocity that sharpens the interface.

Fluid properties are blended:

```
  ρ = α · ρ_water + (1-α) · ρ_air
  ν = α · ν_water + (1-α) · ν_air
```

This is identical to Project 11 and Project 08.

### 4.2 6-DoF Rigid Body Motion — The New Ingredient

A rigid body in 3D space has **six degrees of freedom**:

```
  Translation (3 DoF):          Rotation (3 DoF):
  ┌────────────────────┐       ┌─────────────────────┐
  │  Surge  (x)        │       │  Roll    (about x)  │
  │  Sway   (y)        │       │  Pitch   (about y)  │
  │  Heave  (z)        │       │  Yaw     (about z)  │
  └────────────────────┘       └─────────────────────┘
```

```
                    Yaw ↺ (z-axis)
                        │
                        │     Pitch ↺ (y-axis)
                        │    ╱
                        │   ╱
              ┌─────────┼──╱──────────┐
             ╱│         │ ╱          ╱│
            ╱ │         │╱  BOAT    ╱ │
           ╱  │         ┼─────────╱──→ Surge (x)
          ┌───│─────────│────────┐   │
          │   │         │        │   │
          │   └─────────│────────│───┘
          │  ╱     Roll ↺        │  ╱
          │ ╱      (x-axis)      │ ╱
          │╱                     │╱
          └──────────────────────┘
                  │
                  ↓ Heave (z)
```

**Newton's Second Law (translation):**

```
  F = m · a

  where:
    F = total force on the body (pressure + viscous + gravity)
    m = body mass
    a = acceleration of the centre of mass
```

The fluid exerts pressure forces on every face of the boat patch. Gravity
pulls the boat down. The net force determines the translational acceleration.

**Euler's Equations (rotation):**

```
  τ = I · α + ω × (I · ω)

  where:
    τ = total torque about the centre of mass
    I = moment of inertia tensor (diagonal for symmetric bodies)
    α = angular acceleration
    ω = angular velocity
```

For our small boat, the cross-coupling term `ω × (I · ω)` is very small
since angular velocities are low. But OpenFOAM solves the full equations.

### 4.3 Coupling: Fluid Forces ↔ Body Motion

Each timestep, the solver executes this coupling loop:

```
  ┌─────────────────────────────────────────────────────────────────┐
  │                     PIMPLE outer loop                           │
  │                                                                 │
  │  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐       │
  │  │ Solve p, U   │───→│ Integrate p  │───→│ 6-DoF solver │       │
  │  │ and α on     │    │ and τ over   │    │ F = ma       │       │
  │  │ current mesh │    │ boat patches │    │ τ = Iα       │       │
  │  └──────────────┘    └──────────────┘    └──────┬───────┘       │
  │                                                  │              │
  │                                          new position           │
  │                                          new rotation           │
  │                                                  │              │
  │                                          ┌───────▼──────┐       │
  │                                          │  Move mesh   │       │
  │                                          │  (morph)     │       │
  │                                          └───────┬──────┘       │
  │                                                  │              │
  │                                          back to top            │
  │                                          (next outer corr.)     │
  └─────────────────────────────────────────────────────────────────┘
```

The setting `moveMeshOuterCorrectors yes` in `fvSolution` enables this
re-meshing within the PIMPLE loop, which is essential for accurate coupling.

---

## 5. Dynamic Mesh — How the Mesh Moves

### Mesh Morphing with innerDistance / outerDistance

When the boat moves, the mesh must follow. OpenFOAM uses a **Laplacian
smoothing** approach to redistribute mesh points:

```
  ∇·(γ ∇ cellDisplacement) = 0
```

where γ is a diffusivity field that controls how much each cell moves.

The `innerDistance` and `outerDistance` parameters define three zones:

```
                     outerDistance = 0.08 m
               ┌───────────────────────────────┐
               │                               │
               │   innerDistance = 0.01 m      │
               │   ┌───────────────────┐       │
               │   │                   │       │
               │   │   ┏━━━━━━━━━━━┓   │       │
               │   │   ┃   BOAT    ┃   │       │
               │   │   ┃  (rigid)  ┃   │       │
               │   │   ┗━━━━━━━━━━━┛   │       │
               │   │                   │       │
               │   │  ZONE A: rigid    │       │
               │   │  (moves with body)│       │
               │   └───────────────────┘       │
               │                               │
               │  ZONE B: blended              │
               │  (smooth interpolation)       │
               │                               │
               └───────────────────────────────┘

               ZONE C: fixed (no displacement)
```

| Zone   | Distance from body         | Displacement                  |
|--------|---------------------------|-------------------------------|
| **A**  | < innerDistance (0.01 m)   | Full body displacement (rigid) |
| **B**  | innerDistance → outerDistance | Smoothly blended to zero     |
| **C**  | > outerDistance (0.08 m)   | Zero — mesh stays fixed        |

This means:
- Cells very close to the boat move exactly with it (rigid motion)
- Cells in the transition zone deform smoothly
- Cells far away don't move at all
- The domain boundaries stay fixed

### Why Mesh Quality Can Degrade

As the body moves, cells in Zone B get stretched and compressed. If the body
moves too far in a single timestep, cells can become:
- **Inverted** (negative volume) — solver crashes
- **Highly skewed** — poor accuracy
- **Very thin** — tiny cell volumes → tiny timesteps

This is why:
- Timesteps must be very small (5e-5 s)
- `accelerationRelaxation` and `accelerationDamping` are used to limit
  sudden body motion
- The Courant number limit is lower (maxCo 0.4) than for a static mesh

---

## 6. Boat Geometry and Physical Properties

### STL Geometry

The hull geometry comes from `boat_hull.stl` at the repository root:

```
  File:           boat_hull.stl (binary, 960 triangles)
  Native units:   millimeters (scaled to meters by Allrun)

  Bounding box (mm):
    X: [-77.53, 77.53]     Length ≈ 155 mm
    Y: [-40.85, 40.85]     Width  ≈  82 mm
    Z: [-27.25, 27.25]     Height ≈  54 mm

  Bounding box (m, after scaling):
    X: [-0.0775, 0.0775]   Length ≈ 0.155 m
    Y: [-0.0409, 0.0409]   Width  ≈ 0.082 m
    Z: [-0.0273, 0.0273]   Height ≈ 0.054 m

  Centre:  (0, 0, 0) — symmetric about all three axes
```

```
          Top View (X-Y plane)                Side View (X-Z plane)

         Y                                   Z
         ↑                                   ↑  0.027
         |   ╭───────────────╮               |   ╭──────────────╮
   0.041─|──╱                 ╲──           ─|──╱                ╲
         | ╱     BOAT HULL     ╲             | ╱   BOAT HULL      ╲
         |╱     (top view)      ╲            |╱   (side view)      ╲
  ───────┼───────────────────────→ X    ─────┼──────────────────────→ X
         |╲                     ╱            |╲                    ╱
         | ╲                   ╱             | ╲                  ╱
  -0.041─|──╲                 ╱──     -0.027─|──╲                ╱
         |   ╰───────────────╯               |   ╰──────────────╯
         |                                   |
       -0.078              0.078           -0.078              0.078
```

### Physical Properties (3D-printed PLA toy boat)

| Property              | Value                         | Units         |
|-----------------------|-------------------------------|---------------|
| Mass                  | 0.050                         | kg            |
| Centre of mass        | (0, 0, −0.005)                | m             |
| Ixx (roll)            | 5.0 × 10⁻⁶                   | kg·m²         |
| Iyy (pitch)           | 1.5 × 10⁻⁵                   | kg·m²         |
| Izz (yaw)             | 1.5 × 10⁻⁵                   | kg·m²         |

The centre of mass is 5 mm below the geometric centre. This is physically
reasonable — the bottom of a boat hull is thicker than the top, and the
weight distribution naturally sits low, which improves roll stability.

### Why These Inertia Values?

For a uniform rectangular solid with mass m, length L, width W, height H:

```
  Ixx = m/12 · (W² + H²) = 0.05/12 · (0.082² + 0.054²) ≈ 4.0e-05
  Iyy = m/12 · (L² + H²) = 0.05/12 · (0.155² + 0.054²) ≈ 1.1e-04
  Izz = m/12 · (L² + W²) = 0.05/12 · (0.155² + 0.082²) ≈ 1.3e-04
```

But a boat hull is *not* a solid rectangular block — it's a hollow shell with
most of its mass near the waterline. The actual moments of inertia are much
smaller (roughly 10–20× less) because mass is concentrated near the centre.
Our values (5e-6, 1.5e-5, 1.5e-5) reflect this.

---

## 7. Domain and Mesh Design

### Domain Dimensions

```
           atmosphere (top, z = 0.15 m)
  ┌──────────────────────────────────────────┐
  │                AIR                       │
  │                                          │
  │          ┏━━━━━━━━━┓                     │
  │ inlet    ┃  BOAT   ┃              outlet │
  │──────────┃─────────┃────────────────────→│
  │          ┗━━━━━━━━━┛   waterline z=0     │
  │                                          │
  │               WATER                      │
  │                                          │
  └──────────────────────────────────────────┘
           bottom (z = -0.20 m)

  x: [-0.5, 1.0]   →  1.5 m total  (≈ 10× boat length)
  y: [-0.25, 0.25]  →  0.5 m total  (≈ 6× boat width)
  z: [-0.20, 0.15]  →  0.35 m total

  Water depth: 0.20 m (below waterline)
  Air height:  0.15 m (above waterline)
```

The domain is slightly deeper than Project 11 (0.20 m vs 0.15 m) to give the
boat room to bob up and down without the bottom boundary interfering.

### Background Mesh

```
  Cells:  60 × 20 × 14  =  16,800 background cells
  After snappyHexMesh:    ~100,000–200,000 cells (with refinement)

  Refinement strategy:
  - Level 1: refinementBox around the boat  (±0.15 m in x, ±0.08 m in y/z)
  - Level 2: waterlineBox at z = 0 ± 0.02 m (captures free surface)
  - Level 2-3: on the hull surface itself
  - Level 3: along hull feature edges
```

### Why the Waterline Refinement Box?

The free surface at z ≈ 0 is where all the action happens: waves, splashing,
the hull entering and leaving the water. Fine cells here improve:
- Interface sharpness
- Force calculation accuracy (buoyancy depends on knowing exactly where the
  waterline intersects the hull)
- Wave capture

---

## 8. File-by-File Walkthrough

### 8.1 `0/U`

**Velocity field** — boundary conditions for fluid velocity.

```
  inlet:       fixedValue (0 0 0)          Calm water, no incoming current
  outlet:      pressureInletOutletVelocity  Allows outflow, prevents inflow
  atmosphere:  pressureInletOutletVelocity  Air can flow in/out at top
  bottom:      slip                         Free-slip wall (reduces drag)
  sides:       slip                         Free-slip (symmetric-like)
  boat:        movingWallVelocity           ← THE KEY DIFFERENCE
```

**Why `movingWallVelocity` on the boat?**

In Project 11, the boat was fixed, so we used `noSlip` (velocity = 0 at the
wall). Here, the boat *moves*. The `movingWallVelocity` BC ensures the fluid
velocity at the boat surface matches the hull's velocity. This is a no-slip
condition *in the reference frame of the moving body*:

```
  U_fluid(at wall) = U_body(at wall)
```

If you used `noSlip` instead, the fluid would try to stay at rest while the
mesh moves — the simulation would diverge immediately.

### 8.2 `0/p_rgh`

**Reduced pressure** (p - ρgh).

Using p_rgh instead of p removes the hydrostatic component, which improves
numerical accuracy for free-surface flows.

```
  inlet:       fixedFluxPressure       Pressure adjusts to give correct flux
  outlet:      totalPressure (p0=0)    Reference pressure at outlet
  atmosphere:  totalPressure (p0=0)    Reference pressure at top
  bottom:      fixedFluxPressure       No normal flow through bottom
  sides:       fixedFluxPressure       No normal flow through sides
  boat:        fixedFluxPressure       Pressure adjusts for body forces
```

### 8.3 `0/alpha.water`

**Phase fraction** — identifies where water and air are.

Initially set to 0 everywhere (air). The `setFields` utility fills
everything below z = 0 with α = 1 (water).

```
  inlet:       inletOutlet (inletValue 1)  Water enters from below waterline
  outlet:      inletOutlet (inletValue 0)  Defaults to air if flow reverses
  atmosphere:  inletOutlet (inletValue 0)  Air at top
  bottom:      zeroGradient                No phase gradient at wall
  sides:       zeroGradient                No phase gradient at wall
  boat:        zeroGradient                Phase slides along hull
```

### 8.4 `0/pointDisplacement`

**Point displacement field** — tells the mesh motion solver how much each
mesh point has moved from its original position.

This file is **new in Project 12** — it doesn't exist in Project 11 because
the mesh was static.

```
  class:   pointVectorField    (defined at mesh POINTS, not cell centres)
  object:  pointDisplacement

  boat:        calculated           6-DoF solver provides the displacement
  all others:  fixedValue (0 0 0)   Domain boundaries don't move
```

**Key point:** The `class` must be `pointVectorField` (not `volVectorField`).
Mesh motion operates on *points* (vertices), not cell centres. The `calculated`
type on the boat means "the motion solver will fill in this value."

### 8.5 `constant/transportProperties`

Two-phase fluid properties — identical to Project 11.

```
  Water:  ρ = 998.2 kg/m³,   ν = 1.0 × 10⁻⁶ m²/s
  Air:    ρ = 1.225 kg/m³,   ν = 1.48 × 10⁻⁵ m²/s
  Surface tension: σ = 0.072 N/m
```

### 8.6 `constant/turbulenceProperties`

```
  simulationType  laminar;
```

Laminar flow — appropriate for a small toy boat at low speed. Re ≈ UL/ν =
0.2 × 0.155 / 1e-6 ≈ 31,000 — borderline, but laminar is acceptable for
a learning case.

### 8.7 `constant/g`

Gravity vector — points downward in the z-direction.

```
  dimensions  [0 1 -2 0 0 0 0];
  value       (0 0 -9.81);
```

This is critical for buoyancy. The VOF solver uses gravity to compute
hydrostatic pressure gradients, and the 6-DoF solver uses it for the body's
weight force.

### 8.8 `constant/dynamicMeshDict` — The Key File

This is the **most important file in the project**. It configures the dynamic
mesh motion and the 6-DoF rigid body solver.

```
  dynamicFvMesh       dynamicMotionSolverFvMesh;
  motionSolverLibs    ("libsixDoFRigidBodyMotion.so");
  motionSolver        sixDoFRigidBodyMotion;
```

#### Top-Level Settings

| Setting                | Value                              | Purpose |
|------------------------|------------------------------------|---------|
| `dynamicFvMesh`        | `dynamicMotionSolverFvMesh`       | Use a motion-solver-based dynamic mesh |
| `motionSolverLibs`     | `libsixDoFRigidBodyMotion.so`     | Load the 6-DoF library |
| `motionSolver`         | `sixDoFRigidBodyMotion`           | Use 6-DoF solver as the motion solver |

#### sixDoFRigidBodyMotionCoeffs — Body Properties

| Parameter              | Value                  | Meaning |
|------------------------|------------------------|---------|
| `patches`              | `(boat)`               | Which patches form the rigid body |
| `mass`                 | `0.050`                | Body mass in kg |
| `centreOfMass`         | `(0 0 -0.005)`         | CoM position in m |
| `momentOfInertia`      | `(5e-06 1.5e-05 1.5e-05)` | Ixx, Iyy, Izz in kg·m² |
| `report`               | `on`                   | Print position/rotation each timestep |

#### sixDoFRigidBodyMotionCoeffs — Mesh Motion

| Parameter              | Value   | Meaning |
|------------------------|---------|---------|
| `innerDistance`         | `0.01`  | Cells within 1 cm of body move rigidly |
| `outerDistance`         | `0.08`  | Cells beyond 8 cm don't move at all |

#### sixDoFRigidBodyMotionCoeffs — Solver

```
  solver { type Newmark; }
```

The **Newmark** time integration scheme is used to advance the body's position
and velocity. It is unconditionally stable and second-order accurate. OpenFOAM
also supports `CrankNicolson`, `symplectic`, and others.

#### sixDoFRigidBodyMotionCoeffs — Stability Controls

| Parameter                    | Value  | Purpose |
|------------------------------|--------|---------|
| `accelerationRelaxation`     | `0.7`  | Under-relax acceleration (prevents overshooting) |
| `accelerationDamping`        | `0.9`  | Damp oscillations in acceleration |

These are critical for stability. Without them, the coupling between fluid
forces and body motion can oscillate and diverge ("added mass instability").

#### Constraints — Restricting Degrees of Freedom

```
  fixedLine
  {
      sixDoFRigidBodyMotionConstraint line;
      direction (0 0 1);
  }
```

This constrains the body's translation to a **vertical line** (z-axis only).
The boat can heave (move up/down) but cannot surge (x) or sway (y).

```
  fixedAxis
  {
      sixDoFRigidBodyMotionConstraint axis;
      axis (0 0 1);
  }
```

This constrains rotation about a **fixed axis** — specifically, it *suppresses*
yaw rotation (rotation about z). The boat can still pitch (about y) and roll
(about x).

**Summary of allowed/suppressed DoFs:**

```
  ┌─────────────┬───────────┬─────────────────────────────────┐
  │  Degree     │  Status   │  Reason                         │
  ├─────────────┼───────────┼─────────────────────────────────┤
  │  Surge (x)  │ LOCKED    │  line constraint (z-direction)  │
  │  Sway  (y)  │ LOCKED    │  line constraint (z-direction)  │
  │  Heave (z)  │ FREE  ✓   │  allowed by line direction      │
  │  Roll  (x)  │ FREE  ✓   │  not constrained                │
  │  Pitch (y)  │ FREE  ✓   │  not constrained                │
  │  Yaw   (z)  │ LOCKED    │  axis constraint (z-axis)       │
  └─────────────┴───────────┴─────────────────────────────────┘
```

This is a good starting configuration:
- Heave is the primary motion we want to observe (boat finding its draft)
- Pitch may occur if forces are asymmetric or if waves develop
- Roll is allowed for physical realism
- Surge/sway/yaw are locked for stability in initial testing

### 8.9 `system/blockMeshDict`

Creates the background hexahedral mesh: 60 × 20 × 14 cells.

```
  Domain: [-0.5, 1.0] × [-0.25, 0.25] × [-0.20, 0.15] m
  Cell size: ~25 mm × 25 mm × 25 mm (background)
```

Patches:
```
  inlet (x = -0.5)    — left face
  outlet (x = 1.0)    — right face
  atmosphere (z = 0.15) — top face
  bottom (z = -0.20)  — bottom face
  sides                — front and back (y = ±0.25)
```

### 8.10 `system/controlDict`

| Setting            | Value          | Reason |
|--------------------|----------------|--------|
| `application`      | `interDyMFoam` | Dynamic mesh + VOF |
| `endTime`          | `3`            | 3 seconds for boat to settle |
| `deltaT`           | `5e-05`        | Very small initial timestep |
| `adjustTimeStep`   | `yes`          | Adaptive based on Courant number |
| `maxCo`            | `0.4`          | Lower than typical (0.5) for mesh stability |
| `maxAlphaCo`       | `0.25`         | Lower limit near the interface |
| `maxDeltaT`        | `0.005`        | Cap on timestep growth |
| `writeInterval`    | `0.05`         | Write every 50 ms → 60 frames |

The `libs` section loads the 6-DoF library:
```
  libs ("libsixDoFRigidBodyMotion.so");
```

**Function objects:**
- `sixDoFRigidBodyState` — writes body position, rotation, velocity each step
- `forces` — computes and logs forces/moments on the boat patch

### 8.11 `system/fvSchemes`

Standard interFoam schemes plus mesh motion:

```
  ddtSchemes:           Euler (first-order, robust)
  gradSchemes:          Gauss linear
  divSchemes:           linearUpwind for momentum, vanLeer for alpha
  laplacianSchemes:     Gauss linear corrected (default)
                        + cellDisplacement entry for mesh motion
  fluxRequired:         pcorr, p_rgh, alpha.water
```

The **cellDisplacement** laplacian scheme is new — it controls the mesh motion
equation `∇·(γ ∇ displacement) = 0`.

### 8.12 `system/fvSolution`

Same as Project 11 plus dynamic mesh additions:

```
  New solver:  cellDisplacement  (PCG + DIC preconditioner)
  PIMPLE additions:
    moveMeshOuterCorrectors  yes    ← re-mesh within outer PIMPLE loop
    correctPhi               yes    ← correct fluxes after mesh motion
```

**Why `moveMeshOuterCorrectors yes`?**

With `nOuterCorrectors 2`, the solver does 2 outer PIMPLE iterations. With
this flag on, the mesh is re-positioned at each outer iteration — giving a
tighter coupling between fluid forces and body motion within a single timestep.
Without it, the mesh would only move once per timestep, which can be
insufficiently accurate.

### 8.13 `system/setFieldsDict`

Initialises the water region:

```
  defaultFieldValues:  alpha.water = 0  (air everywhere)
  boxToCell:           everything below z = 0 → alpha.water = 1 (water)
```

The box `(-10 -10 -10) (10 10 0)` is intentionally huge — it fills *all*
cells below the waterline, regardless of domain size.

After `setFields`, the boat hull sits with its centre at the waterline:
- Bottom half submerged in water (α = 1)
- Top half in air (α = 0)

### 8.14 `system/surfaceFeatureExtractDict`

Extracts sharp edges from the STL for better snapping:

```
  boat_hull.stl:
    extractionMethod:  extractFromSurface
    includedAngle:     150°
```

Edges where the surface normal changes by more than 30° (180 - 150) are
marked as "features" and get extra refinement.

### 8.15 `system/snappyHexMeshDict`

Refines and snaps the background mesh to the hull surface.

**Refinement regions:**

| Region           | Level | Purpose                           |
|------------------|-------|-----------------------------------|
| Hull surface     | 2–3   | Sharp hull geometry               |
| Feature edges    | 3     | Crisp edges (bow, stern, keel)    |
| refinementBox    | 1     | Near-body region                  |
| waterlineBox     | 2     | Free-surface resolution           |

```
  locationInMesh: (0.3 0 0.05) — a point in the air, outside the hull
```

This tells snappyHexMesh which side of the hull to keep (the fluid region,
not the interior of the boat).

---

## 9. Running the Simulation

### Quick Start

```bash
cd projects/12_boat_hull_floating
./Allrun
```

### Step-by-Step (manual)

```bash
# 1. Source OpenFOAM environment (if not already done)
source /opt/openfoam11/etc/bashrc    # adjust path for your installation

# 2. Clean any previous run
./Allclean

# 3. Copy and scale the STL
cp ../../boat_hull.stl constant/triSurface/boat_hull.stl
surfaceTransformPoints -scale '(0.001 0.001 0.001)' \
    constant/triSurface/boat_hull.stl \
    constant/triSurface/boat_hull.stl

# 4. Generate background mesh
blockMesh

# 5. Extract surface features
surfaceFeatureExtract

# 6. Snap mesh to hull
snappyHexMesh -overwrite

# 7. Check mesh quality
checkMesh

# 8. Initialise water phase
setFields

# 9. Run the solver
interDyMFoam

# 10. Post-process
paraFoam
```

### Expected Runtime

This case is computationally expensive. Expect:

| Machine                | Estimated runtime |
|------------------------|-------------------|
| Laptop (4 cores)       | 2–6 hours         |
| Workstation (8+ cores) | 1–3 hours         |

For parallel runs:

```bash
# Decompose
decomposePar

# Run on 4 cores
mpirun -np 4 interDyMFoam -parallel

# Reconstruct
reconstructPar
```

You will need a `system/decomposeParDict` for parallel runs (not included by
default — see Exercise 5).

### Monitoring Progress

While the simulation runs, you can watch the 6-DoF output:

```bash
# Watch body position and rotation in real time
tail -f log.interDyMFoam | grep -E "Centre of mass|Orientation"

# Or check the forces file
tail -f postProcessing/forces/0/forces.dat
```

---

## 10. Post-Processing

### 10.1 Tracking Body Motion

The `sixDoFRigidBodyState` function object writes motion data to:

```
  postProcessing/sixDoFRigidBodyState/0/
```

This includes:
- Centre of mass position (x, y, z) over time
- Orientation (Euler angles or quaternion)
- Linear and angular velocities

**Plotting heave over time:**

```bash
# Extract time and z-position from the log
grep "Centre of mass" log.interDyMFoam | awk '{print $NF}' > heave.dat
```

Or in ParaView, load the motion data and plot z-displacement vs time.

### 10.2 Visualising in ParaView

```bash
paraFoam
```

Key visualisations:

1. **Water surface:** Apply a Contour filter on `alpha.water = 0.5`
2. **Hull position:** The boat patch moves with time — scrub the timeline
3. **Pressure on hull:** Colour the boat surface by `p_rgh`
4. **Velocity field:** Clip to the centreline (y = 0) and colour by |U|

**Tip:** Use *File → Save Animation* to create a video of the boat settling.

### 10.3 Force Analysis

The `forces` function object writes to:
```
  postProcessing/forces/0/forces.dat
```

Columns include pressure force, viscous force, and porous force (zero here)
in x, y, z directions. The z-component of the pressure force is the
**buoyancy force** acting on the hull.

At equilibrium:
```
  F_buoyancy (z)  ≈  m × g  =  0.050 × 9.81  =  0.4905 N
```

---

## 11. Expected Results

### Equilibrium Draft Calculation

By Archimedes' principle, the boat floats when:

```
  Weight = Buoyancy
  m · g  = ρ_water · g · V_submerged
  m      = ρ_water · V_submerged

  V_submerged = m / ρ_water = 0.050 / 998.2 ≈ 5.01 × 10⁻⁵ m³
```

For a rough estimate, treating the waterplane area as a rectangle:

```
  A_waterplane ≈ L × W = 0.155 × 0.082 = 0.01271 m²
  draft ≈ V_submerged / A_waterplane = 5.01e-5 / 0.01271 ≈ 0.0039 m ≈ 4 mm
```

So the boat should settle with about **4 mm of hull below the waterline**.

Given the hull height is ~54 mm, only about 7% of the hull is submerged —
this makes sense for a lightweight hollow PLA boat.

### Expected Motion Timeline

```
  Time 0.0 s     Boat at initial position (z_CoM = -0.005 m)
                  Partially submerged (geometric centre at waterline)

  Time 0-0.5 s   Gravity pulls boat down, buoyancy pushes up
                  Oscillatory heave begins (damped by fluid)

  Time 0.5-1.5 s Heave oscillations decay
                  Small pitch oscillations may develop

  Time 1.5-3.0 s Boat reaches steady state
                  z_CoM ≈ -0.005 + adjustment for 4mm draft
                  Pitch ≈ 0° (symmetric hull)
```

```
  Heave (z) vs Time — Expected Behaviour

  z ↑
    │   ╭─╮
    │  ╱   ╲     ╭╮
    │ ╱     ╲   ╱  ╲   ╭╮
    │╱       ╲ ╱    ╲ ╱  ╲─────── equilibrium draft
  ──┼─────────╳──────╳────────→ time
    │          ╲   ╱  ╲╱
    │           ╲ ╱
    │            ╰╯
    │
    │  ← initial drop   ← damped oscillations   ← steady state →
```

### What You Should See

1. **The boat drops slightly** from its initial position
2. **Overshoots** below equilibrium (too much buoyancy → bounces back up)
3. **Damped oscillation** — each bob is smaller than the last
4. **Settles** at the equilibrium draft of ~4 mm
5. **Pitch stays near zero** (symmetric hull in calm water)
6. **Roll stays near zero** (symmetric hull, no lateral forces)

---

## 12. Stability and Troubleshooting

### Why Small Timesteps Matter Even More

With dynamic mesh, three timestep constraints must be satisfied simultaneously:

```
  1. Courant number:    Co = U · Δt / Δx  <  maxCo (0.4)
  2. Interface Courant:  Co_α            <  maxAlphaCo (0.25)
  3. Mesh quality:       No inverted cells after mesh motion
```

Constraint 3 is the new one. If the body moves too far in one step, cells
near it can invert. The `adjustTimeStep` mechanism handles constraints 1 and 2
but does NOT directly monitor mesh quality.

**Rule of thumb:** Start with a very small δt (5e-5) and let the adaptive
timestepping find the right value. Don't increase `maxDeltaT` beyond 0.005.

### Acceleration Relaxation and Damping

These parameters in `dynamicMeshDict` prevent the "added mass instability":

```
  accelerationRelaxation  0.7    (0 = no relaxation, 1 = full relaxation)
  accelerationDamping     0.9    (higher = more damping)
```

**What is added mass instability?**

When the body accelerates, it must also accelerate the surrounding fluid.
This "added mass" creates a feedback loop:

```
  Body accelerates → fluid forces change → body accelerates differently →
  fluid forces change again → oscillation grows → DIVERGENCE
```

Relaxation and damping break this loop by limiting how much the acceleration
can change between iterations.

### Common Problems and Solutions

#### Problem: "Mesh quality check failed" or negative cell volumes

```
  Cause:  Body moved too far in one timestep
  Fix:    Reduce maxCo to 0.2, reduce maxDeltaT to 0.001
          Increase nOuterCorrectors to 3
```

#### Problem: Boat flies off to infinity

```
  Cause:  Incorrect mass, inertia, or constraint setup
  Fix:    Check units (mass in kg, not grams!)
          Verify constraints are correct
          Ensure accelerationRelaxation < 1.0
```

#### Problem: Boat doesn't float (sinks through bottom)

```
  Cause:  Mass too high, or mesh too coarse to resolve pressure forces
  Fix:    Check mass is realistic (0.050 kg for a toy boat)
          Increase mesh refinement around hull
          Verify alpha.water is correctly initialised below waterline
```

#### Problem: Boat doesn't move at all

```
  Cause:  dynamicMeshDict not being read, or wrong patch names
  Fix:    Check log for "sixDoFRigidBodyMotion" messages
          Verify patch name 'boat' matches snappyHexMeshDict
          Make sure libs entry is in controlDict
```

#### Problem: Very slow simulation

```
  Cause:  Timestep is adapting to very small values
  Fix:    Check log for current deltaT
          Coarsen mesh if possible
          Run in parallel (decomposePar + mpirun)
```

#### Problem: Water phase (alpha) becomes noisy near hull

```
  Cause:  Poor mesh quality near moving boundary
  Fix:    Increase nAlphaSubCycles to 2
          Reduce maxAlphaCo to 0.15
          Improve snappyHexMesh settings (more smoothing)
```

---

## 13. Exercises

### Exercise 1: Calm Water Float Test → Current Flow

The default setup starts with zero velocity everywhere (calm water). After the
boat has settled to its equilibrium draft (~1.5 s), try adding a gentle current:

1. Change `0/U` inlet to `fixedValue uniform (0.1 0 0)` — 10 cm/s current
2. Restart from the last written timestep
3. Observe: Does the boat pitch? Does the draft change? Are there waves?

Compare the pressure distribution on the hull with the fixed-hull result from
Project 11.

### Exercise 2: Heavy Boat vs Light Boat

Change the mass in `dynamicMeshDict`:

| Case    | Mass (kg) | Expected draft (mm) | Expected behaviour |
|---------|-----------|---------------------|--------------------|
| Light   | 0.025     | ~2 mm               | Floats higher, faster oscillations |
| Default | 0.050     | ~4 mm               | Baseline |
| Heavy   | 0.100     | ~8 mm               | Sits lower, slower oscillations |
| Too heavy | 0.200   | ~16 mm              | May partially submerge |

Plot heave vs time for each case. How does the oscillation frequency change
with mass?

### Exercise 3: Release from Above the Water

Instead of starting the boat at the waterline, release it from above:

1. In `dynamicMeshDict`, change `centreOfMass` to `(0 0 0.03)` — 30 mm above
   water
2. The boat should fall, splash into the water, and settle

This is a much more dramatic simulation but requires finer timesteps and may
need more mesh motion iterations.

### Exercise 4: Free All Six Degrees of Freedom

Remove the constraints from `dynamicMeshDict`:

```
constraints { }
```

Now the boat can surge, sway, and yaw freely. What happens? Is the simulation
stable? Does the boat drift?

**Warning:** This may be unstable without additional damping or restraints.
Try adding a `linearDamper` restraint to limit runaway drift.

### Exercise 5: Parallel Decomposition

Create `system/decomposeParDict` for parallel execution:

```
numberOfSubdomains  4;
method              scotch;
```

Run with:
```bash
decomposePar
mpirun -np 4 interDyMFoam -parallel
reconstructPar
```

Compare runtime vs serial. Note: parallel dynamic mesh can have load balancing
challenges.

### Exercise 6: Turbulence Model

Switch from laminar to a RANS turbulence model:

1. Change `constant/turbulenceProperties` to use `kOmegaSST`
2. Add `0/k`, `0/omega`, `0/nut` boundary condition files
3. Use appropriate wall functions on the boat patch
4. Re-run and compare forces on the hull

Does turbulence modelling change the equilibrium draft? The oscillation
damping rate?

---

## 14. Cross-References

### Related Notes

| Note  | Topic                              | Relevance |
|-------|------------------------------------|-----------|
| 01    | OpenFOAM directory structure       | Basic case structure applies here |
| 04    | blockMesh                          | Background mesh generation |
| 05    | snappyHexMesh                      | Hull surface meshing |
| 08    | interFoam / VOF method             | Two-phase physics fundamentals |

### Related Projects

| Project | Title                    | Relevance |
|---------|--------------------------|-----------|
| 08      | Dam Break                | VOF basics, free-surface initialisation |
| 11      | Fixed Boat Hull Flow     | **Direct prerequisite** — same geometry, static mesh |
| 12      | **This project**         | Floating hull with 6-DoF motion |

### Progression Path

```
  ┌──────────────────────────────────────────────────────────────────┐
  │                                                                  │
  │  Project 08        Project 11          Project 12                │
  │  ┌──────────┐     ┌──────────────┐    ┌────────────────────┐     │
  │  │Dam Break │ ──→ │Fixed Hull    │ ──→│Floating Hull       │     │
  │  │(VOF)     │     │(VOF + STL)   │    │(VOF + STL + 6-DoF) │     │
  │  └──────────┘     └──────────────┘    └────────────────────┘     │
  │                                                                  │
  │  Learn VOF         Learn STL meshing   Learn dynamic mesh        │
  │  Learn setFields   Learn snappyHex     Learn rigid body motion   │
  │                    Learn interFoam     Learn interDyMFoam        │
  │                                                                  │
  └──────────────────────────────────────────────────────────────────┘
```

### Key OpenFOAM Tutorials to Compare

The official OpenFOAM tutorials contain similar cases:

```
  $FOAM_TUTORIALS/multiphase/interDyMFoam/laminar/floatingObject/
```

This tutorial simulates a floating box (simple geometry) and is an excellent
reference for the `dynamicMeshDict` syntax and solver settings.

---

## Appendix: Complete File Listing

```
projects/12_boat_hull_floating/
├── 0/
│   ├── U                          Velocity BCs (movingWallVelocity on boat)
│   ├── p_rgh                      Reduced pressure BCs
│   ├── alpha.water                Phase fraction BCs
│   └── pointDisplacement          Mesh displacement BCs (pointVectorField)
├── constant/
│   ├── dynamicMeshDict            ★ 6-DoF rigid body + mesh motion config
│   ├── g                          Gravity vector
│   ├── transportProperties        Water/air properties
│   ├── turbulenceProperties       Laminar
│   └── triSurface/                (STL copied here by Allrun)
├── system/
│   ├── blockMeshDict              Background mesh (60×20×14)
│   ├── controlDict                interDyMFoam settings
│   ├── fvSchemes                  Discretisation + mesh motion schemes
│   ├── fvSolution                 Solvers + PIMPLE + mesh motion solver
│   ├── setFieldsDict              Water initialisation (z < 0)
│   ├── snappyHexMeshDict          Hull mesh refinement
│   └── surfaceFeatureExtractDict  Edge extraction
├── Allrun                         Run script (copy STL → mesh → solve)
├── Allclean                       Clean script
└── README.md                      This file
```

---

*This project is part of the OpenFOAM Tutorials series.*
*For questions or improvements, see the repository README.*
