# Tutorial 09: Wind Turbine Blade Simulation with MRF

**Difficulty:** ★★★★☆ Advanced
**Solver:** `simpleFoam` (steady-state, incompressible, turbulent with MRF)
**Turbulence Model:** k-epsilon (RAS)
**Approach:** Multiple Reference Frame (Frozen Rotor)
**Estimated Run Time:** ~5-10 minutes on a modern workstation

---

## Table of Contents

1. [Problem Description](#1-problem-description)
2. [Physical Background — MRF Theory](#2-physical-background--mrf-theory)
3. [MRF vs AMI: Choosing the Right Approach](#3-mrf-vs-ami-choosing-the-right-approach)
4. [Case Structure](#4-case-structure)
5. [Mesh Generation — blockMesh](#5-mesh-generation--blockmesh)
6. [Cell Zone Creation — topoSet](#6-cell-zone-creation--toposet)
7. [MRFProperties — The Key Configuration](#7-mrfproperties--the-key-configuration)
8. [Boundary Conditions](#8-boundary-conditions)
9. [Turbulence Setup](#9-turbulence-setup)
10. [Solver Configuration](#10-solver-configuration)
11. [How to Run](#11-how-to-run)
12. [Expected Results](#12-expected-results)
13. [Post-Processing](#13-post-processing)
14. [Exercises](#14-exercises)
15. [Industrial Applications](#15-industrial-applications)
16. [Troubleshooting](#16-troubleshooting)
17. [References](#17-references)
18. [Cross-References](#18-cross-references)

---

## 1. Problem Description

This tutorial demonstrates the **Multiple Reference Frame (MRF)** technique in OpenFOAM,
a fundamental approach for simulating rotating machinery. Instead of physically rotating
the mesh (which is computationally expensive and requires transient solvers), MRF adds
source terms to the governing equations that approximate the effect of rotation.

We simulate wind flowing past a circular **rotating zone** that represents a simplified
turbine rotor disk. The MRF approach treats this zone as if it were spinning at 100 RPM,
introducing swirl into the flow without any mesh motion.

### Domain Layout

```
                         3.0 m
  ◄─────────────────────────────────────────────────►

  ▲  ┌──────────────────────────────────────────────┐  ▲
  │  │                                              │  │
  │  │                 Wind Direction                │  │
  │  │                  →→→→→→→→→                   │  │
  │  │                                              │  │
  │  │              ╭──────────────╮                 │  │
  │  │             ╱                ╲                │  │
     │  →→→→→→   │    Rotating      │   →→→→→→     │
 2.0 │  →→→→→→   │    MRF Zone      │   →→→→→→     │  2.0 m
  m  │  →→→→→→   │   ω = 100 RPM    │   →→→→→→     │
     │  →→→→→→   │    r = 0.3 m     │   →→→→→→     │
  │  │             ╲                ╱                │  │
  │  │              ╰──────────────╯                 │  │
  │  │            center: (1.0, 0.5)                 │  │
  │  │                                              │  │
  │  │                  →→→→→→→→→                   │  │
  │  │                  Wake Region                  │  │
  ▼  └──────────────────────────────────────────────┘  ▼

     inlet                                     outlet
     U = 10 m/s                                p = 0
```

### Domain Specifications

```
  ┌─────────────────────────────────────────────────────┐
  │  Physical Domain:  3.0 m × 2.0 m × 0.1 m (quasi-2D)│
  │  MRF Zone Center:  (1.0, 0.5, 0.0)                  │
  │  MRF Zone Radius:  0.3 m                             │
  │  MRF Zone Shape:   Cylinder (z-axis aligned)         │
  │  Rotation Speed:   100 RPM = 10.472 rad/s            │
  │  Rotation Axis:    z-axis (0, 0, 1)                  │
  │  Inlet Velocity:   10 m/s (uniform, x-direction)     │
  │  Mesh Resolution:  150 × 100 × 1 cells               │
  │  Reynolds Number:   ~200,000 (based on MRF diameter) │
  └─────────────────────────────────────────────────────┘
```

### Flow Physics Overview

```
  Upstream          MRF Zone                  Downstream
  (Undisturbed)     (Rotation Effects)        (Wake)

     →→→→          ╭───╮  ↗ ↗               →→→↗
     →→→→         ╱  ↻  ╲↗                  →→→→→
     →→→→  ───►  │  ω    │  ───►            →→→↘
     →→→→         ╲     ╱↘                  →→→→→
     →→→→          ╰───╯  ↘ ↘               →→→↗

  Uniform flow   Coriolis + centripetal     Swirling wake
  U = 10 m/s     source terms applied      with velocity
                  inside this region        deficit
```

---

## 2. Physical Background — MRF Theory

### What is MRF (Multiple Reference Frame)?

The **Multiple Reference Frame** method is a steady-state approximation for simulating
rotating machinery. Rather than physically rotating the mesh, MRF modifies the governing
equations inside a designated zone by adding source terms that account for the rotational
effects.

#### The Concept

In a rotating reference frame, an observer spinning with the rotor sees the flow as
steady (or approximately steady). The transformation between the inertial (stationary)
frame and the rotating frame introduces two additional body forces:

```
  ┌─────────────────────────────────────────────────────────────────┐
  │                                                                 │
  │  Navier-Stokes in Rotating Frame:                               │
  │                                                                 │
  │  ∂U_r                                                           │
  │  ──── + (U_r · ∇)U_r = -∇p/ρ + ν∇²U_r                        │
  │   ∂t                                                            │
  │                         - 2ω × U_r        ← Coriolis force     │
  │                         - ω × (ω × r)     ← Centripetal accel. │
  │                                                                 │
  │  where:                                                         │
  │    U_r  = velocity in the rotating frame                        │
  │    ω    = angular velocity vector                               │
  │    r    = position vector from rotation axis                    │
  │                                                                 │
  └─────────────────────────────────────────────────────────────────┘
```

#### Coriolis Force: `-2ω × U_r`

The Coriolis force acts perpendicular to both the rotation axis and the relative velocity.
It is responsible for deflecting the flow and creating the swirling motion that is
characteristic of rotating machinery.

```
  Rotation axis (ω)
        ↑  z
        │
        │    Coriolis force
        │   ──────►
        │  ╱
        │ ╱  Velocity U_r
        │╱──────────►
        O──────────────► x
```

#### Centripetal Acceleration: `-ω × (ω × r)`

The centripetal term always points radially inward toward the rotation axis. It creates
a pressure gradient that balances the centrifugal effect experienced by the fluid in the
rotating frame.

```
                    ω (rotation axis, out of page)
                    ⊙

         centripetal         centripetal
          ◄──────  ●  ──────►
                 center

    Points radially INWARD toward the rotation axis
```

#### Frozen Rotor Approximation

The "frozen rotor" name comes from the fact that the rotor geometry (if present) remains
**frozen** in one angular position. The MRF approach does NOT:

- Move or rotate any mesh cells
- Account for transient rotor-stator interaction
- Capture unsteady wake effects or vortex shedding from blades

The MRF approach DOES:

- Add Coriolis and centripetal source terms inside the rotating zone
- Produce a time-averaged (steady-state) flow field
- Provide a good initial estimate of performance (pressure rise, torque)
- Run much faster than transient sliding-mesh simulations

### When to Use MRF

MRF is appropriate when:
- The flow is approximately steady in the rotating frame
- Rotor-stator spacing is large relative to blade pitch
- You need a quick estimate of machine performance
- You want to initialize a transient AMI simulation

MRF is NOT appropriate when:
- Strong unsteady rotor-stator interaction exists
- You need accurate wake dynamics or vortex shedding
- Blade passing frequency effects are important
- The geometry has large circumferential non-uniformity near the interface

---

## 3. MRF vs AMI: Choosing the Right Approach

OpenFOAM provides two main methods for simulating rotating machinery:

```
  ┌─────────────────────┐              ┌─────────────────────┐
  │                     │              │                     │
  │   MRF (Frozen       │              │   AMI (Arbitrary    │
  │   Rotor)            │              │   Mesh Interface)   │
  │                     │              │                     │
  │  ┌───┐              │              │  ┌───┐ → ┌───┐     │
  │  │ ● │  No mesh     │              │  │ ↻ │   │   │     │
  │  │   │  motion      │              │  │   │   │   │     │
  │  └───┘              │              │  └───┘   └───┘     │
  │                     │              │  Mesh actually      │
  │  Source terms only   │              │  slides/rotates     │
  │                     │              │                     │
  └─────────────────────┘              └─────────────────────┘
```

### Comparison Table

| Feature                  | MRF (This Tutorial)      | AMI (Sliding Mesh)        |
|--------------------------|--------------------------|---------------------------|
| **Type**                 | Steady-state             | Transient                 |
| **Mesh Motion**          | None                     | Actual rotation           |
| **Solver**               | `simpleFoam`             | `pimpleFoam`              |
| **Computational Cost**   | Low (★★☆☆☆)             | High (★★★★★)             |
| **Setup Complexity**     | Moderate                 | Complex (AMI interfaces)  |
| **Physical Accuracy**    | Moderate                 | High                      |
| **Rotor-Stator Effects** | Not captured             | Fully captured            |
| **Blade Wake**           | Time-averaged            | Resolved per time step    |
| **Typical Use**          | Initial design, screening| Detailed analysis         |
| **Run Time**             | Minutes                  | Hours to days             |
| **Interface Treatment**  | cellZone (no interface)  | AMI patch pair            |
| **Multiple Zones**       | Yes (multiple MRF zones) | Yes (multiple AMI zones)  |

### Decision Flowchart

```
  Start
    │
    ▼
  Is the flow approximately ──── Yes ───► Use MRF
  steady in the rotating frame?            (simpleFoam)
    │
    No
    │
    ▼
  Do you need unsteady ──── Yes ───► Use AMI
  rotor-stator effects?               (pimpleFoam)
    │
    No
    │
    ▼
  Is this an initial ──── Yes ───► Use MRF first,
  design study?                     then AMI for final
    │
    No
    │
    ▼
  Use AMI for detailed analysis
```

---

## 4. Case Structure

```
09_wind_turbine_blade/
├── 0/                          ◄── Initial & boundary conditions
│   ├── U                       ◄── Velocity field
│   ├── p                       ◄── Pressure field (kinematic)
│   ├── k                       ◄── Turbulent kinetic energy
│   ├── epsilon                 ◄── Turbulent dissipation rate
│   └── nut                     ◄── Turbulent viscosity
├── constant/                   ◄── Physical properties
│   ├── transportProperties     ◄── Fluid properties (air)
│   ├── turbulenceProperties    ◄── Turbulence model selection
│   ├── MRFProperties           ◄── ★ MRF zone definition ★
│   └── polyMesh/               ◄── Mesh (generated by blockMesh)
├── system/                     ◄── Solver & numerics configuration
│   ├── blockMeshDict           ◄── Background mesh definition
│   ├── topoSetDict             ◄── Cell zone creation for MRF
│   ├── controlDict             ◄── Run control & function objects
│   ├── fvSchemes               ◄── Discretisation schemes
│   └── fvSolution              ◄── Solver & SIMPLE settings
├── Allrun                      ◄── Run script
├── Allclean                    ◄── Clean script
└── README.md                   ◄── This file
```

---

## 5. Mesh Generation — blockMesh

The mesh is a simple rectangular domain created with `blockMesh`. This is intentionally
simplified — a real wind turbine simulation would use `snappyHexMesh` with an actual
blade STL geometry. The purpose here is to demonstrate the MRF concept.

### Domain Geometry

```
  z = 0.1  ┌──────────────────────────────────────────┐
            │                                          │  Patch: frontAndBack
            │  150 × 100 × 1 cells                     │  (empty — 2D)
            │                                          │
  z = 0    └──────────────────────────────────────────┘

            x = 0                                x = 3
            y = 0                                y = 2

            y
            ▲
        2.0 ┤ ╔══════════════════════════════════════╗
            │ ║            topAndBottom               ║
            │ ║                                      ║
            │ ║           ╭──────╮                   ║
        1.0 ┤ ║          │ MRF   │                   ║
            │ ║i          │ zone  │              out  ║
            │ ║n         │r=0.3  │              let  ║
        0.5 ┤ ║l          ╰──────╯                   ║
            │ ║e                                     ║
            │ ║t         center:(1,0.5)              ║
            │ ║            topAndBottom               ║
        0.0 ┤ ╚══════════════════════════════════════╝
            └──┬─────────┬─────────┬─────────┬───────►x
              0.0       1.0       2.0       3.0
```

### blockMeshDict Key Parameters

| Parameter       | Value                          | Purpose                      |
|-----------------|--------------------------------|------------------------------|
| Domain size     | 3.0 × 2.0 × 0.1 m            | Rectangular, quasi-2D        |
| Cell count      | 150 × 100 × 1                 | 15,000 cells total           |
| Cell size (dx)  | 0.02 m                         | Uniform spacing              |
| Cell size (dy)  | 0.02 m                         | Uniform spacing              |
| Grading         | (1 1 1) — uniform              | No stretching                |

### Patch Definitions

| Patch Name     | Type   | Location        | Face Vertices  |
|----------------|--------|-----------------|----------------|
| `inlet`        | patch  | x = 0 plane     | (0 4 7 3)      |
| `outlet`       | patch  | x = 3 plane     | (1 2 6 5)      |
| `topAndBottom`  | patch  | y = 0, y = 2   | Top + bottom   |
| `frontAndBack`  | empty  | z = 0, z = 0.1 | Front + back   |

---

## 6. Cell Zone Creation — topoSet

The `topoSet` utility creates cell sets and cell zones from geometric primitives.
For MRF, we need a **cellZone** (not just a cellSet) because the MRF framework
references zones by name.

### topoSet Workflow

```
  ┌──────────────────────────────────────────────────────────────┐
  │                                                              │
  │  Step 1: cylinderToCell                                      │
  │  ─────────────────────                                      │
  │  Select all cells whose centers fall inside a cylinder:      │
  │                                                              │
  │    • Center axis: from (1.0, 0.5, -1) to (1.0, 0.5, 1)     │
  │    • Radius: 0.3 m                                           │
  │    • Result: cellSet "rotatingZone"                          │
  │                                                              │
  │             ╭───────╮                                        │
  │            ╱         ╲       All cells with centers          │
  │           │  selected │      inside this cylinder            │
  │           │   cells   │      are marked                      │
  │            ╲         ╱                                        │
  │             ╰───────╯                                        │
  │                                                              │
  │  Step 2: setToCellZone                                       │
  │  ────────────────────                                        │
  │  Convert the cellSet to a cellZone:                          │
  │                                                              │
  │    cellSet "rotatingZone" ──► cellZone "rotatingZone"        │
  │                                                              │
  │  The cellZone is stored in the mesh and can be referenced    │
  │  by MRFProperties.                                           │
  │                                                              │
  └──────────────────────────────────────────────────────────────┘
```

### Why Two Steps?

OpenFOAM distinguishes between:

- **cellSet**: A simple collection of cell indices (used for selection/manipulation)
- **cellZone**: A mesh topology entity that solvers can reference by name

MRF requires a `cellZone`, so we first select cells with `cylinderToCell` (creating a
cellSet), then convert it to a cellZone with `setToCellZone`.

### topoSet Source Types

| Source Type         | Description                               |
|---------------------|-------------------------------------------|
| `cylinderToCell`    | Cells inside a cylinder (used here)       |
| `sphereToCell`      | Cells inside a sphere                     |
| `boxToCell`         | Cells inside a bounding box               |
| `surfaceToCell`     | Cells near/inside an STL surface          |
| `zoneToCell`        | Cells from an existing zone               |
| `regionToCell`      | Connected region of cells                 |

---

## 7. MRFProperties — The Key Configuration

The `constant/MRFProperties` file is what makes this an MRF simulation. When
`simpleFoam` reads this file, it adds rotational source terms to the momentum
equation inside the specified cell zone.

### MRFProperties Dictionary Breakdown

```
MRF1                            // Name for this MRF zone (arbitrary)
{
    cellZone    rotatingZone;   // Must match the cellZone from topoSet
    active      yes;            // Enable/disable this MRF zone

    nonRotatingPatches ();      // Patches INSIDE the zone that don't rotate
                                // (e.g., stator blades). Empty = all rotate.

    origin      (1.0 0.5 0);   // Center of rotation
    axis        (0 0 1);       // Rotation axis direction (z-axis here)
    omega       10.472;         // Angular velocity in rad/s
}
```

### Angular Velocity Conversion

```
  RPM to rad/s:   ω = RPM × 2π / 60

  ┌──────────────────────────────────────┐
  │  100 RPM × 2π / 60 = 10.472 rad/s   │
  │   50 RPM × 2π / 60 =  5.236 rad/s   │
  │  200 RPM × 2π / 60 = 20.944 rad/s   │
  │  500 RPM × 2π / 60 = 52.360 rad/s   │
  │ 1000 RPM × 2π / 60 = 104.72 rad/s   │
  │ 1500 RPM × 2π / 60 = 157.08 rad/s   │
  │ 3000 RPM × 2π / 60 = 314.16 rad/s   │
  └──────────────────────────────────────┘
```

### Multiple MRF Zones

You can define multiple MRF zones (e.g., for a multi-stage compressor):

```
// Example: Two counter-rotating rotors
MRF1
{
    cellZone    rotor1Zone;
    origin      (1.0 0.5 0);
    axis        (0 0 1);
    omega       10.472;   // Forward rotation
}

MRF2
{
    cellZone    rotor2Zone;
    origin      (2.0 0.5 0);
    axis        (0 0 1);
    omega       -10.472;  // Counter-rotation
}
```

---

## 8. Boundary Conditions

### Boundary Condition Summary Table

| Field     | Inlet                 | Outlet          | topAndBottom | frontAndBack |
|-----------|-----------------------|-----------------|--------------|--------------|
| **U**     | fixedValue (10 0 0)   | zeroGradient    | slip         | empty        |
| **p**     | zeroGradient          | fixedValue 0    | slip         | empty        |
| **k**     | fixedValue 0.375      | zeroGradient    | slip         | empty        |
| **epsilon** | fixedValue 14.855   | zeroGradient    | slip         | empty        |
| **nut**   | calculated            | calculated      | calculated   | empty        |

### Boundary Condition Rationale

```
  ┌─────────────────────────────────────────────────────────────────┐
  │                                                                 │
  │  INLET (x = 0):                                                │
  │  ├─ U:       fixedValue — prescribe incoming wind speed         │
  │  ├─ p:       zeroGradient — pressure floats at inlet            │
  │  ├─ k:       fixedValue — set turbulence level                  │
  │  ├─ epsilon: fixedValue — set dissipation rate                  │
  │  └─ nut:     calculated — computed from k and epsilon           │
  │                                                                 │
  │  OUTLET (x = 3):                                                │
  │  ├─ U:       zeroGradient — flow exits freely                   │
  │  ├─ p:       fixedValue 0 — reference pressure                  │
  │  ├─ k:       zeroGradient — turbulence exits freely             │
  │  ├─ epsilon: zeroGradient — dissipation exits freely            │
  │  └─ nut:     calculated — computed from k and epsilon           │
  │                                                                 │
  │  TOP & BOTTOM (y = 0, y = 2):                                   │
  │  ├─ All:     slip — free-slip, no friction (far-field approx.)  │
  │  └─ nut:     calculated                                         │
  │                                                                 │
  │  FRONT & BACK (z = 0, z = 0.1):                                 │
  │  └─ All:     empty — 2D simulation (no variation in z)          │
  │                                                                 │
  └─────────────────────────────────────────────────────────────────┘
```

---

## 9. Turbulence Setup

### k-epsilon Model

The standard k-epsilon model is used with two transport equations:

```
  ┌──────────────────────────────────────────────────────────────┐
  │                                                              │
  │  Turbulent kinetic energy (k):                               │
  │                                                              │
  │  ∂k         ∂k        ∂   ⎡        ν_t  ∂k ⎤               │
  │  ── + U_j ─── = ──── ⎢(ν + ───) ──── ⎥ + P_k - ε          │
  │  ∂t        ∂x_j     ∂x_j ⎣        σ_k  ∂x_j⎦              │
  │                                                              │
  │                                                              │
  │  Turbulent dissipation rate (ε):                             │
  │                                                              │
  │  ∂ε         ∂ε        ∂   ⎡        ν_t  ∂ε ⎤    ε          │
  │  ── + U_j ─── = ──── ⎢(ν + ───) ──── ⎥ + ─ (C₁P_k - C₂ε) │
  │  ∂t        ∂x_j     ∂x_j ⎣        σ_ε  ∂x_j⎦    k         │
  │                                                              │
  └──────────────────────────────────────────────────────────────┘
```

### Inlet Turbulence Estimation

```
  Given:
    U_inlet = 10 m/s
    TI = 5% (turbulence intensity — typical for wind conditions)
    L = 0.042 m (turbulent length scale ≈ 0.07 × D_MRF)

  Turbulent kinetic energy:
    k = 1.5 × (U × TI)²
    k = 1.5 × (10 × 0.05)²
    k = 1.5 × 0.25
    k = 0.375 m²/s²

  Turbulent dissipation:
    ε = C_μ^0.75 × k^1.5 / L
    ε = 0.09^0.75 × 0.375^1.5 / 0.042
    ε ≈ 14.855 m²/s³

  Turbulent viscosity:
    ν_t = C_μ × k² / ε
    ν_t = 0.09 × 0.375² / 14.855
    ν_t ≈ 8.52 × 10⁻⁴ m²/s
```

---

## 10. Solver Configuration

### SIMPLE Algorithm

This case uses the **SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations)
algorithm, which is standard for steady-state incompressible flow:

```
  ┌─────────────────────────────────────────────────────┐
  │                                                     │
  │  SIMPLE Iteration Loop:                             │
  │                                                     │
  │  1. Solve momentum (U) with guessed pressure        │
  │  2. Solve pressure correction equation              │
  │  3. Correct velocity with pressure correction       │
  │  4. Solve turbulence equations (k, ε)               │
  │  5. Update turbulent viscosity (ν_t)                │
  │  6. Apply MRF source terms                          │
  │  7. Check convergence → if not converged, go to 1   │
  │                                                     │
  │  ┌───┐  ┌───┐  ┌───┐  ┌─────┐  ┌─────┐            │
  │  │ U │→│ p │→│ U │→│ k,ε │→│ ν_t │→ Repeat     │
  │  └───┘  └───┘  └───┘  └─────┘  └─────┘            │
  │                                                     │
  └─────────────────────────────────────────────────────┘
```

### Linear Solver Settings

| Field     | Solver         | Smoother         | Tolerance | relTol |
|-----------|----------------|------------------|-----------|--------|
| p         | GAMG           | GaussSeidel      | 1e-6      | 0.1    |
| U         | smoothSolver   | symGaussSeidel   | 1e-6      | 0.1    |
| k         | smoothSolver   | symGaussSeidel   | 1e-6      | 0.1    |
| epsilon   | smoothSolver   | symGaussSeidel   | 1e-6      | 0.1    |

### Relaxation Factors

| Field/Equation | Factor | Purpose                                    |
|----------------|--------|--------------------------------------------|
| p (field)      | 0.3    | Pressure under-relaxation (stabilizes)     |
| U (equation)   | 0.7    | Velocity equation relaxation               |
| k (equation)   | 0.7    | Turbulent KE equation relaxation           |
| epsilon (eq.)  | 0.7    | Dissipation equation relaxation            |

### Discretisation Schemes

| Term             | Scheme                            | Properties          |
|------------------|-----------------------------------|---------------------|
| Time             | steadyState                       | No time derivative   |
| Gradient         | Gauss linear                      | 2nd order accurate  |
| div(phi,U)       | bounded Gauss linearUpwind        | 2nd order, bounded  |
| div(phi,k)       | bounded Gauss upwind              | 1st order, stable   |
| div(phi,epsilon) | bounded Gauss upwind              | 1st order, stable   |
| Laplacian        | Gauss linear corrected            | 2nd order           |

---

## 11. How to Run

### Prerequisites

- OpenFOAM v2312 or compatible version installed and sourced
- `simpleFoam`, `blockMesh`, and `topoSet` available in `$PATH`

### Step-by-Step Execution

```bash
# Navigate to the case directory
cd projects/09_wind_turbine_blade

# Option 1: Run everything with the Allrun script
./Allrun

# Option 2: Run each step manually
blockMesh          # Generate the rectangular mesh
topoSet            # Create the cylindrical MRF cell zone
simpleFoam         # Run the solver

# Clean all generated files
./Allclean
```

### What Each Step Does

```
  ┌─────────────┐     ┌─────────────┐     ┌──────────────┐
  │  blockMesh  │────►│   topoSet   │────►│  simpleFoam  │
  │             │     │             │     │              │
  │ Creates the │     │ Selects     │     │ Solves N-S   │
  │ rectangular │     │ cylindrical │     │ with MRF     │
  │ background  │     │ cell zone   │     │ source terms │
  │ mesh        │     │ for MRF     │     │              │
  │             │     │             │     │ Iterates     │
  │ 150×100×1   │     │ r=0.3m at   │     │ until        │
  │ = 15,000    │     │ (1.0,0.5)   │     │ convergence  │
  │ cells       │     │             │     │              │
  └─────────────┘     └─────────────┘     └──────────────┘
```

### Monitoring Convergence

During the simulation, watch the residuals:

```bash
# In a separate terminal, monitor residuals in real-time
foamMonitor -l postProcessing/residuals/0/solverInfo.dat

# Or use gnuplot
gnuplot -e "set logscale y; plot 'postProcessing/residuals/0/solverInfo.dat' \
  using 1:2 title 'Ux' with lines, '' using 1:3 title 'Uy' with lines, \
  '' using 1:4 title 'p' with lines; pause -1"
```

---

## 12. Expected Results

### Flow Field Description

```
  ┌──────────────────────────────────────────────────────────────┐
  │                                                              │
  │  What you should observe:                                    │
  │                                                              │
  │  1. UPSTREAM: Uniform flow at 10 m/s in x-direction          │
  │                                                              │
  │  2. MRF ZONE: Flow is deflected by the Coriolis force:       │
  │     • Velocity vectors show a swirling pattern               │
  │     • Velocity magnitude increases near the zone edge        │
  │     • A tangential velocity component appears                │
  │                                                              │
  │  3. DOWNSTREAM WAKE:                                          │
  │     • Swirling wake extends behind the MRF zone              │
  │     • Velocity deficit (lower speed than freestream)         │
  │     • Gradual recovery toward freestream conditions          │
  │     • Wake asymmetry due to rotation direction               │
  │                                                              │
  │  4. PRESSURE FIELD:                                          │
  │     • High pressure upstream of MRF zone                     │
  │     • Low pressure region inside and downstream              │
  │     • Asymmetric pressure distribution (rotation effect)     │
  │                                                              │
  └──────────────────────────────────────────────────────────────┘
```

### Qualitative Flow Pattern

```
  Inlet                MRF Zone                   Outlet
  ─────                ────────                   ──────

  →→→→→→→→→        ╭───────────╮           →→→→→→→→
  →→→→→→→→→       ╱   ↗  →  ↘  ╲          →→→↗→→→→
  →→→→→→→→→      │  ↑   ↻    ↓  │         →→→→→→→→
  →→→→→→→→→      │  ↖  ←  ↙    │         →→→↘→→→→
  →→→→→→→→→       ╲             ╱          →→→→→→→→
  →→→→→→→→→        ╰───────────╯           →→→→→→→→

  Uniform           Swirling flow           Wake with
  freestream        (rotation effect)       velocity deficit
```

### Convergence Expectations

- Initial residuals may spike as MRF source terms take effect
- Residuals should decrease monotonically after ~50-100 iterations
- Target residuals of 1e-4 should be reached within 500-1000 iterations
- If residuals oscillate, reduce relaxation factors

---

## 13. Post-Processing

### ParaView Visualization

```bash
# Open the case in ParaView
paraFoam

# Or create a .foam file for recent ParaView versions
touch case.foam
paraview --data=case.foam
```

### Recommended Visualizations

1. **Velocity Magnitude Contour**
   - Shows the speed-up around the MRF zone edges
   - Wake velocity deficit is clearly visible

2. **Streamlines**
   - Seed from the inlet
   - Shows how flow is deflected by the rotating zone
   - Reveals the swirling wake structure

3. **Pressure Contour**
   - High pressure upstream, low downstream
   - Asymmetric pattern due to rotation

4. **Velocity Vectors (Glyph filter)**
   - Inside the MRF zone: shows tangential component
   - In the wake: shows swirl decay

5. **Vorticity**
   - Apply the `curl(U)` calculator
   - Strong vorticity inside and at the boundary of the MRF zone

### Command-Line Post-Processing

```bash
# Calculate vorticity field
simpleFoam -postProcess -func "vorticity"

# Calculate Q-criterion
simpleFoam -postProcess -func "Q"

# Sample along a line through the MRF zone center
postProcess -func "
    type            sets;
    libs            (sampling);
    writeControl    writeTime;
    interpolationScheme cellPoint;
    setFormat       raw;
    sets
    (
        centerline
        {
            type    lineUniform;
            axis    x;
            start   (0 0.5 0.05);
            end     (3 0.5 0.05);
            nPoints 300;
        }
    );
    fields (U p k);
"
```

---

## 14. Exercises

### Exercise 1: Vary the Rotation Speed

Modify `constant/MRFProperties` to test different RPMs:

```
omega   5.236;    // 50 RPM — mild rotation
omega   10.472;   // 100 RPM — baseline (current)
omega   20.944;   // 200 RPM — strong rotation
omega   52.360;   // 500 RPM — very strong rotation (may need lower relaxation)
```

Observe how the wake swirl and velocity deficit change with RPM.

### Exercise 2: Change the MRF Zone Size

Modify `system/topoSetDict` to change the cylinder radius:

```
radius  0.15;    // Smaller rotor zone
radius  0.30;    // Baseline (current)
radius  0.50;    // Larger rotor zone
radius  0.80;    // Very large zone (check domain size!)
```

Re-run `topoSet` and `simpleFoam` after changing. How does the zone size affect the wake?

### Exercise 3: Move the MRF Zone

Change the MRF zone position to see how upstream/downstream distance affects results:

```
// In topoSetDict:
p1 (0.5 0.5 -1);  p2 (0.5 0.5 1);   // Closer to inlet
p1 (1.5 0.5 -1);  p2 (1.5 0.5 1);   // Default position
p1 (2.0 0.5 -1);  p2 (2.0 0.5 1);   // Closer to outlet

// Don't forget to update MRFProperties origin to match!
```

### Exercise 4: Add a Second MRF Zone (Counter-Rotating)

Add a second rotor zone downstream to simulate a counter-rotating turbine:

```
// In topoSetDict — add a second cylinder:
{
    name    rotatingZone2;
    type    cellSet;
    action  new;
    source  cylinderToCell;
    sourceInfo { p1 (2.0 0.5 -1); p2 (2.0 0.5 1); radius 0.3; }
}
{
    name    rotatingZone2;
    type    cellZoneSet;
    action  new;
    source  setToCellZone;
    sourceInfo { set rotatingZone2; }
}

// In MRFProperties — add counter-rotating zone:
MRF2
{
    cellZone    rotatingZone2;
    active      yes;
    nonRotatingPatches ();
    origin      (2.0 0.5 0);
    axis        (0 0 1);
    omega       -10.472;  // Counter-rotation (negative)
}
```

### Exercise 5: Refine the Mesh

Test mesh sensitivity by modifying `system/blockMeshDict`:

```
// Coarse:  75 × 50 × 1 = 3,750 cells
// Medium: 150 × 100 × 1 = 15,000 cells (current)
// Fine:   300 × 200 × 1 = 60,000 cells
// Very fine: 450 × 300 × 1 = 135,000 cells
```

Compare velocity profiles through the MRF zone center for each mesh.

### Exercise 6: Try k-omega SST Turbulence Model

Replace k-epsilon with k-omega SST:

1. Change `constant/turbulenceProperties`:
   ```
   RASModel    kOmegaSST;
   ```
2. Replace `0/epsilon` with `0/omega`
3. Adjust inlet values for omega:
   ```
   omega = epsilon / (Cmu * k) = 14.855 / (0.09 * 0.375) = 440.1 1/s
   ```

### Exercise 7: Transition to AMI (Advanced)

For a more accurate simulation, convert this case to use the AMI sliding mesh approach:

1. Create a cylindrical mesh region with `snappyHexMesh` or `cfMesh`
2. Define AMI patch pairs at the interface
3. Switch to `pimpleFoam` (transient solver)
4. Add `dynamicMeshDict` with rotation
5. Use much smaller time steps (Courant number < 1)

This is significantly more complex but physically more accurate.

---

## 15. Industrial Applications

The MRF technique demonstrated here is widely used in industry:

```
  ┌──────────────────────────────────────────────────────────┐
  │                                                          │
  │  WIND TURBINES                                           │
  │  ┌───┐                                                   │
  │  │ ╱ │  • Rotor performance prediction                   │
  │  │╱  │  • Annual energy production estimates              │
  │  │╲  │  • Wake interaction studies                        │
  │  │ ╲ │  • Yaw misalignment effects                       │
  │  └───┘                                                   │
  │                                                          │
  │  CENTRIFUGAL PUMPS & FANS                                │
  │  ╭───╮                                                   │
  │  │ ↻ │  • Impeller performance curves                    │
  │  ╰───╯  • Head-flow characteristics                      │
  │         • Efficiency optimization                        │
  │                                                          │
  │  AXIAL COMPRESSORS                                       │
  │  ═══╪═══  • Stage performance matching                   │
  │     │     • Surge/stall prediction                       │
  │           • Multi-stage interactions                     │
  │                                                          │
  │  MIXING VESSELS                                          │
  │  ┌─┬─┐                                                   │
  │  │ │ │  • Impeller mixing efficiency                     │
  │  │ ↻ │  • Power consumption prediction                   │
  │  └───┘  • Blending time estimates                        │
  │                                                          │
  │  MARINE PROPELLERS                                       │
  │  ───╤───  • Thrust and torque prediction                 │
  │     │     • Cavitation onset                             │
  │           • Hull-propeller interaction                   │
  │                                                          │
  │  HELICOPTER ROTORS                                       │
  │  ──┼──    • Hover performance                            │
  │    │      • Forward flight approximation                 │
  │           • Ground effect studies                        │
  │                                                          │
  └──────────────────────────────────────────────────────────┘
```

---

## 16. Troubleshooting

### Common Issues

| Problem                          | Cause                          | Solution                       |
|----------------------------------|--------------------------------|--------------------------------|
| `cellZone not found`             | topoSet not run or failed      | Run `topoSet` before solver    |
| Divergence at startup            | RPM too high, relaxation too high | Reduce omega, reduce relaxation |
| Residuals plateau high           | Mesh too coarse                 | Refine mesh near MRF zone      |
| `empty` patch error              | Wrong face assignment           | Check blockMeshDict faces       |
| No rotation effect visible       | MRF zone too small or omega=0  | Check MRFProperties settings    |
| Oscillating residuals            | CFL-like instability           | Reduce p relaxation to 0.2     |

### Debugging Tips

```bash
# Check if the cellZone was created correctly
checkMesh -allGeometry -allTopology

# Verify cellZone exists in the mesh
foamToVTK -cellZone rotatingZone

# Check MRF is being read by the solver
simpleFoam | grep -i "MRF"

# Visualize just the MRF zone cells
paraFoam  # Then use "Extract Block" filter on the rotatingZone
```

---

## 17. References

1. **OpenFOAM User Guide — MRF**
   - Section on Multiple Reference Frame zones
   - `$WM_PROJECT_DIR/doc/` or https://www.openfoam.com/documentation/guides

2. **OpenFOAM Tutorial: Mixer (MRF)**
   - `$FOAM_TUTORIALS/incompressible/simpleFoam/mixerVessel2D`
   - Standard MRF tutorial case included with OpenFOAM

3. **Luo, J.Y., Issa, R.I., Gosman, A.D. (1994)**
   - "Prediction of Impeller-Induced Flows in Mixing Vessels Using
     Multiple Frames of Reference"
   - Institution of Chemical Engineers Symposium Series, 136, pp. 549-556

4. **Tabor, G., Gosman, A.D., Issa, R.I. (1996)**
   - "Numerical Simulation of the Flow in a Mixing Vessel Stirred by
     a Rushton Turbine"
   - Proc. 1st European-Africa Conf. on Wind Engineering

5. **Jasak, H. (1996)**
   - "Error Analysis and Estimation for the Finite Volume Method with
     Applications to Fluid Flows"
   - PhD Thesis, Imperial College London

6. **OpenFOAM Source Code**
   - `src/finiteVolume/cfdTools/general/MRF/` — MRF implementation
   - `src/finiteVolume/cfdTools/general/MRF/MRFZone.C` — Source term calculation

---

## 18. Cross-References

### Related Tutorials in This Repository

- **Tutorial 01** (if exists): Basic mesh generation with blockMesh
- **Tutorial 02** (if exists): simpleFoam basics
- See `notes/` directory for additional learning materials

### Related OpenFOAM Utilities

| Utility          | Purpose                                      |
|------------------|----------------------------------------------|
| `blockMesh`      | Structured mesh generation                    |
| `topoSet`        | Cell/face/point set and zone manipulation     |
| `checkMesh`      | Mesh quality analysis                         |
| `simpleFoam`     | Steady-state incompressible solver            |
| `pimpleFoam`     | Transient incompressible solver (for AMI)     |
| `foamToVTK`      | Convert results to VTK format                 |
| `postProcess`    | Run function objects on existing results      |
| `paraFoam`       | Launch ParaView with OpenFOAM reader          |

### OpenFOAM Keyword Reference

| Keyword in This Case  | File                    | Purpose                    |
|-----------------------|-------------------------|----------------------------|
| `MRF1`               | MRFProperties           | MRF zone definition        |
| `cellZone`           | MRFProperties           | Zone name reference        |
| `omega`              | MRFProperties           | Angular velocity (rad/s)   |
| `cylinderToCell`     | topoSetDict             | Cylinder cell selection    |
| `setToCellZone`      | topoSetDict             | Convert set to zone        |
| `steadyState`        | fvSchemes               | No time derivative         |
| `SIMPLE`             | fvSolution              | Pressure-velocity coupling |
| `GAMG`               | fvSolution              | Geometric agglomeration    |
| `linearUpwind`       | fvSchemes               | 2nd order bounded scheme   |

---

*This tutorial is part of the OpenFoam-Tutorials repository.*
*For questions or improvements, please open an issue or pull request.*
