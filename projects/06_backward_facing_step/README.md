# Backward-Facing Step Flow

> **Level:** Intermediate · **Solver:** `simpleFoam` · **Flow:** Steady-state, incompressible, turbulent

---

## Table of Contents

1. [Problem Description](#problem-description)
2. [Physics & Theory](#physics--theory)
3. [Case Structure](#case-structure)
4. [Mesh Generation](#mesh-generation)
5. [Boundary Conditions](#boundary-conditions)
6. [Turbulence Setup](#turbulence-setup)
7. [Solver Configuration](#solver-configuration)
8. [How to Run](#how-to-run)
9. [Expected Results](#expected-results)
10. [Post-Processing Guide](#post-processing-guide)
11. [Exercises](#exercises)
12. [References](#references)

---

## Problem Description

This case simulates steady-state turbulent flow over a **backward-facing step** — one of
the most important and widely studied validation benchmarks in computational fluid
dynamics. A sudden expansion in a channel creates a region of separated flow downstream
of the step, followed by reattachment further along the bottom wall.

```
  Inlet U∞                             Reattachment point
    →→→→                                      ↓
  ┌──────────┐────────────────────────────────────────────────────┐
  │          │    Recirculation         →→→→→→→→→→→→→→→→→→→→→→→→│
  │  →→→→→→  │    ╭──────────╮         →→→→→→→→→→→→→→→→→→→→→→→→│
  │  →→→→→→  │    │  primary │  ←←←←   →→→→→→→→→→→→→→→→→→→→→→→→│
  │  →→→→→→  │    │  eddy   │          →→→→→→→→→→→→→→→→→→→→→→→→│
  │  →→→→→→  ├────╰──────────╯─────────┬──→→→→→→→→→→→→→→→→→→→→→│
  │          │ Step    Xr ≈ 6H         │                        │
  └──────────┘                          └────────────────────────┘
  ← 4H →      ←          30H                                  →
```

The backward-facing step is deceptively simple in geometry but extraordinarily rich in
flow physics. It has served as a primary test case for turbulence model validation for
over four decades and continues to challenge RANS models, LES approaches, and DNS
simulations alike.

### Why This Case Matters

- **Separated-flow benchmark** — The reattachment length Xr/H is a universally
  accepted metric for model validation. Experimental data from Driver & Seegmiller
  (1985), Kim et al. (1980), and Armaly et al. (1983) provide gold-standard comparisons.
- **Turbulence model stress test** — Separation and reattachment involve adverse
  pressure gradients, streamline curvature, and anisotropic turbulence that expose
  weaknesses in common RANS closures.
- **Mesh sensitivity showcase** — The recirculation zone is extremely sensitive to
  near-wall resolution, making this an ideal case for learning about y+ requirements,
  wall functions, and grid convergence studies.
- **Practical relevance** — Sudden expansions appear in combustors, diffusers, HVAC
  ducts, electronic cooling channels, and many other engineering applications.

### Key Parameters

| Parameter              | Symbol     | Value                         |
|------------------------|------------|-------------------------------|
| Step height            | H          | 0.0127 m (12.7 mm)           |
| Inlet velocity         | U∞         | 10 m/s                        |
| Kinematic viscosity    | ν          | 1.5 × 10⁻⁵ m²/s (air)       |
| Reynolds number        | Re_H       | ≈ 8,467 (based on H)         |
| Expansion ratio        | ER         | 2:1                           |
| Upstream length        | L_up       | 4H = 0.0508 m                |
| Downstream length      | L_down     | 30H = 0.381 m                |
| Channel height (inlet) | h          | H (above step)                |
| Channel height (outlet)| 2H         | 0.0254 m                     |
| Expected reattachment  | Xr/H       | ≈ 6.0 ± 0.5                  |

> **Note on Reynolds number:** Different researchers define Re using different
> reference lengths (step height H, hydraulic diameter, or upstream channel height).
> Re_H ≈ 8,467 based on step height is common; using the full downstream channel
> height gives $Re \approx 16{,}933$. Some references cite $Re$ based on upstream centreline
> velocity and channel height which can give $Re \approx 36{,}000$.

---

## Physics & Theory

### Governing Equations

For steady-state, incompressible flow, the Reynolds-Averaged Navier-Stokes (RANS)
equations govern the mean flow:

**Continuity:**

```
∂Ūᵢ/∂xᵢ = 0
```

**Momentum:**

```
Ūⱼ ∂Ūᵢ/∂xⱼ = -1/ρ ∂p̄/∂xᵢ + ∂/∂xⱼ [ν(∂Ūᵢ/∂xⱼ + ∂Ūⱼ/∂xᵢ) - u'ᵢu'ⱼ]
```

The Reynolds stress tensor `u'ᵢu'ⱼ` is modelled using the Boussinesq hypothesis:

```
-u'ᵢu'ⱼ = νₜ(∂Ūᵢ/∂xⱼ + ∂Ūⱼ/∂xᵢ) - 2/3 k δᵢⱼ
```

where νₜ is the turbulent eddy viscosity and k is the turbulent kinetic energy.

### Flow Separation and Reattachment

When the flow encounters the sudden expansion at the step edge, the following physical
processes occur:

1. **Separation** — The flow cannot follow the abrupt geometry change. A free shear
   layer detaches from the step edge, with a strong inflection point in the velocity
   profile.

2. **Primary recirculation zone** — Below the shear layer and behind the step, a
   large clockwise (in the standard orientation) vortex forms. This is the primary
   recirculation bubble, bounded by the step face, the bottom wall, and the shear
   layer.

3. **Shear layer development** — The free shear layer grows through Kelvin-Helmholtz
   instability, entraining fluid from both the main flow above and the recirculation
   zone below. Turbulent mixing in this layer controls the reattachment process.

4. **Reattachment** — At a distance Xr downstream of the step, the shear layer
   impinges on the bottom wall. The reattachment point is defined as the location
   where the time-averaged wall shear stress changes sign (from negative to positive).

5. **Recovery region** — Downstream of reattachment, the flow gradually recovers
   toward a fully developed channel flow profile. This recovery region extends for
   approximately 30–50H.

```
  Flow structure detail:
                                                                  
                 Free shear layer                                 
            ╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱                          
           ╱  Mixing & entrainment     ╱                          
          ╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱                           
  ┌────┐                                                          
  │    │  ←←←←←←←←←  (reverse flow)                              
  │    │  ╭────────────────╮                                      
  │    │  │   Primary      │                                      
  │    │  │   recirculation│                                      
  │    │  │   vortex       │     Reattachment                     
  │    ├──╰────────────────╯─────── ★ ────── Recovery →→→→        
  │    │  Step                      Xr                            
  │    │                                                          
```

### The Reattachment Length (Xr/H)

The normalised reattachment length Xr/H is the single most important validation
metric for backward-facing step simulations:

| Study                       | Re        | Xr/H   | Method          |
|-----------------------------|-----------|--------|-----------------|
| Driver & Seegmiller (1985)  | 37,400    | 6.26   | Experiment      |
| Kim et al. (1980)           | 132,000   | 7.0    | Experiment      |
| Armaly et al. (1983)        | 8,000     | 6.0    | Experiment      |
| Adams & Johnston (1988)     | 38,000    | 6.4    | Experiment      |
| Menter (1994) k-ω SST      | 37,400    | 6.1    | RANS CFD        |
| Standard k-ε               | Various   | 4.5–5  | RANS CFD        |

> **Key insight:** Standard k-ε consistently under-predicts the reattachment length
> by 15–25%, while k-ω SST typically predicts Xr/H within 5–10% of experimental
> values. This is one of the primary reasons k-ω SST became the industry standard
> for separated flows.

### Why k-omega SST for This Problem

The k-ω SST (Shear Stress Transport) model by Menter (1994) was specifically developed
to handle the weaknesses of existing two-equation models in separated flows:

1. **Near-wall behaviour** — Uses k-ω formulation in the inner boundary layer,
   avoiding the need for damping functions and providing superior performance in
   adverse pressure gradient regions.

2. **Free-stream robustness** — Switches to k-ε behaviour in the outer region and
   free shear layers, avoiding the sensitivity of k-ω to free-stream turbulence
   levels.

3. **Shear stress limiter** — The SST modification limits the eddy viscosity using
   Bradshaw's assumption that the principal turbulent shear stress is proportional
   to the turbulent kinetic energy: `τ = ρa₁k` where a₁ = 0.31. This prevents
   over-prediction of eddy viscosity in adverse pressure gradient flows, which is
   the primary failure mode of standard k-ε in separated flows.

4. **Blending function** — A smooth blending function F₁ transitions between the
   k-ω and k-ε formulations based on wall distance, ensuring the right model is
   active in each region:

```
  F₁ = tanh(arg₁⁴)
  
  arg₁ = min[max(√k/(β*ωy), 500ν/(y²ω)), 4ρσ_ω₂k/(CD_kω y²)]
```

The transport equations for k and ω in the SST model are:

**Turbulent kinetic energy (k):**
```
∂(ρk)/∂t + ∂(ρUⱼk)/∂xⱼ = P̃k - β*ρkω + ∂/∂xⱼ[(μ + σk μt) ∂k/∂xⱼ]
```

**Specific dissipation rate (ω):**
```
∂(ρω)/∂t + ∂(ρUⱼω)/∂xⱼ = αS² - βρω² + ∂/∂xⱼ[(μ + σω μt) ∂ω/∂xⱼ]
                             + 2(1-F₁)ρσ_ω₂ (1/ω) ∂k/∂xⱼ ∂ω/∂xⱼ
```

### Turbulent Inlet Conditions

The inlet turbulence quantities are estimated from the turbulence intensity (I) and
a mixing length scale (l):

```
  k = 3/2 (U∞ · I)²  =  3/2 × (10 × 0.01)²  =  0.015  →  rounded to 0.24
                          (using ~4% intensity for developed channel flow)
  
  ω = k^0.5 / (C_μ^0.25 · l)  ≈  440 s⁻¹
  
  where C_μ = 0.09, l ≈ 0.07 × h_channel
```

> The exact inlet values have moderate influence on the reattachment length. A
> sensitivity study varying I from 1% to 10% typically shifts Xr/H by ±0.5.

---

## Case Structure

```
06_backward_facing_step/
├── 0/                          ← Initial & boundary conditions
│   ├── U                       ← Velocity field
│   ├── p                       ← Pressure field (kinematic)
│   ├── k                       ← Turbulent kinetic energy
│   ├── omega                   ← Specific dissipation rate
│   └── nut                     ← Turbulent viscosity
├── constant/                   ← Physical properties
│   ├── transportProperties     ← Fluid properties (ν for air)
│   ├── turbulenceProperties    ← Turbulence model selection
│   └── polyMesh/               ← Mesh (generated by blockMesh)
├── system/                     ← Simulation controls
│   ├── blockMeshDict           ← Mesh definition
│   ├── controlDict             ← Time control & function objects
│   ├── fvSchemes               ← Discretisation schemes
│   └── fvSolution              ← Linear solvers & SIMPLE
├── Allrun                      ← Run script
├── Allclean                    ← Clean script
└── README.md                   ← This file
```

### File Purposes

| File                   | Purpose                                                    |
|------------------------|------------------------------------------------------------|
| `0/U`                  | Velocity: inlet profile, no-slip walls, empty for 2D      |
| `0/p`                  | Pressure: zero-gradient inlet, fixed outlet                |
| `0/k`                  | TKE: inlet value, wall functions on walls                  |
| `0/omega`              | Omega: inlet value, wall functions on walls                |
| `0/nut`                | Eddy viscosity: calculated at inlet/outlet, wall functions |
| `constant/transport..` | Newtonian fluid, ν = 1.5e-5 m²/s (air at 20°C)           |
| `constant/turbulence.` | RAS → k-omega SST selection                               |
| `system/blockMeshDict` | Multi-block structured mesh definition                     |
| `system/controlDict`   | 2000 iterations, write every 200, function objects         |
| `system/fvSchemes`     | Steady-state, linear upwind for U, upwind for k/omega     |
| `system/fvSolution`    | SIMPLE, GAMG for p, relaxation factors                     |

---

## Mesh Generation

### Geometry Definition

The mesh is created using `blockMesh` with three structured hexahedral blocks:

```
  Vertex numbering (back plane, z = -0.5H):

    v3 ─────────────── v4 ─────────────────────────── v5
    │                  │                               │
    │                  │                               │
    │    Block 0       │        Block 1                │
    │   (upstream)     │    (downstream top)           │
    │   80 × 40        │       300 × 40                │
    │                  │                               │
    v0 ─────────────── v1                              v2
                       │                               │
                       │        Block 2                │
                       │    (downstream bottom)        │
                       │       300 × 40                │
                       │                               │
                       v6 ────────────────────────── v7

  Coordinates (in units of H, convertToMeters = 0.0127):
    v0 = (-4, 0)    v1 = (0, 0)     v2 = (30, 0)
    v3 = (-4, 1)    v4 = (0, 1)     v5 = (30, 1)
    v6 = (0, -1)    v7 = (30, -1)
```

### Block Details

| Block | Region              | Size   | Cells     | Grading                     |
|-------|---------------------|--------|-----------|-----------------------------|
| 0     | Upstream channel    | 4H×1H  | 80 × 40  | x: 0.25 (refine at step)   |
| 1     | Downstream top      | 30H×1H | 300 × 40 | x: 4.0 (expand downstream) |
| 2     | Downstream bottom   | 30H×1H | 300 × 40 | x: 4.0 (expand downstream) |

**Total cells:** approximately 27,200 (2D, single cell in z-direction)

### Grading Strategy

The mesh grading is designed to resolve critical flow features:

- **Near the step edge (x-direction):** Cells are compressed toward the step corner
  where separation occurs. The upstream block uses x-grading ratio 0.25 (cells get
  smaller toward the step). The downstream blocks use ratio 4.0 (fine near step,
  expanding downstream).

- **Near walls (y-direction):** Each block uses symmetric bi-directional grading
  with ratio 4:1 toward the walls, concentrating cells in the boundary layers.

- **Near the step face:** The junction between blocks 1 and 2 at y = 0 places the
  finest y-cells at the shear layer location — exactly where resolution is most
  critical.

### y+ Considerations

For wall-function approaches (as used here), the first cell centre should be in the
log-layer region:

```
  Target y+ ≈ 30-100 for wall functions
  
  y₁ = y+ · ν / u_τ
  
  where u_τ = √(τ_w/ρ) ≈ U∞ × √(Cf/2)
  
  For flat plate at Re_x ~ 10⁵: Cf ≈ 0.005
  u_τ ≈ 10 × √(0.0025) ≈ 0.5 m/s
  y₁ = 30 × 1.5e-5 / 0.5 ≈ 0.9 mm
```

The mesh grading places the first cell at approximately 1 mm from the walls, giving
y+ ≈ 30–50, appropriate for wall-function treatment. For resolved boundary layers
(y+ < 1), the mesh would need ~10× more cells in the wall-normal direction.

---

## Boundary Conditions

### Patch Overview

```
                         topWall
  ┌──────────────────────────────────────────────────────────┐
  │                                                          │
  │                                                          │
  inlet                                                  outlet
  │                                                          │
  │          step (top + vertical)                           │
  ├──────────┐                                               │
  │          │                                               │
  └──────────┘               bottomWall                      ┘
              
              frontAndBack: both z-planes (empty for 2D)
```

### Velocity (U)

| Patch        | Type        | Value           | Notes                        |
|--------------|-------------|-----------------|------------------------------|
| inlet        | fixedValue  | (10 0 0) m/s    | Uniform inflow               |
| outlet       | zeroGradient| —               | Convective outflow           |
| topWall      | noSlip      | (0 0 0) m/s     | No-slip wall                 |
| bottomWall   | noSlip      | (0 0 0) m/s     | No-slip wall                 |
| step         | noSlip      | (0 0 0) m/s     | Step top face + vertical     |
| frontAndBack | empty       | —               | 2D constraint                |

### Pressure (p)

| Patch        | Type        | Value           | Notes                        |
|--------------|-------------|-----------------|------------------------------|
| inlet        | zeroGradient| —               | Pressure floats at inlet     |
| outlet       | fixedValue  | 0 m²/s²         | Reference pressure           |
| topWall      | zeroGradient| —               | Zero normal gradient         |
| bottomWall   | zeroGradient| —               | Zero normal gradient         |
| step         | zeroGradient| —               | Zero normal gradient         |
| frontAndBack | empty       | —               | 2D constraint                |

### Turbulent Kinetic Energy (k)

| Patch        | Type              | Value       | Notes                      |
|--------------|-------------------|-------------|----------------------------|
| inlet        | fixedValue        | 0.24 m²/s²  | ~4% turbulence intensity   |
| outlet       | zeroGradient      | —           | Convective outflow         |
| topWall      | kqRWallFunction   | 0.24        | Wall function for k        |
| bottomWall   | kqRWallFunction   | 0.24        | Wall function for k        |
| step         | kqRWallFunction   | 0.24        | Wall function for k        |
| frontAndBack | empty             | —           | 2D constraint              |

### Specific Dissipation Rate (omega)

| Patch        | Type                | Value       | Notes                    |
|--------------|---------------------|-------------|--------------------------|
| inlet        | fixedValue          | 440 s⁻¹     | From mixing length       |
| outlet       | zeroGradient        | —           | Convective outflow       |
| topWall      | omegaWallFunction   | 440         | Automatic wall treatment |
| bottomWall   | omegaWallFunction   | 440         | Automatic wall treatment |
| step         | omegaWallFunction   | 440         | Automatic wall treatment |
| frontAndBack | empty               | —           | 2D constraint            |

### Turbulent Viscosity (nut)

| Patch        | Type              | Value     | Notes                       |
|--------------|-------------------|-----------|-----------------------------|
| inlet        | calculated        | 0         | Computed from k and omega   |
| outlet       | calculated        | 0         | Computed from k and omega   |
| topWall      | nutkWallFunction  | 0         | Standard wall function      |
| bottomWall   | nutkWallFunction  | 0         | Standard wall function      |
| step         | nutkWallFunction  | 0         | Standard wall function      |
| frontAndBack | empty             | —         | 2D constraint               |

### Wall Function Summary

The three wall function types work together to provide consistent near-wall modelling:

- **`kqRWallFunction`** — Provides a zero-gradient (Neumann) condition for k at the
  wall, which is the correct physical behaviour since k has a non-zero value at the
  wall in wall-function approaches.

- **`omegaWallFunction`** — Sets omega based on the analytical near-wall solution,
  blending between viscous and log-layer values. This is the most critical wall
  function for the k-ω SST model.

- **`nutkWallFunction`** — Computes νₜ from the log-law, providing the correct eddy
  viscosity profile near the wall. Automatically adapts based on the local y+ value.

---

## Turbulence Setup

### Model Selection

The turbulence model is configured in `constant/turbulenceProperties`:

```
simulationType  RAS;

RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
```

Setting `printCoeffs on` writes the model coefficients to the log file at startup,
allowing verification of default values:

| Coefficient | Value  | Description                    |
|-------------|--------|--------------------------------|
| αk1         | 0.85   | σk in k-ω region               |
| αk2         | 1.0    | σk in k-ε region               |
| αω1         | 0.5    | σω in k-ω region               |
| αω2         | 0.856  | σω in k-ε region               |
| β1          | 0.075  | β in k-ω region                |
| β2          | 0.0828 | β in k-ε region                |
| βStar       | 0.09   | β* (= Cμ)                      |
| γ1          | 5/9    | γ in k-ω region                |
| γ2          | 0.44   | γ in k-ε region                |
| a1          | 0.31   | Bradshaw's constant            |
| b1          | 1.0    | Blending parameter             |
| c1          | 10.0   | Blending parameter             |

### Fluid Properties

Configured in `constant/transportProperties`:

```
transportModel  Newtonian;
nu              nu [0 2 -1 0 0 0 0] 1.5e-05;
```

This corresponds to air at approximately 20°C and atmospheric pressure:

| Property              | Value              |
|-----------------------|--------------------|
| Kinematic viscosity ν | 1.5 × 10⁻⁵ m²/s  |
| Density ρ (implicit)  | ~1.2 kg/m³        |
| Temperature           | ~20°C              |
| Fluid                 | Air                |

---

## Solver Configuration

### SIMPLE Algorithm

This case uses the **SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations)
algorithm for steady-state pressure-velocity coupling. The `consistent yes` option
enables the SIMPLEC variant for improved convergence:

```
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;

    residualControl
    {
        p               1e-4;
        U               1e-5;
        k               1e-5;
        omega           1e-5;
    }
}
```

### Linear Solvers

| Field | Solver       | Smoother       | Tolerance | Notes                      |
|-------|-------------|----------------|-----------|----------------------------|
| p     | GAMG        | GaussSeidel    | 1e-6      | Multigrid for pressure     |
| U     | smoothSolver| symGaussSeidel | 1e-7      | Symmetric Gauss-Seidel     |
| k     | smoothSolver| symGaussSeidel | 1e-7      | Same as U                  |
| omega | smoothSolver| symGaussSeidel | 1e-7      | Same as U                  |

> **GAMG** (Geometric-Algebraic Multi-Grid) is highly efficient for the pressure
> equation because it handles the elliptic nature of the Poisson-type equation well.
> For convection-dominated equations (U, k, ω), smoothSolver is more appropriate.

### Relaxation Factors

| Field/Equation | Factor | Notes                                           |
|----------------|--------|-------------------------------------------------|
| p              | 0.3    | Conservative for pressure stability             |
| U              | 0.7    | Standard for momentum                           |
| k              | 0.7    | Standard for turbulence                         |
| omega          | 0.7    | Standard for turbulence                         |

Lower relaxation factors improve stability but slow convergence. If the solution
diverges, try reducing p to 0.2 and U/k/omega to 0.5.

### Discretisation Schemes

| Term           | Scheme                            | Notes                       |
|----------------|-----------------------------------|-----------------------------|
| ddt            | steadyState                       | No time derivative          |
| grad           | Gauss linear                      | Second-order gradients      |
| div(phi,U)     | bounded Gauss linearUpwind grad(U)| Second-order, bounded       |
| div(phi,k)     | bounded Gauss upwind              | First-order for stability   |
| div(phi,omega) | bounded Gauss upwind              | First-order for stability   |
| laplacian      | Gauss linear corrected            | Second-order with correction|

> **Why linearUpwind for U but upwind for k/omega?** The velocity field benefits
> from second-order accuracy to resolve the separation point and shear layer
> accurately. The turbulence equations are more prone to unboundedness, so first-order
> upwind provides the necessary stability. This is a common and well-tested combination
> for industrial RANS simulations.

---

## How to Run

### Prerequisites

- OpenFOAM v2312 or compatible version installed and sourced
- Bash shell environment

### Running the Case

**Option 1: Using the Allrun script (recommended)**

```bash
cd projects/06_backward_facing_step
./Allrun
```

This will:
1. Generate the mesh with `blockMesh`
2. Check mesh quality with `checkMesh`
3. Run the solver with `simpleFoam`

**Option 2: Step by step**

```bash
cd projects/06_backward_facing_step

# Generate mesh
blockMesh

# Verify mesh quality
checkMesh

# Run solver
simpleFoam
```

**Option 3: Parallel execution (for larger meshes)**

```bash
# Decompose
decomposePar

# Run in parallel (e.g., 4 cores)
mpirun -np 4 simpleFoam -parallel

# Reconstruct
reconstructPar
```

### Cleaning Up

```bash
./Allclean
```

### Monitoring Convergence

While the solver is running, monitor residuals in another terminal:

```bash
# Watch residuals live
foamMonitor -l postProcessing/residuals/0/residuals.dat

# Or use gnuplot
gnuplot -e "set logscale y; plot 'log.simpleFoam' using 1:2 with lines"

# Or simply grep residuals
grep "Solving for Ux" log.simpleFoam | tail -20
```

---

## Expected Results

### Convergence Behaviour

The simulation should converge within **800–1500 iterations**. Typical convergence
pattern:

| Phase          | Iterations | Behaviour                                       |
|----------------|------------|-------------------------------------------------|
| Initial        | 0–100      | Rapid residual drop, flow establishes           |
| Transition     | 100–400    | Residuals may stall or oscillate slightly       |
| Convergence    | 400–1000   | Steady exponential residual decay               |
| Converged      | 800–1500   | All residuals below thresholds                  |

### Reattachment Length

The primary validation metric is the reattachment length:

```
  Expected:  Xr/H ≈ 5.5–6.5  (k-ω SST prediction)
  
  Experimental reference:
    Driver & Seegmiller (1985):  Xr/H = 6.26 ± 0.10
    Armaly et al. (1983):        Xr/H ≈ 6.0 (at similar Re)
```

### Flow Field Features

After convergence, you should observe:

1. **Recirculation bubble** — A closed streamline pattern behind the step, visible
   in velocity vectors or streamlines.

2. **Negative velocity region** — The x-velocity is negative (reverse flow) in the
   recirculation zone.

3. **Pressure recovery** — Pressure increases through the expansion and reaches a
   new equilibrium downstream of reattachment.

4. **Turbulent kinetic energy peak** — Maximum k occurs in the free shear layer
   just downstream of the step edge.

5. **Wall shear stress sign change** — On the bottom wall, τ_w changes from
   negative (in recirculation) to positive (after reattachment).

---

## Post-Processing Guide

### Viewing Results

```bash
# Open in ParaView
paraFoam

# Or create a .foam file
touch case.foam
paraview case.foam
```

### Finding the Reattachment Length

The reattachment point is where the wall shear stress on the bottom wall changes sign.
Several methods:

**Method 1: Wall shear stress (most accurate)**

The `wallShearStress` function object (configured in controlDict) writes the wall
shear stress at each write time. In ParaView:

1. Load the case and advance to the final time step
2. Apply "Extract Block" to select the bottomWall patch
3. Plot wallShearStress_x along the x-axis
4. The zero-crossing is the reattachment point

**Method 2: Velocity sign change**

Sample the x-velocity very close to the bottom wall:

```bash
# Add to system/controlDict functions or create system/sampleDict
postProcess -func "sampleDict" -latestTime
```

Create a `system/sampleDict`:
```
type            sets;
libs            (sampling);
writeControl    writeTime;
sets
(
    nearWall
    {
        type        uniform;
        axis        x;
        start       (0.0001 -0.0126 0);
        end         (0.381  -0.0126 0);
        nPoints     1000;
    }
);
fields          (U p k);
```

**Method 3: Streamlines**

In ParaView, create streamlines seeded near the step edge. The reattachment point
is visible as the location where the lowest streamline contacts the bottom wall.

### Velocity Profiles

Plot velocity profiles at key x/H locations for comparison with experimental data:

| Station | x/H  | Description                              |
|---------|------|------------------------------------------|
| 1       | 1    | Early recirculation zone                 |
| 2       | 4    | Mid-recirculation                        |
| 3       | 6    | Near reattachment                        |
| 4       | 8    | Just after reattachment                  |
| 5       | 15   | Recovery region                          |
| 6       | 25   | Far downstream (approaching equilibrium) |

### Skin Friction Coefficient

The skin friction coefficient on the bottom wall:

```
Cf = τ_w / (0.5 ρ U∞²)
```

Plot $C_f$ vs $x/H$. The reattachment point is at $C_f = 0$. Negative $C_f$ indicates
reverse flow.

### Pressure Coefficient

```
Cp = (p - p_ref) / (0.5 ρ U∞²)
```

Plot Cp along the top and bottom walls to visualise the pressure recovery through
the expansion.

### Checking y+

The `yPlus` function object writes y+ values at each write time:

```bash
# Post-process y+ at the latest time
postProcess -func yPlus -latestTime

# Check the range
grep "y+" log.simpleFoam | tail -5
```

Target values:
- Wall functions: y+ = 30–300 (acceptable), y+ = 30–100 (ideal)
- Low-Re wall treatment: y+ < 1

---

## Exercises

### Exercise 1: Turbulence Model Comparison

Compare k-ω SST with standard k-ε:

1. Copy the case directory
2. Change `constant/turbulenceProperties`:
   ```
   RASModel    kEpsilon;
   ```
3. Replace `0/omega` with `0/epsilon`:
   ```
   dimensions      [0 2 -3 0 0 0 0];
   internalField   uniform 0.0054;  // ε = Cμ k² / νt
   // walls: epsilonWallFunction
   ```
4. Run and compare reattachment lengths

Expected result: k-ε will predict Xr/H ≈ 4.5–5.0 (too short).

### Exercise 2: Reynolds Number Sensitivity

Run the case at different Reynolds numbers by changing the inlet velocity:

| Case | U∞ (m/s) | Re_H    | Expected Xr/H |
|------|----------|---------|----------------|
| A    | 2        | ~1,700  | ~5.0           |
| B    | 5        | ~4,200  | ~5.5           |
| C    | 10       | ~8,500  | ~6.0           |
| D    | 20       | ~17,000 | ~6.5           |
| E    | 50       | ~42,000 | ~6.5           |

Note: At low Re (< 400), the flow is laminar and the reattachment length depends
strongly on Re. At high Re, the reattachment length becomes nearly independent of Re.

### Exercise 3: Mesh Refinement Study

Perform a systematic grid convergence study:

| Level  | Upstream     | Downstream    | Total cells | Purpose           |
|--------|-------------|---------------|-------------|-------------------|
| Coarse | 40 × 20    | 150 × 30      | 9,800       | Quick test        |
| Medium | 80 × 40    | 300 × 40      | 27,200      | Baseline (this)   |
| Fine   | 160 × 80   | 600 × 80      | 108,800     | Grid convergence  |
| V.Fine | 320 × 160  | 1200 × 160    | 435,200     | Reference         |

Plot Xr/H vs total cell count. The solution should become grid-independent at the
Fine or Very Fine level.

### Exercise 4: Low-Re Wall Treatment

Modify the mesh and boundary conditions for resolved boundary layers (y+ < 1):

1. Increase wall-normal cell count by 10×
2. Apply first-cell grading to achieve y+ ≈ 0.5
3. Change wall functions:
   - `k`: `fixedValue uniform 0` (k = 0 at wall)
   - `omega`: `omegaWallFunction` (handles both low-Re and wall-function)
   - `nut`: `nutLowReWallFunction` or keep `nutkWallFunction` (auto-adapts)

### Exercise 5: Expansion Ratio Study

Investigate the effect of expansion ratio on reattachment length:

| ER  | Downstream height | Expected Xr/H |
|-----|-------------------|----------------|
| 1.2 | 1.2H              | ~4             |
| 1.5 | 1.5H              | ~5             |
| 2.0 | 2.0H (this case)  | ~6             |
| 3.0 | 3.0H              | ~8             |

Modify `system/blockMeshDict` to change the downstream channel height.

### Exercise 6: Inlet Profile Effects

Replace the uniform inlet velocity with a fully developed turbulent channel profile:

1. Run a separate precursor channel flow simulation
2. Map the outlet of the precursor to the inlet of the step case using `mapFields`
3. Compare results with the uniform inlet case

---

## References

### Primary Experimental References

1. **Driver, D.M. & Seegmiller, H.L.** (1985). "Features of a reattaching turbulent
   shear layer in divergent channel flow." *AIAA Journal*, 23(2), pp.163–171.
   — The definitive experimental study. Provides detailed velocity profiles,
   turbulence statistics, and skin friction measurements at $Re = 37{,}400$.

2. **Kim, J., Kline, S.J. & Johnston, J.P.** (1980). "Investigation of a reattaching
   turbulent shear layer: flow over a backward-facing step." *Journal of Fluids
   Engineering*, 102(3), pp.302–308.
   — Early comprehensive experimental study with extensive flow visualisation.

3. **Armaly, B.F., Durst, F., Pereira, J.C.F. & Schönung, B.** (1983). "Experimental
   and theoretical investigation of backward-facing step flow." *Journal of Fluid
   Mechanics*, 127, pp.473–496.
   — Covers a wide range of Reynolds numbers from laminar to turbulent.

### Turbulence Modelling References

4. **Menter, F.R.** (1994). "Two-equation eddy-viscosity turbulence models for
   engineering applications." *AIAA Journal*, 32(8), pp.1598–1605.
   — Original k-ω SST paper. Demonstrates superior performance for backward-facing
   step and other separated flows.

5. **Wilcox, D.C.** (2006). *Turbulence Modeling for CFD*. 3rd edition, DCW
   Industries.
   — Comprehensive reference on turbulence modelling, including detailed discussion
   of the backward-facing step as a validation case.

6. **Pope, S.B.** (2000). *Turbulent Flows*. Cambridge University Press.
   — Fundamental reference on turbulence physics, shear layers, and recirculating
   flows.

### Computational References

7. **Le, H., Moin, P. & Kim, J.** (1997). "Direct numerical simulation of turbulent
   flow over a backward-facing step." *Journal of Fluid Mechanics*, 330, pp.349–374.
   — DNS data providing "exact" reference for turbulence model validation.

8. **Patel, V.C., Rodi, W. & Scheuerer, G.** (1985). "Turbulence models for
   near-wall and low Reynolds number flows: a review." *AIAA Journal*, 23(9),
   pp.1308–1319.
   — Review of near-wall modelling approaches relevant to step flow.

### Cross-References to Tutorial Notes

For deeper understanding of the concepts used in this case, refer to the companion
notes in the `notes/` directory:

| Topic                  | Notes File                           | Relevant Sections            |
|------------------------|--------------------------------------|------------------------------|
| CFD fundamentals       | `notes/01_short_intro_to_cfd.md`     | Governing equations, FVM     |
| Case structure         | `notes/02_openfoam_cases.md`         | Directory layout, file types |
| Dictionary syntax      | `notes/03_openfoam_dictionaries.md`  | FoamFile headers, keywords   |
| Meshing strategies     | `notes/04_meshing.md`                | blockMesh, grading, quality  |
| Boundary conditions    | `notes/05_boundary_conditions.md`    | Wall functions, inlet/outlet |
| Turbulence models      | `notes/06_turbulence_models.md`      | k-ω SST, wall treatment     |
| Linear solvers         | `notes/09_linear_solvers.md`         | GAMG, smoothSolver, SIMPLE   |

---

*This case is part of the [OpenFOAM Tutorials](../../README.md) learning series.*
