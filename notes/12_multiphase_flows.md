# Multiphase Flow Modeling in OpenFOAM

## Lecture 7 — Models in OpenFOAM: The Multiphase Part

---

> **Cross-references:**
> - [01 — Short Intro to CFD](01_short_intro_to_cfd.md)
> - [03 — OpenFOAM Dictionaries](03_openfoam_dictionaries.md)
> - [04 — Meshing](04_meshing.md)
> - [05 — Boundary Conditions](05_boundary_conditions.md)
> - [06 — Turbulence Models](06_turbulence_models.md)
> - [08 — CFL Number](08_cfl_number.md)
> - [10 — icoFoam Solver Analysis](10_icofoam_solver_analysis.md)

---

## 1. What Is Multiphase Flow?

A **multiphase flow** is any flow in which two or more distinct phases — gas,
liquid, or solid — coexist and interact. The key word is *distinct*: each phase
has its own physical properties (density, viscosity), and there is a recognizable
interface or boundary between them.

### Types of Phase Combinations

| Combination     | Example                                  |
|-----------------|------------------------------------------|
| Gas – Liquid    | Ocean waves, rain, boiling, spray nozzles |
| Liquid – Liquid | Oil-water separation, emulsions           |
| Gas – Solid     | Fluidized beds, pneumatic conveying       |
| Liquid – Solid  | Sediment transport, slurry flow           |
| Gas – Liquid – Solid | Flotation cells, three-phase reactors |

### Real-World Applications

Multiphase flows are everywhere in engineering and nature:

- **Naval / marine:** Ship hulls moving through water with a free surface
  (see `projects/11_boat_hull_fixed/` and `projects/12_boat_hull_floating/`)
- **Civil / hydraulic:** Dam breaks, spillway flows, river flooding
  (see `projects/08_dam_break/`)
- **Chemical engineering:** Bubble columns, stirred reactors, distillation
- **Energy:** Steam generation in boilers, cavitation in turbines
- **Automotive:** Fuel injection sprays, water management on windshields
- **Environmental:** Oil spills, sediment transport in rivers

### Why Multiphase Flow Is Challenging

```
  ┌──────────────────────────────────────────────────────────────────────┐
  │              WHY IS MULTIPHASE CFD SO HARD?                         │
  ├──────────────────────────────────────────────────────────────────────┤
  │                                                                      │
  │  1. INTERFACE TRACKING                                               │
  │     The boundary between phases moves, deforms, breaks up,          │
  │     and merges — all within the simulation                           │
  │                                                                      │
  │  2. PROPERTY DISCONTINUITIES                                         │
  │     Density jumps by 1000× across a water-air interface             │
  │     (ρ_water = 1000 kg/m³,  ρ_air = 1 kg/m³)                       │
  │                                                                      │
  │  3. SURFACE TENSION                                                  │
  │     A force that acts only at the interface — not in the bulk        │
  │     Introduces capillary pressure jumps and parasitic currents       │
  │                                                                      │
  │  4. NUMERICAL STABILITY                                              │
  │     Large density ratios cause stiff systems                        │
  │     Interface smearing degrades accuracy if schemes are wrong        │
  │                                                                      │
  │  5. TIME STEPPING                                                    │
  │     The interface has its own CFL constraint (maxAlphaCo)           │
  │     Often much stricter than the flow CFL                            │
  │                                                                      │
  └──────────────────────────────────────────────────────────────────────┘
```

> **⚠️ Key insight:** Single-phase CFD is already complex (see
> [06 — Turbulence Models](06_turbulence_models.md)). Multiphase adds an
> entirely new layer: tracking where each phase *is* in space and time, and
> correctly handling the physics at the interface.

---

## 2. Multiphase Modeling Approaches

There is no single "best" method for all multiphase flows. The choice depends on
the physics: Is the interface sharp and well-defined? Are the dispersed entities
(bubbles, droplets, particles) small or large? How many are there?

### 2.1 Euler-Euler Methods

Both phases are treated as interpenetrating continua. Each cell may contain a
mixture of phases, described by a **volume fraction** field.

#### Volume of Fluid (VOF)

The most popular method in OpenFOAM for flows with a clear interface (free
surfaces, sloshing, waves). A single scalar field α tracks which phase occupies
each cell.

- **Strengths:** Inherently conserves mass; handles topology changes (breaking,
  merging) naturally; relatively simple to implement
- **Weakness:** Interface can smear over several cells without special treatment

#### Mixture Model

Treats the mixture as a single fluid with averaged properties. The phases share
a velocity field but can have different volume fractions. Simpler and cheaper
than full Euler-Euler, but less accurate for flows with strong phase interaction.

#### Full Euler-Euler (Two-Fluid Model)

Each phase has its own momentum equation. Phases interact through drag,
lift, virtual mass, and turbulent dispersion forces. Used for dispersed flows
(many bubbles or droplets) where resolving each interface is impractical.

### 2.2 Euler-Lagrange Methods

The continuous phase (e.g., gas) is solved on the Eulerian grid, while individual
particles, droplets, or bubbles are tracked as Lagrangian entities moving through
the domain. Each particle has its own position, velocity, and properties.

- **Strengths:** Natural for sprays, particle-laden flows; each particle can
  carry detailed information (size, temperature, species)
- **Weakness:** Expensive when particle count is very large (> 10⁶); coupling
  between phases can be complex

### 2.3 Level-Set Method

Uses a signed distance function $\phi$ to represent the interface ($\phi = 0$ at the
interface, $\phi > 0$ in one phase, $\phi < 0$ in the other). The interface is always
sharp, but mass conservation is not guaranteed without special corrections.
Less common in OpenFOAM than VOF.

### Comparison: When to Use Each Approach

| Criterion                  | VOF             | Euler-Euler (two-fluid) | Euler-Lagrange   | Level-Set       |
|----------------------------|-----------------|-------------------------|------------------|-----------------|
| Interface type             | Sharp, resolved | Dispersed (many bubbles)| Dispersed particles| Sharp, resolved |
| Mass conservation          | ✓ Excellent     | ✓ Good                  | ✓ Good           | ✗ Needs fixing  |
| Topology changes           | ✓ Automatic     | N/A                     | ✗ Difficult      | ✓ Automatic     |
| Computational cost         | Moderate        | High                    | Varies           | Moderate        |
| Number of dispersed bodies | Few (interfaces)| Many (10³–10⁶)         | Many (10³–10⁶)  | Few (interfaces)|
| OpenFOAM support           | ★★★ Excellent   | ★★ Good                 | ★★ Good          | ★ Limited       |
| Typical application        | Free surface    | Bubble columns          | Sprays, particles| Droplet dynamics|

> **💡 Tip:** For most engineering free-surface problems (waves, sloshing, dam
> breaks, ship hulls), **VOF with interFoam** is the go-to method in OpenFOAM.
> All three projects in this repository use it.

---

## 3. Volume of Fluid (VOF) Method — Deep Dive

The VOF method is the backbone of multiphase simulation in OpenFOAM. This section
explains how it works in detail.

### 3.1 The α (Alpha) Field

The volume fraction **α** is a scalar field defined on every cell:

```
  α = 1    →  Cell is 100% phase 1  (e.g., water)
  α = 0    →  Cell is 100% phase 2  (e.g., air)
  0 < α < 1  →  Cell contains the interface
```

Visually, here is what α looks like across a free surface:

```
  ┌──────────────────────────────────────────────────────────────────┐
  │                                                                  │
  │   AIR   α = 0     α = 0     α = 0     α = 0     α = 0          │
  │        ┌─────┐   ┌─────┐   ┌─────┐   ┌─────┐   ┌─────┐        │
  │        │     │   │     │   │     │   │     │   │     │        │
  │        │  0  │   │  0  │   │  0  │   │  0  │   │  0  │        │
  │        └─────┘   └─────┘   └─────┘   └─────┘   └─────┘        │
  │  ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─  INTERFACE  ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─   │
  │        ┌─────┐   ┌─────┐   ┌─────┐   ┌─────┐   ┌─────┐        │
  │        │     │   │     │   │     │   │     │   │     │        │
  │        │ 0.3 │   │ 0.6 │   │ 0.8 │   │ 0.5 │   │ 0.2 │        │
  │        └─────┘   └─────┘   └─────┘   └─────┘   └─────┘        │
  │        ┌─────┐   ┌─────┐   ┌─────┐   ┌─────┐   ┌─────┐        │
  │        │     │   │     │   │     │   │     │   │     │        │
  │        │  1  │   │  1  │   │  1  │   │  1  │   │  1  │        │
  │        └─────┘   └─────┘   └─────┘   └─────┘   └─────┘        │
  │  WATER  α = 1     α = 1     α = 1     α = 1     α = 1          │
  │                                                                  │
  └──────────────────────────────────────────────────────────────────┘
```

> **📌 Note:** The interface is never infinitely sharp in VOF — it always
> occupies at least 1–3 cells. The goal of interface sharpening schemes is
> to keep this transition as narrow as possible.

### 3.2 The Transport Equation

The α field is advected by the flow velocity **U**:

```
    ∂α
   ──── + ∇·(α U) = 0
    ∂t

  In words: the rate of change of α in a cell equals the net flux
  of α through its faces, carried by the velocity field.
```

This is a pure advection equation — there is no diffusion term. However,
numerical diffusion from the discretization can smear the interface, which is
why special treatment is needed.

### 3.3 Interface Sharpening — The MULES Algorithm

OpenFOAM uses the **MULES** (Multidimensional Universal Limiter with Explicit
Solution) algorithm to keep the interface sharp. MULES adds an artificial
**compression term** to the α transport equation:

```
    ∂α
   ──── + ∇·(α U) + ∇·(α(1−α) U_c) = 0
    ∂t
                      ↑
                      Compression velocity
                      Only active at the interface (where 0 < α < 1)
                      Pushes α toward 0 or 1
```

The compression velocity U_c acts perpendicular to the interface, squeezing
the transition region. The strength of compression is controlled by the
`cAlpha` parameter (set in `fvSolution`):

| `cAlpha` value | Effect                                       |
|----------------|----------------------------------------------|
| 0              | No compression (interface will smear)        |
| 1              | Standard compression (recommended default)   |
| 2              | Aggressive compression (very sharp, may oscillate) |

### 3.4 Property Interpolation

With α known in every cell, the mixture properties are computed as
volume-weighted averages:

```
    Density:     ρ  = α · ρ₁  +  (1 − α) · ρ₂
    Viscosity:   ν  = α · ν₁  +  (1 − α) · ν₂

    Example (water–air):
    ───────────────────
    Cell with α = 1.0:  ρ = 1.0 × 1000 + 0.0 × 1   = 1000 kg/m³  (pure water)
    Cell with α = 0.0:  ρ = 0.0 × 1000 + 1.0 × 1   = 1    kg/m³  (pure air)
    Cell with α = 0.5:  ρ = 0.5 × 1000 + 0.5 × 1   = 500.5 kg/m³ (interface)
```

> **⚠️ Warning:** The 1000:1 density ratio between water and air is one of the
> reasons multiphase simulations are more numerically challenging than
> single-phase ones. The momentum equation sees a thousand-fold jump in
> inertia across just a few cells.

---

## 4. OpenFOAM Multiphase Solvers

OpenFOAM provides a family of solvers for different multiphase scenarios. The
most commonly used is `interFoam`, but understanding the full family helps you
choose the right tool for your problem.

### Solver Overview

```
  ┌──────────────────────────────────────────────────────────────────────┐
  │                  OPENFOAM MULTIPHASE SOLVER FAMILY                  │
  ├──────────────────────────────────────────────────────────────────────┤
  │                                                                      │
  │  interFoam                                                           │
  │  ├── The workhorse: 2 incompressible, isothermal, immiscible fluids │
  │  ├── VOF with MULES interface sharpening                            │
  │  └── Used in: projects/08_dam_break/                                │
  │              projects/11_boat_hull_fixed/                            │
  │                                                                      │
  │  interDyMFoam                                                        │
  │  ├── interFoam + dynamic mesh capability                            │
  │  ├── Mesh moves/deforms each time step                              │
  │  └── Used in: projects/12_boat_hull_floating/ (6-DoF body motion)  │
  │                                                                      │
  │  multiphaseInterFoam                                                 │
  │  ├── More than 2 immiscible phases                                  │
  │  └── Each phase pair has its own surface tension                    │
  │                                                                      │
  │  interPhaseChangeFoam                                                │
  │  ├── VOF + phase change (cavitation, boiling, condensation)         │
  │  └── Mass transfer between phases                                   │
  │                                                                      │
  │  reactingMultiphaseEulerFoam                                         │
  │  ├── Full Euler-Euler multi-phase with reactions                    │
  │  ├── Heat transfer, mass transfer, chemical reactions               │
  │  └── The "Swiss Army knife" — complex but powerful                  │
  │                                                                      │
  └──────────────────────────────────────────────────────────────────────┘
```

### Solver Comparison Table

| Solver                        | Phases | Mesh    | Phase Change | Compressible | Typical Use Case                     |
|-------------------------------|--------|---------|--------------|--------------|--------------------------------------|
| `interFoam`                   | 2      | Static  | No           | No           | Dam break, sloshing, waves           |
| `interDyMFoam`                | 2      | Dynamic | No           | No           | Floating bodies, wave-structure      |
| `multiphaseInterFoam`         | N      | Static  | No           | No           | Oil-water-gas separation             |
| `interPhaseChangeFoam`        | 2      | Static  | Yes          | No           | Cavitation, boiling                  |
| `compressibleInterFoam`       | 2      | Static  | No           | Yes          | High-speed gas-liquid                |
| `reactingMultiphaseEulerFoam` | N      | Static  | Yes          | Yes          | Bubble columns, fluidized beds       |

> **💡 Tip:** Start with `interFoam`. It handles the vast majority of
> incompressible free-surface problems. Only move to more complex solvers when
> you need a feature that `interFoam` lacks (dynamic mesh, phase change, etc.).

---

## 5. Setting Up a Multiphase Case in OpenFOAM

This section walks through every file you need to configure for a VOF simulation
with `interFoam`. We use water and air as the two phases — the most common
combination. Examples are drawn from `projects/08_dam_break/`.

### 5.1 Transport Properties

**File:** `constant/transportProperties`

This file defines the two phases and their physical properties:

```c
// From: projects/08_dam_break/constant/transportProperties

phases (water air);     // Phase names — must match alpha field name

water
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1e-06;   // Kinematic viscosity [m²/s]
    rho             rho [ 1 -3 0 0 0 0 0 ] 1000;    // Density [kg/m³]
}

air
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1.48e-05; // Kinematic viscosity [m²/s]
    rho             rho [ 1 -3 0 0 0 0 0 ] 1;       // Density [kg/m³]
}

sigma           sigma [ 1 0 -2 0 0 0 0 ] 0.07;      // Surface tension [N/m]
                                                      // (water-air at ~20°C ≈ 0.072 N/m)
```

> **📌 Note:** The phase names (`water`, `air`) determine the name of the volume
> fraction field: `alpha.water`. You'll see this in the `0/` directory and in
> `setFieldsDict`.

### 5.2 Gravity

**File:** `constant/g`

Multiphase flows need gravity — it drives buoyancy and determines which phase
floats and which sinks.

```c
// From: projects/08_dam_break/constant/g

dimensions      [0 1 -2 0 0 0 0];   // m/s²
value           (0 -9.81 0);         // Gravity pointing downward (-y)
```

> **⚠️ Warning:** Forgetting to set gravity (or setting it to zero) in a
> multiphase simulation is a common mistake. Without gravity, the phases
> have no reason to separate, and a dam break becomes a pressure-driven
> expansion with no falling water.

### 5.3 Initializing the α Field with setFieldsDict

**File:** `system/setFieldsDict`

Before running the solver, you must define *where* each phase is. The
`setFields` utility reads this dictionary and sets the `alpha.water` field
accordingly.

```c
// From: projects/08_dam_break/system/setFieldsDict

defaultFieldValues
(
    volScalarFieldValue alpha.water 0    // Default: everything is air (α = 0)
);

regions
(
    boxToCell
    {
        box (0 0 -1) (0.1461 0.292 1);  // Rectangular region for water column
        fieldValues
        (
            volScalarFieldValue alpha.water 1  // Inside box: water (α = 1)
        );
    }
);
```

```
  Visualization of setFields result:
  ┌─────────────────────────────────────────────┐
  │                                             │
  │  α = 0 (air)                                │
  │                                             │
  │  ┌──────────┐                               │
  │  │          │                               │
  │  │  α = 1   │  ← Water column               │
  │  │ (water)  │     defined by box region      │
  │  │          │                               │
  │  │          │                               │
  ├──┴──────────┴───────────────────────────────┤
  │  WALL (bottom)                              │
  └─────────────────────────────────────────────┘
```

The workflow:

1. Define initial `alpha.water` in `0/alpha.water` (usually all zeros)
2. Run `setFields` to overwrite α based on `setFieldsDict`
3. Run `interFoam`

> **⚠️ Warning — Critical Mistake:** If you forget to run `setFields` before
> `interFoam`, the entire domain will be a single phase (all air or all water
> depending on your default). The simulation will run but produce meaningless
> results. Always verify your initial α field before solving.

### 5.4 Boundary Conditions for alpha.water

**File:** `0/alpha.water`

```c
// From: projects/08_dam_break/0/alpha.water

dimensions      [0 0 0 0 0 0 0];    // Dimensionless (volume fraction)
internalField   uniform 0;           // Will be overwritten by setFields

boundaryField
{
    walls
    {
        type            zeroGradient;  // No flux of α through walls
    }
    atmosphere
    {
        type            inletOutlet;   // Allows air to enter/leave at top
        inletValue      uniform 0;     // If flow enters, it brings air (α = 0)
        value           uniform 0;
    }
    defaultFaces
    {
        type            empty;         // 2D case — front/back faces
    }
}
```

Common boundary condition choices for `alpha.water`:

| Boundary Type     | BC for α                   | Explanation                              |
|-------------------|----------------------------|------------------------------------------|
| Solid wall        | `zeroGradient`             | No phase flux through walls              |
| Open atmosphere   | `inletOutlet` (value = 0)  | Air can enter; water can leave           |
| Inlet (water)     | `fixedValue` (value = 1)   | Pure water enters                        |
| Inlet (air)       | `fixedValue` (value = 0)   | Pure air enters                          |
| Outlet            | `inletOutlet` (value = 0)  | Default to air if backflow occurs        |
| 2D front/back     | `empty`                    | Standard for 2D cases                    |
| Symmetry plane    | `symmetryPlane`            | Mirror condition                         |

### 5.5 controlDict Settings for Multiphase

**File:** `system/controlDict`

Multiphase simulations almost always use **adaptive time stepping** to keep the
interface CFL number under control:

```c
// From: projects/08_dam_break/system/controlDict

application     interFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1;                 // 1 second of physical time

deltaT          0.001;             // Initial time step (will be adjusted)

adjustTimeStep  yes;               // ← KEY: Enable adaptive time stepping
maxCo           0.5;               // Maximum flow Courant number
maxAlphaCo      0.5;               // Maximum interface Courant number
maxDeltaT       0.01;              // Maximum allowed time step

writeControl    adjustableRunTime;
writeInterval   0.05;              // Write every 0.05 seconds
purgeWrite      0;

functions
{
    // Optional: track interface position, forces, etc.
}
```

> **💡 Tip:** Note the **two** Courant number limits: `maxCo` for the flow and
> `maxAlphaCo` for the interface. The solver uses whichever gives the smaller
> time step. See [08 — CFL Number](08_cfl_number.md) for a thorough explanation.

### 5.6 Numerical Schemes (fvSchemes)

**File:** `system/fvSchemes`

The divergence schemes for the α transport equation are critical for interface
sharpness:

```c
// From: projects/08_dam_break/system/fvSchemes

ddtSchemes
{
    default         Euler;            // First-order time (standard for VOF)
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    // Momentum equation
    div(rhoPhi,U)   Gauss linearUpwind grad(U);

    // Alpha transport — THIS IS THE CRITICAL PART
    div(phi,alpha)  Gauss vanLeer;                    // Bounded, TVD scheme
    div(phirb,alpha) Gauss interfaceCompression;      // Compression term

    // Turbulence (if used)
    div(phi,k)      Gauss upwind;
    div(phi,omega)  Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
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
```

The two key divergence entries for VOF:

| Entry                     | Scheme                     | Purpose                        |
|---------------------------|----------------------------|--------------------------------|
| `div(phi,alpha)`          | `Gauss vanLeer`            | Advection of α (bounded, TVD) |
| `div(phirb,alpha)`        | `Gauss interfaceCompression`| Interface compression term    |

> **⚠️ Warning — Interface Smearing:** Using `Gauss upwind` or `Gauss linear`
> for `div(phi,alpha)` will smear the interface over many cells. Always use a
> bounded scheme like `vanLeer`, `vanAlbada`, or `MUSCL` for the alpha transport.

### 5.7 Solver Settings (fvSolution)

**File:** `system/fvSolution`

```c
// From: projects/08_dam_break/system/fvSolution

solvers
{
    "alpha.water.*"
    {
        // MULES settings for alpha
        nAlphaCorr      2;          // Number of alpha correction loops
        nAlphaSubCycles 1;          // Sub-cycling: solve α multiple times per Δt
                                     // Increase to 2-4 if interface is unstable
        cAlpha          1;          // Compression factor (0=none, 1=standard, 2=aggressive)
    }

    pcorr
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-5;
        relTol          0;
    }

    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-7;
        relTol          0.01;
    }

    p_rghFinal
    {
        $p_rgh;
        relTol          0;         // Tighter tolerance for final corrector
    }

    "U.*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-6;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor   no;         // Often "no" for multiphase
    nOuterCorrectors    1;          // 1 = PISO mode (standard for interFoam)
    nCorrectors         3;          // Pressure-velocity corrections
    nNonOrthogonalCorrectors 1;     // For non-orthogonal meshes
}
```

> **📌 Note:** In multiphase VOF, the pressure variable is `p_rgh` (pressure
> minus the hydrostatic component ρgh), not plain `p`. This removes the large
> hydrostatic gradient from the pressure equation, improving numerical accuracy.

Key parameters for the α solver:

| Parameter         | Typical Value | Effect                                        |
|-------------------|---------------|-----------------------------------------------|
| `nAlphaCorr`      | 1–3           | More corrections → sharper interface, slower   |
| `nAlphaSubCycles` | 1–4           | Sub-cycles within each time step for stability |
| `cAlpha`          | 1             | Interface compression strength                 |

---

## 6. Surface Tension

### 6.1 The Continuum Surface Force (CSF) Model

In VOF, the interface is not a discrete boundary — it's a smooth transition over
a few cells. Surface tension is modeled as a **body force** distributed over this
transition region, using the CSF model by Brackbill et al. (1992):

```
    F_σ = σ · κ · ∇α

    where:
      σ  = surface tension coefficient [N/m]
      κ  = interface curvature = −∇·(∇α / |∇α|)
      ∇α = gradient of α (points from light to heavy phase)
```

The force is concentrated where ∇α is large (i.e., at the interface) and
vanishes in the bulk of each phase.

### 6.2 The Sigma Parameter

Surface tension is set in `constant/transportProperties`:

```c
sigma           sigma [ 1 0 -2 0 0 0 0 ] 0.07;   // Surface tension [N/m]
```

Common values at ~20°C:

| Interface       | σ (N/m)  |
|-----------------|----------|
| Water – Air     | 0.072    |
| Mercury – Air   | 0.47     |
| Oil – Air       | 0.02–0.03|
| Oil – Water     | 0.02–0.05|
| Molten steel – Gas | 1.5–1.8 |

### 6.3 When Does Surface Tension Matter?

The **Weber number** determines the relative importance of inertia vs. surface
tension:

```
           ρ · U² · L
    We  =  ───────────
               σ

    We >> 1  →  Inertia dominates, surface tension is minor
                (e.g., large dam break, ship hull flows)
    We ~ 1   →  Both matter (e.g., jet breakup, droplet impact)
    We << 1  →  Surface tension dominates (e.g., capillary flows, micro-droplets)
```

| Flow Type              | Typical We | Surface Tension Effect       |
|------------------------|------------|------------------------------|
| Dam break              | > 1000     | Negligible                   |
| Ship hull waves        | > 100      | Minor                        |
| Spray droplet breakup  | 10–100     | Important                    |
| Dripping faucet        | 1–10       | Dominant                     |
| Microfluidics          | < 1        | Dominant                     |

> **⚠️ Warning — Parasitic Currents:** The CSF model can generate spurious
> velocities ("parasitic currents") near the interface, especially in cases
> where surface tension dominates. These are numerical artifacts caused by
> errors in curvature estimation. Finer meshes and higher-order gradient
> schemes help reduce them.

---

## 7. The CFL Condition for Multiphase Flows

Time stepping in multiphase simulations requires special attention because the
interface transport has its own stability constraint. See
[08 — CFL Number](08_cfl_number.md) for the general theory.

### 7.1 Two Courant Numbers

In a multiphase simulation with `interFoam`, OpenFOAM reports **two** Courant
numbers every time step:

```
  Courant Number mean: 0.0152 max: 0.423
  Interface Courant Number mean: 0.00243 max: 0.387
```

```
  ┌──────────────────────────────────────────────────────────────────┐
  │            TWO CFL CONSTRAINTS IN MULTIPHASE FLOWS              │
  ├──────────────────────────────────────────────────────────────────┤
  │                                                                  │
  │  Flow CFL:                                                       │
  │                                                                  │
  │         |U| · Δt                                                 │
  │  Co  = ──────────   controlled by  maxCo                        │
  │           Δx                                                     │
  │                                                                  │
  │  Interface CFL:                                                  │
  │                                                                  │
  │           |U_interface| · Δt                                     │
  │  Co_α  = ──────────────────  controlled by  maxAlphaCo          │
  │                Δx                                                │
  │                                                                  │
  │  The solver uses:  Δt = min(Δt_Co, Δt_CoAlpha, maxDeltaT)      │
  │                                                                  │
  └──────────────────────────────────────────────────────────────────┘
```

### 7.2 Why the Interface Needs Tighter Control

The α transport equation is **explicit** — solved with MULES outside the
pressure-velocity coupling. Explicit methods are conditionally stable: if
information travels farther than one cell per time step, the solution becomes
unstable and the interface breaks apart.

Recommended settings:

| Setting        | Conservative | Standard | Aggressive (risky) |
|----------------|--------------|----------|--------------------|
| `maxCo`        | 0.3          | 0.5      | 0.9                |
| `maxAlphaCo`   | 0.3          | 0.5      | 0.9                |
| `maxDeltaT`    | 0.001        | 0.01     | 0.1                |

> **💡 Tip:** Always start conservative (maxCo = maxAlphaCo = 0.5). Only relax
> these if the simulation is stable and you need faster turnaround. If the
> interface starts looking ragged or the solver diverges, tighten them.

---

## 8. Real Case Studies from This Repository

### 8.1 Dam Break — `projects/08_dam_break/`

```
  ┌─────────────────────────────────────────────────────────────┐
  │                     DAM BREAK SETUP                         │
  │                                                             │
  │           atmosphere (open top)                             │
  │  ┌───────────────────────────────────────────────────────┐  │
  │  │                                                       │  │
  │  │  AIR (α = 0)                                          │  │
  │  │                                                       │  │
  │  │  ┌──────────┐                                         │  │
  │  │  │ WATER    │                                         │  │
  │  │  │ (α = 1)  │              ┌──┐                       │  │
  │  │  │          │              │  │ ← Obstacle             │  │
  │  │  │          │              │  │                        │  │
  │  └──┴──────────┴──────────────┴──┴───────────────────────┘  │
  │           walls (bottom + sides)                            │
  └─────────────────────────────────────────────────────────────┘
```

**What it demonstrates:**
- Basic VOF setup with `interFoam`
- Gravity-driven flow (water column collapses under gravity)
- Laminar flow (no turbulence model needed)
- `setFieldsDict` to initialize the water column
- Adaptive time stepping with `maxCo` and `maxAlphaCo`

**Key settings:**
- Solver: `interFoam`
- Phases: water ($\rho$ = 1000, $\nu$ = 1e-06) and air ($\rho$ = 1, $\nu$ = 1.48e-05)
- Time: adaptive, maxCo = 0.5, maxAlphaCo = 0.5
- Mesh: structured (blockMesh)

**What the student learns:** How to set up a complete multiphase case from
scratch, configure all the necessary files, and visualize the free surface
evolution in ParaView.

### 8.2 Boat Hull (Fixed) — `projects/11_boat_hull_fixed/`

**What it demonstrates:**
- VOF with a complex geometry (boat hull STL)
- `snappyHexMesh` for meshing around a hull
  (see [04 — Meshing](04_meshing.md) for snappyHexMesh details)
- Free surface interaction with a solid body
- Steady-state-like multiphase (hull fixed in a flow)

**Key settings:**
- Solver: `interFoam`
- Geometry: imported from `boat_hull.stl` via snappyHexMesh
- Mesh: snappyHexMesh with surface refinement near hull and waterline
- Phases: water and air (same properties as dam break)
- Turbulence: likely k-ω SST (see [06 — Turbulence Models](06_turbulence_models.md))

**What the student learns:** How to combine VOF multiphase with complex geometry
meshing — the realistic workflow for naval CFD. The interaction between mesh
refinement at the free surface and VOF accuracy is a key lesson.

### 8.3 Boat Hull (Floating) — `projects/12_boat_hull_floating/`

**What it demonstrates:**
- Dynamic mesh with VOF: `interDyMFoam`
- 6 degrees of freedom (6-DoF) rigid body motion
- The hull moves in response to hydrodynamic forces (heave, pitch, etc.)
- Mesh deforms to follow the body motion

**Key settings:**
- Solver: `interDyMFoam` (interFoam + dynamic mesh)
- 6-DoF motion: `dynamicMeshDict` specifies the rigid body solver
- Motion constraints: may restrict to heave + pitch only
- Mesh motion: typically `displacementLaplacian` or similar

**What the student learns:** The most advanced case in the series — how to
couple fluid forces with body motion in a multiphase context. This is the
basis for seakeeping, wave energy, and offshore engineering simulations.

> **📌 Note:** These three projects form a progression of increasing complexity:
>
> | Project                 | Solver         | Geometry   | Mesh       | Motion  |
> |-------------------------|----------------|------------|------------|---------|
> | `08_dam_break`          | `interFoam`    | Simple box | blockMesh  | None    |
> | `11_boat_hull_fixed`    | `interFoam`    | STL hull   | snappyHex  | None    |
> | `12_boat_hull_floating` | `interDyMFoam` | STL hull   | snappyHex  | 6-DoF   |

---

## 9. Common Pitfalls & Troubleshooting

### Pitfall 1: Interface Smearing

> **Problem:** The interface between phases spreads over many cells, making
> the free surface look like a fuzzy gradient instead of a sharp boundary.
>
> **Cause:** Wrong divergence scheme for `div(phi,alpha)` — typically using
> `Gauss upwind` or `Gauss linear` instead of a bounded scheme.
>
> **Fix:** Use `Gauss vanLeer` for `div(phi,alpha)` and
> `Gauss interfaceCompression` for `div(phirb,alpha)`. Also check that
> `cAlpha` is set to 1 (not 0) in `fvSolution`.

### Pitfall 2: Forgetting to Run setFields

> **Problem:** The simulation runs but the entire domain is one phase — no
> interface, no two-phase behavior.
>
> **Cause:** You defined `setFieldsDict` but never ran the `setFields` utility
> before starting the solver.
>
> **Fix:** Always run `setFields` after `blockMesh` (or `snappyHexMesh`) and
> before `interFoam`. Check `0/alpha.water` in ParaView to verify the initial
> condition looks correct.

### Pitfall 3: CFL Too High for Interface Sharpness

> **Problem:** The interface becomes ragged, develops spurious droplets, or the
> solver diverges.
>
> **Cause:** `maxCo` or `maxAlphaCo` set too high, allowing the interface to
> jump multiple cells per time step.
>
> **Fix:** Reduce both `maxCo` and `maxAlphaCo` to 0.5 or lower. Increase
> `nAlphaSubCycles` to 2 for additional stability.

### Pitfall 4: Wrong Boundary Conditions for Alpha

> **Problem:** Water appears or disappears at boundaries unexpectedly, or the
> simulation crashes near open boundaries.
>
> **Cause:** Using `fixedValue 1` (water) on a boundary that should allow air
> to enter, or using `zeroGradient` on an open atmosphere boundary.
>
> **Fix:** Use `inletOutlet` with `inletValue 0` on atmosphere/outlet patches.
> Use `zeroGradient` on solid walls. See the BC table in § 5.4.

### Pitfall 5: Density Ratio Instability

> **Problem:** Solver diverges or produces oscillations, especially with very
> large density ratios (e.g., water/air = 1000:1, liquid metal/gas > 7000:1).
>
> **Cause:** The momentum equation becomes very stiff when adjacent cells have
> vastly different densities.
>
> **Fix:** Tighten time stepping (`maxCo` = 0.3), increase `nAlphaSubCycles`,
> ensure under-relaxation is appropriate. For extreme density ratios, consider
> using the `GAMG` preconditioner for the pressure equation.

### Pitfall 6: Parasitic Currents from Surface Tension

> **Problem:** Unphysical velocities appear near the interface even in cases
> that should be static (e.g., a stationary droplet).
>
> **Cause:** The CSF model computes interface curvature from the α gradient,
> which introduces numerical errors that manifest as spurious velocities.
>
> **Fix:** Refine the mesh near the interface, use higher-order gradient schemes,
> and ensure the α field is sharp (`cAlpha` = 1). For very surface-tension-
> dominated flows, consider using `isoAdvector` (geometric VOF) if available
> in your OpenFOAM version.

### Pitfall 7: Gravity Direction Wrong or Missing

> **Problem:** Water and air do not separate properly, or water "floats" upward.
>
> **Cause:** The `constant/g` file has the wrong sign or magnitude, or is missing.
>
> **Fix:** Check that `g` points in the correct direction with the correct
> magnitude: `(0 -9.81 0)` for gravity in the −y direction.

### Troubleshooting Checklist

```
  ┌──────────────────────────────────────────────────────────────────┐
  │              MULTIPHASE DEBUGGING CHECKLIST                     │
  ├──────────────────────────────────────────────────────────────────┤
  │                                                                  │
  │  □  Did you run setFields?                                      │
  │  □  Is alpha.water correct? (Check in ParaView at t = 0)       │
  │  □  Is gravity set correctly in constant/g?                     │
  │  □  Are phase densities and viscosities physically correct?     │
  │  □  Is adjustTimeStep = yes?                                    │
  │  □  Are maxCo and maxAlphaCo ≤ 0.5?                            │
  │  □  Is div(phi,alpha) using vanLeer (not upwind)?              │
  │  □  Is div(phirb,alpha) using interfaceCompression?            │
  │  □  Is cAlpha = 1 in fvSolution?                               │
  │  □  Are BCs for alpha correct? (zeroGradient on walls,         │
  │       inletOutlet on atmosphere/outlets)                        │
  │  □  Is p_rgh used as the pressure variable (not p)?            │
  │  □  Does checkMesh report acceptable quality?                  │
  │                                                                  │
  └──────────────────────────────────────────────────────────────────┘
```

---

## 10. Quick Reference

```
  ╔══════════════════════════════════════════════════════════════════════╗
  ║                  MULTIPHASE FLOW CHEAT SHEET                       ║
  ╠══════════════════════════════════════════════════════════════════════╣
  ║                                                                      ║
  ║  DEFAULT SOLVER:  interFoam                                          ║
  ║    Two incompressible, immiscible fluids with VOF                    ║
  ║                                                                      ║
  ║  NEED DYNAMIC MESH?  →  interDyMFoam                                ║
  ║  NEED PHASE CHANGE?  →  interPhaseChangeFoam                        ║
  ║  NEED > 2 PHASES?    →  multiphaseInterFoam                         ║
  ║                                                                      ║
  ║  KEY FILES:                                                          ║
  ║    constant/transportProperties  → phases, ρ, ν, σ                  ║
  ║    constant/g                    → gravity vector                    ║
  ║    system/setFieldsDict          → initial α distribution           ║
  ║    0/alpha.water                 → α boundary conditions            ║
  ║                                                                      ║
  ║  CRITICAL SCHEMES:                                                   ║
  ║    div(phi,alpha)    →  Gauss vanLeer                               ║
  ║    div(phirb,alpha)  →  Gauss interfaceCompression                  ║
  ║                                                                      ║
  ║  TIME STEPPING:                                                      ║
  ║    adjustTimeStep  yes                                               ║
  ║    maxCo           0.5     (flow Courant)                            ║
  ║    maxAlphaCo      0.5     (interface Courant)                       ║
  ║                                                                      ║
  ║  MULES SETTINGS:                                                     ║
  ║    cAlpha           1       (compression, 0=off, 1=standard)        ║
  ║    nAlphaCorr       2       (correction loops)                      ║
  ║    nAlphaSubCycles  1       (sub-cycling, increase for stability)   ║
  ║                                                                      ║
  ║  PRESSURE VARIABLE:  p_rgh  (not p!)                                ║
  ║                                                                      ║
  ║  WORKFLOW:                                                           ║
  ║    blockMesh → snappyHexMesh (if complex) → setFields → interFoam  ║
  ║                                                                      ║
  ╚══════════════════════════════════════════════════════════════════════╝
```

### Solver Selection Flowchart

```
  START: Is the flow multiphase?
    │
    ├── NO  → Use single-phase solver
    │         (simpleFoam, pimpleFoam, icoFoam, etc.)
    │
    └── YES → How many phases?
              │
              ├── 2 phases (most common)
              │     │
              │     ├── Incompressible, no phase change?
              │     │     │
              │     │     ├── Static mesh  →  interFoam  ★
              │     │     └── Moving body  →  interDyMFoam
              │     │
              │     ├── Phase change (cavitation/boiling)?
              │     │     └── interPhaseChangeFoam
              │     │
              │     └── Compressible?
              │           └── compressibleInterFoam
              │
              └── > 2 phases
                    │
                    ├── Immiscible, no reactions
                    │     └── multiphaseInterFoam
                    │
                    └── With reactions / heat transfer
                          └── reactingMultiphaseEulerFoam
```

---

## 11. Where to Go Next

| Topic | Note |
|-------|------|
| CFD fundamentals | [01 — Short Intro to CFD](01_short_intro_to_cfd.md) |
| Case structure for multiphase | [02 — OpenFOAM Cases](02_openfoam_cases.md) |
| Dictionary files (fvSchemes, fvSolution) | [03 — OpenFOAM Dictionaries](03_openfoam_dictionaries.md) |
| Meshing (blockMesh, snappyHexMesh) | [04 — Meshing](04_meshing.md) |
| Boundary conditions for α, U, p | [05 — Boundary Conditions](05_boundary_conditions.md) |
| Turbulence in multiphase flows | [06 — Turbulence Models](06_turbulence_models.md) |
| Parallel runs for large cases | [07 — Parallelization](07_parallelization.md) |
| CFL number and time stepping | [08 — CFL Number](08_cfl_number.md) |
| Pressure solvers (p_rgh) | [09 — Linear Solvers](09_linear_solvers.md) |
| Solver internals (PISO algorithm) | [10 — icoFoam Solver Analysis](10_icofoam_solver_analysis.md) |
| **Hands-on:** Dam break (VOF basics) | `projects/08_dam_break/` |
| **Hands-on:** Fixed hull (VOF + snappyHex) | `projects/11_boat_hull_fixed/` |
| **Hands-on:** Floating hull (dynamic mesh + 6-DoF) | `projects/12_boat_hull_floating/` |

---

*Last updated: 2025. Part of the [OpenFOAM-Tutorials](../README.md) repository.*
