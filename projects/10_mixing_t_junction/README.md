# Mixing T-Junction — Hot/Cold Water Thermal Mixing

| Detail | Value |
|---|---|
| **Difficulty** | ⭐⭐⭐ Intermediate–Advanced |
| **Solver** | `buoyantSimpleFoam` (steady-state, buoyant, compressible with energy) |
| **Algorithm** | SIMPLE |
| **Turbulence** | k-$\varepsilon$ (RAS) |
| **Mesh** | `blockMesh` (multi-block structured T-shape) |
| **Physics** | Forced convection, scalar transport, thermal mixing |
| **Key Learning** | Temperature field, thermal stratification, multi-block meshing |
| **OpenFOAM Version** | 6+ / v2012+ |

---

## Table of Contents

1. [Problem Description](#problem-description)
2. [Physics and Theory](#physics-and-theory)
3. [Practical Relevance](#practical-relevance)
4. [Case Structure](#case-structure)
5. [Geometry and Mesh](#geometry-and-mesh)
6. [Boundary Conditions](#boundary-conditions)
7. [Fluid and Thermophysical Properties](#fluid-and-thermophysical-properties)
8. [Turbulence Setup](#turbulence-setup)
9. [Numerical Schemes and Solver Settings](#numerical-schemes-and-solver-settings)
10. [How to Run](#how-to-run)
11. [Expected Results](#expected-results)
12. [Post-Processing Guide](#post-processing-guide)
13. [Exercises and Experiments](#exercises-and-experiments)
14. [References](#references)
15. [Related Notes in This Repository](#related-notes-in-this-repository)

---

## Problem Description

This tutorial simulates the mixing of **hot and cold water streams** in a
T-junction pipe. Hot water enters through the main horizontal pipe, and cold
water enters from a vertical branch pipe. The two streams meet at the junction
and mix as they flow toward the outlet.

```
                     Cold water (300K)
                          ↓ ↓ ↓
                     branchInlet
                    ┌─────┴───┴─────┐
                    │               │
                    │  Branch Pipe  │  0.05m wide
                    │   (0.5 m/s)   │  0.1m tall
  ┌─────────────────┴───────────────┴──────────────────────────────┐
  │                 │               │                              │
  │  Hot water      │   Junction    │       Mixed flow region      │
  │  (350K)         │   ~~~T~~~     │   ~~~Temperature mixing~~~   │→→→ outlet
  │  →→→ (1.0 m/s) │   ~~~T~~~     │   Thermal stratification     │    (mixed T)
  │                 │               │                              │
  └─────────────────┴───────────────┴──────────────────────────────┘
  mainInlet         x=0.125  x=0.175                          x=0.5
  (x=0)
```

### Key Features of This Flow

```
  Side view of the T-junction flow physics:
  ==========================================

         Cold (300K)
           ↓ ↓ ↓
      ┌────┴───┴────┐
      │  ↓   ↓   ↓  │
      │  ↓   ↓   ↓  │
  ────┴──↓───↓───↓──┴────────────────────────
  →→→→→→→↘   ↓   ↙→→→→→→→→→→→→→→→→→→→→→→→→→
  →→→ Hot ↘  ↓  ↙  ~~~Mixed~~~   →→→→→→→→→→→
  →→→→→→→→→↘ ↓ ↙→→→→→→→→→→→→→→→→→→→→→→→→→→→
  ──────────────────────────────────────────── 
             ◜◝
          Recirculation
            zone

  Temperature stratification downstream:

  Top:     ≈300K  ░░░░░░░░░░  (cold stream hugs top wall)
  Middle:  ≈325K  ▒▒▒▒▒▒▒▒▒▒  (mixing zone)
  Bottom:  ≈350K  ▓▓▓▓▓▓▓▓▓▓  (hot stream along bottom)

  Far downstream → fully mixed ≈ 330K
```

- **Thermal mixing**: Two streams at different temperatures interact and
  exchange heat through convective transport and turbulent diffusion
- **Thermal stratification**: In the near-junction region, the cold stream
  rides over the hot stream due to density differences and momentum
- **Recirculation zone**: A small recirculation develops at the junction
  where the branch flow deflects the main flow
- **Temperature gradients**: Sharp thermal gradients near the junction
  that smooth out further downstream

### Why Study T-Junction Mixing?

T-junctions are among the most common piping geometries in engineering.
Understanding how temperature, velocity, and turbulence interact at the
junction is essential for:

- Predicting **outlet temperature** (energy balance)
- Assessing **thermal fatigue** risk from temperature fluctuations
- Designing **mixing devices** for chemical processes
- Evaluating **HVAC system** performance

---

## Physics and Theory

### Governing Equations

`buoyantSimpleFoam` solves the steady-state Reynolds-Averaged Navier-Stokes
(RANS) equations with the energy equation for a compressible fluid:

**Continuity (mass conservation):**

```
∂ρ/∂t + ∇·(ρU) = 0
```

**Momentum (RANS):**

```
∇·(ρUU) = -∇p + ∇·(μ_eff ∇U) + ρg

where μ_eff = μ + μ_t  (molecular + turbulent viscosity)
```

**Energy (enthalpy form):**

```
∇·(ρUh) = ∇·(α_eff ∇h) + ∂p/∂t

where α_eff = μ/Pr + μ_t/Pr_t  (effective thermal diffusivity)
      h = Cp × T               (sensible enthalpy)
```

**Buoyancy coupling:**

```
p_rgh = p - ρg·r

where p_rgh is the pressure minus the hydrostatic component
      g     is the gravitational acceleration vector
      r     is the position vector
```

### Reynolds Number

For the main pipe:

```
Re_main = U_main × D_h / ν
        = 1.0 × 0.1 / 1e-6
        = 100,000

where U_main = 1.0 m/s       (main inlet velocity)
      D_h    = 0.1 m         (hydraulic diameter ≈ pipe height for 2D)
      ν      = 1e-6 m²/s     (kinematic viscosity of water)
```

For the branch pipe:

```
Re_branch = U_branch × D_h / ν
          = 0.5 × 0.05 / 1e-6
          = 25,000

where U_branch = 0.5 m/s     (branch inlet velocity)
      D_h      = 0.05 m      (branch hydraulic diameter)
```

Both Reynolds numbers are well above the laminar-turbulent transition
(Re $\approx$ 2300 for pipe flow), justifying the use of a turbulence model.

### Thermal Mixing Theory

**Energy balance at the junction:**

```
T_mixed = (ṁ_main × Cp × T_main + ṁ_branch × Cp × T_branch) /
          (ṁ_main × Cp + ṁ_branch × Cp)

For equal Cp:
T_mixed = (ṁ_main × T_main + ṁ_branch × T_branch) / (ṁ_main + ṁ_branch)
```

**Mass flow rates (per unit depth for 2D):**

```
ṁ_main   = ρ × U_main   × A_main   = ρ × 1.0 × 0.1   = 0.1ρ
ṁ_branch = ρ × U_branch × A_branch = ρ × 0.5 × 0.05  = 0.025ρ

Flow ratio R = ṁ_branch / ṁ_main = 0.025 / 0.1 = 0.25
```

**Ideal mixed temperature:**

```
T_mixed = (0.1 × 350 + 0.025 × 300) / (0.1 + 0.025)
        = (35 + 7.5) / 0.125
        = 42.5 / 0.125
        = 340 K
```

### Mixing Quality Metrics

**Mixing coefficient (0 = unmixed, 1 = perfectly mixed):**

```
η_mix = 1 - σ_T / σ_T,unmixed

where σ_T         = standard deviation of T across a cross-section
      σ_T,unmixed = std dev at the junction entrance (maximum non-uniformity)
```

**Temperature variance decay:**

```
σ²_T(x) ∝ exp(-x / L_mix)

where L_mix is the characteristic mixing length
```

Good mixing requires L_mix $\approx$ 10–20 hydraulic diameters downstream
of the junction for turbulent flow.

### Turbulent Prandtl Number

The turbulent thermal diffusivity is modelled as:

```
α_t = ν_t / Pr_t

where Pr_t ≈ 0.85  (turbulent Prandtl number, default in OpenFOAM)
      ν_t           (turbulent kinematic viscosity from k-ε model)
```

This means turbulent heat transport is enhanced relative to momentum
transport. The k-$\varepsilon$ model provides $\nu$_t, and `alphat` provides the
wall-bounded turbulent thermal diffusivity.

---

## Practical Relevance

### Nuclear Piping — Thermal Fatigue

T-junctions in nuclear power plant piping carry coolant at different
temperatures. Temperature fluctuations at the junction wall cause cyclic
thermal stresses, leading to **thermal fatigue** and potential pipe failure.
The Civaux-1 incident (1998, France) involved cracking in a residual heat
removal line T-junction, motivating extensive research into thermal
mixing and striping phenomena.

```
  Nuclear thermal fatigue mechanism:
  ┌──────────────────────────────────────────┐
  │  Pipe wall                               │
  │  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ │
  │  ───── T fluctuations ──────────────     │
  │  ░hot░ ▒cold▒ ░hot░ ▒cold▒ ░hot░        │  ← thermal striping
  │  ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ │
  │  Pipe wall                               │
  │                                          │
  │  Cyclic ΔT → thermal stress → fatigue    │
  │  crack initiation → propagation → leak   │
  └──────────────────────────────────────────┘
```

### HVAC Systems

In heating, ventilation, and air conditioning systems, hot and cold air
streams mix at T-junctions in ductwork. Understanding mixing quality
and temperature uniformity at downstream locations is critical for
occupant comfort and energy efficiency.

### Chemical Process Mixing

Chemical reactors and mixing vessels often use T-junction configurations
to combine reactants. The mixing quality directly affects reaction
rates, selectivity, and product quality. Poor mixing can lead to
hot spots, side reactions, and safety hazards.

---

## Case Structure

```
09_mixing_t_junction/
├── 0/                              ← Initial & boundary conditions
│   ├── U                           ← Velocity field
│   ├── p                           ← Pressure (total)
│   ├── p_rgh                       ← Pressure (minus hydrostatic)
│   ├── T                           ← Temperature field ★ KEY FIELD
│   ├── k                           ← Turbulent kinetic energy
│   ├── epsilon                     ← Turbulent dissipation rate
│   ├── nut                         ← Turbulent viscosity
│   └── alphat                      ← Turbulent thermal diffusivity
├── constant/                       ← Physical properties
│   ├── g                           ← Gravitational acceleration
│   ├── thermophysicalProperties    ← Thermo model & fluid properties
│   ├── transportProperties         ← Transport model
│   ├── turbulenceProperties        ← Turbulence model selection
│   └── polyMesh/                   ← Mesh (generated by blockMesh)
├── system/                         ← Solver & numerical settings
│   ├── blockMeshDict               ← Mesh definition (T-shape)
│   ├── controlDict                 ← Run control
│   ├── fvSchemes                   ← Discretization schemes
│   └── fvSolution                  ← Linear solvers & SIMPLE
├── Allrun                          ← Run script
├── Allclean                        ← Clean script
└── README.md                       ← This file
```

### What Makes This Case Special

Compared to simpler tutorials (e.g., lid-driven cavity, backward-facing
step), this case introduces:

1. **Temperature as an active scalar**: The energy equation is solved,
   and temperature affects fluid density (through the equation of state)
2. **Buoyancy**: `buoyantSimpleFoam` couples the energy and momentum
   equations through the Boussinesq-like body force term
3. **Multi-block mesh**: The T-shape requires 4 separate hex blocks
   that share common faces — a key skill for complex geometries
4. **Turbulent heat transfer**: The `alphat` field and turbulent Prandtl
   number model provide turbulent thermal diffusivity

---

## Geometry and Mesh

### Domain Geometry

```
  y
  ↑
  |
  0.2  ────── v16────────v17 ──────────────────────
              │  branch   │
              │  inlet    │
              │           │
              │  Block 3  │  branch pipe
              │  (branch) │  width:  0.05m (x: 0.125 → 0.175)
              │           │  height: 0.1m  (y: 0.1 → 0.2)
  0.1  v8────v9──────────v10────────────────────v11
       │      │           │                      │
       │  Blk │   Blk 1   │      Block 2        │
       │  0   │ (junction) │   (right section)   │
       │(left)│           │                      │
       │      │           │                      │
  0.0  v0────v1──────────v2────────────────────v3──→ x
       |      |           |                      |
       0    0.125       0.175                   0.5
```

### Vertex Numbering

The mesh uses 20 vertices (10 on the front face z=0, 10 on the back face
z=0.01):

| Vertex | x | y | z | Description |
|---|---|---|---|---|
| 0 | 0.0 | 0.0 | 0.0 | Main pipe, bottom-left, front |
| 1 | 0.125 | 0.0 | 0.0 | Junction left edge, bottom, front |
| 2 | 0.175 | 0.0 | 0.0 | Junction right edge, bottom, front |
| 3 | 0.5 | 0.0 | 0.0 | Main pipe, bottom-right, front |
| 4–7 | — | — | 0.01 | Back-face mirrors of vertices 0–3 |
| 8 | 0.0 | 0.1 | 0.0 | Main pipe, top-left, front |
| 9 | 0.125 | 0.1 | 0.0 | Junction top-left / branch bottom-left |
| 10 | 0.175 | 0.1 | 0.0 | Junction top-right / branch bottom-right |
| 11 | 0.5 | 0.1 | 0.0 | Main pipe, top-right, front |
| 12–15 | — | — | 0.01 | Back-face mirrors of vertices 8–11 |
| 16 | 0.125 | 0.2 | 0.0 | Branch pipe, top-left, front |
| 17 | 0.175 | 0.2 | 0.0 | Branch pipe, top-right, front |
| 18–19 | — | — | 0.01 | Back-face mirrors of vertices 16–17 |

### Block Decomposition

```
  ┌────────────────────────────────────────────────────────┐
  │  Block Layout (4 hex blocks forming a T-shape)         │
  │                                                        │
  │           ┌───────────┐                                │
  │           │  Block 3  │  20 × 40 × 1 cells            │
  │           │  (branch) │  branch pipe                   │
  │  ┌────────┼───────────┼────────────────────────┐       │
  │  │Block 0 │  Block 1  │       Block 2          │       │
  │  │ (left) │(junction) │      (right)           │       │
  │  │50×40×1 │ 20×40×1   │    130×40×1            │       │
  │  └────────┼───────────┼────────────────────────┘       │
  │           └───────────┘                                │
  └────────────────────────────────────────────────────────┘
```

| Block | Region | Cells (x $\times$ y $\times$ z) | Dimensions (m) |
|---|---|---|---|
| 0 | Main pipe, left of junction | 50 $\times$ 40 $\times$ 1 | 0.125 $\times$ 0.1 $\times$ 0.01 |
| 1 | Main pipe, junction zone | 20 $\times$ 40 $\times$ 1 | 0.05 $\times$ 0.1 $\times$ 0.01 |
| 2 | Main pipe, right of junction | 130 $\times$ 40 $\times$ 1 | 0.325 $\times$ 0.1 $\times$ 0.01 |
| 3 | Branch pipe | 20 $\times$ 40 $\times$ 1 | 0.05 $\times$ 0.1 $\times$ 0.01 |

**Total cells:** 50$\times$40 + 20$\times$40 + 130$\times$40 + 20$\times$40 = 2000 + 800 + 5200 + 800
= **8,800 cells**

### Multi-Block Connectivity

The key to a multi-block mesh is **shared vertices** at block interfaces.
OpenFOAM automatically merges faces that share the same vertex coordinates:

```
  Block 0 right face  ←→  Block 1 left face   (vertices 1,9,5,13)
  Block 1 right face  ←→  Block 2 left face   (vertices 2,10,6,14)
  Block 1 top face    ←→  Block 3 bottom face (vertices 9,10,13,14)
```

> **Note:** The `frontAndBack` faces are set to `empty` type, making this
> effectively a 2D simulation in the x-y plane.

---

## Boundary Conditions

### Patch Summary

| Patch | Type | Location | Description |
|---|---|---|---|
| `mainInlet` | `patch` | x = 0 face | Hot water inlet (350K) |
| `branchInlet` | `patch` | y = 0.2 face (branch top) | Cold water inlet (300K) |
| `outlet` | `patch` | x = 0.5 face | Pipe exit |
| `walls` | `wall` | All solid surfaces | No-slip, adiabatic walls |
| `frontAndBack` | `empty` | z = 0 and z = 0.01 faces | 2D constraint |

### Velocity (`U`)

| Patch | Type | Value | Notes |
|---|---|---|---|
| `mainInlet` | `fixedValue` | `(1 0 0)` | 1 m/s in x-direction |
| `branchInlet` | `fixedValue` | `(0 -0.5 0)` | 0.5 m/s downward (-y) |
| `outlet` | `zeroGradient` | — | Fully developed outflow |
| `walls` | `noSlip` | `(0 0 0)` | No-slip condition |
| `frontAndBack` | `empty` | — | 2D constraint |

### Pressure (`p_rgh`) — Buoyant Pressure

| Patch | Type | Value | Notes |
|---|---|---|---|
| `mainInlet` | `zeroGradient` | — | Pressure adjusts to match flow |
| `branchInlet` | `zeroGradient` | — | Pressure adjusts to match flow |
| `outlet` | `fixedValue` | `0` | Reference pressure at outlet |
| `walls` | `fixedFluxPressure` | `0` | Consistent with body forces |
| `frontAndBack` | `empty` | — | 2D constraint |

> **Why `fixedFluxPressure`?** This BC adjusts the pressure gradient at
> walls to be consistent with the specified velocity (no-slip) and body
> force (gravity). It is required when buoyancy is active.

### Pressure (`p`) — Total Pressure

| Patch | Type | Value | Notes |
|---|---|---|---|
| All patches | `calculated` | `0` | Derived from `p_rgh` + $\rho$gh |
| `frontAndBack` | `empty` | — | 2D constraint |

### Temperature (`T`) — ★ The Key Field

| Patch | Type | Value | Notes |
|---|---|---|---|
| `mainInlet` | `fixedValue` | `350` K | Hot stream |
| `branchInlet` | `fixedValue` | `300` K | Cold stream |
| `outlet` | `zeroGradient` | — | Temperature exits freely |
| `walls` | `zeroGradient` | — | Adiabatic (no heat flux) |
| `frontAndBack` | `empty` | — | 2D constraint |

> **Why adiabatic walls?** We assume the pipe is perfectly insulated, so
> all temperature changes come from mixing, not heat loss through walls.
> This isolates the mixing physics for study.

### Turbulent Kinetic Energy (`k`)

| Patch | Type | Value | Notes |
|---|---|---|---|
| `mainInlet` | `fixedValue` | `0.06` m$^2$/s$^2$ | TI $\approx$ 3.5% at 1.0 m/s |
| `branchInlet` | `fixedValue` | `0.015` m$^2$/s$^2$ | TI $\approx$ 3.5% at 0.5 m/s |
| `outlet` | `zeroGradient` | — | Free outflow |
| `walls` | `kqRWallFunction` | `0.06` | Wall function for k |
| `frontAndBack` | `empty` | — | 2D constraint |

**Turbulence intensity calculation:**

```
k = 1.5 × (U × TI)²

Main:   k = 1.5 × (1.0 × 0.05)²  ≈ 0.00375  (using TI=5%)
        k = 1.5 × (1.0 × 0.2)²   ≈ 0.06     (using TI=20% — as specified)
Branch: k = 1.5 × (0.5 × 0.2)²   ≈ 0.015    (using TI=20%)
```

> **Note:** Higher turbulence intensity (20%) is used here to represent
> developed pipe flow with upstream disturbances.

### Turbulent Dissipation Rate (`epsilon`)

| Patch | Type | Value | Notes |
|---|---|---|---|
| `mainInlet` | `fixedValue` | `0.0495` m$^2$/s$^3$ | Estimated from k and L |
| `branchInlet` | `fixedValue` | `0.00877` m$^2$/s$^3$ | Estimated from k and L |
| `outlet` | `zeroGradient` | — | Free outflow |
| `walls` | `epsilonWallFunction` | `0.0495` | Wall function for $\varepsilon$ |
| `frontAndBack` | `empty` | — | 2D constraint |

**Dissipation rate estimation:**

```
ε = Cμ^0.75 × k^1.5 / l

where Cμ = 0.09
      l  = 0.07 × D_h  (mixing length ≈ 7% of hydraulic diameter)

Main:   l = 0.07 × 0.1 = 0.007 m
        ε = 0.09^0.75 × 0.06^1.5 / 0.007 ≈ 0.0495 m²/s³

Branch: l = 0.07 × 0.05 = 0.0035 m
        ε = 0.09^0.75 × 0.015^1.5 / 0.0035 ≈ 0.00877 m²/s³
```

### Turbulent Viscosity (`nut`)

| Patch | Type | Value | Notes |
|---|---|---|---|
| `mainInlet` | `calculated` | `0` | Computed from k and $\varepsilon$ |
| `branchInlet` | `calculated` | `0` | Computed from k and $\varepsilon$ |
| `outlet` | `calculated` | `0` | Computed from k and $\varepsilon$ |
| `walls` | `nutkWallFunction` | `0` | Standard wall function for $\nu$_t |
| `frontAndBack` | `empty` | — | 2D constraint |

### Turbulent Thermal Diffusivity (`alphat`)

| Patch | Type | Value | Notes |
|---|---|---|---|
| `mainInlet` | `calculated` | `0` | Computed from $\nu$_t and Pr_t |
| `branchInlet` | `calculated` | `0` | Computed from $\nu$_t and Pr_t |
| `outlet` | `calculated` | `0` | Computed from $\nu$_t and Pr_t |
| `walls` | `compressible::alphatWallFunction` | `0` | Wall function for $\alpha$_t |
| `frontAndBack` | `empty` | — | 2D constraint |

> **Why `compressible::alphatWallFunction`?** Since `buoyantSimpleFoam`
> uses a compressible formulation, the turbulent thermal diffusivity BC
> must also be from the compressible namespace.

---

## Fluid and Thermophysical Properties

### Thermophysical Model

The `thermophysicalProperties` dictionary configures the thermodynamic
model used by `buoyantSimpleFoam`:

| Setting | Value | Description |
|---|---|---|
| `type` | `heRhoThermo` | Enthalpy-based, density-calculating thermo |
| `mixture` | `pureMixture` | Single-component fluid |
| `transport` | `const` | Constant transport properties |
| `thermo` | `hConst` | Constant specific heat capacity |
| `equationOfState` | `perfectGas` | Ideal gas law ($\rho$ = p/(RT)) |
| `specie` | `specie` | Basic species model |
| `energy` | `sensibleEnthalpy` | Energy variable is sensible enthalpy |

> **Note on equation of state:** We use `perfectGas` here for simplicity,
> which means density varies with temperature and pressure. For water at
> these conditions, `incompressiblePerfectGas` or a `Boussinesq` model
> would be more physically accurate, but `perfectGas` is sufficient for
> demonstrating thermal mixing physics.

### Mixture Properties

| Property | Symbol | Value | Units |
|---|---|---|---|
| Molecular weight | M | 18.0 | g/mol |
| Specific heat capacity | Cp | 4182 | J/(kg·K) |
| Heat of formation | Hf | 0 | J/kg |
| Dynamic viscosity | $\mu$ | 1$\times$10$^{-3}$ | Pa·s |
| Prandtl number | Pr | 7.0 | — |

### Transport Properties

| Property | Symbol | Value | Dimensions |
|---|---|---|---|
| Transport model | — | Newtonian | — |
| Kinematic viscosity | $\nu$ | 1$\times$10$^{-6}$ | [0 2 -1 0 0 0 0] m$^2$/s |

### Gravity

| Component | Value | Units |
|---|---|---|
| g_x | 0 | m/s$^2$ |
| g_y | -9.81 | m/s$^2$ |
| g_z | 0 | m/s$^2$ |

---

## Turbulence Setup

### Model Selection

| Setting | Value |
|---|---|
| Simulation type | RAS |
| RAS model | kEpsilon |
| Turbulence | on |
| Print coefficients | on |

### k-$\varepsilon$ Model Constants (Defaults)

| Constant | Value | Description |
|---|---|---|
| C$\mu$ | 0.09 | Eddy viscosity coefficient |
| C1 | 1.44 | Production coefficient |
| C2 | 1.92 | Destruction coefficient |
| $\sigma$_k | 1.0 | Turbulent Prandtl number for k |
| $\sigma$_$\varepsilon$ | 1.3 | Turbulent Prandtl number for $\varepsilon$ |

### Wall Treatment

The k-$\varepsilon$ model uses **standard wall functions** for near-wall treatment:

- `kqRWallFunction` — applies the equilibrium assumption for k at walls
- `epsilonWallFunction` — computes $\varepsilon$ from the equilibrium log-law
- `nutkWallFunction` — computes $\nu$_t from k and the wall distance
- `compressible::alphatWallFunction` — wall function for $\alpha$_t

> **Important:** Wall functions require the first cell centre to be in
> the log-law region (30 < y$^+$ < 300). With this mesh resolution, y$^+$
> values may be approximate. For production work, verify y$^+$ values
> and consider using `kOmegaSST` with wall-resolved meshes.

---

## Numerical Schemes and Solver Settings

### controlDict — Run Control

| Parameter | Value | Description |
|---|---|---|
| `application` | `buoyantSimpleFoam` | Steady-state buoyant solver |
| `startFrom` | `startTime` | Start from t=0 |
| `endTime` | `2000` | Maximum 2000 iterations |
| `deltaT` | `1` | Pseudo-time step = 1 |
| `writeControl` | `timeStep` | Write at fixed iteration intervals |
| `writeInterval` | `200` | Write every 200 iterations |
| `purgeWrite` | `3` | Keep only last 3 time directories |

### fvSchemes — Discretization

| Category | Scheme | Notes |
|---|---|---|
| `ddtSchemes` | `steadyState` | No time derivative (SIMPLE) |
| `gradSchemes` | `Gauss linear` | Second-order gradient |
| `div(phi,U)` | `bounded Gauss linearUpwind` | Blended upwind for velocity |
| `div(phi,h)` | `bounded Gauss linearUpwind` | Blended upwind for enthalpy |
| `div(phi,k)` | `bounded Gauss upwind` | First-order for stability |
| `div(phi,epsilon)` | `bounded Gauss upwind` | First-order for stability |
| `laplacianSchemes` | `Gauss linear corrected` | Second-order Laplacian |
| `interpolationSchemes` | `linear` | Linear interpolation |
| `snGradSchemes` | `corrected` | Non-orthogonal correction |

> **Why `bounded` schemes?** The `bounded` keyword ensures that the
> divergence schemes account for the continuity error, preventing
> unbounded behaviour in steady-state simulations.

### fvSolution — Linear Solvers and SIMPLE

**Linear solvers:**

| Field | Solver | Preconditioner | Tolerance | Relative Tolerance |
|---|---|---|---|---|
| `p_rgh` | `GAMG` | `GaussSeidel` | 1$\times$10$^{-7}$ | 0.01 |
| `U`, `h`, `k`, `epsilon` | `PBiCGStab` | `DILU` | 1$\times$10$^{-7}$ | 0.1 |

**SIMPLE algorithm settings:**

| Parameter | Value |
|---|---|
| Non-orthogonal correctors | 0 |
| Reference pressure cell | 0 |
| Reference pressure value | 0 |

**Under-relaxation factors:**

| Field/Equation | Factor | Notes |
|---|---|---|
| `p_rgh` | 0.3 | Conservative for pressure stability |
| `U` | 0.7 | Standard for SIMPLE |
| `h` | 0.7 | Enthalpy (temperature) |
| `k` | 0.7 | Turbulent kinetic energy |
| `epsilon` | 0.7 | Dissipation rate |

**Residual control:**

All fields converge to 1$\times$10$^{-4}$. The solver stops when all residuals drop
below these thresholds, even before reaching `endTime`.

---

## How to Run

### Prerequisites

- OpenFOAM 6+ or OpenFOAM v2012+ installed and sourced
- `buoyantSimpleFoam` available in your path

### Step-by-Step

**1. Navigate to the case directory:**

```bash
cd projects/09_mixing_t_junction
```

**2. Run the complete simulation:**

```bash
./Allrun
```

This executes:
1. `blockMesh` — generates the T-junction mesh
2. `checkMesh` — validates mesh quality
3. `buoyantSimpleFoam` — runs the steady-state simulation

**3. Or run each step manually:**

```bash
# Generate the mesh
blockMesh

# Verify mesh quality
checkMesh

# Run the solver
buoyantSimpleFoam
```

**4. Monitor convergence:**

```bash
# Watch residuals in real-time
tail -f log.buoyantSimpleFoam | grep "Solving for"

# Or use foamMonitor (if available)
foamMonitor -l postProcessing/residuals/0/residuals.dat
```

**5. Post-process:**

```bash
# Open in ParaView
paraFoam

# Or convert to VTK format
foamToVTK
```

**6. Clean up and start fresh:**

```bash
./Allclean
```

### Expected Runtime

On a modern single-core machine, expect:
- `blockMesh`: < 1 second
- `buoyantSimpleFoam`: 1–5 minutes (depending on convergence)
- Total: ~2–5 minutes

---

## Expected Results

### Convergence

The simulation should converge within 500–1500 iterations. Residuals
should drop monotonically for all fields. If residuals oscillate or
diverge, reduce the under-relaxation factors (especially `p_rgh`).

```
  Typical residual history:

  Residual
  1e+0  ─┐
         │╲
  1e-1  ─│ ╲  ← U, h
         │  ╲╲
  1e-2  ─│   ╲╲╲  ← k, epsilon
         │    ╲╲╲╲
  1e-3  ─│     ╲╲╲╲╲
         │      ╲╲╲╲╲  ← p_rgh
  1e-4  ─│───────╲╲╲╲╲──── convergence threshold
         │        ~~~~
  1e-5  ─┘
         0    500   1000   1500   2000
                  Iteration
```

### Temperature Field

```
  Temperature distribution at convergence:

  Branch (300K)
       ↓
  ┌────┴────┐
  │░░░░░░░░░│  ← 300K (cold)
  │░░░░░░░░░│
  ┼─────────┼────────────────────────────────┐
  │▓▓▓▓░░░░░│░░▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒│
  │▓▓▓▓▓▒▒░░│▒▒▒▒▒▒▓▓▓▓▓▓▓▓▓▓▓▓▒▒▒▒▒▒▒▒▒▒▒│ → outlet
  │▓▓▓▓▓▓▓▓▒│▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓│
  └─────────┴────────────────────────────────┘
  mainInlet                              outlet
  (350K)

  ░ = 300–315K   ▒ = 315–335K   ▓ = 335–350K
```

**Key observations:**

1. **At the junction:** Sharp temperature interface where hot and cold
   streams meet. The cold stream is deflected by the main flow.

2. **Downstream mixing zone:** Temperature gradients gradually smooth
   out due to turbulent diffusion. The cold stream initially hugs the
   top wall (carried by its downward momentum and main flow).

3. **Far downstream:** Temperature approaches the theoretical mixed
   value of approximately **340K** (mass-flow-weighted average).

4. **Thermal stratification:** Near the junction, a clear temperature
   gradient exists from bottom (hot) to top (cold). This stratification
   persists for several pipe diameters downstream.

### Velocity Field

- The main flow accelerates slightly at the junction as the branch flow
  adds mass
- A small recirculation zone may form just downstream of the branch
  entry on the upper wall
- The combined flow develops toward a parabolic-like profile downstream

### Bulk Quantities

| Quantity | Expected Value |
|---|---|
| Outlet bulk temperature | ~340K |
| Mass-averaged velocity at outlet | ~1.125 m/s |
| Maximum velocity at junction | ~1.3–1.5 m/s |
| Temperature range at outlet | 330–345K (not fully mixed) |

---

## Post-Processing Guide

### Temperature Contours (ParaView)

1. Open the case: `paraFoam`
2. Select the `T` field
3. Apply a colour map (e.g., "Cool to Warm" with range 300–350K)
4. This reveals the mixing pattern and thermal stratification

### Velocity Vectors at Junction

1. In ParaView, apply the **Glyph** filter to `U`
2. Set glyph type to "Arrow"
3. Scale by magnitude
4. Focus on the junction region to see flow deflection

### Temperature Profile at Outlet

Extract the temperature distribution across the outlet face to quantify
mixing quality:

```bash
# Using postProcess utility
postProcess -func "sampleDict" -latestTime
```

Create a `system/sampleDict` file to sample along lines:

```
type            sets;
libs            ("libsampling.so");
writeControl    writeTime;
sets
(
    outletLine
    {
        type    uniform;
        axis    y;
        start   (0.499 0.0 0.005);
        end     (0.499 0.1 0.005);
        nPoints 100;
    }
    junctionLine
    {
        type    uniform;
        axis    y;
        start   (0.2 0.0 0.005);
        end     (0.2 0.1 0.005);
        nPoints 100;
    }
);
fields (T U);
```

### Mixing Quality Assessment

Calculate the mixing coefficient at different downstream positions:

```
η_mix(x) = 1 - σ_T(x) / (T_hot - T_cold) × 2

where σ_T(x) is the standard deviation of temperature at position x

η_mix → 0:  completely unmixed (separate hot/cold streams)
η_mix → 1:  perfectly mixed (uniform temperature)
```

Typical mixing evolution:

```
  η_mix
  1.0  ─┐                              ●───── fully mixed
        │                          ●●
  0.8  ─│                      ●●
        │                   ●●
  0.6  ─│                ●●
        │             ●●
  0.4  ─│          ●●
        │       ●●
  0.2  ─│    ●●
        │ ●●
  0.0  ─●─────────────────────────────────
        junction  5D   10D   15D   20D
              Distance downstream
```

### Residual Monitoring

```bash
# Plot residuals using gnuplot
foamLog log.buoyantSimpleFoam
cd logs/
gnuplot
> set logscale y
> plot "p_rgh_0" w l, "Ux_0" w l, "h_0" w l, "k_0" w l, "epsilon_0" w l
```

---

## Exercises and Experiments

### Exercise 1: Vary the Flow Ratio

Change the branch inlet velocity to study different mixing scenarios:

| Case | U_branch (m/s) | Flow ratio R | Expected T_out (K) |
|---|---|---|---|
| Baseline | 0.5 | 0.25 | 340 |
| Low branch | 0.2 | 0.10 | 345 |
| Equal flows | 1.0 | 0.50 | 333 |
| High branch | 2.0 | 1.00 | 325 |

Modify `0/U`:
```
branchInlet
{
    type    fixedValue;
    value   uniform (0 -1.0 0);   // change velocity here
}
```

**Question:** How does the flow ratio affect mixing length and thermal
stratification?

### Exercise 2: Change Temperature Difference

Keep velocities fixed but vary the temperature difference:

| Case | T_main (K) | T_branch (K) | $\Delta T$ (K) |
|---|---|---|---|
| Small $\Delta T$ | 320 | 300 | 20 |
| Baseline | 350 | 300 | 50 |
| Large $\Delta T$ | 400 | 300 | 100 |

**Question:** Does the buoyancy effect become more significant with
larger temperature differences? Look for changes in the flow pattern.

### Exercise 3: Switch to k-$\omega$ SST Model

Replace the k-$\varepsilon$ model with k-$\omega$ SST for better near-wall behaviour:

1. Modify `constant/turbulenceProperties`:
   ```
   RASModel    kOmegaSST;
   ```

2. Replace `0/epsilon` with `0/omega`:
   ```
   omega = epsilon / (Cmu × k)
   ```

3. Update wall functions accordingly

**Question:** Does the turbulence model significantly affect the
predicted temperature field?

### Exercise 4: Add Heat Loss Through Walls

Change the wall thermal BC from adiabatic to a fixed heat flux or
convective condition:

```
walls
{
    type    fixedGradient;
    gradient uniform 1000;    // heat flux in W/m²
}
```

**Question:** How does wall heat loss affect the downstream temperature
profile?

### Exercise 5: Transient Simulation (LES for Thermal Striping)

For advanced users, convert to a transient LES simulation to capture
thermal striping (high-frequency temperature fluctuations at the wall):

1. Change solver to `buoyantPimpleFoam`
2. Switch turbulence to LES with Smagorinsky or WALE model
3. Refine the mesh significantly (wall-resolved, y$^+$ < 1)
4. Set appropriate time step (CFL < 0.5)
5. Run for sufficient flow-through times

This is directly relevant to nuclear thermal fatigue assessment where
the frequency and amplitude of temperature fluctuations determine
fatigue life.

### Exercise 6: 3D Simulation

Extend the 2D case to 3D by:

1. Changing `frontAndBack` from `empty` to `symmetryPlane` or extending
   the mesh in the z-direction
2. Using circular pipe cross-sections with `snappyHexMesh`
3. Comparing 2D vs 3D mixing patterns

---

## References

1. **OpenFOAM User Guide** — buoyantSimpleFoam solver documentation
   https://www.openfoam.com/documentation/guides/latest/doc/

2. **Thermal mixing in T-junctions** — Westin, J. et al. (2008),
   "Experiments and Unsteady CFD-Calculations of Thermal Mixing in a
   T-Junction," OECD/NEA/CSNI Workshop on Benchmarking of CFD Codes
   for Application to Nuclear Reactor Safety (CFD4NRS).

3. **Nuclear thermal fatigue** — Braillard, O., Edelin, D. (2000),
   "Advanced experimental tools designed for the assessment of the
   thermal load applied to the mixing tee and nozzle geometries,"
   ASME PVP Conference.

4. **Civaux-1 incident** — Chapuliot, S. et al. (2005), "Hydro-thermal-
   mechanical analysis of thermal fatigue in a mixing tee," Nuclear
   Engineering and Design, Vol. 235, pp. 575–596.

5. **k-$\varepsilon$ turbulence model** — Launder, B.E. and Spalding, D.B. (1974),
   "The numerical computation of turbulent flows," Computer Methods in
   Applied Mechanics and Engineering, Vol. 3, pp. 269–289.

6. **Turbulent mixing** — Dimotakis, P.E. (2005), "Turbulent Mixing,"
   Annual Review of Fluid Mechanics, Vol. 37, pp. 329–356.

7. **OpenFOAM thermophysical models** — Greenshields, C.J. (2021),
   OpenFOAM User Guide, Chapter 7: Thermophysical Models.

8. **Thermal striping in nuclear systems** — IAEA-TECDOC-1627 (2009),
   "Pressurized Thermal Shock in Nuclear Power Plants: Good Practices
   for Assessment."

---

## Related Notes in This Repository

| Note | Topic | Relevance |
|---|---|---|
| [01 — Introduction to CFD](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md) | Fundamental CFD concepts | Governing equations overview |
| [02 — OpenFOAM Cases](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md) | Case directory structure | Understanding 0/, constant/, system/ |
| [03 — OpenFOAM Dictionaries](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md) | Dictionary file format | FoamFile headers and syntax |
| [04 — Meshing](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/04_meshing.md) | Meshing theory | Multi-block mesh concepts |
| [05 — Boundary Conditions](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md) | BC types and usage | Wall functions, inlet/outlet BCs |
| [06 — Schemes and Solvers](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_schemes_and_solvers.md) | Numerical methods | SIMPLE algorithm, discretization |
| [08 — CFL Number](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md) | Time stepping | Relevant for transient extensions |

> **Further reading:**
> - [Meshing Notes](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/04_meshing.md) for multi-block mesh theory
> - [Boundary Conditions](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md) for wall functions and thermal BCs
> - [Schemes and Solvers](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_schemes_and_solvers.md) for SIMPLE algorithm details

---

*This tutorial is part of the [OpenFOAM Tutorials](https://github.com/djeada/OpenFoam-Tutorials/blob/main/README.md) collection.*
