# Viscosity and Viscosity Models in OpenFOAM

[All notes](../README.md#-notes-overview) · [Previous: Note 10](10_icofoam_solver_analysis.md) · [Next: Note 12](12_multiphase_flows.md)

> **Version scope:** Classic OpenFOAM Foundation 6, matching the project launchers.
> Other releases/forks may need different dictionaries or executables; see
> [version guidance](../README.md#version-and-example-scope) and the
> [current reference cases](02_openfoam_cases.md#current-reference-cases).

<details>
<summary>Contents</summary>

  - [Dynamic vs Kinematic Viscosity](#dynamic-vs-kinematic-viscosity)
  - [Why OpenFOAM Uses Kinematic Viscosity](#why-openfoam-uses-kinematic-viscosity)
- [Newtonian vs Non-Newtonian Fluids](#newtonian-vs-non-newtonian-fluids)
  - [Newtonian Fluids](#newtonian-fluids)
  - [Non-Newtonian Fluids](#non-newtonian-fluids)
  - [Shear Stress vs Shear Rate — Visual Comparison](#shear-stress-vs-shear-rate--visual-comparison)
  - [Bingham (Yield-Stress) Fluids](#bingham-yield-stress-fluids)
- [Viscosity Models Available in OpenFOAM](#viscosity-models-available-in-openfoam)
  - [Newtonian](#newtonian)
  - [Solver Support and Coefficient Units](#solver-support-and-coefficient-units)
  - [3.2 Power Law](#32-power-law)
  - [Cross Power Law](#cross-power-law)
  - [Bird-Carreau](#bird-carreau)
  - [Herschel-Bulkley](#herschel-bulkley)
  - [Model Comparison Summary](#model-comparison-summary)
- [Configuring Viscosity in OpenFOAM](#configuring-viscosity-in-openfoam)
  - [The `transportProperties` Dictionary](#the-transportproperties-dictionary)
  - [Structure of the Dictionary](#structure-of-the-dictionary)
  - [Complete Working Examples](#complete-working-examples)
- [Effective Viscosity in Turbulent Flows](#effective-viscosity-in-turbulent-flows)
  - [How It Works in OpenFOAM](#how-it-works-in-openfoam)
  - [Turbulent Viscosity Ratio](#turbulent-viscosity-ratio)
- [Temperature-Dependent Viscosity](#temperature-dependent-viscosity)
  - [When Temperature Matters](#when-temperature-matters)
- [Common Pitfalls & Troubleshooting](#common-pitfalls--troubleshooting)
  - [Pitfall #1: Dynamic vs Kinematic Viscosity Mix-Up](#pitfall-1-dynamic-vs-kinematic-viscosity-mix-up)
  - [Pitfall #2: Wrong Dimensions in the Dictionary](#pitfall-2-wrong-dimensions-in-the-dictionary)
  - [Pitfall #3: Non-Newtonian Models Diverging](#pitfall-3-non-newtonian-models-diverging)
  - [Pitfall #4: Missing `transportProperties` Keywords](#pitfall-4-missing-transportproperties-keywords)
  - [Pitfall #5: Non-Newtonian with Turbulence](#pitfall-5-non-newtonian-with-turbulence)
- [Practical Tips](#practical-tips)
  - [Checking Available Models](#checking-available-models)
  - [Inspecting the Viscosity Field](#inspecting-the-viscosity-field)
  - [Dimensional Consistency Check](#dimensional-consistency-check)
- [Quick Reference](#quick-reference)
  - [Model Selection Quick Guide](#model-selection-quick-guide)
- [Where to Go Next](#where-to-go-next)

</details>

Viscosity is a fluid's **resistance to deformation** under shear stress. Imagine two parallel plates with fluid between them — the top plate moves while the bottom plate stays fixed. The force needed to maintain that motion depends on the fluid's viscosity.

```
        Velocity u
        ──────────►
  ┌──────────────────────────────────────┐  ◄── Moving plate
  │  → → → → → → → → → → → → → → → →     │
  │  → → → → → → → → → → → → → → →       │
  │  → → → → → → → → → → → → →           │      Shear stress
  │  → → → → → → → → → → →               │      τ = μ (du/dy)
  │  → → → → → → → → →                   │
  │  → → → → → → →                       │
  │  → → → → →                           │
  │  → → →                               │
  │  →                                   │
  ├──────────────────────────────────────┤  ◄── Fixed plate
  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
```

For steady, fully developed flow with constant viscosity and no streamwise
pressure gradient, the velocity profile is linear for a Newtonian fluid — this is **Couette flow**, the simplest example of viscous shear.

### Dynamic vs Kinematic Viscosity

| Symbol | Name                 | Definition | SI Units     | Physical Meaning                         |
|--------|----------------------|------------|--------------|------------------------------------------|
| $\mu$      | Dynamic viscosity    | —          | $\mathrm{Pa\,s}$ (kg/m·s)| Shear stress per unit shear rate         |
| $\nu$      | Kinematic viscosity  | $\nu = \mu / \rho$  | $\mathrm{m}^2/\mathrm{s}$         | Momentum diffusivity                     |

**Dynamic viscosity ($\mu$)** relates shear stress to shear rate directly:

$$
\tau = \mu \cdot \frac{du}{dy}
$$

where:

- $\tau$ = shear stress [Pa]
- $\mu$ = dynamic viscosity [$\mathrm{Pa\,s}$]
- $du/dy$ = velocity gradient [1/s]

**Kinematic viscosity ($\nu$)** is dynamic viscosity divided by density. It appears naturally in the incompressible Navier-Stokes equations:

$$
\frac{\partial u}{\partial t} + (u \cdot \nabla)u = -\frac{1}{\rho}\nabla p + \nu \nabla^{2}u
$$

Here $\nu = \mu/\rho$ controls momentum diffusion. This simplified Laplacian
form assumes constant viscosity and incompressibility. Variable-viscosity
generalized Newtonian flow requires the divergence of the viscous stress tensor,
not simply replacing a constant in $\nu\nabla^2 u$.

### Why OpenFOAM Uses Kinematic Viscosity

Most OpenFOAM incompressible solvers (`icoFoam`, `simpleFoam`, `pimpleFoam`) divide the entire momentum equation by density $\rho$. This means:

- Pressure becomes $p/\rho$ (kinematic pressure, units $\mathrm{m}^2/\mathrm{s}^2$)
- Viscosity becomes $\nu = \mu/\rho$ (kinematic viscosity, units $\mathrm{m}^2/\mathrm{s}$)

> **⚠️ Warning:** When you see `nu` in OpenFOAM dictionaries, it is **always**
> kinematic viscosity $\nu$ ($\mathrm{m}^2/\mathrm{s}$), not dynamic viscosity $\mu$ ($\mathrm{Pa\,s}$).
> If a data sheet gives you $\mu$, divide by $\rho$ first!

**Typical values:**

| Fluid          | $\mu$ ($\mathrm{Pa\,s}$)          | $\rho$ ($\mathrm{kg/m}^3$)  | $\nu$ ($\mathrm{m}^2/\mathrm{s}$)           |
|----------------|--------------------|------------|---------------------|
| Air (20 °C)    | $1.825 \times 10^{-5}$      | 1.205      | $1.516 \times 10^{-5}$       |
| Water (20 °C)  | $1.002 \times 10^{-3}$      | 998.2      | $1.004 \times 10^{-6}$       |
| Engine oil     | 0.1 – 0.3         | 880        | $1.1 \times 10^{-4}$ – $3.4 \times 10^{-4}$ |
| Honey          | 2 – 10            | 1420       | $1.4 \times 10^{-3}$ – $7 \times 10^{-3}$   |
| Blood (37 °C)  | $3$–$4 \times 10^{-3}$      | 1060       | $2.8$–$3.8 \times 10^{-6}$   |

## Newtonian vs Non-Newtonian Fluids

The distinction comes down to one question: **does viscosity change with shear rate?**

### Newtonian Fluids

For a Newtonian fluid, viscosity is constant regardless of how fast the fluid is being sheared:

$$
\tau = \mu \cdot \dot{\gamma} \quad \text{(linear relationship)}
$$

where $\dot{\gamma} = du/dy$ = shear rate [1/s]

Most common engineering fluids are Newtonian: **water, air, most gases, light oils, and simple organic solvents**.

### Non-Newtonian Fluids

For non-Newtonian fluids, the **apparent viscosity** $\eta$ changes with shear rate $\dot{\gamma}$. There are several categories:

| Type                | Behavior                          | Examples                          |
|---------------------|-----------------------------------|-----------------------------------|
| Shear-thinning      | $\eta$ decreases as $\dot{\gamma}$ increases       | Blood, paint, polymer melts, ketchup |
| Shear-thickening    | $\eta$ increases as $\dot{\gamma}$ increases       | Cornstarch in water, wet sand     |
| Bingham / Yield     | No flow below critical stress $\tau_0$  | Toothpaste, concrete, drilling mud|
| Thixotropic         | $\eta$ decreases over time at constant $\dot{\gamma}$ | Yogurt, some gels              |
| Rheopectic          | $\eta$ increases over time at constant $\dot{\gamma}$  | Gypsum paste                   |

> These examples cover time-independent generalized Newtonian models.
> Thixotropic/rheopectic fluids need a model that evolves material structure or
> history; check the available models in your release before choosing or extending one.

### Shear Stress vs Shear Rate — Visual Comparison

```
Shear
Stress
τ
▲
│                                              ╱  Shear-thickening
│                                           ╱╱   (dilatant)
│                                        ╱╱
│                                     ╱╱
│                                  ╱╱
│                               ╱╱
│                            ╱╱
│                         ╱╱                      Newtonian
│                      ╱╱                     ╱
│                   ╱╱                    ╱╱
│                ╱╱                   ╱╱
│             ╱╱                  ╱╱
│          ╱╱                 ╱╱
│       ╱╱                ╱╱
│    ╱╱               ╱╱─────────────────  Shear-thinning
│ ╱╱            ╱╱╱╱                        (pseudoplastic)
│╱         ╱╱╱╱
│     ╱╱╱╱
┼────────────────────────────────────────────────► γ̇  (Shear Rate)
```

```
  Apparent
  Viscosity
  η ▲
    │
    │ ╲                                        Shear-thinning
    │  ╲╲
    │    ╲╲╲
    │       ╲╲╲╲
    │           ╲╲╲╲╲╲╲╲────────────────────
    │
    │ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─  Newtonian (constant)
    │
    │                          ╱╱╱╱╱╱╱╱╱╱╱╱╱
    │                   ╱╱╱╱╱╱
    │            ╱╱╱╱╱╱╱                       Shear-thickening
    │     ╱╱╱╱╱╱
    │╱╱╱╱╱
    ┼──────────────────────────────────────────► γ̇  (Shear Rate)
```

### Bingham (Yield-Stress) Fluids

These fluids behave as a **solid** below a critical shear stress $\tau_0$
and flow as a fluid above it:

```
  τ ▲
    │                                    ╱╱  Herschel-Bulkley
    │                                 ╱╱     (yield + shear-thinning)
    │                              ╱╱
    │                           ╱╱
    │                        ╱╱    ╱╱╱╱╱  Bingham plastic
    │                     ╱╱  ╱╱╱╱╱       (yield + linear)
    │                  ╱╱╱╱╱╱
    │               ╱╱╱╱╱
    │            ╱╱╱╱
    │         ╱╱╱
    │      ╱╱╱
    │   ╱╱╱
    ├───╱─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─  τ₀  (yield stress)
    │  │
    │  │  ◄── No flow below τ₀
    ┼──┼───────────────────────────────────► γ̇
```

Real-world examples: toothpaste stays on the brush ($\tau < \tau_0$) but flows when you squeeze ($\tau > \tau_0$). Concrete doesn't flow under its own weight until vibrated.

## Viscosity Models Available in OpenFOAM

OpenFOAM provides several viscosity models through the `transportModel`
keyword in `constant/transportProperties`. The solver then computes the
apparent kinematic viscosity field `nu` at each cell and time step.

### Newtonian

The simplest model — constant viscosity everywhere.

**Equation:**

$$
\nu = \text{constant}
$$

**Parameters:**

| Parameter | Symbol | Meaning             | Units |
|-----------|--------|---------------------|-------|
| `nu`      | $\nu$      | Kinematic viscosity | $\mathrm{m}^2/\mathrm{s}$  |

**Configuration:**

```c
// constant/transportProperties — Newtonian fluid

FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 1e-06;  // Water at 20°C
```

> **💡 Tip:** The `[0 2 -1 0 0 0 0]` is OpenFOAM's dimension set:
> `[kg m s K mol A cd]`. So `[0 2 -1 0 0 0 0]` means $\mathrm{m}^2/\mathrm{s}$ — exactly
> the units of kinematic viscosity. See [OpenFOAM Dictionaries](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md)
> for details on dimension sets.

### Solver Support and Coefficient Units

These dictionaries describe the classic Foundation 6 incompressible transport
model interface. Solvers such as `simpleFoam`/`pimpleFoam` construct a transport
model; stock `icoFoam` 6 reads a constant `nu` and does **not** activate a
non-Newtonian model merely because `transportModel` is present. Use a solver
that instantiates and updates the selected model. See the [v6 rheology guide](https://doc.cfd.direct/openfoam/user-guide-v6/transport-rheology).

In dimensional continuum notation, $\nu=k_{phys}\dot\gamma^{n-1}$ requires
$k_{phys}$ in $\mathrm{m}^2\mathrm{s}^{n-2}$. The v6 dictionaries instead
assign `k` viscosity dimensions and apply the power to a shear rate normalized
by a reference of 1 s⁻¹. Keep their `[0 2 -1 0 0 0 0]` dimension set; do not
replace it with the continuum coefficient’s n-dependent dimensions. `k` in
Bird–Carreau denotes a time constant and is a different parameter.

### 3.2 Power Law

Models shear-thinning or shear-thickening behavior with a simple power
function.

**Equation:**

$$
\nu(\dot{\gamma}) = k \cdot \dot{\gamma}^{n-1}
$$

Clamped: $\nu_{\min} \leq \nu \leq \nu_{\max}$

**Parameters:**

| Parameter | Symbol | Meaning                        | Units      |
|-----------|--------|--------------------------------|------------|
| `nuMax`   | $\nu_{\max}$  | Upper viscosity limit          | $\mathrm{m}^2/\mathrm{s}$       |
| `nuMin`   | $\nu_{\min}$  | Lower viscosity limit          | $\mathrm{m}^2/\mathrm{s}$       |
| `k`       | k      | Consistency index              | $\mathrm{m}^2/\mathrm{s}$ in the v6 dictionary; see coefficient convention above     |
| `n`       | n      | Power-law index                | dimensionless |

- **$n < 1$** → shear-thinning (pseudoplastic)
- **$n = 1$** → Newtonian
- **$n > 1$** → shear-thickening (dilatant)

**Configuration:**

```c
// constant/transportProperties — Power Law fluid (polymer melt)

transportModel  powerLaw;

powerLawCoeffs
{
    nuMax       [0 2 -1 0 0 0 0] 1e-02;
    nuMin       [0 2 -1 0 0 0 0] 1e-06;
    k           [0 2 -1 0 0 0 0] 0.01;
    n           [0 0  0 0 0 0 0] 0.5;    // Shear-thinning
}
```

> **⚠️ Warning:** Always set `nuMin` and `nuMax` to prevent numerical
> blow-up. Without bounds, very low or very high shear rates produce
> extreme viscosity values that crash the solver.

### Cross Power Law

The Cross model handles the transition between a zero-shear-rate plateau ($\eta_0$) and an infinite-shear-rate plateau ($\eta_\infty$) — more physically realistic than the plain power law.

**Equation:**

$$
\nu(\dot{\gamma}) = \nu_\infty + \frac{\nu_0 - \nu_\infty}{1 + (m \cdot \dot{\gamma})^n}
$$

**Parameters:**

| Parameter | Symbol | Meaning                             | Units         |
|-----------|--------|-------------------------------------|---------------|
| `nu0`     | $\nu_0$     | Zero-shear-rate kinematic viscosity | $\mathrm{m}^2/\mathrm{s}$          |
| `nuInf`   | $\nu_\infty$     | Infinite-shear-rate viscosity       | $\mathrm{m}^2/\mathrm{s}$          |
| `m`       | m      | Consistency time constant           | s             |
| `n`       | n      | Power-law index                     | dimensionless |

```
  ν ▲
    │
 ν₀ ├─────╲
    │       ╲╲
    │         ╲╲╲           Cross model gives smooth
    │            ╲╲╲        S-shaped transition on
    │               ╲╲╲╲    a log-log plot
    │                   ╲╲╲╲
 ν∞ ├─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─╲╲───────────
    │
    ┼─────────────────────────────────────► γ̇  (log scale)
```

**Configuration:**

```c
// constant/transportProperties — Cross Power Law (blood flow)

transportModel  CrossPowerLaw;

CrossPowerLawCoeffs
{
    nu0         [0 2 -1 0 0 0 0] 3.3e-06;   // Zero-shear viscosity
    nuInf       [0 2 -1 0 0 0 0] 3.0e-06;   // Infinite-shear viscosity
    m           [0 0  1 0 0 0 0] 3.736;      // Time constant
    n           [0 0  0 0 0 0 0] 0.357;      // Power index
}
```

### Bird-Carreau

Similar to the Cross model but uses a different mathematical form. Very popular for polymer solutions and biological fluids.

**Equation:**

$$
\nu(\dot{\gamma}) = \nu_\infty + (\nu_0 - \nu_\infty) \cdot \left[1 + (k \cdot \dot{\gamma})^{2}\right]^{(n-1)/2}
$$

**Parameters:**

| Parameter | Symbol | Meaning                             | Units         |
|-----------|--------|-------------------------------------|---------------|
| `nu0`     | $\nu_0$     | Zero-shear-rate kinematic viscosity | $\mathrm{m}^2/\mathrm{s}$          |
| `nuInf`   | $\nu_\infty$     | Infinite-shear-rate viscosity       | $\mathrm{m}^2/\mathrm{s}$          |
| `k`       | k      | Relaxation time                     | s             |
| `n`       | n      | Power-law index (< 1 for thinning)  | dimensionless |

**Configuration:**

```c
// constant/transportProperties — Bird-Carreau (blood flow simulation)

transportModel  BirdCarreau;

BirdCarreauCoeffs
{
    nu0         [0 2 -1 0 0 0 0] 3.3e-06;   // Zero-shear viscosity
    nuInf       [0 2 -1 0 0 0 0] 3.0e-06;   // Infinite-shear viscosity
    k           [0 0  1 0 0 0 0] 8.2;        // Relaxation time
    n           [0 0  0 0 0 0 0] 0.2;        // Power index
}
```

> **📌 Note:** Bird-Carreau and Cross models produce similar results for
> many fluids. Bird-Carreau is more common in biomedical CFD literature,
> while Cross is often preferred in polymer engineering. Choose based on
> available experimental data for your fluid.

### Herschel-Bulkley

Models **yield-stress fluids** — materials that behave as solids below a threshold stress and flow like a power-law fluid above it.

**Equation:**

$$
\tau = \tau_y + K \cdot \dot{\gamma}^n \quad (\text{yielded material; physical stress in Pa})
$$

$$
\nu(\dot{\gamma}) = \min(\nu_0, \tau_0/\dot{\gamma} + k \cdot \dot{\gamma}^{n-1})
$$

**Parameters:**

| Parameter | Symbol | Meaning                        | Units         |
|-----------|--------|--------------------------------|---------------|
| `nu0`     | $\nu_0$     | Limiting viscosity (yield reg.)| $\mathrm{m}^2/\mathrm{s}$          |
| `tau0`    | $\tau_0$     | Yield stress                   | $\mathrm{m}^2/\mathrm{s}^2$  *      |
| `k`       | k      | Consistency index              | $\mathrm{m}^2/\mathrm{s}$ in the v6 dictionary; see coefficient convention above        |
| `n`       | n      | Power-law index                | dimensionless |

> \* Because OpenFOAM works with kinematic quantities (divided by $\rho$),
> `tau0` has units of $\mathrm{m}^2/\mathrm{s}^2$ (kinematic stress), not Pa:
> `tau0 = tau_y / rho`. In the second equation, $\tau_0$ denotes this
> kinematic value, and `k` follows the normalized-shear-rate convention above.
> The finite `nu0` cap regularizes the ideal yield law: below yield it permits
> slow viscous creep rather than an exactly rigid unyielded solid.

**Configuration:**

```c
// constant/transportProperties — Herschel-Bulkley (concrete slurry)

transportModel  HerschelBulkley;

HerschelBulkleyCoeffs
{
    nu0         [0 2 -1 0 0 0 0] 1e+03;     // Large "solid" viscosity
    tau0        [0 2 -2 0 0 0 0] 1.0;       // Yield stress (kinematic)
    k           [0 2 -1 0 0 0 0] 0.01;      // Consistency index
    n           [0 0  0 0 0 0 0] 0.5;       // Power index
}
```

**Special cases of Herschel-Bulkley:**

- **$n = 1, \tau_0 > 0$** → Bingham plastic
- **$n = 1, \tau_0 = 0$** → Newtonian
- **$n < 1, \tau_0 = 0$** → Power-law shear-thinning
- **$n < 1, \tau_0 > 0$** → General Herschel-Bulkley

### Model Comparison Summary

```
┌────────────────────────────────────────────────────────────────────────┐
│              VISCOSITY MODEL DECISION TREE                             │
├──────────────────────────────────────────────────────��────────────────┤
│                                                                        │
│  Does the fluid have a yield stress?                                   │
│      │                                                                 │
│      ├── YES ──► HerschelBulkley                                       │
│      │           (set τ₀ > 0, choose n for post-yield behavior)        │
│      │                                                                 │
│      └── NO ──► Is viscosity constant?                                 │
│                  │                                                     │
│                  ├── YES ──► Newtonian                                 │
│                  │                                                     │
│                  └── NO ──► Do you have plateau data (ν₀, ν∞)?         │
│                              │                                         │
│                              ├── YES ──► BirdCarreau or CrossPowerLaw  │
│                              │                                         │
│                              └── NO ──► powerLaw                       │
│                                         (simpler, fewer parameters)    │
│                                                                        │
└────────────────────────────────────────────────────────────────────────┘
```

## Configuring Viscosity in OpenFOAM

### The `transportProperties` Dictionary

The viscosity model is configured in `constant/transportProperties`.
This file is read at the start of the simulation and (for non-Newtonian
models) the viscosity field is updated every time step.

```
  Case Directory
  ├── 0/
  │   ├── U
  │   ├── p
  │   └── ...
  ├── constant/
  │   ├── transportProperties    ◄── Viscosity model lives here
  │   ├── turbulenceProperties
  │   └── polyMesh/
  └── system/
      ├── controlDict
      ├── fvSchemes
      └── fvSolution
```

See [OpenFOAM Cases](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md) for the full case directory structure.

### Structure of the Dictionary

Every `transportProperties` file follows this pattern:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

// 1. Select the model
transportModel  <modelName>;

// 2. For Newtonian, just set nu:
nu              [0 2 -1 0 0 0 0] <value>;

// 3. For non-Newtonian, add a coefficients sub-dictionary:
<modelName>Coeffs
{
    // model-specific parameters with dimensions
}
```

### Complete Working Examples

#### Example 1: Water (Newtonian)

A common water-like example. The cavity project uses the same Newtonian syntax
but sets `nu = 0.001` to obtain `Re = 100`; this `1e-06` example is not its value.

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 1e-06;
```

#### Example 2: Blood Flow (Bird-Carreau)

Blood is shear-thinning — red blood cells aggregate at low shear rates,
increasing viscosity. At higher shear rates the aggregates break apart.

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

transportModel  BirdCarreau;

BirdCarreauCoeffs
{
    nu0         [0 2 -1 0 0 0 0] 3.3e-06;
    nuInf       [0 2 -1 0 0 0 0] 3.0e-06;
    k           [0 0  1 0 0 0 0] 8.2;
    n           [0 0  0 0 0 0 0] 0.2;
}
```

> **Teaching values:** The Cross/Bird–Carreau coefficients above illustrate
> dictionary syntax and weak shear thinning. They are not a validated blood
> rheology fit. Density, temperature, composition, and shear-rate range must
> match the measured material data; fit both plateaus and the transition range
> before using a model quantitatively.

#### Example 3: Polymer Melt (Power Law)

A molten polymer with strong shear-thinning behavior:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

transportModel  powerLaw;

powerLawCoeffs
{
    nuMax       [0 2 -1 0 0 0 0] 1e+00;
    nuMin       [0 2 -1 0 0 0 0] 1e-06;
    k           [0 2 -1 0 0 0 0] 0.1;
    n           [0 0  0 0 0 0 0] 0.3;
}
```

#### Example 4: Concrete Slurry (Herschel-Bulkley)

Fresh concrete is a yield-stress fluid — it holds its shape until a critical stress is exceeded, then flows with shear-thinning behavior.

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

transportModel  HerschelBulkley;

HerschelBulkleyCoeffs
{
    nu0         [0 2 -1 0 0 0 0] 100;
    tau0        [0 2 -2 0 0 0 0] 0.5;
    k           [0 2 -1 0 0 0 0] 0.005;
    n           [0 0  0 0 0 0 0] 0.65;
}
```

## Effective Viscosity in Turbulent Flows

In turbulent simulations, the **total effective viscosity** is the sum of molecular (laminar) and turbulent contributions:

$$
\nu_{\text{eff}} = \nu + \nu_t
$$

where:

- $\nu$ = molecular kinematic viscosity (from transportProperties)
- $\nu_t$ = turbulent (eddy) viscosity (computed by turbulence model)
- $\nu_{\text{eff}}$ = total effective viscosity (used in momentum equation)

### How It Works in OpenFOAM

```
┌─────────────────────┐     ┌──────────────────────────┐
│  transportProperties│     │  turbulenceProperties    │
│                     │     │                          │
│  transportModel     │     │  simulationType  RAS;    │
│  Newtonian;         │     │  RAS { model kEpsilon; } │
│  nu  1e-06;         │     │                          │
└────────┬────────────┘     └────────────┬─────────────┘
         │                               │
         ▼                               ▼
    ┌──────────┐                   ┌────────────┐
    │  ν field │                   │  ν_t field │
    │  (nu)    │                   │  (nut)     │
    └────┬─────┘                   └─────┬──────┘
         │                               │
         └───────────┬───────────────────┘
                     ▼
              ┌──────────────┐
              │  ν_eff field │
              │  (nuEff)     │
              │  = ν + ν_t   │
              └──────────────┘
```

**Fields in OpenFOAM:**

| Field   | File | Meaning                        | Source                    |
|---------|------|--------------------------------|---------------------------|
| `nu`    | —    | Molecular kinematic viscosity  | `transportProperties`     |
| `nut`   | `0/nut` | Turbulent kinematic viscosity | Turbulence model ($k$-$\varepsilon$, $k$-$\omega$, etc.) |
| `nuEff` | —    | Effective total viscosity       | $\nu + \nu_t$ (computed internally) |

> **📌 Note:** In laminar simulations (`simulationType laminar;`), $\nu_t = 0$
> everywhere, so $\nu_{\text{eff}} = \nu$. The molecular viscosity from `transportProperties`
> is the only viscosity that matters.

For detailed coverage of how turbulence models compute $\nu_t$, see [Turbulence Models](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md).

### Turbulent Viscosity Ratio

A useful diagnostic: the ratio $\nu_t / \nu$ tells you how much the turbulence dominates over molecular diffusion:

| $\nu_t / \nu$    | Regime                                         |
|------------|-------------------------------------------------|
| 0          | Laminar (no turbulence model or laminar sim)    |
| 1 – 10     | Low-turbulence regions (near walls, wake edges) |
| 10 – 1000  | Fully turbulent regions                         |
| > 1000     | Possibly unphysical — check your setup          |

## Temperature-Dependent Viscosity

For many real fluids, viscosity is strongly temperature-dependent. OpenFOAM's standard `transportProperties` models assume **isothermal** (constant temperature) conditions.

### When Temperature Matters

| Fluid    | $\mu$ at 20 °C ($\mathrm{Pa\,s}$) | $\mu$ at 80 °C ($\mathrm{Pa\,s}$) | Change   |
|----------|--------------------|--------------------|----------|
| Water    | $1.002 \times 10^{-3}$      | $0.354 \times 10^{-3}$      | −65 %    |
| Engine oil | 0.2              | 0.02               | −90 %    |
| Air      | $1.825 \times 10^{-5}$      | $2.09 \times 10^{-5}$       | +15 %    |

For **non-isothermal** simulations, you need one of:

I. **Sutherland's law** (gases): built into compressible solvers like `buoyantSimpleFoam` and `rhoCentralFoam`

$$
\mu(T) = \frac{A_s \cdot T^{3/2}}{T + T_s}
$$

For air: $A_s = 1.458 \times 10^{-6}$ kg/(m·s·K$^{0.5}$), $T_s = 110.4$ K

II. **Custom `coded` function objects** or **custom transport models** for liquid-phase temperature dependence

III. **Polynomial fits** via the `polynomial` transport model in some OpenFOAM distributions

> **💡 Tip:** If temperature variation is small (< 10–15 °C), a constant
> viscosity at the mean temperature is usually a reasonable approximation.
> For larger variations, consider the `buoyant*` family of solvers that
> handle variable properties natively.

## Common Pitfalls & Troubleshooting

### Pitfall #1: Dynamic vs Kinematic Viscosity Mix-Up

The most common mistake. Material data sheets typically report **dynamic viscosity $\mu$** in $\mathrm{Pa\,s}$ or cP, but OpenFOAM expects **kinematic viscosity $\nu$** in $\mathrm{m}^2/\mathrm{s}$.

$$
\nu = \mu / \rho
$$

Example: Water at 20 °C

- $\mu = 1.002 \times 10^{-3}$ $\mathrm{Pa\,s}$
- $\rho = 998.2$ $\mathrm{kg/m}^3$
- $\nu = 1.002 \times 10^{-3} / 998.2 = 1.004 \times 10^{-6}$ $\mathrm{m}^2/\mathrm{s}$ ✓

Common mistake: using $\mu$ = 1e-3 directly as nu → 1000× too high!

| Symptom                                   | Likely Cause                     | Fix                              |
|-------------------------------------------|----------------------------------|----------------------------------|
| Re is 1000× too low                       | Used $\mu$ instead of $\nu$              | Divide by density                |
| Flow looks completely laminar when it shouldn't | Viscosity too high          | Check units — should be $\mathrm{m}^2/\mathrm{s}$     |
| Wildly wrong drag/lift coefficients       | Wrong viscosity magnitude        | Verify $\nu$ against known values    |

### Pitfall #2: Wrong Dimensions in the Dictionary

OpenFOAM checks dimensions at runtime. A mismatch causes an immediate crash with a dimension error.

```
  // WRONG — Pa·s dimensions for kinematic viscosity
  nu    [1 -1 -1 0 0 0 0] 1e-06;    // ← This is [kg/(m·s)] = Pa·s

  // CORRECT — m²/s dimensions
  nu    [0 2 -1 0 0 0 0] 1e-06;     // ← This is [m²/s]
```

**Common dimension sets:**

| Quantity                | Dimensions `[kg m s K mol A cd]` |
|-------------------------|-----------------------------------|
| Kinematic viscosity $\nu$   | `[0 2 -1 0 0 0 0]`               |
| Dynamic viscosity $\mu$     | `[1 -1 -1 0 0 0 0]`              |
| Kinematic stress $\tau/\rho$    | `[0 2 -2 0 0 0 0]`               |
| Time constant           | `[0 0  1 0 0 0 0]`               |
| Dimensionless           | `[0 0  0 0 0 0 0]`               |

### Pitfall #3: Non-Newtonian Models Diverging

Non-Newtonian viscosity models compute $\nu$ as a function of the local
shear rate. Extreme values of shear rate can cause:

- **Very high viscosity** in near-stagnant regions → stiff system, slow convergence
- **Very low viscosity** at high-shear regions → instability, oscillation

**Solutions:**

1. Always set `nuMin` and `nuMax` bounds (power-law model)
2. For Herschel-Bulkley, use a finite `nu0` (not infinity) to regularize the yield region
3. Use under-relaxation in `fvSolution` for non-Newtonian cases
4. Start with a converged Newtonian solution, then switch to non-Newtonian

### Pitfall #4: Missing `transportProperties` Keywords

```
  // WRONG — missing transportModel keyword
  nu    [0 2 -1 0 0 0 0] 1e-06;

  // CORRECT
  transportModel  Newtonian;
  nu              [0 2 -1 0 0 0 0] 1e-06;
```

If the `transportModel` keyword is missing, the solver may use a default
or crash with an unhelpful error message.

### Pitfall #5: Non-Newtonian with Turbulence

When combining non-Newtonian models with turbulence:

- The molecular viscosity $\nu$ varies spatially → $\nu_{\text{eff}} = \nu(\dot{\gamma}) + \nu_t$
- This can cause convergence issues, especially with RANS models
- Standard wall functions assume Newtonian behavior near the wall
- Consider using low-Re turbulence models or resolving the viscous sublayer

> **⚠️ Warning:** Most standard wall functions in OpenFOAM are designed
> for Newtonian fluids. Non-Newtonian turbulent simulations near walls
> require careful validation. See [Boundary Conditions](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md)
> for wall function details.

## Practical Tips

### Checking Available Models

List all available transport models in your OpenFOAM installation:

```bash
# Set transportModel to a nonsense value and run — OpenFOAM will list valid options
# In transportProperties:
#   transportModel  banana;
# Then run the solver:
simpleFoam -help
# In a disposable copy of a compatible simpleFoam case, an invalid model
# selection will report available types at startup; icoFoam 6 ignores it.
# simpleFoam > log.model-check 2>&1
# sed -n '1,80p' log.model-check
# Output will include something like:
#   Valid transportModel types:
#       BirdCarreau CrossPowerLaw HerschelBulkley Newtonian powerLaw
```

### Inspecting the Viscosity Field

For non-Newtonian models, you can write the viscosity field to inspect it:

```c
// In system/controlDict, add to the functions section:

functions
{
    writeNu
    {
        type            writeObjects;
        libs            ("libutilityFunctionObjects.so");
        objects         (nu);
        writeControl    writeTime;
    }
}
```

Then visualize in ParaView to check for unphysical values.

### Dimensional Consistency Check

Before running, mentally verify your Reynolds number:

$$
Re = U \cdot L / \nu
$$

Hypothetical water-filled cavity (not project 01), $U = 1$ m/s, $L = 0.1$ m, $\nu = 10^{-6}$ $\mathrm{m}^2/\mathrm{s}$
$Re = 1 \times 0.1 / 10^{-6} = 100{,}000$ (turbulent!)

If that's not what you intended, check your $\nu$ value.

## Quick Reference

```
╔══════════════════════════════════════════════════════════════════════════════╗
║                    VISCOSITY MODELS — QUICK REFERENCE                        ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  Model              Equation                     Key Parameters              ║
║  ─────              ────────                     ──────────────              ║
║  Newtonian          ν = const                    nu                          ║
║  powerLaw           ν = k · γ̇ⁿ⁻¹               k, n, nuMin, nuMax            ║
║  CrossPowerLaw      ν = ν∞ + (ν₀−ν∞)/(1+(mγ̇)ⁿ)  nu0, nuInf, m, n             ║
║  BirdCarreau        ν = ν∞ + (ν₀−ν∞)[1+(kγ̇)²]^…  nu0, nuInf, k, n            ║
║  HerschelBulkley    τ = τ₀ + k · γ̇ⁿ            nu0, tau0, k, n               ║
║                                                                              ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  Typical ν values:                                                           ║
║  ─────────────────                                                           ║
║  Air (20 °C)        1.516 × 10⁻⁵ m²/s                                        ║
║  Water (20 °C)      1.004 × 10⁻⁶ m²/s                                        ║
║  Blood (37 °C)      ~3 × 10⁻⁶ m²/s  (shear-rate dependent)                   ║
║  Engine oil         ~1 × 10⁻⁴ m²/s                                           ║
║                                                                              ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  Remember:  OpenFOAM uses KINEMATIC viscosity  ν = μ/ρ   [m²/s]              ║
║             Effective:  ν_eff = ν + ν_t   (turbulent flows)                  ║
║             Configure in:  constant/transportProperties                      ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

### Model Selection Quick Guide

| Fluid Type              | Recommended Model    | Typical n  | Example Fluid        |
|-------------------------|----------------------|------------|----------------------|
| Simple liquids / gases  | `Newtonian`          | —          | Water, air, oils     |
| Weak shear-thinning     | `CrossPowerLaw`      | 0.3 – 0.8 | Blood, dilute polymer|
| Strong shear-thinning   | `powerLaw`           | 0.1 – 0.5 | Polymer melts        |
| Biological fluids       | `BirdCarreau`        | 0.2 – 0.6 | Blood, synovial fluid|
| Yield-stress fluids     | `HerschelBulkley`    | 0.3 – 1.0 | Concrete, mud, paste |

## Where to Go Next

| Topic                              | Note                                         |
|------------------------------------|----------------------------------------------|
| CFD fundamentals                   | `notes/01_short_intro_to_cfd.md`             |
| OpenFOAM case structure            | `notes/02_openfoam_cases.md`                 |
| Dictionary syntax & dimensions     | `notes/03_openfoam_dictionaries.md`          |
| Meshing                            | `notes/04_meshing.md`                        |
| Boundary conditions & wall funcs   | `notes/05_boundary_conditions.md`            |
| Turbulence models & $\nu_t$            | `notes/06_turbulence_models.md`          |
| Parallel computing                 | `notes/07_parallelization.md`                |
| CFL number & time stepping         | `notes/08_cfl_number.md`                     |
| Linear solvers                     | `notes/09_linear_solvers.md`                 |
| icoFoam solver analysis            | `notes/10_icofoam_solver_analysis.md`        |
| **Hands-on:** Lid-driven cavity    | `projects/01_lid_driven_cavity/`             |
