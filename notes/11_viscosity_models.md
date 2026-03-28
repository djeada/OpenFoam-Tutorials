# Viscosity and Viscosity Models in OpenFOAM

> **Lecture 7 — Part 1: Models in OpenFOAM — Viscosity**
>
> **Prerequisites:** [Short Intro to CFD](01_short_intro_to_cfd.md),
> [OpenFOAM Dictionaries](03_openfoam_dictionaries.md),
> [Turbulence Models](06_turbulence_models.md)
>
> **Key files:** `constant/transportProperties`

---

## 1. What Is Viscosity?

Viscosity is a fluid's **resistance to deformation** under shear stress.
Imagine two parallel plates with fluid between them — the top plate moves
while the bottom plate stays fixed. The force needed to maintain that
motion depends on the fluid's viscosity.

```
        Velocity u
        ──────────►
  ┌──────────────────────────────────────┐  ◄── Moving plate
  │  → → → → → → → → → → → → → → → →  │
  │  → → → → → → → → → → → → → → →    │
  │  → → → → → → → → → → → → →        │      Shear stress
  │  → → → → → → → → → → →            │      τ = μ (du/dy)
  │  → → → → → → → → →                │
  │  → → → → → → →                    │
  │  → → → → →                        │
  │  → → →                            │
  │  →                                 │
  ├──────────────────────────────────────┤  ◄── Fixed plate
  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
```

The velocity profile is linear for a Newtonian fluid — this is
**Couette flow**, the simplest example of viscous shear.

### Dynamic vs Kinematic Viscosity

| Symbol | Name                 | Definition | SI Units     | Physical Meaning                         |
|--------|----------------------|------------|--------------|------------------------------------------|
| μ      | Dynamic viscosity    | —          | Pa·s (kg/m·s)| Shear stress per unit shear rate         |
| ν      | Kinematic viscosity  | ν = μ / ρ  | m²/s         | Momentum diffusivity                     |

**Dynamic viscosity (μ)** relates shear stress to shear rate directly:

```
  τ = μ · (du/dy)

  where:
    τ   = shear stress        [Pa]
    μ   = dynamic viscosity   [Pa·s]
    du/dy = velocity gradient  [1/s]
```

**Kinematic viscosity (ν)** is dynamic viscosity divided by density.
It appears naturally in the incompressible Navier-Stokes equations:

```
  ∂u/∂t + (u · ∇)u = −(1/ρ)∇p + ν ∇²u

  Here ν = μ/ρ controls how fast momentum diffuses through the fluid.
```

### Why OpenFOAM Uses Kinematic Viscosity

Most OpenFOAM incompressible solvers (`icoFoam`, `simpleFoam`, `pimpleFoam`)
divide the entire momentum equation by density ρ. This means:

- Pressure becomes p/ρ (kinematic pressure, units m²/s²)
- Viscosity becomes ν = μ/ρ (kinematic viscosity, units m²/s)

> **⚠️ Warning:** When you see `nu` in OpenFOAM dictionaries, it is **always**
> kinematic viscosity ν (m²/s), not dynamic viscosity μ (Pa·s).
> If a data sheet gives you μ, divide by ρ first!

**Typical values:**

| Fluid          | μ (Pa·s)          | ρ (kg/m³)  | ν (m²/s)           |
|----------------|--------------------|------------|---------------------|
| Air (20 °C)    | 1.825 × 10⁻⁵      | 1.205      | 1.516 × 10⁻⁵       |
| Water (20 °C)  | 1.002 × 10⁻³      | 998.2      | 1.004 × 10⁻⁶       |
| Engine oil     | 0.1 – 0.3         | 880        | 1.1 × 10⁻⁴ – 3.4 × 10⁻⁴ |
| Honey          | 2 – 10            | 1420       | 1.4 × 10⁻³ – 7 × 10⁻³   |
| Blood (37 °C)  | 3 – 4 × 10⁻³      | 1060       | 2.8 – 3.8 × 10⁻⁶   |

---

## 2. Newtonian vs Non-Newtonian Fluids

The distinction comes down to one question: **does viscosity change
with shear rate?**

### Newtonian Fluids

For a Newtonian fluid, viscosity is constant regardless of how fast the
fluid is being sheared:

```
  τ = μ · γ̇       (linear relationship)

  where γ̇ = du/dy = shear rate [1/s]
```

Most common engineering fluids are Newtonian: **water, air, most gases,
light oils, and simple organic solvents**.

### Non-Newtonian Fluids

For non-Newtonian fluids, the **apparent viscosity** η changes with
shear rate γ̇. There are several categories:

| Type                | Behavior                          | Examples                          |
|---------------------|-----------------------------------|-----------------------------------|
| Shear-thinning      | η decreases as γ̇ increases       | Blood, paint, polymer melts, ketchup |
| Shear-thickening    | η increases as γ̇ increases       | Cornstarch in water, wet sand     |
| Bingham / Yield     | No flow below critical stress τ₀  | Toothpaste, concrete, drilling mud|
| Thixotropic         | η decreases over time at constant γ̇ | Yogurt, some gels              |
| Rheopectic          | η increases over time at constant γ̇  | Gypsum paste                   |

> **📌 Note:** OpenFOAM's built-in transport models handle the first three
> categories (shear-thinning, shear-thickening, yield-stress). Thixotropic
> and rheopectic behavior requires custom implementations.

### Shear Stress vs Shear Rate — Visual Comparison

```
  Shear
  Stress
  τ ▲
    │                                      ╱  Shear-thickening
    │                                   ╱╱    (dilatant)
    │                                ╱╱
    │                             ╱╱
    │                          ╱╱
    │                       ╱╱       ╱╱╱╱╱╱╱╱╱╱  Newtonian
    │                    ╱╱     ╱╱╱╱╱
    │                 ╱╱   ╱╱╱╱
    │              ╱╱ ╱╱╱╱
    │           ╱╱╱╱╱
    │        ╱╱╱╱              ╱──────────────  Shear-thinning
    │     ╱╱╱           ╱╱╱╱╱╱                  (pseudoplastic)
    │  ╱╱╱       ╱╱╱╱╱╱╱
    │╱╱╱  ╱╱╱╱╱╱╱
    │╱╱╱╱╱
    ┼──────────────────────────────────────────► γ̇  (Shear Rate)
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

These fluids behave as a **solid** below a critical shear stress τ₀
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

Real-world examples: toothpaste stays on the brush (τ < τ₀) but flows
when you squeeze (τ > τ₀). Concrete doesn't flow under its own weight
until vibrated.

---

## 3. Viscosity Models Available in OpenFOAM

OpenFOAM provides several viscosity models through the `transportModel`
keyword in `constant/transportProperties`. The solver then computes the
apparent kinematic viscosity field `nu` at each cell and time step.

### 3.1 Newtonian

The simplest model — constant viscosity everywhere.

**Equation:**

```
  ν = constant
```

**Parameters:**

| Parameter | Symbol | Meaning             | Units |
|-----------|--------|---------------------|-------|
| `nu`      | ν      | Kinematic viscosity | m²/s  |

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
> `[kg m s K mol A cd]`. So `[0 2 -1 0 0 0 0]` means m²/s — exactly
> the units of kinematic viscosity. See [OpenFOAM Dictionaries](03_openfoam_dictionaries.md)
> for details on dimension sets.

### 3.2 Power Law

Models shear-thinning or shear-thickening behavior with a simple power
function.

**Equation:**

```
  ν(γ̇) = k · γ̇ⁿ⁻¹

  Clamped:  nuMin ≤ ν ≤ nuMax
```

**Parameters:**

| Parameter | Symbol | Meaning                        | Units      |
|-----------|--------|--------------------------------|------------|
| `nuMax`   | ν_max  | Upper viscosity limit          | m²/s       |
| `nuMin`   | ν_min  | Lower viscosity limit          | m²/s       |
| `k`       | k      | Consistency index              | m²/s^n     |
| `n`       | n      | Power-law index                | dimensionless |

- **n < 1** → shear-thinning (pseudoplastic)
- **n = 1** → Newtonian
- **n > 1** → shear-thickening (dilatant)

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

### 3.3 Cross Power Law

The Cross model handles the transition between a zero-shear-rate plateau
(η₀) and an infinite-shear-rate plateau (η∞) — more physically realistic
than the plain power law.

**Equation:**

```
  ν(γ̇) = ν∞ + (ν₀ − ν∞) / (1 + (m · γ̇)ⁿ)
```

**Parameters:**

| Parameter | Symbol | Meaning                             | Units         |
|-----------|--------|-------------------------------------|---------------|
| `nu0`     | ν₀     | Zero-shear-rate kinematic viscosity | m²/s          |
| `nuInf`   | ν∞     | Infinite-shear-rate viscosity       | m²/s          |
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

### 3.4 Bird-Carreau

Similar to the Cross model but uses a different mathematical form.
Very popular for polymer solutions and biological fluids.

**Equation:**

```
  ν(γ̇) = ν∞ + (ν₀ − ν∞) · [1 + (k · γ̇)²]^((n−1)/2)
```

**Parameters:**

| Parameter | Symbol | Meaning                             | Units         |
|-----------|--------|-------------------------------------|---------------|
| `nu0`     | ν₀     | Zero-shear-rate kinematic viscosity | m²/s          |
| `nuInf`   | ν∞     | Infinite-shear-rate viscosity       | m²/s          |
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

### 3.5 Herschel-Bulkley

Models **yield-stress fluids** — materials that behave as solids below
a threshold stress and flow like a power-law fluid above it.

**Equation:**

```
  τ = τ₀ + k · γ̇ⁿ          (for τ > τ₀)

  ν(γ̇) = min(ν₀, τ₀/γ̇ + k · γ̇ⁿ⁻¹)
```

**Parameters:**

| Parameter | Symbol | Meaning                        | Units         |
|-----------|--------|--------------------------------|---------------|
| `nu0`     | ν₀     | Limiting viscosity (yield reg.)| m²/s          |
| `tau0`    | τ₀     | Yield stress                   | m²/s²  *      |
| `k`       | k      | Consistency index              | m²/s^n        |
| `n`       | n      | Power-law index                | dimensionless |

> \* Because OpenFOAM works with kinematic quantities (divided by ρ),
> `tau0` has units of m²/s² (kinematic stress), not Pa.

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
- **n = 1, τ₀ > 0** → Bingham plastic
- **n = 1, τ₀ = 0** → Newtonian
- **n < 1, τ₀ = 0** → Power-law shear-thinning
- **n < 1, τ₀ > 0** → General Herschel-Bulkley

### 3.6 Model Comparison Summary

```
┌───────────────────────────────────────────────────────────────────────┐
│              VISCOSITY MODEL DECISION TREE                           │
├──────────────────────────────────────────────────────��────────────────┤
│                                                                       │
│  Does the fluid have a yield stress?                                 │
│      │                                                                │
│      ├── YES ──► HerschelBulkley                                     │
│      │           (set τ₀ > 0, choose n for post-yield behavior)      │
│      │                                                                │
│      └── NO ──► Is viscosity constant?                               │
│                  │                                                     │
│                  ├── YES ──► Newtonian                                │
│                  │                                                     │
│                  └── NO ──► Do you have plateau data (ν₀, ν∞)?       │
│                              │                                        │
│                              ├── YES ──► BirdCarreau or CrossPowerLaw│
│                              │                                        │
│                              └── NO ──► powerLaw                     │
│                                         (simpler, fewer parameters)   │
│                                                                       │
└───────────────────────────────────────────────────────────────────────┘
```

---

## 4. Configuring Viscosity in OpenFOAM

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

See [OpenFOAM Cases](02_openfoam_cases.md) for the full case directory structure.

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

The simplest and most common configuration. Used in
`projects/01_lid_driven_cavity/` and many tutorial cases.

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

> **💡 Tip:** For blood flow at physiological conditions (37 °C, ρ ≈ 1060 kg/m³),
> typical values are ν₀ ≈ 3.3 × 10⁻⁶ m²/s and ν∞ ≈ 3.0 × 10⁻⁶ m²/s.
> The shear-thinning effect is moderate compared to polymer melts.

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

Fresh concrete is a yield-stress fluid — it holds its shape until a
critical stress is exceeded, then flows with shear-thinning behavior.

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

---

## 5. Effective Viscosity in Turbulent Flows

In turbulent simulations, the **total effective viscosity** is the sum
of molecular (laminar) and turbulent contributions:

```
  ν_eff = ν + ν_t

  where:
    ν     = molecular kinematic viscosity  (from transportProperties)
    ν_t   = turbulent (eddy) viscosity     (computed by turbulence model)
    ν_eff = total effective viscosity       (used in momentum equation)
```

### How It Works in OpenFOAM

```
┌─────────────────────┐     ┌─────────────────────────┐
│  transportProperties│     │  turbulenceProperties    │
│                     │     │                          │
│  transportModel     │     │  simulationType  RAS;    │
│  Newtonian;         │     │  RAS { model kEpsilon; } │
│  nu  1e-06;         │     │                          │
└────────┬────────────┘     └────────────┬─────────────┘
         │                               │
         ▼                               ▼
    ┌─────────┐                    ┌───────────┐
    │  ν field│                    │  ν_t field │
    │  (nu)   │                    │  (nut)     │
    └────┬────┘                    └─────┬─────┘
         │                               │
         └───────────┬───────────────────┘
                     ▼
              ┌──────────────┐
              │  ν_eff field │
              │  (nuEff)     │
              │  = ν + ν_t   │
              └──────────────┘
```

**Key fields in OpenFOAM:**

| Field   | File | Meaning                        | Source                    |
|---------|------|--------------------------------|---------------------------|
| `nu`    | —    | Molecular kinematic viscosity  | `transportProperties`     |
| `nut`   | `0/nut` | Turbulent kinematic viscosity | Turbulence model (k-ε, k-ω, etc.) |
| `nuEff` | —    | Effective total viscosity       | ν + ν_t (computed internally) |

> **📌 Note:** In laminar simulations (`simulationType laminar;`), ν_t = 0
> everywhere, so ν_eff = ν. The molecular viscosity from `transportProperties`
> is the only viscosity that matters.

For detailed coverage of how turbulence models compute ν_t, see
[Turbulence Models](06_turbulence_models.md).

### Turbulent Viscosity Ratio

A useful diagnostic: the ratio ν_t / ν tells you how much the turbulence
dominates over molecular diffusion:

| ν_t / ν    | Regime                                         |
|------------|-------------------------------------------------|
| 0          | Laminar (no turbulence model or laminar sim)    |
| 1 – 10     | Low-turbulence regions (near walls, wake edges) |
| 10 – 1000  | Fully turbulent regions                         |
| > 1000     | Possibly unphysical — check your setup          |

---

## 6. Temperature-Dependent Viscosity

For many real fluids, viscosity is strongly temperature-dependent.
OpenFOAM's standard `transportProperties` models assume **isothermal**
(constant temperature) conditions.

### When Temperature Matters

| Fluid    | μ at 20 °C (Pa·s) | μ at 80 °C (Pa·s) | Change   |
|----------|--------------------|--------------------|----------|
| Water    | 1.002 × 10⁻³      | 0.354 × 10⁻³      | −65 %    |
| Engine oil | 0.2              | 0.02               | −90 %    |
| Air      | 1.825 × 10⁻⁵      | 2.09 × 10⁻⁵       | +15 %    |

For **non-isothermal** simulations, you need one of:

1. **Sutherland's law** (gases): built into compressible solvers like
   `buoyantSimpleFoam` and `rhoCentralFoam`

```
  μ(T) = A_s · T^(3/2) / (T + T_s)

  For air:  A_s = 1.458 × 10⁻⁶ kg/(m·s·K^0.5)
            T_s = 110.4 K
```

2. **Custom `coded` function objects** or **custom transport models**
   for liquid-phase temperature dependence

3. **Polynomial fits** via the `polynomial` transport model in some
   OpenFOAM distributions

> **💡 Tip:** If temperature variation is small (< 10–15 °C), a constant
> viscosity at the mean temperature is usually a reasonable approximation.
> For larger variations, consider the `buoyant*` family of solvers that
> handle variable properties natively.

---

## 7. Common Pitfalls & Troubleshooting

### Pitfall #1: Dynamic vs Kinematic Viscosity Mix-Up

The most common mistake. Material data sheets typically report **dynamic
viscosity μ** in Pa·s or cP, but OpenFOAM expects **kinematic viscosity
ν** in m²/s.

```
  ν = μ / ρ

  Example: Water at 20 °C
    μ = 1.002 × 10⁻³ Pa·s
    ρ = 998.2 kg/m³
    ν = 1.002e-3 / 998.2 = 1.004 × 10⁻⁶ m²/s  ✓

  Common mistake: using μ = 1e-3 directly as nu → 1000× too high!
```

| Symptom                                   | Likely Cause                     | Fix                              |
|-------------------------------------------|----------------------------------|----------------------------------|
| Re is 1000× too low                       | Used μ instead of ν              | Divide by density                |
| Flow looks completely laminar when it shouldn't | Viscosity too high          | Check units — should be m²/s     |
| Wildly wrong drag/lift coefficients       | Wrong viscosity magnitude        | Verify ν against known values    |

### Pitfall #2: Wrong Dimensions in the Dictionary

OpenFOAM checks dimensions at runtime. A mismatch causes an immediate
crash with a dimension error.

```
  // WRONG — Pa·s dimensions for kinematic viscosity
  nu    [1 -1 -1 0 0 0 0] 1e-06;    // ← This is [kg/(m·s)] = Pa·s

  // CORRECT — m²/s dimensions
  nu    [0 2 -1 0 0 0 0] 1e-06;     // ← This is [m²/s]
```

**Common dimension sets:**

| Quantity                | Dimensions `[kg m s K mol A cd]` |
|-------------------------|-----------------------------------|
| Kinematic viscosity ν   | `[0 2 -1 0 0 0 0]`               |
| Dynamic viscosity μ     | `[1 -1 -1 0 0 0 0]`              |
| Kinematic stress τ/ρ    | `[0 2 -2 0 0 0 0]`               |
| Time constant           | `[0 0  1 0 0 0 0]`               |
| Dimensionless           | `[0 0  0 0 0 0 0]`               |

### Pitfall #3: Non-Newtonian Models Diverging

Non-Newtonian viscosity models compute ν as a function of the local
shear rate. Extreme values of shear rate can cause:

- **Very high viscosity** in near-stagnant regions → stiff system, slow convergence
- **Very low viscosity** at high-shear regions → instability, oscillation

**Solutions:**

1. Always set `nuMin` and `nuMax` bounds (power-law model)
2. For Herschel-Bulkley, use a finite `nu0` (not infinity) to regularize
   the yield region
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

- The molecular viscosity ν varies spatially → ν_eff = ν(γ̇) + ν_t
- This can cause convergence issues, especially with RANS models
- Standard wall functions assume Newtonian behavior near the wall
- Consider using low-Re turbulence models or resolving the viscous sublayer

> **⚠️ Warning:** Most standard wall functions in OpenFOAM are designed
> for Newtonian fluids. Non-Newtonian turbulent simulations near walls
> require careful validation. See [Boundary Conditions](05_boundary_conditions.md)
> for wall function details.

---

## 8. Practical Tips

### Checking Available Models

List all available transport models in your OpenFOAM installation:

```bash
# Set transportModel to a nonsense value and run — OpenFOAM will list valid options
# In transportProperties:
#   transportModel  banana;
# Then run the solver:
icoFoam 2>&1 | head -30
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

```
  Re = U · L / ν

  Example: Lid-driven cavity, U = 1 m/s, L = 0.1 m, ν = 1e-6 m²/s
  Re = 1 × 0.1 / 1e-6 = 100,000  (turbulent!)

  If that's not what you intended, check your ν value.
```

---

## 9. Quick Reference

```
╔══════════════════════════════════════════════════════════════════════════════╗
║                    VISCOSITY MODELS — QUICK REFERENCE                      ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                            ║
║  Model              Equation                     Key Parameters            ║
║  ─────              ────────                     ──────────────            ║
║  Newtonian          ν = const                    nu                        ║
║  powerLaw           ν = k · γ̇ⁿ⁻¹               k, n, nuMin, nuMax        ║
║  CrossPowerLaw      ν = ν∞ + (ν₀−ν∞)/(1+(mγ̇)ⁿ)  nu0, nuInf, m, n       ║
║  BirdCarreau        ν = ν∞ + (ν₀−ν∞)[1+(kγ̇)²]^…  nu0, nuInf, k, n      ║
║  HerschelBulkley    τ = τ₀ + k · γ̇ⁿ            nu0, tau0, k, n          ║
║                                                                            ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                            ║
║  Typical ν values:                                                         ║
║  ─────────────────                                                         ║
║  Air (20 °C)        1.516 × 10⁻⁵ m²/s                                    ║
║  Water (20 °C)      1.004 × 10⁻⁶ m²/s                                    ║
║  Blood (37 °C)      ~3 × 10⁻⁶ m²/s  (shear-rate dependent)              ║
║  Engine oil         ~1 × 10⁻⁴ m²/s                                       ║
║                                                                            ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                            ║
║  Remember:  OpenFOAM uses KINEMATIC viscosity  ν = μ/ρ   [m²/s]          ║
║             Effective:  ν_eff = ν + ν_t   (turbulent flows)               ║
║             Configure in:  constant/transportProperties                    ║
║                                                                            ║
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

---

## Where to Go Next

| Topic                              | Note                                         |
|------------------------------------|----------------------------------------------|
| CFD fundamentals                   | `notes/01_short_intro_to_cfd.md`             |
| OpenFOAM case structure            | `notes/02_openfoam_cases.md`                 |
| Dictionary syntax & dimensions     | `notes/03_openfoam_dictionaries.md`          |
| Meshing                            | `notes/04_meshing.md`                        |
| Boundary conditions & wall funcs   | `notes/05_boundary_conditions.md`            |
| Turbulence models & ν_t            | `notes/06_turbulence_models.md`              |
| Parallel computing                 | `notes/07_parallelization.md`                |
| CFL number & time stepping         | `notes/08_cfl_number.md`                     |
| Linear solvers                     | `notes/09_linear_solvers.md`                 |
| icoFoam solver analysis            | `notes/10_icofoam_solver_analysis.md`        |
| **Hands-on:** Lid-driven cavity    | `projects/01_lid_driven_cavity/`             |

---

*Last updated: 2025. Part of the [OpenFOAM-Tutorials](../README.md) repository.*
