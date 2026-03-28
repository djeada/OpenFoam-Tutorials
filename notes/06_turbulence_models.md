# Turbulence and Turbulence Modeling in OpenFOAM

## A Comprehensive Guide for CFD Practitioners

---

## 1. What Is Turbulence?

Turbulence is the chaotic, seemingly random motion of fluid particles. Almost every
flow in engineering practice is turbulent — from air flowing over a car to water in a
pipe. Understanding and modeling turbulence is the single most important challenge
in computational fluid dynamics (CFD).

### The Reynolds Number

The **Reynolds number** (Re) determines whether a flow is laminar or turbulent:

$$Re = \frac{\rho \cdot U \cdot L}{\mu} = \frac{U \cdot L}{\nu}$$

where:
- $\rho$ = fluid density [kg/m³]
- $U$ = characteristic velocity [m/s]
- $L$ = characteristic length [m]
- $\mu$ = dynamic viscosity [Pa·s]
- $\nu$ = kinematic viscosity [m²/s]

| Flow Type     | Typical Re Range | Example                          |
|---------------|------------------|----------------------------------|
| Laminar       | $Re < 2,300$       | Honey pouring, microfluidics     |
| Transitional  | $2,300 < Re < 4,000$ | Pipe flow near critical $Re$    |
| Turbulent     | $Re > 4,000$       | Most engineering flows           |
| Highly turbulent | $Re > 1,000,000$ | Aircraft in flight, rivers     |

> **💡 Tip:** Our lid-driven cavity project (`projects/01_lid_driven_cavity/`) uses
> `ν = 0.001 m²/s` with a lid velocity of 1 m/s and cavity size ~1 m, giving
> $Re \approx 1,000$ — firmly laminar, which is why it uses `icoFoam` with no turbulence
> model. The NACA airfoil project uses `ν = 1e-05 m²/s`, yielding $Re \approx 100,000$ —
> turbulent, requiring a $k$-$\varepsilon$ model with `simpleFoam`.

### Laminar-to-Turbulent Transition

Consider flow over a flat plate. Near the leading edge the flow is smooth and
orderly (laminar). As the boundary layer grows, instabilities amplify until the
flow becomes fully turbulent:

```
  Free stream  →→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→→
                                                              
  ═══════╗            ~~~~~~~~~~~~        ≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋
         ║  Laminar    ~~~~~~~~~~~~ Trans   ≋≋≋≋ Turbulent ≋≋≋≋
         ║  (smooth)   ~~~~~~~~~~~~        ≋≋≋≋ (chaotic)  ≋≋≋≋
  ═══════╝            ~~~~~~~~~~~~        ≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋≋
  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░  WALL
  
  ←── δ thin ──→←─ δ grows ──→←──── δ thick, irregular ────→
  
        Re_x < 5×10⁵              Re_x > 5×10⁵
```

Key features of turbulent flow:
- **Irregular**: Velocity fluctuates randomly in time and space
- **Diffusive**: Greatly enhances mixing of momentum, heat, and mass
- **3D and rotational**: Always three-dimensional, always has vorticity
- **Dissipative**: Kinetic energy is continuously converted to heat
- **Multi-scale**: Contains eddies ranging from the flow domain size down to tiny Kolmogorov scales

---

## 2. The Energy Cascade

One of the most fundamental concepts in turbulence is the **energy cascade**,
described by Richardson (1922) and formalized by Kolmogorov (1941).

```
  ┌─────────────────────────────────────────────────────────────────┐
  │                                                                 │
  │   LARGE EDDIES  ◉◉◉◉◉◉◉◉◉◉◉◉         Scale ~ L (domain size) │
  │   (Energy-containing range)            Production of energy     │
  │          │                                                      │
  │          │  Inertial transfer (no dissipation)                  │
  │          ▼                                                      │
  │   MEDIUM EDDIES  ○○○○○○○○○○○           Scale ~ intermediate    │
  │   (Inertial subrange)                  E(κ) ~ κ^(-5/3)        │
  │          │                                                      │
  │          │  Inertial transfer continues                         │
  │          ▼                                                      │
  │   SMALL EDDIES   ∘∘∘∘∘∘∘∘∘∘∘∘∘∘∘∘     Scale ~ η (Kolmogorov) │
  │   (Dissipation range)                  Viscous dissipation → 🔥│
  │                                                                 │
  │   Energy IN (production) ═══════════► Energy OUT (dissipation)  │
  └─────────────────────────────────────────────────────────────────┘
```

### Kolmogorov Scales

The smallest turbulent scales are set by viscosity and the dissipation rate $\varepsilon$:

Length scale:

$$\eta = (\nu^{3}/\varepsilon)^{1/4}$$

Time scale:

$$\tau = (\nu/\varepsilon)^{1/2}$$

Velocity scale:

$$v = (\nu \cdot \varepsilon)^{1/4}$$

The ratio of the largest to smallest scales grows with Reynolds number:

$$L / \eta \sim Re^{3/4}$$

This is why we cannot resolve all scales directly — even a moderate $Re = 10,000$
flow would need approximately `(10,000)^(9/4) ≈ 10^9.75` grid points for a 3D DNS.
For $Re = 1,000,000$ (a car on a highway), you would need $\sim 10^{27}$ cells. This is why
turbulence **modeling** exists.

> **⚠️ Key insight:** The entire field of turbulence modeling exists because we
> cannot afford to resolve every eddy. Each modeling approach makes different
> trade-offs between accuracy and computational cost.

---

## 3. Turbulence Modeling Hierarchy

```
  ┌──────────────────────────────────────────────────────────────────────┐
  │                                                                      │
  │           DNS  (Direct Numerical Simulation)                         │
  │           ● Resolves ALL scales from L down to η                     │
  │           ● No modeling error — only discretization error             │
  │           ● Cost ~ Re³                                               │
  │           ▲                                                          │
  │           │ Increasing cost & accuracy                               │
  │           │                                                          │
  │           LES  (Large Eddy Simulation)                               │
  │           ● Resolves large eddies, models small eddies               │
  │           ● Sub-grid scale (SGS) model for unresolved scales         │
  │           ● Cost: 10-100× RANS                                       │
  │           ▲                                                          │
  │           │                                                          │
  │           DES / Hybrid  (Detached Eddy Simulation)                   │
  │           ● RANS near walls, LES in free stream                      │
  │           ● Practical compromise for many flows                      │
  │           ▲                                                          │
  │           │                                                          │
  │           RANS  (Reynolds-Averaged Navier-Stokes)                    │
  │           ● Models ALL turbulent scales                              │
  │           ● Solves time-averaged equations                           │
  │           ● k-ε, k-ω, SST, Spalart-Allmaras, RSM                   │
  │           ▲                                                          │
  │           │                                                          │
  │           Laminar  (No turbulence model)                             │
  │           ● Valid only for low-Re flows                              │
  │           ● Example: our lid-driven cavity (Re ≈ 1000)              │
  │                                                                      │
  └──────────────────────────────────────────────────────────────────────┘
    Low cost ◄────────────────────────────────────────────► High cost
    Low fidelity ◄────────────────────────────────────────► High fidelity
```

### Which Approach Should You Use?

| Approach | Grid Points (3D) | Time Steps | Typical Use Case |
|----------|-------------------|------------|------------------|
| Laminar  | $10^{3} - 10^{5}$         | $\sim 10^{2}$       | Microfluidics, creeping flows |
| RANS     | $10^{5} - 10^{7}$         | Steady or $\sim 10^{3}$ | Most engineering design |
| DES      | $10^{6} - 10^{8}$         | $\sim 10^{4}$       | Separated flows, bluff bodies |
| LES      | $10^{7} - 10^{9}$         | $\sim 10^{4} - 10^{5}$  | Combustion, acoustics, mixing |
| DNS      | $10^{9} - 10^{12}$        | $\sim 10^{5} - 10^{6}$  | Fundamental research only |

---

## 4. RANS Models — Deep Dive

RANS (Reynolds-Averaged Navier-Stokes) models decompose every flow variable into
a **mean** part and a **fluctuating** part:

$$u(x,t) = \bar{U}(x) + u'(x,t)$$

Averaging the Navier-Stokes equations produces the **Reynolds stress tensor**
`-ρ⟨u'ᵢu'ⱼ⟩` — six unknown quantities that must be modeled. This is the famous
**closure problem** of turbulence.

### 4.1 Spalart-Allmaras (SA)

**Equations solved:** 1 (modified turbulent viscosity $\tilde{\nu}$)

The Spalart-Allmaras model solves a single transport equation for a quantity
related to the turbulent (eddy) viscosity. It was designed specifically for
aerospace boundary layer flows.

- **Strengths**: Simple, stable, low memory, good for attached boundary layers
- **Weaknesses**: Poor for free shear flows, jet spreading, recirculation zones
- **Best for**: External aerodynamics, wing/fuselage flows, turbomachinery blades
- **OpenFOAM name**: `SpalartAllmaras`
- **Required fields**: `nuTilda`, `nut`

### 4.2 $k$-$\varepsilon$ Standard

**Equations solved:** 2 (turbulent kinetic energy $k$, dissipation rate $\varepsilon$)

The standard $k$-$\varepsilon$ model by Launder & Spalding (1974) is the most widely used
turbulence model in industrial CFD. It computes the eddy viscosity as:

$$\nu_t = C_\mu \cdot k^{2} / \varepsilon$$

Standard coefficients:
- $C_\mu = 0.09$
- $C_1 = 1.44$
- $C_2 = 1.92$
- $\sigma_k = 1.0$
- $\sigma_\varepsilon = 1.3$

> **📌 Note:** These are exactly the coefficients used in our NACA airfoil
> project — see `projects/03_naca_airfoil_analysis/constant/turbulenceProperties`.

- **Strengths**: Robust, well-validated, fast convergence, enormous body of literature
- **Weaknesses**: Poor in adverse pressure gradients, over-predicts spreading in round jets, requires wall functions (cannot integrate to wall)
- **Best for**: Internal flows, mixing, general-purpose industrial CFD
- **OpenFOAM name**: `kEpsilon`
- **Required fields**: `k`, `epsilon`, `nut`

### 4.3 $k$-$\varepsilon$ Realizable

**Equations solved:** 2 ($k$, $\varepsilon$)

The realizable $k$-$\varepsilon$ model (Shih et al., 1995) modifies the $\varepsilon$ equation and makes
$C_\mu$ a function of the local flow conditions rather than a constant. This ensures
the "realizability" constraint (normal Reynolds stresses must be positive).

- **Strengths**: Better for rotation, boundary layers with strong adverse pressure gradients, separation, recirculation
- **Weaknesses**: Slightly less robust than standard $k$-$\varepsilon$
- **Best for**: Flows with strong streamline curvature, vortices, rotation
- **OpenFOAM name**: `realizableKE`
- **Required fields**: `k`, `epsilon`, `nut`

### 4.4 $k$-$\omega$ (Wilcox)

**Equations solved:** 2 ($k$, specific dissipation rate $\omega$)

The $k$-$\omega$ model solves for $\omega = \varepsilon/(C_\mu \cdot k)$ instead of $\varepsilon$. This reformulation provides
superior performance near solid walls.

- **Strengths**: Excellent near-wall behavior, no wall functions needed (can integrate to $y^{+} \approx 1$), good for transitional flows
- **Weaknesses**: Sensitive to freestream $\omega$ values, can give different results depending on freestream BC
- **Best for**: Boundary layer flows, low-Re turbulence
- **OpenFOAM name**: `kOmega`
- **Required fields**: `k`, `omega`, `nut`

### 4.5 $k$-$\omega$ SST (Menter)

**Equations solved:** 2 ($k$, $\omega$) with a blending function

The **Shear Stress Transport (SST)** model by Menter (1994) is widely considered
the best general-purpose RANS model. It blends $k$-$\omega$ (near walls) with $k$-$\varepsilon$
(in the free stream), getting the best of both:

```
  ┌─────────────────────────────────────────────────┐
  │                                                 │
  │   Free stream:   k-ε behavior (via ω)          │
  │         ▲        (robust, insensitive to BC)    │
  │         │                                       │
  │   Blending zone: F₁ blending function           │
  │         │        (smooth transition)            │
  │         ▼                                       │
  │   Near wall:     k-ω behavior                   │
  │                  (accurate boundary layer)      │
  │                                                 │
  │   ═══════════════════════════════ WALL           │
  └─────────────────────────────────────────────────┘
```

Additionally, the SST model includes a **shear stress limiter** that prevents
over-prediction of eddy viscosity in adverse pressure gradient flows, which
plagues both standard $k$-$\varepsilon$ and $k$-$\omega$.

- **Strengths**: Best all-around RANS model, excellent for adverse pressure gradients, separation prediction, no freestream sensitivity
- **Weaknesses**: Slightly more expensive than $k$-$\varepsilon$, still RANS limitations for massive separation
- **Best for**: Almost everything — default choice for RANS
- **OpenFOAM name**: `kOmegaSST`
- **Required fields**: `k`, `omega`, `nut`

### 4.6 Reynolds Stress Models (RSM)

**Equations solved:** 7 (six Reynolds stress components + $\varepsilon$ or $\omega$)

Instead of assuming isotropic eddy viscosity, RSM solves transport equations
for each component of the Reynolds stress tensor directly.

- **Strengths**: Captures anisotropy, swirling flows, strong curvature effects
- **Weaknesses**: Expensive, harder to converge, numerically stiff, many model constants
- **Best for**: Cyclone separators, swirl combustors, strongly anisotropic flows
- **OpenFOAM name**: `LRR` (Launder-Reece-Rodi) or `SSG`
- **Required fields**: `R`, `epsilon` (or `omega`), `nut`

### RANS Model Comparison Table

| Model | Eqns | Best For | Wall Treatment | Relative Cost | Accuracy |
|-------|------|----------|----------------|---------------|----------|
| Spalart-Allmaras | 1 | Aerospace external aero | Low-Re or wall fn | ★☆☆☆☆ | Attached BL: ★★★★☆ |
| $k$-$\varepsilon$ standard | 2 | General industrial | Wall functions | ★★☆☆☆ | General: ★★★☆☆ |
| $k$-$\varepsilon$ realizable | 2 | Rotation, separation | Wall functions | ★★☆☆☆ | Separation: ★★★★☆ |
| $k$-$\omega$ Wilcox | 2 | Near-wall flows | Resolves BL | ★★☆☆☆ | Near-wall: ★★★★☆ |
| **$k$-$\omega$ SST** | **2** | **Almost everything** | **Both options** | **★★★☆☆** | **General: ★★★★☆** |
| RSM (LRR/SSG) | 7 | Anisotropic flows | Wall functions | ★★★★☆ | Anisotropy: ★★★★★ |

> **💡 Recommendation:** If you are unsure which model to use, **start with $k$-$\omega$ SST**.
> It is the best general-purpose RANS model and the default recommendation by most
> CFD experts. Only switch to another model if you have a specific reason.

---

## 5. Wall Treatment — The Most Critical Detail

The near-wall region is where turbulence modeling matters most. Getting wall
treatment wrong is the #1 source of error in RANS simulations.

### The Boundary Layer Structure

```
  Distance from wall (y)
  │
  │                                    ╱ Outer layer
  │                                 ╱    (wake region)
  │                              ╱
  │                           ╱
  │         ╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱          Log-law layer
  │       ╱    u⁺ = (1/κ)·ln(y⁺) + B    (y⁺ = 30-300)
  │     ╱      where κ ≈ 0.41, B ≈ 5.2
  │   ╱
  │  ┃         Buffer layer               (y⁺ = 5-30)
  │  ┃         (transition — neither law works well)
  │  ┃
  │  ┃         Viscous sublayer            (y⁺ < 5)
  │  ┃         u⁺ = y⁺  (linear profile)
  │  ┃
  ══════════════════════════════════════════════ WALL  (y⁺ = 0)
```

### What Is $y^{+}$?

The dimensionless wall distance $y^{+}$ is defined as:

$$y^{+} = \frac{y \cdot u_\tau}{\nu}$$

where $u_\tau = \sqrt{\tau_w / \rho}$ is the friction velocity, $\tau_w$ = wall shear stress, $y$ = distance from wall to cell center, $\nu$ = kinematic viscosity

### $y^{+}$ Requirements by Approach

| Wall Treatment | $y^{+}$ of First Cell | Description |
|----------------|-------------------|-------------|
| **Wall functions** | $30 < y^{+} < 300$ | Bridges the viscous sublayer with empirical functions. Cheaper mesh. Used with standard $k$-$\varepsilon$. |
| **Low-Re / Resolve** | $y^{+} \approx 1$ | Resolves the entire boundary layer. More cells needed. Used with $k$-$\omega$, SST. |
| **Enhanced wall treatment** | $y^{+} \approx 1$ (ideal), tolerates higher | Blends both approaches. Most flexible. |

> **⚠️ Warning — The Buffer Layer Trap:** Never place your first cell center in the
> buffer layer ($5 < y^{+} < 30$). Neither wall functions nor direct resolution work well
> there. Either aim for $y^{+} \approx 1$ (resolve) or $y^{+} \approx 30\text{-}100$ (wall functions). The buffer
> layer is a no-man's-land.

### Estimating First Cell Height from $y^{+}$

To design your mesh, estimate the first cell height $\Delta y$:

Step 1: Estimate $Re_L = U \cdot L / \nu$

Step 2: Estimate skin friction coefficient:

$$c_f \approx 0.058 \cdot Re_L^{-0.2} \quad \text{(turbulent flat plate)}$$

Step 3: Estimate wall shear stress:

$$\tau_w = 0.5 \cdot c_f \cdot \rho \cdot U^{2}$$

Step 4: Friction velocity:

$$u_\tau = \sqrt{\tau_w / \rho}$$

Step 5: First cell height:

$$\Delta y = y^{+}_{target} \cdot \nu / u_\tau$$

**Example for our NACA airfoil project:**

$$U = 1 \text{ m/s}, \quad L = 1 \text{ m}, \quad \nu = 1\text{e-05 m}^{2}\text{/s}$$

$$Re = 1 \times 1 / 1\text{e-05} = 100{,}000$$

$$c_f \approx 0.058 \times (100{,}000)^{-0.2} = 0.058 \times 0.1 = 0.0058$$

$$\tau_w \approx 0.5 \times 0.0058 \times 1.225 \times 1^{2} \approx 0.00355 \text{ Pa}$$

$$u_\tau = \sqrt{0.00355 / 1.225} \approx 0.0538 \text{ m/s}$$

For $y^{+} = 30$ (wall functions): $\Delta y \approx 30 \times 1\text{e-05} / 0.0538 \approx 0.0056$ m

For $y^{+} = 1$ (resolve BL): $\Delta y \approx 1 \times 1\text{e-05} / 0.0538 \approx 0.000186$ m

### Wall Functions in OpenFOAM

OpenFOAM provides several wall function boundary conditions:

| Field | Wall Function BC | Low-Re / Resolve BC | Description |
|-------|-----------------|---------------------|-------------|
| `nut` | `nutkWallFunction` | `nutLowReWallFunction` | Turbulent viscosity at wall |
| `k` | `kqRWallFunction` | `fixedValue` (= 0) | Turbulent kinetic energy |
| `epsilon` | `epsilonWallFunction` | `epsilonLowReWallFunction` | Dissipation rate |
| `omega` | `omegaWallFunction` | `fixedValue` (large) | Specific dissipation |

> **📌 Note:** Our NACA airfoil project uses wall functions (`nutkWallFunction`,
> `kqRWallFunction`, `epsilonWallFunction`) — consistent with the $k$-$\varepsilon$ model
> and a mesh designed for $y^{+} \approx 30\text{-}100$.

---

## 6. Large Eddy Simulation (LES)

LES resolves the **large, energy-containing eddies** directly and models only the
small-scale (sub-grid) eddies. This provides much higher fidelity than RANS for
flows with large-scale unsteadiness.

### The Filtering Concept

```
  ┌──────────────────────────────────────────────────────────────┐
  │  Full turbulent field     →    Filtered (resolved) field     │
  │                                + Sub-grid (modeled) field    │
  │                                                              │
  │  u(x,t) = ū(x,t) + u'(x,t)                                │
  │            ▲               ▲                                 │
  │            │               │                                 │
  │         Resolved        Modeled by                           │
  │         on the grid     SGS model                            │
  │                                                              │
  │  Filter width Δ ≈ grid spacing (implicit filtering)          │
  └──────────────────────────────────────────────────────────────┘
```

Unlike RANS, the filter in LES is **spatial**, not temporal. The resolved field
is still time-dependent, so LES captures transient vortex dynamics, shedding,
and large-scale mixing.

### Common LES Sub-Grid Scale Models

| SGS Model | OpenFOAM Name | Description | Best For |
|-----------|---------------|-------------|----------|
| Smagorinsky | `Smagorinsky` | Simplest. Constant coefficient Cs ≈ 0.1-0.2. Too dissipative near walls. | Quick estimates, free shear |
| Dynamic Smagorinsky | `dynamicKEqn` | Coefficient computed dynamically from resolved field. Much better. | General LES |
| WALE | `WALE` | Wall-Adapting Local Eddy-viscosity. Correct near-wall scaling without dynamic procedure. | Wall-bounded LES |
| k-equation SGS | `kEqn` | Solves a transport equation for sub-grid k. | Complex flows |
| Sigma | `sigma` | Good near-wall behavior. Less dissipative than Smagorinsky. | Transitional flows |

### Grid Requirements for LES

LES requires **much finer grids** than RANS, especially in the wall-normal and
streamwise directions:

| Direction | RANS (wall functions) | LES |
|-----------|----------------------|-----|
| Wall-normal ($y^{+}_1$) | 30-100 | 1-2 |
| Streamwise ($\Delta x^{+}$) | ~1000 | 50-100 |
| Spanwise ($\Delta z^{+}$) | N/A (often 2D) | 15-40 |

> **⚠️ Warning:** Running LES on a RANS-quality mesh gives results that are
> **worse** than RANS — you get the cost of LES with the accuracy of neither.
> Always verify that your mesh resolves at least 80% of the turbulent kinetic
> energy. If your SGS model is doing most of the work, your "LES" is really
> just expensive RANS.

### When Is LES Worth the Cost?

- Flows dominated by **large-scale unsteadiness** (bluff body shedding, combustion)
- When **acoustic predictions** are needed
- When RANS **fails consistently** (massive separation, strong mixing)
- When you need **time-resolved** statistics, not just means

---

## 7. Direct Numerical Simulation (DNS)

DNS resolves **every** turbulent scale — from the largest eddies down to the
Kolmogorov microscale. No turbulence model is used at all.

### What DNS Resolves

```
  Energy spectrum E(κ)
  │
  │ ╲
  │   ╲
  │     ╲          DNS resolves
  │       ╲        ← ALL of this →
  │         ╲╲
  │           ╲╲╲
  │              ╲╲╲╲╲
  │                   ╲╲╲╲╲╲╲╲╲
  └──────────────────────────────── Wave number κ
    Large eddies ───────────► Small eddies
```

### Computational Cost

The number of grid points required scales as:

$$N \sim Re^{9/4} \quad \text{(in 3D)}$$

| Re | Grid Points | Feasibility |
|----|-------------|-------------|
| 100 | $\sim 10^{4}$ | Laptop |
| 1,000 | $\sim 10^{7}$ | Workstation |
| 10,000 | $\sim 10^{10}$ | Supercomputer |
| 100,000 | $\sim 10^{13}$ | Beyond current capability |
| 1,000,000 | $\sim 10^{16}$ | Science fiction |

DNS is a **research tool**, not an engineering tool. It is used to:
- Validate RANS and LES models
- Study fundamental turbulence physics
- Generate reference data for model development

> **💡 Perspective:** The highest-Re DNS ever performed (as of ~2023) reached
> $Re_\tau \approx 10,000$ in a channel flow, requiring billions of grid points and months
> on top-tier supercomputers. A typical aircraft operates at $Re \approx 10^{7}\text{-}10^{8}$.
> DNS of full aircraft is decades away.

---

## 8. Setting Up Turbulence in OpenFOAM

This section uses **actual code from our projects** to show exactly how
turbulence is configured.

### 8.1 The turbulenceProperties Dictionary

This file lives in `constant/turbulenceProperties` and controls which turbulence
model is active.

**From `projects/03_naca_airfoil_analysis/constant/turbulenceProperties`:**

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  RAS;

RAS
{
    RASModel        kEpsilon;

    turbulence      on;

    printCoeffs     on;

    kEpsilonCoeffs
    {
        Cmu             0.09;
        C1              1.44;
        C2              1.92;
        sigmak          1.0;
        sigmaEps        1.3;
    }
}
```

#### simulationType Options

| Value | Meaning | When to Use |
|-------|---------|-------------|
| `laminar` | No turbulence model | Low-Re flows (Re < ~2000). Our lid-driven cavity uses icoFoam which defaults to laminar. |
| `RAS` | Reynolds-Averaged (RANS) | Most engineering simulations. Our NACA airfoil project. |
| `LES` | Large Eddy Simulation | High-fidelity transient simulations. |

#### Available RAS Models in OpenFOAM

| OpenFOAM Name | Model | Equations |
|---------------|-------|-----------|
| `kEpsilon` | Standard $k$-$\varepsilon$ | $k$, $\varepsilon$ |
| `realizableKE` | Realizable $k$-$\varepsilon$ | $k$, $\varepsilon$ |
| `RNGkEpsilon` | RNG $k$-$\varepsilon$ | $k$, $\varepsilon$ |
| `kOmega` | Wilcox $k$-$\omega$ | $k$, $\omega$ |
| `kOmegaSST` | Menter SST | $k$, $\omega$ |
| `SpalartAllmaras` | Spalart-Allmaras | $\tilde{\nu}$ |
| `LRR` | Launder-Reece-Rodi RSM | $R_{ij}$, $\varepsilon$ |
| `SSG` | Speziale-Sarkar-Gatski RSM | $R_{ij}$, $\varepsilon$ |
| `LaunderSharmaKE` | Low-Re $k$-$\varepsilon$ | $k$, $\varepsilon$ |
| `LamBremhorstKE` | Low-Re $k$-$\varepsilon$ variant | $k$, $\varepsilon$ |
| `v2f` | $v^{2}$-$f$ model | $k$, $\varepsilon$, $v^{2}$, $f$ |
| `kOmegaSSTLM` | SST with Langtry-Menter transition | $k$, $\omega$, $\gamma$, $Re_\theta$ |

### 8.2 Turbulent Field Initialization (Boundary & Initial Conditions)

Every turbulence model requires specific field files in the `0/` directory.

**From `projects/03_naca_airfoil_analysis/0/k` (turbulent kinetic energy):**

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      k;
}

dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0.01;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 0.01;
    }

    outlet
    {
        type            zeroGradient;
    }

    walls
    {
        type            kqRWallFunction;
        value           uniform 0.01;
    }

    airfoil
    {
        type            kqRWallFunction;
        value           uniform 0.01;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

**From `projects/03_naca_airfoil_analysis/0/epsilon` (dissipation rate):**

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      epsilon;
}

dimensions      [0 2 -3 0 0 0 0];
internalField   uniform 0.01;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 0.01;
    }

    outlet
    {
        type            zeroGradient;
    }

    walls
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }

    airfoil
    {
        type            epsilonWallFunction;
        value           uniform 0.01;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

**From `projects/03_naca_airfoil_analysis/0/nut` (turbulent viscosity):**

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      nut;
}

dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            calculated;
        value           uniform 0;
    }

    outlet
    {
        type            calculated;
        value           uniform 0;
    }

    walls
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    airfoil
    {
        type            nutkWallFunction;
        value           uniform 0;
    }

    frontAndBack
    {
        type            empty;
    }
}
```

### 8.3 Transport Properties

The kinematic viscosity $\nu$ is critical for determining the Reynolds number
and must be consistent with your turbulence model assumptions.

**From `projects/03_naca_airfoil_analysis/constant/transportProperties`:**

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

transportModel  Newtonian;

nu              nu [ 0 2 -1 0 0 0 0 ] 1e-05;
```

Compare with the **laminar** lid-driven cavity case:

**From `projects/01_lid_driven_cavity/constant/transportProperties`:**

```c
transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 0.001;  // Higher viscosity → lower Re → laminar
```

> **💡 Key difference:** The NACA airfoil uses $\nu = 1\text{e-05}$ (air-like), giving $Re \approx 100,000$
> (turbulent). The cavity uses $\nu = 0.001$, giving $Re \approx 1,000$ (laminar). The viscosity
> choice drives the entire turbulence modeling strategy.

### 8.4 Solver and Numerical Scheme Requirements

Turbulence models require matching solver settings in `system/fvSolution` and
discretization schemes in `system/fvSchemes`.

**From `projects/03_naca_airfoil_analysis/system/fvSchemes` (div schemes for $k$-$\varepsilon$):**

```c
divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);
    div((nuEff*dev(T(grad(U))))) Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div(phi,epsilon) Gauss upwind;
    div(phi,k) Gauss upwind;
}
```

**From `projects/03_naca_airfoil_analysis/system/fvSolution` (turbulence solvers & relaxation):**

```c
solvers
{
    k
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

    epsilon
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent yes;
    residualControl
    {
        p       1e-3;
        U       1e-4;
        "(k|epsilon)" 1e-4;
    }
}

relaxationFactors
{
    equations
    {
        U       0.7;
        k       0.7;
        epsilon 0.7;
    }
}
```

### 8.5 Calculating Initial Values of $k$ and $\varepsilon$

The inlet values of $k$ and $\varepsilon$ should not be arbitrary. Calculate them from
the **turbulence intensity** ($I$) and a **length scale** ($l$):

$$k = 1.5 \cdot (U \cdot I)^{2}$$

$$\varepsilon = C_\mu^{3/4} \cdot k^{3/2} / l$$

where:
- $I$ = turbulence intensity (typically 0.01 to 0.10 for external flows)
- $l$ = turbulence length scale (often $0.07 \times$ hydraulic diameter)
- $C_\mu = 0.09$

**Example for our NACA airfoil ($U = 1$ m/s, $I = 5\%$, $l = 0.07$ m):**

$$k = 1.5 \times (1.0 \times 0.05)^{2} = 1.5 \times 0.0025 = 0.00375 \text{ m}^{2}\text{/s}^{2}$$

$$\varepsilon = 0.09^{0.75} \times 0.00375^{1.5} / 0.07 = 0.0302 \times 0.000230 / 0.07 \approx 0.0001 \text{ m}^{2}\text{/s}^{3}$$

> **⚠️ Caution:** The values `k = 0.01` and `ε = 0.01` used in our project files
> are reasonable initial estimates but not precisely calculated from turbulence
> intensity. For production simulations, always calculate proper values.

---

## 9. Model Selection Flowchart

Use this decision tree to choose the right turbulence modeling approach:

```
  Is the flow turbulent? (Re > ~2000-4000)
  │
  ├─ NO ──→  Use LAMINAR (simulationType laminar)
  │          Example: lid-driven cavity at Re=1000 → icoFoam
  │
  └─ YES ─→  Do you need time-resolved details?
              │
              ├─ NO ──→  Use RANS (simulationType RAS)
              │          │
              │          ├─ Is the geometry simple, attached flow?
              │          │   └─ YES → k-ε standard or SA
              │          │
              │          ├─ Flow separation or adverse pressure gradient?
              │          │   └─ YES → k-ω SST  ← (BEST DEFAULT CHOICE)
              │          │
              │          ├─ Strong rotation or swirl?
              │          │   └─ YES → k-ε realizable or RSM
              │          │
              │          └─ Strongly anisotropic turbulence?
              │              └─ YES → RSM (LRR or SSG)
              │
              └─ YES ─→  Can you afford fine mesh + many time steps?
                          │
                          ├─ NO ──→  Use DES / Hybrid RANS-LES
                          │          (RANS near walls, LES away)
                          │
                          └─ YES ─→  Is the Re moderate (< ~50,000)?
                                      │
                                      ├─ YES → LES (Smagorinsky, WALE)
                                      │
                                      └─ NO ──→  LES (very expensive)
                                                  Consider DES instead
```

> **💡 The golden rule:** When in doubt, use **$k$-$\omega$ SST** for RANS. It works
> well for the widest range of flows and is the default recommendation in
> most OpenFOAM tutorials and industrial practice.

---

## 10. Multiphase Flow and Turbulence

Multiphase flow modeling adds another layer of complexity to turbulence. When
multiple fluid phases (e.g., water and air) interact, the interface dynamics
couple strongly with turbulence.

### Volume-of-Fluid (VOF) Method

The VOF method tracks the interface between immiscible fluids using a scalar
volume fraction field $\alpha$:

```
  ┌──────────────────────────────────────────────────────┐
  │                                                      │
  │   α = 0          α = 0          α = 0                │
  │   (pure air)     (pure air)     (pure air)           │
  │                                                      │
  │ ─ ─ ─ ╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╲╲╲╲╲╲╲╲╲ ─ ─ ─ ─ ─       │
  │       ╱  0<α<1  interface  0<α<1  ╲                  │
  │ ─ ─ ╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╱╲╲ ─ ─ ─ ─       │
  │                                                      │
  │   α = 1          α = 1          α = 1                │
  │   (pure water)   (pure water)   (pure water)         │
  │                                                      │
  └──────────────────────────────────────────────────────┘
  
  Fluid properties are interpolated:
    ρ = α·ρ_water + (1-α)·ρ_air
    μ = α·μ_water + (1-α)·μ_air
```

The VOF transport equation is:

$$\frac{\partial \alpha}{\partial t} + \nabla \cdot (\alpha \cdot U) = 0$$

### interFoam Solver

OpenFOAM's `interFoam` is a VOF-based solver for two incompressible, isothermal
immiscible fluids. It uses the MULES (Multidimensional Universal Limiter for
Explicit Solution) algorithm for interface sharpness.

**Example `system/controlDict` for interFoam:**

```c
application     interFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         10;

deltaT          0.001;
writeControl    timeStep;
writeInterval   100;

purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
```

**Example `constant/transportProperties` for multiphase:**

```c
phases
(
    water
    {
        transportModel  Newtonian;
        nu              1e-06;
        rho             1000;
    }
    air
    {
        transportModel  Newtonian;
        nu              1.5e-05;
        rho             1;
    }
);

sigma           0.07;    // Surface tension [N/m]
```

### Turbulence in Multiphase Flows

Turbulence modeling in multiphase flows requires special care:

- The turbulence model applies to the **mixture**, not individual phases
- Interface dynamics can generate or suppress turbulence
- Most RANS models work with interFoam ($k$-$\varepsilon$, $k$-$\omega$ SST)
- LES of multiphase flows is an active research area
- Surface tension and buoyancy effects interact with turbulence

---

## 11. Practical Guidelines and Common Mistakes

### The Top 10 Turbulence Modeling Mistakes

| # | Mistake | Consequence | Fix |
|---|---------|-------------|-----|
| 1 | Wrong $y^{+}$ for wall treatment | Wildly incorrect wall shear stress | Check $y^{+}$ after meshing, before trusting results |
| 2 | $y^{+}$ in buffer layer (5-30) | Neither wall functions nor resolve works | Redesign mesh for $y^{+} \approx 1$ or $y^{+} \approx 30\text{-}100$ |
| 3 | Arbitrary $k$, $\varepsilon$ inlet values | Wrong turbulence levels propagate | Calculate from turbulence intensity and length scale |
| 4 | 2D LES | LES requires 3D — turbulence is inherently 3D | Always run LES in 3D with spanwise extent |
| 5 | LES on RANS mesh | Resolves neither large nor small eddies | Generate LES-quality mesh first |
| 6 | Steady-state for massively separated flow | RANS cannot capture unsteady shedding | Use transient RANS (URANS) or DES/LES |
| 7 | Ignoring convergence of turbulence eqns | $k$ and $\varepsilon$ not converged → wrong $\nu_t$ | Monitor residuals for ALL equations |
| 8 | Using $k$-$\varepsilon$ for adverse pressure gradient | Over-predicts attached flow, misses separation | Use $k$-$\omega$ SST |
| 9 | Wrong relaxation factors | Divergence or extremely slow convergence | Start with 0.3-0.7, reduce if diverging |
| 10 | Not checking turbulence ratio $\nu_t/\nu$ | Unrealistically high values indicate problems | Should be $O(10^{1})$-$O(10^{3})$ for most flows |

### Convergence Tips for Turbulence Simulations

```
  Relaxation factors (SIMPLE/SIMPLEC):
  ┌─────────────────────────────────────┐
  │  Field     │ Aggressive │  Safe     │
  ├─────────────────────────────────────┤
  │  p         │ 0.3        │  0.1-0.2  │
  │  U         │ 0.7        │  0.5      │
  │  k         │ 0.7        │  0.3-0.5  │
  │  epsilon   │ 0.7        │  0.3-0.5  │
  │  omega     │ 0.7        │  0.3-0.5  │
  └─────────────────────────────────────┘
```

> **💡 Tip:** If your simulation diverges, try these steps in order:
> 1. Reduce relaxation factors (especially for $k$/$\varepsilon$/$\omega$)
> 2. Switch turbulence convection to `upwind` (first-order) temporarily
> 3. Initialize with a potential flow solution (`potentialFoam`)
> 4. Start with a coarser mesh, then map the solution to a finer mesh

### Checking Your Results

Always verify these quantities after a turbulence simulation:

1. **$y^{+}$ distribution**: Run `yPlus` utility — values should be appropriate for your wall treatment
2. **Residuals**: All equations should be converged (typically $< 10^{-4}$)
3. **Turbulence ratio** ($\nu_t/\nu$): Should be reasonable (not $> 10^{5}$)
4. **Mass balance**: Check that mass is conserved
5. **Grid independence**: Run at least 2-3 mesh levels to check convergence

---

## 12. Quick Reference Tables

### OpenFOAM Solver Selection for Turbulence

| Solver | Flow Type | Turbulence | Phase |
|--------|-----------|------------|-------|
| `icoFoam` | Transient, incompressible | Laminar only | Single |
| `pisoFoam` | Transient, incompressible | RANS or LES | Single |
| `simpleFoam` | Steady-state, incompressible | RANS | Single |
| `pimpleFoam` | Transient, incompressible | RANS or LES | Single |
| `interFoam` | Transient, incompressible | RANS or LES | Two-phase VOF |
| `rhoSimpleFoam` | Steady-state, compressible | RANS | Single |
| `rhoPimpleFoam` | Transient, compressible | RANS or LES | Single |

### Complete Field File Checklist by Model

| Model | k | epsilon | omega | nut | nuTilda | R |
|-------|---|---------|-------|-----|---------|---|
| kEpsilon | ✅ | ✅ | ❌ | ✅ | ❌ | ❌ |
| realizableKE | ✅ | ✅ | ❌ | ✅ | ❌ | ❌ |
| kOmega | ✅ | ❌ | ✅ | ✅ | ❌ | ❌ |
| kOmegaSST | ✅ | ❌ | ✅ | ✅ | ❌ | ❌ |
| SpalartAllmaras | ❌ | ❌ | ❌ | ✅ | ✅ | ❌ |
| LRR / SSG | ❌ | ✅ | ❌ | ✅ | ❌ | ✅ |

### Dimensions Reference

| Field | Dimensions `[kg m s K mol A cd]` | Unit |
|-------|----------------------------------|------|
| k | `[0 2 -2 0 0 0 0]` | $\text{m}^{2}\text{/s}^{2}$ |
| epsilon | `[0 2 -3 0 0 0 0]` | $\text{m}^{2}\text{/s}^{3}$ |
| omega | `[0 0 -1 0 0 0 0]` | 1/s |
| nut | `[0 2 -1 0 0 0 0]` | $\text{m}^{2}\text{/s}$ |
| nuTilda | `[0 2 -1 0 0 0 0]` | $\text{m}^{2}\text{/s}$ |

### Switching from $k$-$\varepsilon$ to $k$-$\omega$ SST (Migration Checklist)

If you want to switch from $k$-$\varepsilon$ (as used in our NACA project) to $k$-$\omega$ SST:

```
  1. constant/turbulenceProperties:
     - Change RASModel from "kEpsilon" to "kOmegaSST"

  2. Delete:  0/epsilon
     Create:  0/omega    (dimensions [0 0 -1 0 0 0 0])

  3. Update wall BCs:
     - nut:     nutkWallFunction → nutUSpaldingWallFunction (optional, works either way)
     - k:       kqRWallFunction  → kqRWallFunction (no change needed)
     - omega:   new file with omegaWallFunction at walls

  4. system/fvSchemes:
     - Change div(phi,epsilon) to div(phi,omega)

  5. system/fvSolution:
     - Change epsilon solver to omega solver
     - Update residualControl to use "(k|omega)" instead of "(k|epsilon)"

  6. Refine mesh near walls:
     - k-ω SST works best with y⁺ ≈ 1 (can also use wall functions)
```

---

## 13. Summary

```
  ╔═══════════════════════════════════════════════════════════════════╗
  ║                  TURBULENCE MODELING CHEAT SHEET                 ║
  ╠═══════════════════════════════════════════════════════════════════╣
  ║                                                                   ║
  ║  Re < 2000     →  Laminar (no model)                             ║
  ║  Re > 4000     →  Need turbulence model                          ║
  ║                                                                   ║
  ║  DEFAULT CHOICE: k-ω SST  (kOmegaSST)                           ║
  ║                                                                   ║
  ║  Attached BL, aerospace  →  Spalart-Allmaras                     ║
  ║  General industrial      →  k-ε standard                         ║
  ║  Separation, pressure gr →  k-ω SST  ★                          ║
  ║  Strong swirl / rotation →  RSM (LRR/SSG)                       ║
  ║  Transient, high detail  →  LES (WALE or dynamic)               ║
  ║  Fundamental research    →  DNS (if Re allows)                   ║
  ║                                                                   ║
  ║  WALL TREATMENT:                                                  ║
  ║    Wall functions  →  y⁺ ≈ 30-100  (cheaper mesh)               ║
  ║    Resolve BL      →  y⁺ ≈ 1      (more accurate)               ║
  ║    NEVER y⁺ = 5-30 (buffer layer)                                ║
  ║                                                                   ║
  ║  INITIAL CONDITIONS:                                              ║
  ║    k = 1.5·(U·I)²           I = turbulence intensity             ║
  ║    ε = Cμ^0.75·k^1.5/l      l = turbulence length scale         ║
  ║                                                                   ║
  ╚═══════════════════════════════════════════════════════════════════╝
```

---

*This document is part of the [OpenFOAM Tutorials](../README.md) series.
See `projects/03_naca_airfoil_analysis/` for a complete working example of
turbulent flow simulation using the $k$-$\varepsilon$ model with simpleFoam.*
