# 05 — Heated Pipe Flow (Conjugate Heat Transfer)

| Detail | Value |
|---|---|
| **Difficulty** | ⭐⭐ Intermediate |
| **Solver** | `buoyantSimpleFoam` (steady-state, buoyant, turbulent, compressible with energy equation) |
| **Algorithm** | SIMPLE |
| **Turbulence Model** | k-epsilon (RAS) |
| **Reynolds Number** | Re $\approx$ 3333 (based on hydraulic diameter) |
| **Prandtl Number** | Pr = 0.71 (air) |
| **OpenFOAM Version** | 6+ |

---

## Table of Contents

1. [Problem Description](#1-problem-description)
2. [Physics and Theory](#2-physics-and-theory)
3. [Case Structure](#3-case-structure)
4. [Mesh Details](#4-mesh-details)
5. [Boundary Conditions](#5-boundary-conditions)
6. [Solver Configuration](#6-solver-configuration)
7. [How to Run](#7-how-to-run)
8. [Expected Results](#8-expected-results)
9. [Post-Processing Guide](#9-post-processing-guide)
10. [Exercises](#10-exercises)
11. [Troubleshooting](#11-troubleshooting)
12. [References and Further Reading](#12-references-and-further-reading)

---

## 1. Problem Description

This simulation models **forced convection heat transfer** inside a heated rectangular duct.
Cold air enters at 300 K and is heated by the hot walls (350 K) as it flows through the
2-meter-long channel. This is a classic problem in thermal-fluid engineering, illustrating
how thermal boundary layers develop and how heat is transferred from a wall to a moving
fluid stream.

### Why This Problem Matters

Heat exchangers, cooling channels, HVAC ducts, and electronics cooling systems all involve
fluid flowing past heated or cooled surfaces. Understanding how temperature profiles
develop, how heat flux varies along the length, and how turbulence affects mixing is
fundamental to designing efficient thermal systems.

### Geometry and Setup

```
                      T_wall = 350 K (hot wall)
  ┌─────────────────────────────────────────────────────────────────────────────┐
  │░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░│
  ├─────────────────────────────────────────────────────────────────────────────┤
  │                                                                             │
  │   T_in = 300 K  →→→→→→→→  heated fluid  →→→→→→→→  T_out > 300 K  →→→     │
  │   (cold air)               (warming up)              (heated air)           │
  │                                                                             │
  │   U_in = 0.5 m/s           Thermal BL                Fully/partially       │
  │   (uniform)                 develops                  developed             │
  │                                                                             │
  ├─────────────────────────────────────────────────────────────────────────────┤
  │░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░│
  └─────────────────────────────────────────────────────────────────────────────┘
                      T_wall = 350 K (hot wall)

  ◄────────────────────── L = 2.0 m ──────────────────────►

                         H = 0.1 m (duct height)
```

### Temperature Profile Development Along the Duct

```
  Inlet (x=0)          Mid-duct (x=1m)         Outlet (x=2m)
  ─────────────         ─────────────           ─────────────
  │           │         │    ╱──╲    │          │   ╱────╲   │
  │           │         │  ╱    ╲  │          │  ╱      ╲  │
  │  uniform  │   →→→   │╱  warm  ╲│    →→→   │╱  warmer ╲│
  │  T = 300K │         │╲  core  ╱│          │╲  core   ╱│
  │           │         │  ╲    ╱  │          │  ╲      ╱  │
  │           │         │    ╲──╱    │          │   ╲────╱   │
  ─────────────         ─────────────           ─────────────
  T_wall=350K           T_wall=350K             T_wall=350K

  Thin thermal BL       Growing thermal BL      Thick thermal BL
  High wall heat flux    Moderate heat flux      Lower heat flux
```

---

## 2. Physics and Theory

### 2.1 Governing Equations

This problem solves the **Reynolds-Averaged Navier-Stokes (RANS) equations** coupled with
the **energy equation** for a compressible (buoyant) flow:

**Continuity (mass conservation):**

$$\partial \rho / \partial t + \nabla \cdot (\rho U) = 0$$

**Momentum (Navier-Stokes):**

$$\partial (\rho U) / \partial t + \nabla \cdot (\rho U U) = -\nabla p + \nabla \cdot (\mu_{eff} \nabla U) + \rho g$$

**Energy (enthalpy form):**

$$\partial (\rho h) / \partial t + \nabla \cdot (\rho U h) = \nabla \cdot (\alpha_{eff} \nabla h) + \partial p / \partial t$$

where:
- `ρ` = density (calculated from ideal gas law: $\rho$ = p/(RT))
- `U` = velocity vector
- `p` = pressure
- `h` = specific enthalpy
- `μ_eff` = effective (molecular + turbulent) viscosity
- `α_eff` = effective thermal diffusivity
- `g` = gravitational acceleration

### 2.2 Key Dimensionless Numbers

#### Reynolds Number (Re)

The Reynolds number characterizes the ratio of inertial to viscous forces:

$$Re = (U \times D_h) / \nu$$

For this case:
- U = 0.5 m/s (inlet velocity)
- $D_h = 2 \times H \times W/(H+W) \approx 2 \times 0.1 \times 0.01/(0.1+0.01) \approx 0.018$ m (hydraulic diameter for 2D duct $\approx 2H = 0.2$ m)
- $\nu = 1.5 \times 10^{-5}$ m$^2$/s (kinematic viscosity of air at 300K)

For the 2D approximation (infinite width): $D_h \approx 2H = 0.2$ m

$$Re = (0.5 \times 0.2) / (1.5 \times 10^{-5}) \approx 6667$$

This is in the **transitional-to-turbulent** regime ($Re > 2300$ for pipe flow), justifying
the use of a turbulence model.

#### Prandtl Number (Pr)

The Prandtl number relates momentum diffusivity to thermal diffusivity:

$$Pr = (\mu \times Cp) / k_{thermal} = \nu / \alpha$$

For air: **Pr $\approx$ 0.71**

This means the thermal boundary layer is slightly thicker than the velocity boundary
layer — characteristic of gases.

#### Nusselt Number (Nu)

The Nusselt number quantifies the enhancement of heat transfer due to convection over
pure conduction:

$$Nu = (h_{conv} \times D_h) / k_{thermal}$$

where:
- `h_conv` = convective heat transfer coefficient [W/(m$^2$·K)]
- `k_thermal` = fluid thermal conductivity [W/(m·K)]

**Expected correlations for fully-developed turbulent pipe flow:**

Dittus-Boelter correlation (heating):
$$Nu = 0.023 \times Re^{0.8} \times Pr^{0.4}$$

For our case: $Nu \approx 0.023 \times 6667^{0.8} \times 0.71^{0.4} \approx$ **24.8**

Gnielinski correlation (more accurate for transitional flow):
$$Nu = (f/8)(Re - 1000)Pr / [1 + 12.7(f/8)^{0.5} (Pr^{2/3} - 1)]$$

where `f` is the Darcy friction factor from the Moody chart.

### 2.3 Thermal Boundary Layer

When a fluid at uniform temperature enters a heated duct, a **thermal boundary layer**
develops at the walls:

1. **Entrance region (thermally developing):** The thermal boundary layer grows from
   the wall inward. Heat flux is highest here because the temperature gradient at the
   wall is steepest.

2. **Fully-developed region:** The thermal boundary layer fills the entire duct cross-
   section. The shape of the dimensionless temperature profile no longer changes with
   axial distance. The Nusselt number reaches a constant value.

The **thermal entry length** is approximately:

$$L_{thermal} \approx 0.05 \times Re \times Pr \times D_h$$

For our case: L_thermal $\approx$ 0.05 $\times$ 6667 $\times$ 0.71 $\times$ 0.2 $\approx$ **47 m**

Since our duct is only 2 m long, we are in the **thermally developing** region, where
the local Nusselt number varies along the duct length and is higher near the inlet.

### 2.4 Forced vs Natural Convection

This simulation includes gravity and buoyancy effects (via `buoyantSimpleFoam`), but
the dominant heat transfer mechanism is **forced convection**. The relative importance is
characterized by the Richardson number:

$$Ri = Gr / Re^2 = (g \times \beta \times \Delta T \times D_h^3) / (\nu^2 \times Re^2)$$

- If Ri << 1: forced convection dominates (our case)
- If Ri >> 1: natural convection dominates
- If Ri $\approx$ 1: mixed convection

### 2.5 Why buoyantSimpleFoam?

`buoyantSimpleFoam` is chosen because:

1. **Energy equation:** Unlike `simpleFoam`, it solves for temperature/enthalpy, which
   is essential for heat transfer problems.
2. **Buoyancy coupling:** Temperature differences cause density variations, which drive
   buoyant flows. Even when forced convection dominates, including buoyancy makes the
   simulation more physical.
3. **Steady-state SIMPLE:** Appropriate for steady forced convection problems where we
   seek the converged steady-state solution.
4. **Turbulence:** Supports RANS models (k-epsilon) for turbulent flow.
5. **Compressible formulation:** Uses the ideal gas law to couple density to temperature,
   which is more physical than assuming constant density when temperature varies.

---

## 3. Case Structure

```
05_heated_pipe_flow/
├── README.md                              ← This documentation file
├── Allrun                                 ← Script to run the entire simulation
├── Allclean                               ← Script to clean generated files
│
├── 0/                                     ← Initial & boundary conditions (t = 0)
│   ├── U                                  ← Velocity field [m/s]
│   ├── p                                  ← Total pressure [Pa]
│   ├── p_rgh                              ← Pressure minus hydrostatic [Pa]
│   ├── T                                  ← Temperature field [K]
│   ├── k                                  ← Turbulent kinetic energy [m²/s²]
│   ├── epsilon                            ← Turbulent dissipation rate [m²/s³]
│   ├── nut                                ← Turbulent kinematic viscosity [m²/s]
│   └── alphat                             ← Turbulent thermal diffusivity [kg/(m·s)]
│
├── constant/                              ← Time-invariant properties
│   ├── g                                  ← Gravitational acceleration
│   ├── thermophysicalProperties           ← Thermal properties (air as ideal gas)
│   ├── transportProperties                ← Kinematic transport properties
│   ├── turbulenceProperties               ← Turbulence model selection (k-epsilon)
│   └── polyMesh/                          ← Mesh (generated by blockMesh)
│
└── system/                                ← Solver & numerical settings
    ├── blockMeshDict                      ← Mesh generation dictionary
    ├── controlDict                        ← Simulation control parameters
    ├── fvSchemes                           ← Discretization schemes
    └── fvSolution                          ← Linear solvers & algorithm settings
```

### File Descriptions

| File | Purpose | Key Settings |
|---|---|---|
| `0/U` | Velocity boundary conditions | Inlet: 0.5 m/s, Walls: no-slip |
| `0/p_rgh` | Pressure (minus hydrostatic) | Used by buoyant solver instead of `p` |
| `0/p` | Total pressure | Calculated from `p_rgh` + $\rho$gh |
| `0/T` | Temperature field | Inlet: 300K, Walls: 350K |
| `0/k` | Turbulent kinetic energy | Wall functions at walls |
| `0/epsilon` | Turbulent dissipation | Wall functions at walls |
| `0/nut` | Turbulent viscosity | Wall functions at walls |
| `0/alphat` | Turbulent thermal diffusivity | alphatWallFunction at walls |
| `constant/g` | Gravity vector | (0, -9.81, 0) m/s$^2$ |
| `constant/thermophysicalProperties` | Gas properties | Air: Cp=1005, Pr=0.71 |
| `constant/transportProperties` | Kinematic viscosity | $\nu$ = 1.5$\times$10$^{-5}$ m$^2$/s |
| `constant/turbulenceProperties` | Turbulence model | k-epsilon (RAS) |
| `system/blockMeshDict` | Mesh definition | 100$\times$20$\times$1 cells, 2m$\times$0.1m$\times$0.01m |
| `system/controlDict` | Run control | 1000 iterations, write every 100 |
| `system/fvSchemes` | Discretization | Upwind for convection, linear for diffusion |
| `system/fvSolution` | Solvers & relaxation | GAMG for pressure, smoothSolver for rest |

---

## 4. Mesh Details

### Domain Geometry

The domain is a **2D rectangular duct** (approximated by a thin 3D slab with `empty`
front/back boundaries):

```
  ┌──────────────────────────────────────────────────────────┐
  │                                                          │ y = 0.1 m
  │              100 cells (x-direction)                     │
  │         ←──────────────────────────→                     │
  │                                                          │ 20 cells
  │                                                          │ (y-direction)
  └──────────────────────────────────────────────────────────┘ y = 0
  x = 0                                                x = 2 m

  z-direction: 1 cell (0.01 m thick) — empty boundaries
```

| Parameter | Value |
|---|---|
| Length (x) | 2.0 m |
| Height (y) | 0.1 m |
| Depth (z) | 0.01 m (single cell, 2D approximation) |
| Cells in x | 100 |
| Cells in y | 20 |
| Cells in z | 1 |
| **Total cells** | **2000** |
| Cell aspect ratio | 1:1 ($\Delta x$ = 0.02 m, $\Delta y$ = 0.005 m → AR = 4:1) |
| Grading | Uniform (1:1:1) |

### Mesh Quality Notes

- The mesh is intentionally coarse for fast turnaround and learning purposes.
- **No grading** is applied — for production simulations, you would want wall-normal
  clustering (grading in y-direction) to resolve the boundary layer properly.
- The first cell height of $\Delta y$ = 0.005 m gives an approximate y+ value suitable for
  wall functions (y+ $\approx$ 30–100 range).
- For higher accuracy, increase y-direction cells to 40–60 and add grading toward walls.

---

## 5. Boundary Conditions

### 5.1 Patch Summary

| Patch | Type | Location | Description |
|---|---|---|---|
| `inlet` | `patch` | x = 0 plane | Left face — air enters here |
| `outlet` | `patch` | x = 2 m plane | Right face — air exits here |
| `topWall` | `wall` | y = 0.1 m plane | Upper heated wall |
| `bottomWall` | `wall` | y = 0 plane | Lower heated wall |
| `frontAndBack` | `empty` | z = 0, z = 0.01 m | 2D approximation (no flux in z) |

### 5.2 Boundary Conditions Table (All 8 Variables)

| Variable | `inlet` | `outlet` | `topWall` | `bottomWall` | `frontAndBack` |
|---|---|---|---|---|---|
| **U** | fixedValue (0.5 0 0) | zeroGradient | noSlip | noSlip | empty |
| **p_rgh** | zeroGradient | fixedValue 0 | fixedFluxPressure | fixedFluxPressure | empty |
| **p** | calculated | calculated | calculated | calculated | empty |
| **T** | fixedValue 300 | inletOutlet 300 | fixedValue 350 | fixedValue 350 | empty |
| **k** | fixedValue 0.1 | zeroGradient | kqRWallFunction 0.1 | kqRWallFunction 0.1 | empty |
| **epsilon** | fixedValue 0.01 | zeroGradient | epsilonWallFunction 0.01 | epsilonWallFunction 0.01 | empty |
| **nut** | calculated 0 | calculated 0 | nutkWallFunction 0 | nutkWallFunction 0 | empty |
| **alphat** | calculated 0 | calculated 0 | alphatWallFunction 0 | alphatWallFunction 0 | empty |

### 5.3 Boundary Condition Explanations

#### Velocity (U)
- **Inlet:** Uniform velocity of 0.5 m/s in the x-direction. This represents a plug
  flow profile entering the duct.
- **Outlet:** `zeroGradient` — the velocity profile is extrapolated from the interior,
  appropriate for fully or partially developed flow exiting the domain.
- **Walls:** `noSlip` — the fluid velocity is zero at the wall surface, enforcing the
  viscous no-slip condition.

#### Pressure (p_rgh)
- `p_rgh` is the "reduced pressure," defined as `p_rgh = p - ρgh`. This formulation
  removes the hydrostatic pressure from the solution, improving numerical stability for
  buoyant flows.
- **Inlet:** `zeroGradient` — pressure adjusts freely to accommodate the imposed velocity.
- **Outlet:** `fixedValue 0` — sets the reference pressure at the outlet.
- **Walls:** `fixedFluxPressure` — adjusts the pressure gradient at walls to satisfy the
  momentum equation (accounts for gravity body force).

#### Temperature (T)
- **Inlet:** `fixedValue 300 K` — cold air enters the duct.
- **Outlet:** `inletOutlet` — behaves as `zeroGradient` for outflow but switches to
  `fixedValue 300 K` if reverse flow occurs (numerical safety).
- **Walls:** `fixedValue 350 K` — the walls are held at a constant elevated temperature,
  representing a heated surface (e.g., electrically heated or in contact with a hot medium).

#### Turbulence Variables (k, epsilon, nut, alphat)
- **Inlet:** Fixed turbulence levels representing low-intensity incoming turbulence.
  The turbulence intensity is approximately I = sqrt(2k/3) / U $\approx$ sqrt(2$\times$0.1/3) / 0.5 $\approx$ 52%.
  For lower turbulence intensity, reduce k accordingly.
- **Outlet:** `zeroGradient` — turbulence is convected out of the domain.
- **Walls:** Wall functions bridge the viscous sublayer, allowing coarser meshes:
  - `kqRWallFunction` — provides the near-wall treatment for k
  - `epsilonWallFunction` — calculates epsilon at the wall-adjacent cell
  - `nutkWallFunction` — calculates turbulent viscosity from k at near-wall cells
  - `compressible::alphatWallFunction` — wall function for turbulent thermal diffusivity
    with turbulent Prandtl number Prt = 0.85

> **Further reading on boundary conditions:**
> - [Boundary Conditions Guide](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md)
> - [Turbulence Models](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md)

---

## 6. Solver Configuration

### 6.1 Discretization Schemes (fvSchemes)

| Category | Scheme | Rationale |
|---|---|---|
| `ddtSchemes` | `steadyState` | No time derivatives — we seek the steady-state solution |
| `gradSchemes` | `Gauss linear` | Second-order gradient evaluation using cell-face interpolation |
| `div(phi,U)` | `bounded Gauss linearUpwind` | Second-order upwind for momentum — good accuracy with stability |
| `div(phi,h)` | `bounded Gauss upwind` | First-order upwind for enthalpy — robust for energy equation |
| `div(phi,k)` | `bounded Gauss upwind` | First-order upwind for turbulence — robustness over accuracy |
| `div(phi,epsilon)` | `bounded Gauss upwind` | First-order upwind for dissipation |
| `laplacianSchemes` | `Gauss linear corrected` | Second-order diffusion with non-orthogonal correction |
| `snGradSchemes` | `corrected` | Non-orthogonal correction for surface-normal gradients |

The `bounded` keyword ensures boundedness of the convection schemes, preventing
non-physical oscillations — critical for stability in compressible SIMPLE solvers.

### 6.2 Linear Solvers and Algorithm (fvSolution)

**Pressure solver (p_rgh):**
- **GAMG** (Geometric-Algebraic Multi-Grid) with Gauss-Seidel smoother
- Tolerance: $1 \times 10^{-6}$, relative tolerance: 0.01
- GAMG is optimal for pressure equations due to their elliptic nature

**Momentum, energy, turbulence (U, h, k, epsilon):**
- **smoothSolver** with symmetric Gauss-Seidel (symGaussSeidel) smoother
- Tolerance: $1 \times 10^{-6}$, relative tolerance: 0.1
- Efficient for hyperbolic/parabolic transport equations

**SIMPLE algorithm settings:**
- No non-orthogonal corrections needed (orthogonal mesh)
- Relaxation factors balance convergence speed vs. stability:
  - p_rgh: 0.3 (conservative — pressure is sensitive)
  - U: 0.7, h: 0.7, k: 0.7, epsilon: 0.7 (moderate)
  - rho: 1.0 (no under-relaxation for density)

**Convergence criteria:**
- All residuals must fall below $1 \times 10^{-4}$

> **Further reading on linear solvers:**
> - [Linear Solvers and Convergence](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/09_linear_solvers.md)

### 6.3 Run Control (controlDict)

| Parameter | Value | Description |
|---|---|---|
| `application` | `buoyantSimpleFoam` | Steady-state buoyant flow solver with energy |
| `endTime` | 1000 | Maximum number of SIMPLE iterations |
| `deltaT` | 1 | Pseudo-time step (iteration counter for steady-state) |
| `writeInterval` | 100 | Save results every 100 iterations |
| `purgeWrite` | 3 | Keep only the last 3 saved time directories |
| `writeFormat` | ascii | Human-readable output |

The `fieldAverage` function object computes time-averaged fields for `U` and `T`, which
can be useful for analyzing the converged flow patterns.

---

## 7. How to Run

### Prerequisites

- OpenFOAM 6 (or compatible version) installed and sourced
- ParaView (optional, for visualization)

### Step-by-Step

#### Option A: Using the Allrun Script

```bash
# Navigate to the case directory
cd projects/05_heated_pipe_flow

# Make the scripts executable (if not already)
chmod +x Allrun Allclean

# Run the entire simulation
./Allrun
```

#### Option B: Manual Execution

```bash
# Navigate to the case directory
cd projects/05_heated_pipe_flow

# Step 1: Generate the mesh
blockMesh

# Step 2: Verify the mesh (optional but recommended)
checkMesh

# Step 3: Run the solver
buoyantSimpleFoam

# Step 4: Visualize results
paraFoam
```

#### Cleaning the Case

To remove all generated data and start fresh:

```bash
./Allclean
```

### Monitoring Convergence

While the solver is running, you can monitor residuals in real-time:

```bash
# In a separate terminal, plot residuals (requires gnuplot)
foamMonitor -l postProcessing/residuals/0/residuals.dat

# Or simply watch the log file
tail -f log.buoyantSimpleFoam
```

You should see residuals for `Ux`, `Uy`, `h`, `k`, `epsilon`, and `p_rgh` decreasing
over iterations. The simulation converges when all residuals fall below the target
values ($1 \times 10^{-4}$) or when the solution no longer changes.

---

## 8. Expected Results

### 8.1 Convergence Behavior

- Initial residuals may be high (especially for pressure and energy)
- Expect convergence within **300–800 iterations** for this case
- The energy equation (`h`) typically converges slower than momentum
- Residuals should reach $O(10^{-4})$ or lower

### 8.2 Temperature Field

The temperature field shows the development of thermal boundary layers:

```
  Temperature Profile at Different Axial Locations:
  
  y                    y                    y
  ↑  T=350K            ↑  T=350K            ↑  T=350K
  │──┐                 │──╮                 │───╮
  │  │                 │   ╲                │    ╲
  │  │                 │    │               │     │
  │  │  T≈300K         │    │  T≈310K       │     │  T≈325K
  │  │  (uniform)      │    │  (warming)    │     │  (warmer)
  │  │                 │   ╱                │    ╱
  │──┘                 │──╯                 │───╯
  │  T=350K            │  T=350K            │  T=350K
  └──────→ T           └──────→ T           └──────→ T
   x = 0 (inlet)       x = 1 m (mid)        x = 2 m (outlet)
```

Key observations:
- Near the **inlet**, the core temperature remains close to 300 K with thin hot layers
  near the walls.
- Moving **downstream**, the heated layer grows inward from both walls.
- At the **outlet**, the bulk fluid temperature has increased significantly.
- The temperature profile is symmetric about the duct centerline (y = 0.05 m).

### 8.3 Velocity Field

- The velocity develops from a uniform inlet profile toward a turbulent-like profile.
- Near walls, velocity decreases sharply to zero (no-slip condition).
- The centerline velocity is slightly higher than the mean to satisfy continuity.
- Buoyancy effects may cause slight asymmetry (hot fluid rises), but forced convection
  dominates in this configuration.

### 8.4 Heat Transfer Analysis

**Bulk temperature rise:**
The energy balance gives the expected outlet bulk temperature:

$$T_{out,bulk} \approx T_{in} + (2 \times h_{conv} \times L \times \Delta T) / (\dot{m} \times Cp)$$

where ṁ is the mass flow rate. The actual outlet temperature will depend on the
converged heat transfer coefficient.

**Wall heat flux:**
The wall heat flux is highest near the inlet (where the temperature gradient is steepest)
and decreases along the duct as the fluid warms up and the driving temperature difference
between wall and fluid core diminishes.

### 8.5 Turbulence Fields

- `k` is highest near the walls where velocity gradients are large
- `epsilon` peaks very close to the walls
- `nut` (turbulent viscosity) reaches maximum values slightly away from the wall
- Turbulence enhances mixing, which increases heat transfer compared to laminar flow

---

## 9. Post-Processing Guide

### 9.1 Opening in ParaView

```bash
# From the case directory
paraFoam

# Or create a .foam file and open manually
touch heated_pipe.foam
paraview heated_pipe.foam
```

### 9.2 Temperature Contour Plot

1. Open the case in ParaView
2. Click **Apply** to load the data
3. In the **Properties** panel, select the last time step
4. From the variable dropdown (top toolbar), select **T**
5. The color map will show temperature distribution from 300 K (blue) to 350 K (red)
6. Adjust the color scale: **Edit → Rescale to Data Range**

### 9.3 Velocity Vector Plot

1. Apply a **Glyph** filter (Filters → Alphabetical → Glyph)
2. Set Glyph Type to **Arrow**
3. Orient by: **U**
4. Scale by: **U** magnitude
5. Apply to see velocity vectors showing the flow pattern

### 9.4 Temperature and Velocity Line Plots

To extract profiles at specific locations:

1. Apply a **Plot Over Line** filter
2. Set the line endpoints:
   - For a **cross-section at x = 1 m**: Point1 = (1, 0, 0.005), Point2 = (1, 0.1, 0.005)
   - For the **centerline**: Point1 = (0, 0.05, 0.005), Point2 = (2, 0.05, 0.005)
3. Select the variable(s) to plot (T, U magnitude, k, etc.)
4. Click **Apply** — a line chart will appear

### 9.5 Calculating Nusselt Number

You can post-process the wall heat flux to calculate the local Nusselt number:

```bash
# Using OpenFOAM's wallHeatFlux utility (if available)
buoyantSimpleFoam -postProcess -func wallHeatFlux

# Or use the postProcess utility
postProcess -func "wallHeatFlux" -latestTime
```

Then in ParaView or with a script:

$$Nu_{local}(x) = q_{wall}(x) \times D_h / (k_{thermal} \times (T_{wall} - T_{bulk}(x)))$$

### 9.6 Extracting Data with OpenFOAM Utilities

```bash
# Sample along a line (e.g., centerline temperature)
postProcess -func "singleGraph"

# Write wall heat flux
postProcess -func "wallHeatFlux" -latestTime

# Calculate y+ at walls
postProcess -func "yPlus" -latestTime
```

---

## 10. Exercises

These exercises progressively explore the physics and numerics of heated duct flow.
They are ordered from simple parameter changes to more significant modifications.

### Exercise 1: Vary Wall Temperature

**Goal:** Understand the effect of driving temperature difference on heat transfer.

- Change `T_wall` from 350 K to 400 K, 500 K, and 600 K
- Edit both wall patches in `0/T`
- Compare temperature profiles at the outlet cross-section
- **Question:** How does the bulk outlet temperature change? Is it linear with $\Delta T$?

### Exercise 2: Change Inlet Velocity (Reynolds Number)

**Goal:** Explore the relationship between flow rate and heat transfer.

- Try U_in = 0.1, 0.5, 1.0, and 2.0 m/s (edit `0/U`)
- Remember to adjust `k` and `epsilon` inlet values for consistency:
  - k $\approx$ 1.5 $\times$ (U $\times$ I)$^2$ where I is turbulence intensity (e.g., 5%)
  - epsilon $\approx$ C_$\mu$^0.75 $\times$ k^1.5 / l, where l is the turbulence length scale
- Plot Nu vs Re and compare with Dittus-Boelter correlation
- **Question:** At what velocity does the flow transition from laminar-like to turbulent?

### Exercise 3: Refine the Mesh

**Goal:** Perform a mesh sensitivity study.

- Change the mesh resolution in `system/blockMeshDict`:
  - Coarse: 50$\times$10$\times$1 (500 cells)
  - Medium: 100$\times$20$\times$1 (2000 cells) ← current
  - Fine: 200$\times$40$\times$1 (8000 cells)
  - Very fine: 400$\times$80$\times$1 (32000 cells)
- Add y-direction grading for better wall resolution:
  ```
  simpleGrading (1 ((0.5 0.5 0.25)(0.5 0.5 4)) 1)
  ```
- Compare outlet temperature profiles across meshes
- **Question:** At what resolution do results become mesh-independent?

### Exercise 4: Add Wall Grading

**Goal:** Improve wall resolution for better heat transfer prediction.

- Modify `blockMeshDict` to use multi-grading in the y-direction:
  ```
  simpleGrading (1 ((0.3 0.4 0.2)(0.4 0.2 1)(0.3 0.4 5)) 1)
  ```
- Check y+ values using `postProcess -func yPlus`
- Target y+ $\approx$ 30–100 for wall functions, or y+ < 1 for resolving the viscous sublayer
- **Question:** How much does the Nusselt number change with better wall resolution?

### Exercise 5: Compare Turbulence Models

**Goal:** Assess the sensitivity of heat transfer to turbulence modeling.

- Try different turbulence models by editing `constant/turbulenceProperties`:
  - `kEpsilon` (current)
  - `kOmegaSST` (requires `0/omega` instead of `0/epsilon`)
  - `RNGkEpsilon`
  - `realizableKE`
- Remember to update wall functions accordingly
- **Question:** Which model gives results closest to empirical correlations?

> See [Turbulence Models](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md) for detailed guidance.

### Exercise 6: Non-Uniform Wall Temperature

**Goal:** Simulate a more realistic heating scenario.

- Replace the uniform wall temperature with a temperature profile using `codedFixedValue`:
  ```cpp
  topWall
  {
      type            codedFixedValue;
      value           uniform 300;
      name            linearWallTemp;
      code
      #{
          const vectorField& Cf = patch().Cf();
          scalarField& field = *this;
          forAll(Cf, i)
          {
              scalar x = Cf[i].x();
              field[i] = 300 + 25*x;  // Linear ramp from 300K to 350K
          }
      #};
  }
  ```
- **Question:** How does the heat flux distribution change compared to uniform wall T?

### Exercise 7: Constant Heat Flux Wall

**Goal:** Compare constant-temperature vs constant-heat-flux boundary conditions.

- Change wall BC from `fixedValue` (constant T) to `fixedGradient` (constant heat flux):
  ```cpp
  topWall
  {
      type            fixedGradient;
      gradient        uniform 500;  // dT/dn at the wall
  }
  ```
- **Question:** How does the wall temperature distribution differ? Which BC gives higher
  average heat transfer?

### Exercise 8: Longer Duct

**Goal:** Approach fully-developed conditions.

- Extend the duct to 10 m (change `blockMeshDict` vertex coordinates)
- Increase x-direction cells to 500 (maintaining $\Delta x$ $\approx$ 0.02 m)
- Track the bulk temperature and Nusselt number along the length
- **Question:** Does the Nusselt number reach the fully-developed value?
  Compare with the theoretical value for constant wall temperature.

### Exercise 9: Parallel Execution

**Goal:** Learn domain decomposition and parallel simulation.

- Create `system/decomposeParDict` for 4 processors:
  ```cpp
  numberOfSubdomains 4;
  method scotch;
  ```
- Run in parallel:
  ```bash
  decomposePar
  mpirun -np 4 buoyantSimpleFoam -parallel
  reconstructPar
  ```
- **Question:** Is the result identical to the serial run? What is the speedup?

> See [Parallelization Guide](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/07_parallelization.md) for details.

### Exercise 10: Validate Against Analytical Solutions

**Goal:** Build confidence in the simulation by comparing with known results.

- For **laminar flow** (reduce Re by lowering velocity to ~0.01 m/s and turning off
  turbulence: `simulationType laminar`):
  - Compare the velocity profile with the analytical parabolic Poiseuille flow profile
  - Compare the Nusselt number with $Nu = 7.54$ (constant wall T, parallel plates, fully developed)
- For **turbulent flow**:
  - Compare with the Dittus-Boelter and Gnielinski correlations
  - Plot Nu vs x/D_h and compare with correlations for the thermal entry region

---

## 11. Troubleshooting

### Common Issues

| Problem | Likely Cause | Solution |
|---|---|---|
| Solver crashes immediately | Missing field files | Verify all 8 field files exist in `0/` |
| "No mesh" error | blockMesh not run | Run `blockMesh` before the solver |
| Divergence (residuals blow up) | Relaxation too aggressive | Reduce relaxation factors (p_rgh: 0.2, U: 0.5) |
| Very slow convergence | Under-relaxation too conservative | Increase relaxation factors slightly |
| Temperature exceeds 350 K | Numerical oscillations | Use `bounded` upwind schemes for energy |
| Negative k or epsilon | Turbulence instability | Ensure `bounded` in div schemes; check inlet values |
| Non-physical reverse flow | Pressure BC issue | Ensure outlet `p_rgh` is `fixedValue` |
| y+ too low for wall functions | Mesh too fine near walls | Either coarsen wall cells or switch to low-Re model |

### Checking Mesh Quality

```bash
checkMesh
```

Look for:
- **Max aspect ratio** < 100 (ideally < 10)
- **Max non-orthogonality** < 70° (ideally < 40°)
- **Max skewness** < 4 (ideally < 1)
- **Mesh OK** message at the end

### Checking Convergence

After the solver finishes, examine the log file:

```bash
# Check final residuals
grep "Solving for Ux" log.buoyantSimpleFoam | tail -5
grep "Solving for h"  log.buoyantSimpleFoam | tail -5
grep "Solving for p_rgh" log.buoyantSimpleFoam | tail -5
```

If residuals stagnate above the target, try:
1. Increasing `endTime` (more iterations)
2. Tightening relaxation factors
3. Refining the mesh
4. Switching to second-order schemes gradually

---

## 12. References and Further Reading

### Theory & Correlations

- **Incropera, F.P. & DeWitt, D.P.** — *Fundamentals of Heat and Mass Transfer*, 7th Ed.
  Chapter 8: Internal Flow (pipe flow heat transfer, Nusselt correlations)
- **Cengel, Y.A. & Cimbala, J.M.** — *Fluid Mechanics: Fundamentals and Applications*
- **Bejan, A.** — *Convection Heat Transfer*, 4th Ed.

### OpenFOAM Documentation

- [OpenFOAM User Guide — buoyantSimpleFoam](https://www.openfoam.com/documentation/guides/latest/doc/guide-applications-solvers-heat-transfer-buoyantSimpleFoam.html)
- [OpenFOAM Wiki — Heat Transfer](https://openfoamwiki.net/index.php/Main_Page)
- [OpenFOAM Tutorials — buoyantCavity](https://github.com/OpenFOAM/OpenFOAM-6/tree/master/tutorials/heatTransfer/buoyantSimpleFoam)

### Related Notes in This Repository

> These notes provide the theoretical foundation for concepts used in this project:

| Note | Relevance |
|---|---|
| [01 — Introduction to CFD](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md) | Governing equations (continuity, momentum, energy) |
| [02 — OpenFOAM Cases](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md) | Case directory structure and file organization |
| [03 — OpenFOAM Dictionaries](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md) | Dictionary format, FoamFile headers, syntax |
| [04 — Meshing](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/04_meshing.md) | blockMesh usage, mesh quality, grading |
| [05 — Boundary Conditions](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md) | BC types used in this project (fixedValue, wall functions) |
| [06 — Turbulence Models](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md) | k-epsilon model theory, wall functions, model selection |
| [07 — Parallelization](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/07_parallelization.md) | Running cases in parallel (Exercise 9) |
| [08 — CFL Number](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md) | Stability and time stepping (relevant for transient extensions) |
| [09 — Linear Solvers](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/09_linear_solvers.md) | GAMG, smoothSolver, convergence, residuals |

### Related Projects

| Project | What It Teaches | How It Connects |
|---|---|---|
| [01 — Lid-Driven Cavity](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/01_lid_driven_cavity) | Basic CFD setup, icoFoam | Foundation for understanding BCs and mesh |
| [04 — NACA Airfoil](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis) | External aerodynamics, turbulence | Same k-epsilon model, different application |

---

*This project is part of the [OpenFOAM Tutorials](https://github.com/djeada/OpenFoam-Tutorials/blob/main/README.md) repository — a
hands-on learning path from CFD fundamentals to advanced simulations.*
