# OpenFOAM Dictionaries — The Complete Reference

Everything in OpenFOAM is controlled by plain-text **dictionary** files. There is no GUI, no
binary project file, and no database — just human-readable text organized into directories.
This design makes OpenFOAM cases easy to version-control, script, and reproduce.

This note is the deep-dive reference for every dictionary you will encounter in a typical
simulation. For boundary condition theory see [05_boundary_conditions.md](05_boundary_conditions.md),
for turbulence model selection see [06_turbulence_models.md](06_turbulence_models.md), and for
solver math see [09_linear_solvers.md](09_linear_solvers.md).

---

## 1 — Case Directory Structure

Every OpenFOAM case follows a strict three-directory layout:

```
 ┌─────────────────────────────────────────────────────────────────────────┐
 │                       OpenFOAM Case Directory                          │
 ├─────────────────────────────────────────────────────────────────────────┤
 │                                                                        │
 │   YourProject/                                                         │
 │   ├── 0/                    ◄── Initial & boundary conditions          │
 │   │   ├── U                     (one file per field variable)          │
 │   │   ├── p                                                            │
 │   │   ├── k                     ← only for turbulent cases             │
 │   │   ├── epsilon               ← only for turbulent cases             │
 │   │   └── nut                   ← only for turbulent cases             │
 │   │                                                                    │
 │   ├── constant/             ◄── Physical properties + mesh             │
 │   │   ├── transportProperties                                          │
 │   │   ├── turbulenceProperties                                         │
 │   │   ├── triSurface/           ← STL geometry files (if any)          │
 │   │   └── polyMesh/             ← mesh data (generated)                │
 │   │       ├── blockMeshDict     ← only for blockMesh-based meshes      │
 │   │       ├── points                                                   │
 │   │       ├── faces                                                    │
 │   │       ├── owner                                                    │
 │   │       └── neighbour                                                │
 │   │                                                                    │
 │   └── system/               ◄── Simulation controls                    │
 │       ├── controlDict           ← master control file                  │
 │       ├── fvSchemes             ← discretization schemes               │
 │       ├── fvSolution            ← linear solvers & algorithms          │
 │       ├── decomposeParDict      ← parallel decomposition               │
 │       └── snappyHexMeshDict     ← complex meshing (if used)            │
 │                                                                        │
 └─────────────────────────────────────────────────────────────────────────┘
```

| Directory   | Purpose | Typical files |
|-------------|---------|---------------|
| `0/`        | Initial and boundary conditions for every solved field | `U`, `p`, `k`, `epsilon`, `nut`, `omega` |
| `constant/` | Physical properties of the fluid and the mesh data | `transportProperties`, `turbulenceProperties`, `polyMesh/` |
| `system/`   | Numerical controls: solvers, schemes, time stepping, I/O | `controlDict`, `fvSchemes`, `fvSolution`, `decomposeParDict` |

> **Tip:** The laminar lid-driven cavity (`projects/01_lid_driven_cavity`) only has `U` and `p`
> in `0/`, while the turbulent NACA airfoil (`projects/03_naca_airfoil_analysis`) also has `k`,
> `epsilon`, and `nut`. The files you need in `0/` depend on which solver and turbulence model
> you are running.

---

## 2 — The FoamFile Header

**Every** OpenFOAM dictionary starts with a mandatory header block that tells the framework
how to parse the file:

```
 ┌───────────────────────────────────────────────────────────────────────┐
 │  FoamFile                                                            │
 │  {                                                                   │
 │      version     2.0;          ◄── File format version               │
 │      format      ascii;        ◄── ascii  or  binary                 │
 │      class       dictionary;   ◄── C++ class used to read this file  │
 │      location    "system";     ◄── (optional) directory hint         │
 │      object      controlDict;  ◄── Name that matches the filename   │
 │  }                                                                   │
 └───────────────────────────────────────────────────────────────────────┘
```

### Header Field Reference

| Field      | Required? | Description | Common values |
|------------|-----------|-------------|---------------|
| `version`  | Yes | File format version | `2.0` (always) |
| `format`   | Yes | Data encoding | `ascii`, `binary` |
| `class`    | Yes | C++ class that reads the file | See table below |
| `location` | No  | Hint for the directory the file lives in | `"system"`, `"constant"`, `"0"` |
| `object`   | Yes | Logical name — must match the filename | `controlDict`, `U`, `p`, etc. |

### Common `class` Values

| Class | Used for | Example file |
|-------|----------|-------------|
| `dictionary` | General key-value settings | `controlDict`, `fvSchemes`, `transportProperties` |
| `volScalarField` | Scalar field (one value per cell) | `p`, `k`, `epsilon`, `nut` |
| `volVectorField` | Vector field (three components per cell) | `U` |
| `volSymmTensorField` | Symmetric tensor field | `R` (Reynolds stress) |
| `surfaceScalarField` | Scalar on face centres | `phi` (flux field) |
| `pointScalarField` | Scalar at mesh points | Point-based displacement |
| `uniformDimensionedVectorField` | Single global vector | `g` (gravity) |

> **Tip:** If you see `class volScalarField;` you know the file describes a per-cell scalar
> field and must contain `dimensions`, `internalField`, and `boundaryField` entries.

**Real example** — from `projects/01_lid_driven_cavity/0/p`:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;   // ← tells OpenFOAM this is a scalar field
    object      p;
}
```

**Real example** — from `projects/01_lid_driven_cavity/0/U`:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;   // ← tells OpenFOAM this is a vector field
    object      U;
}
```

---

## 3 — Dictionary Syntax

OpenFOAM dictionaries use a C++-inspired syntax. Understanding these building blocks lets
you read and write any file in the case.

### 3.1 — Key-Value Pairs

The most basic element. A keyword followed by a value terminated by a semicolon:

```c
application     icoFoam;        // string value
startTime       0;              // scalar value
runTimeModifiable true;         // boolean (switch)
```

### 3.2 — Sub-Dictionaries (Nested Braces)

Group related keywords inside braces:

```c
PISO
{
    nCorrectors                 2;
    nNonOrthogonalCorrectors    0;
}
```

Sub-dictionaries can nest to arbitrary depth (e.g. `boundaryField → inlet → { ... }`).

### 3.3 — Lists and Arrays

Round parentheses define lists:

```c
internalField   uniform (0 0 0);              // inline vector
vertices
(
    (0 0 0)    // vertex 0
    (1 0 0)    // vertex 1
    (1 1 0)    // vertex 2
);
```

### 3.4 — Regex Matching for Boundary Names

Boundary names can be POSIX regular expressions in double quotes:

```c
boundaryField
{
    ".*Wall"          // matches movingWall, fixedWalls, etc.
    {
        type    noSlip;
    }
    "(inlet|outlet)"  // matches inlet OR outlet
    {
        type    zeroGradient;
    }
}
```

### 3.5 — Include Directives and Macros

| Directive | Purpose | Example |
|-----------|---------|---------|
| `#include "filename"` | Include another dictionary file | `#include "initialConditions"` |
| `#includeIfPresent "f"` | Include if file exists, skip otherwise | `#includeIfPresent "debug"` |
| `#calc "expr"` | Evaluate C++ expression at parse time | `#calc "0.1 * 0.05"` |
| `$keyword` | Variable expansion (reuse a value) | `$nu;` |
| `$:path.to.key` | Scoped variable lookup | `$:SIMPLE.nNonOrthogonalCorrectors;` |

### 3.6 — Dimensional Units

OpenFOAM enforces physical dimensions at runtime using a 7-element array:

```
    [ kg   m   s   K   mol   A   cd ]
      ──   ─   ─   ─   ───   ─   ──
      Mass Len Time Temp Mol  Amp Lum
```

This is covered in depth in the next section.

---

## 4 — The Dimensions System

OpenFOAM tracks physical dimensions for every field and constant using the SI base units.
If you accidentally assign wrong dimensions, the solver **will** crash with a dimensional
inconsistency error — this is a powerful safety net.

```
 ┌──────────────────────────────────────────────────────────────────────┐
 │         The 7 Base Dimensions in OpenFOAM                           │
 │                                                                     │
 │  Index:   [  0      1      2      3      4      5      6   ]       │
 │  Unit:    [  kg     m      s      K      mol    A      cd  ]       │
 │  Name:    [ Mass  Length  Time  Temp   Moles  Amps  Candela]       │
 │                                                                     │
 │  Example — velocity:       [ 0  1  -1  0  0  0  0 ]               │
 │            means:  m^1 · s^(-1)  →  m/s   ✓                        │
 │                                                                     │
 │  Example — kinematic pressure: [ 0  2  -2  0  0  0  0 ]           │
 │            means:  m^2 · s^(-2)  →  m²/s²  ✓                       │
 └──────────────────────────────────────────────────────────────────────┘
```

### Common Dimension Sets

| Quantity | Dimensions | SI Unit | OpenFOAM notation |
|----------|-----------|---------|-------------------|
| Velocity | $m/s$ | $m \cdot s^{-1}$ | `[0 1 -1 0 0 0 0]` |
| Kinematic pressure ($p/\rho$) | $m^{2}/s^{2}$ | $m^{2} \cdot s^{-2}$ | `[0 2 -2 0 0 0 0]` |
| Absolute pressure | $kg/(m \cdot s^{2})$ | Pa | `[1 -1 -2 0 0 0 0]` |
| Kinematic viscosity ($\nu$) | $m^{2}/s$ | $m^{2} \cdot s^{-1}$ | `[0 2 -1 0 0 0 0]` |
| Dynamic viscosity ($\mu$) | $kg/(m \cdot s)$ | $Pa \cdot s$ | `[1 -1 -1 0 0 0 0]` |
| Density | $kg/m^{3}$ | $kg \cdot m^{-3}$ | `[1 -3 0 0 0 0 0]` |
| Turbulent kinetic energy ($k$) | $m^{2}/s^{2}$ | $m^{2} \cdot s^{-2}$ | `[0 2 -2 0 0 0 0]` |
| Turbulent dissipation ($\varepsilon$) | $m^{2}/s^{3}$ | $m^{2} \cdot s^{-3}$ | `[0 2 -3 0 0 0 0]` |
| Specific turbulent dissipation ($\omega$) | $1/s$ | $s^{-1}$ | `[0 0 -1 0 0 0 0]` |
| Turbulent viscosity ($\nu_t$) | $m^{2}/s$ | $m^{2} \cdot s^{-1}$ | `[0 2 -1 0 0 0 0]` |
| Temperature | K | K | `[0 0 0 1 0 0 0]` |
| Dimensionless | — | — | `[0 0 0 0 0 0 0]` |
| Force per unit volume | N/m³ | kg·m⁻²·s⁻² | `[1 -2 -2 0 0 0 0]` |
| Gravity | m/s² | m·s⁻² | `[0 1 -2 0 0 0 0]` |

### Reading Dimensions from Project Files

From `projects/01_lid_driven_cavity/constant/transportProperties`:

```c
nu              [0 2 -1 0 0 0 0] 0.001;
//               │ │  │
//               │ │  └─ s^(-1)  ──┐
//               │ └─── m^2  ──────┤──▶  m²/s  (kinematic viscosity)
//               └───── kg^0 ──────┘
```

From `projects/03_naca_airfoil_analysis/0/epsilon`:

```c
dimensions      [0 2 -3 0 0 0 0];
//               │ │  │
//               │ │  └─ s^(-3)  ──┐
//               │ └─── m^2  ──────┤──▶  m²/s³  (turbulent dissipation rate)
//               └───── kg^0 ──────┘
```

> **Tip:** A quick sanity check — kinematic quantities (used by incompressible solvers like
> `icoFoam` and `simpleFoam`) always have `kg` exponent = 0, because density has been divided
> out. If you see `[1 ...]` as the first element, you are in absolute (compressible) units.

---

## 5 — `controlDict` — Master Simulation Control

The `controlDict` is the **first file the solver reads**. It determines which application runs,
how time advances, and what gets written to disk.

### 5.1 — Lid-Driven Cavity (Transient / icoFoam)

From `projects/01_lid_driven_cavity/system/controlDict`:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}

application     icoFoam;

startFrom       startTime;
startTime       0;

stopAt          endTime;
endTime         0.5;

deltaT          0.001;

writeControl    timeStep;
writeInterval   1;

purgeWrite      0;

writeFormat     ascii;
writePrecision  6;
writeCompression off;

timeFormat      general;
timePrecision   6;

runTimeModifiable true;
```

### 5.2 — NACA Airfoil (Steady-State / simpleFoam)

From `projects/03_naca_airfoil_analysis/system/controlDict`:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     simpleFoam;

startFrom       startTime;
startTime       0;

stopAt          endTime;
endTime         1000;

deltaT          0.1;

writeControl    timeStep;
writeInterval   1;

purgeWrite      0;

writeFormat     ascii;
writePrecision  6;
writeCompression off;

timeFormat      general;
timePrecision   6;

runTimeModifiable true;
```

### 5.3 — Transient vs Steady-State Comparison

```
 ┌──────────────────────────────────┬──────────────────────────────────────┐
 │      TRANSIENT  (icoFoam)        │      STEADY-STATE  (simpleFoam)      │
 ├──────────────────────────────────┼──────────────────────────────────────┤
 │  deltaT = real physical Δt       │  deltaT = pseudo-time step           │
 │  endTime = physical end time     │  endTime = iteration count           │
 │  ddtSchemes → Euler / backward   │  ddtSchemes → steadyState            │
 │  Algorithm: PISO or PIMPLE       │  Algorithm: SIMPLE                   │
 │  CFL condition must be met       │  No CFL constraint                   │
 │  Example: endTime 0.5 (seconds)  │  Example: endTime 1000 (iterations)  │
 └──────────────────────────────────┴──────────────────────────────────────┘
```

### 5.4 — controlDict Keyword Reference

| Keyword | Type | Description | Common values |
|---------|------|-------------|---------------|
| `application` | word | Solver to run | `icoFoam`, `simpleFoam`, `pisoFoam`, `pimpleFoam`, `rhoSimpleFoam` |
| `startFrom` | word | Where to start | `firstTime`, `startTime`, `latestTime` |
| `startTime` | scalar | Start time value (used when `startFrom startTime`) | `0` |
| `stopAt` | word | When to stop | `endTime`, `writeNow`, `noWriteNow`, `nextWrite` |
| `endTime` | scalar | End time / iteration count | `0.5`, `1000` |
| `deltaT` | scalar | Time step | `0.001`, `1` |
| `writeControl` | word | What triggers writing output | `timeStep`, `runTime`, `adjustableRunTime`, `clockTime`, `cpuTime` |
| `writeInterval` | scalar | Interval between writes (in units of writeControl) | `1`, `0.05`, `100` |
| `purgeWrite` | int | Keep only N latest time dirs (0 = keep all) | `0`, `5`, `10` |
| `writeFormat` | word | Output data format | `ascii`, `binary` |
| `writePrecision` | int | Significant digits in output | `6`, `8`, `12` |
| `writeCompression` | switch | Compress output files | `on`, `off` |
| `timeFormat` | word | Time directory naming format | `fixed`, `scientific`, `general` |
| `timePrecision` | int | Digits in time directory names | `6` |
| `runTimeModifiable` | switch | Re-read controlDict during run | `true`, `false` |
| `adjustTimeStep` | switch | Auto-adjust deltaT to meet maxCo | `yes`, `no` |
| `maxCo` | scalar | Maximum Courant number (with adjustTimeStep) | `0.5`, `1.0` |
| `functions` | dict | Runtime post-processing (probes, forces, etc.) | Sub-dictionary |

### 5.5 — How controlDict Drives the Solver

```
 ┌─────────────┐
 │  Read        │
 │  controlDict │
 └──────┬───────┘
        │
        ▼
 ┌─────────────────┐     startFrom / startTime
 │  Set t = t_start │◄────────────────────────
 └──────┬──────────┘
        │
        ▼
 ┌─────────────────────────────────────────────────┐
 │              TIME LOOP  (while t < endTime)      │
 │  ┌───────────────────────────────────────────┐   │
 │  │  t = t + deltaT                           │   │
 │  │  Solve momentum, pressure, turbulence ... │   │
 │  │                                           │   │
 │  │  if (writeControl condition met)          │   │
 │  │      Write fields to disk ──────────────────▶ 0.001/  0.002/ ...
 │  │                                           │   │
 │  │  if (runTimeModifiable)                   │   │
 │  │      Re-read controlDict                  │   │
 │  └───────────────────────────────────────────┘   │
 └──────────────────┬──────────────────────────────┘
                    │
                    ▼
             ┌──────────┐
             │   Done   │
             └──────────┘
```

---

## 6 — `fvSchemes` — Discretization Schemes

The `fvSchemes` dictionary tells OpenFOAM **how** to convert the continuous PDEs into discrete
algebraic equations. Each mathematical operator (time derivative, gradient, divergence,
Laplacian) needs a numerical scheme.

### 6.1 — What Each Scheme Category Represents

```
 ┌───────────────────────────────────────────────────────────────────────┐
 │                   PDE Terms  →  fvSchemes Categories                 │
 │                                                                      │
 │   ∂U/∂t          ──────────────▶  ddtSchemes        (time deriv.)   │
 │                                                                      │
 │   ∇p, ∇U         ──────────────▶  gradSchemes       (gradient)      │
 │                                                                      │
 │   ∇·(φU)         ──────────────▶  divSchemes        (divergence)    │
 │                                                                      │
 │   ∇·(ν∇U)       ──────────────▶  laplacianSchemes  (Laplacian)     │
 │                                                                      │
 │   face values     ──────────────▶  interpolationSchemes              │
 │                                                                      │
 │   ∂/∂n (face     ──────────────▶  snGradSchemes     (surface-       │
 │    normal grad)                                       normal grad)   │
 └───────────────────────────────────────────────────────────────────────┘
```

### 6.2 — Lid-Driven Cavity fvSchemes (Transient)

From `projects/01_lid_driven_cavity/system/fvSchemes`:

```c
ddtSchemes
{
    default         Euler;
}

gradSchemes
{
    default         Gauss linear;
    grad(p)         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);
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

### 6.3 — NACA Airfoil fvSchemes (Steady-State + Turbulence)

From `projects/03_naca_airfoil_analysis/system/fvSchemes`:

```c
ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);
    div((nuEff*dev(T(grad(U))))) Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div(phi,epsilon) Gauss upwind;
    div(phi,k) Gauss upwind;
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

fluxRequired
{
    default         no;
    p;
}
```

> **Key difference:** The transient case uses `Euler` for `ddtSchemes`, while the steady-state
> case uses `steadyState`. The NACA case also needs extra `divSchemes` entries for the turbulence
> transport equations (`div(phi,k)`, `div(phi,epsilon)`) and the effective viscosity stress
> tensor terms.

### 6.4 — ddtSchemes (Time Derivative $\partial/\partial t$)

| Scheme | Order | Type | Description |
|--------|-------|------|-------------|
| `steadyState` | — | Steady | Sets time derivative to zero; for SIMPLE-based solvers |
| `Euler` | 1st | Implicit | First-order backward Euler; stable, diffusive |
| `backward` | 2nd | Implicit | Second-order backward differencing; less diffusive |
| `CrankNicolson <α>` | 2nd | Implicit | Blends Euler ($\alpha=0$) and Crank-Nicolson ($\alpha=1$); $\alpha=0.9$ typical |
| `localEuler` | 1st | Local | Local time stepping for steady-state acceleration |

### 6.5 — gradSchemes (Gradient $\nabla$)

| Scheme | Description | Typical use |
|--------|-------------|-------------|
| `Gauss linear` | Green-Gauss with linear interpolation to faces | Default, good for most cases |
| `leastSquares` | Least-squares gradient; better on irregular meshes | Unstructured meshes |
| `cellLimited Gauss linear <k>` | Limits gradient to prevent overshoots; k=1 is full limiting | High-gradient regions |
| `faceLimited Gauss linear <k>` | Limits at faces rather than cells | Alternative limiting |

### 6.6 — divSchemes (Divergence $\nabla \cdot$)

This is the most important category for **convection** (advection) terms. The choice
of divergence scheme controls numerical diffusion and stability.

| Scheme | Order | Bounded? | Description |
|--------|-------|----------|-------------|
| `Gauss upwind` | 1st | Yes | Most stable, most diffusive; good start for turbulence |
| `Gauss linear` | 2nd | No | Central differencing; accurate but can oscillate |
| `Gauss linearUpwind grad(φ)` | 2nd | Yes | Upwind-biased with gradient correction; good default |
| `Gauss limitedLinear <k>` | 2nd | TVD | Blends linear and upwind; k=1 is most limited |
| `Gauss vanLeer` | 2nd | TVD | Van Leer limiter; good balance of accuracy/stability |
| `Gauss MUSCL` | 2nd | TVD | Monotonic Upstream-centred Scheme; sharp gradients |
| `Gauss QUICK` | 3rd | No | Quadratic upstream interpolation; only for structured meshes |

> **Tip:** Start with `Gauss upwind` for turbulence quantities (`k`, `epsilon`) and
> `Gauss linearUpwind grad(U)` for velocity. Only move to higher-order schemes once the
> simulation is stable.

### 6.7 — laplacianSchemes (Laplacian $\nabla^{2}$)

General format: `Gauss <interpolation> <snGrad>`

| Combination | Description | When to use |
|------------|-------------|-------------|
| `Gauss linear corrected` | Second-order, explicit non-orthogonal correction | Default for good-quality meshes |
| `Gauss linear uncorrected` | No correction; first-order on non-orthogonal meshes | Very orthogonal meshes only |
| `Gauss linear limited <k>` | Blend of corrected/uncorrected; k=0.5 typical | Moderate non-orthogonality |
| `Gauss linear orthogonal` | Alias for uncorrected; assumes perfect orthogonality | Cartesian/hex meshes |

### 6.8 — interpolationSchemes and snGradSchemes

| Category | Default | Description |
|----------|---------|-------------|
| `interpolationSchemes` | `linear` | How cell-centre values are interpolated to face centres |
| `snGradSchemes` | `corrected` | Surface-normal gradient at faces; same options as laplacian snGrad |

### 6.9 — Accuracy vs Stability Trade-Off

```
     Stability                                        Accuracy
     ◄──────────────────────────────────────────────────────────►
     ┌──────────┐  ┌───────────┐  ┌──────────────┐  ┌─────────┐
     │  upwind  │  │ vanLeer   │  │ linearUpwind │  │ linear  │
     │  (1st)   │  │ (2nd TVD) │  │ (2nd)        │  │ (2nd)   │
     └──────────┘  └───────────┘  └──────────────┘  └─────────┘
     Most stable                                     Most accurate
     Most diffusive                                  Can oscillate
```

---

## 7 — `fvSolution` — Solvers and Algorithms

The `fvSolution` dictionary controls **how** the discretized equations are solved: which linear
solver to use for each field, convergence tolerances, and the pressure-velocity coupling
algorithm.

### 7.1 — Lid-Driven Cavity fvSolution

From `projects/01_lid_driven_cavity/system/fvSolution`:

```c
solvers
{
    p
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0.05;
    }

    pFinal
    {
        $p;
        relTol          0;
    }

    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
}

PISO
{
    nCorrectors     2;
    nNonOrthogonalCorrectors 0;
    pRefCell        0;
    pRefValue       0;
}
```

### 7.2 — NACA Airfoil fvSolution

From `projects/03_naca_airfoil_analysis/system/fvSolution`:

```c
solvers
{
    p
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0.05;
    }

    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

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
    fields
    {
        p       0.3;
    }
    equations
    {
        U       0.7;
        k       0.7;
        epsilon 0.7;
    }
}
```

### 7.3 — Linear Solver Selection Guide

```
              Which linear solver should I use?
              ─────────────────────────────────
                         │
                         ▼
                 Is the matrix symmetric?
                 (pressure usually is)
                    /            \
                  YES             NO
                  /                \
                 ▼                  ▼
          ┌───────────┐     ┌──────────────┐
          │    PCG     │     │  Is the case  │
          │  + DIC     │     │   very large? │
          └───────────┘     └──┬─────────┬──┘
                               /          \
                             YES           NO
                             /              \
                            ▼                ▼
                    ┌───────────┐    ┌──────────────┐
                    │   GAMG    │    │ smoothSolver  │
                    │ (multigrid│    │ + symGaussSeidel
                    │  fastest  │    │  or GaussSeidel│
                    │  for large│    └──────────────┘
                    │  cases)   │
                    └───────────┘
```

### 7.4 — Linear Solver Reference

| Solver | Full Name | Matrix Type | Best for |
|--------|-----------|-------------|----------|
| `PCG` | Preconditioned Conjugate Gradient | Symmetric | Pressure (incompressible) |
| `PBiCGStab` | Preconditioned Bi-Conjugate Gradient Stabilized | Asymmetric | U, k, $\varepsilon$, $\omega$ (general purpose) |
| `smoothSolver` | Iterative smoother | Any | U, k, $\varepsilon$ — simple and robust |
| `GAMG` | Geometric-Algebraic Multi-Grid | Any | Large cases; fastest for pressure |
| `diagonal` | Diagonal solver | Diagonal | Trivial systems only |
| `PBiCG` | Preconditioned Bi-Conjugate Gradient (legacy) | Asymmetric | Older synonym; prefer PBiCGStab |

### 7.5 — Preconditioner Reference

| Preconditioner | Full Name | Use with | Notes |
|----------------|-----------|----------|-------|
| `DIC` | Diagonal Incomplete Cholesky | `PCG` (symmetric) | Default for pressure |
| `FDIC` | Faster DIC | `PCG` (symmetric) | Less memory, similar performance |
| `DILU` | Diagonal Incomplete LU | `PBiCGStab` (asymmetric) | Default for velocity, turbulence |
| `GAMG` | Multigrid as preconditioner | Any | Heavyweight but powerful |
| `none` | No preconditioning | Any | Only for testing |

### 7.6 — Tolerance Settings

| Keyword | Description | Typical value |
|---------|-------------|---------------|
| `tolerance` | Absolute residual to stop iterating | `1e-06` (p), `1e-05` (U, k, $\varepsilon$) |
| `relTol` | Relative drop from initial residual | `0.05`–`0.1` (SIMPLE), `0` (final corrector in PISO) |
| `minIter` | Minimum iterations before checking convergence | `1` |
| `maxIter` | Maximum iterations per time step | `1000` (default) |

> **Tip:** In PISO/PIMPLE, the final corrector step should use `relTol 0;` so the solver
> iterates until the absolute `tolerance` is met. This is why the cavity case defines `pFinal`
> with `$p;` (inherits all settings from `p`) but overrides `relTol 0;`.

### 7.7 — Pressure-Velocity Coupling Algorithms

| Algorithm | Full Name | Solver type | How it works |
|-----------|-----------|-------------|--------------|
| `SIMPLE` | Semi-Implicit Method for Pressure-Linked Equations | Steady-state | One pressure correction per iteration; needs relaxation |
| `PISO` | Pressure-Implicit with Splitting of Operators | Transient | Multiple pressure correctors per time step; no relaxation needed |
| `PIMPLE` | Merged PISO-SIMPLE | Both | Outer SIMPLE loops + inner PISO correctors; most flexible |

### 7.8 — Algorithm Keywords

| Keyword | Used in | Description | Typical value |
|---------|---------|-------------|---------------|
| `nCorrectors` | PISO/PIMPLE | Number of pressure corrector steps | `2`–`4` |
| `nNonOrthogonalCorrectors` | All | Extra corrections for non-orthogonal meshes | `0`–`3` |
| `nOuterCorrectors` | PIMPLE | Number of outer (SIMPLE-like) loops | `1`–`50` |
| `pRefCell` | PISO/PIMPLE | Reference cell for pressure (closed domains) | `0` |
| `pRefValue` | PISO/PIMPLE | Reference pressure value | `0` |
| `consistent` | SIMPLE | Use SIMPLEC variant (less relaxation needed) | `yes`, `no` |
| `residualControl` | SIMPLE | Stop iterating when residuals drop below thresholds | Sub-dictionary |

### 7.9 — Relaxation Factors

Relaxation is essential for **SIMPLE** to prevent divergence. Values between 0 and 1 control
how aggressively the solution updates each iteration.

| Factor range | Behaviour |
|-------------|-----------|
| 0.0–0.3 | Very conservative, slow convergence, high stability |
| 0.3–0.5 | Typical for pressure |
| 0.5–0.8 | Typical for velocity and turbulence |
| 0.8–1.0 | Aggressive; faster but may diverge |

The NACA airfoil case uses:

```c
relaxationFactors
{
    fields
    {
        p       0.3;       // conservative for pressure
    }
    equations
    {
        U       0.7;       // moderate for velocity
        k       0.7;       // moderate for turbulence
        epsilon 0.7;       // moderate for turbulence
    }
}
```

> **Tip:** If your steady-state simulation diverges, first try lowering the relaxation factors
> (e.g., p → 0.2, U → 0.5). If it converges but slowly, try raising them. With
> `consistent yes;` (SIMPLEC), you can often use higher p relaxation (0.5–0.7).

---

## 8 — `transportProperties` — Fluid Properties

This dictionary defines the transport model and fluid properties, primarily the kinematic
viscosity for incompressible solvers.

### Lid-Driven Cavity (Laminar)

From `projects/01_lid_driven_cavity/constant/transportProperties`:

```c
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}

transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 0.001;
```

### NACA Airfoil (Turbulent)

From `projects/03_naca_airfoil_analysis/constant/transportProperties`:

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

### Transport Model Reference

| Model | Description | When to use |
|-------|-------------|-------------|
| `Newtonian` | Constant viscosity | Most engineering flows (water, air) |
| `CrossPowerLaw` | Shear-thinning/thickening power law | Blood flow, polymers |
| `BirdCarreau` | Bird-Carreau shear-thinning model | More accurate non-Newtonian |
| `HerschelBulkley` | Yield-stress fluid with power law | Cement, toothpaste, mud |
| `powerLaw` | Simple power-law viscosity | Basic non-Newtonian |

> **Note:** The viscosity value differs greatly between the two projects:
> - Cavity: $\nu = 0.001$ $m^{2}/s$ $\rightarrow$ $Re = UL/\nu = 1 \cdot 0.1/0.001 = 100$ (laminar)
> - Airfoil: $\nu = 1 \times 10^{-5}$ $m^{2}/s$ $\rightarrow$ $Re \approx 1 \cdot 1/10^{-5} = 100{,}000$ (turbulent — needs $k$-$\varepsilon$ model)

---

## 9 — `turbulenceProperties` — Turbulence Model Selection

This dictionary selects the turbulence modelling approach. It is **always required**, even for
laminar cases (where you set `simulationType laminar;`).

### NACA Airfoil (k-epsilon)

From `projects/03_naca_airfoil_analysis/constant/turbulenceProperties`:

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

### Simulation Type Reference

| `simulationType` | Description | Required extra files |
|-------------------|-------------|---------------------|
| `laminar` | No turbulence model | None |
| `RAS` | Reynolds-Averaged Simulation (RANS) | `k`, `epsilon` or `k`, `omega`, `nut` in `0/` |
| `LES` | Large Eddy Simulation | Sub-grid scale fields in `0/` |
| `DES` | Detached Eddy Simulation | Hybrid RANS/LES fields |

### Common RAS Models

| Model | Fields needed | Best for |
|-------|--------------|----------|
| `kEpsilon` | k, epsilon, nut | General purpose; industrial flows |
| `kOmega` | k, omega, nut | Near-wall flows; better separation |
| `kOmegaSST` | k, omega, nut | Best all-round RANS model; recommended default |
| `SpalartAllmaras` | nuTilda, nut | Aerospace; one-equation model |
| `realizableKE` | k, epsilon, nut | Improved k-ε for jets, separation |
| `LRR` | R, epsilon, nut | Reynolds Stress Model; anisotropic flows |

> **Tip:** For most external aerodynamics cases, `kOmegaSST` is the recommended starting point.
> The NACA project uses `kEpsilon` for simplicity, but for production work you would switch.
> See [06_turbulence_models.md](06_turbulence_models.md) for detailed theory.

---

## 10 — Boundary Condition Dictionaries (`0/` Directory)

Each file in `0/` defines one field variable. The structure is always:

```
 ┌───────────────────────────────────────────────────────────────────┐
 │  FoamFile { ... }                                                │
 │                                                                  │
 │  dimensions      [  ...  ];      ◄── physical dimensions         │
 │                                                                  │
 │  internalField   uniform <val>;  ◄── initial value everywhere    │
 │                                                                  │
 │  boundaryField                                                   │
 │  {                                                               │
 │      patchName1                                                  │
 │      {                                                           │
 │          type    <bcType>;       ◄── boundary condition type     │
 │          value   uniform <val>;  ◄── (if required by type)       │
 │      }                                                           │
 │      patchName2 { ... }                                          │
 │  }                                                               │
 └───────────────────────────────────────────────────────────────────┘
```

### 10.1 — Laminar Case: Lid-Driven Cavity

**Velocity** — from `projects/01_lid_driven_cavity/0/U`:

```c
dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    movingWall
    {
        type            fixedValue;
        value           uniform (1 0 0);
    }

    fixedWalls
    {
        type            noSlip;
    }

    frontAndBack
    {
        type            empty;
    }

    inlet
    {
        type            fixedValue;
        value           uniform (1 0 0.5);
    }

    outlet
    {
        type            zeroGradient;
    }
}
```

**Pressure** — from `projects/01_lid_driven_cavity/0/p`:

```c
dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    movingWall
    {
        type            zeroGradient;
    }

    fixedWalls
    {
        type            zeroGradient;
    }

    frontAndBack
    {
        type            empty;
    }

    inlet
    {
        type            totalPressure;
        p0              uniform 1e5;
    }

    outlet
    {
        type            fixedValue;
        value           uniform 0;
    }
}
```

### 10.2 — Turbulent Case: NACA Airfoil

The turbulent case needs five field files. Here are the key entries:

**Velocity (U)** — from `projects/03_naca_airfoil_analysis/0/U`:

```c
dimensions      [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);

boundaryField
{
    inlet       { type fixedValue;    value uniform (1 0 0); }
    outlet      { type zeroGradient; }
    walls       { type noSlip; }
    airfoil     { type noSlip; }
    frontAndBack { type empty; }
}
```

**Turbulent kinetic energy (k)** — from `projects/03_naca_airfoil_analysis/0/k`:

```c
dimensions      [0 2 -2 0 0 0 0];
internalField   uniform 0.01;

boundaryField
{
    inlet       { type fixedValue;        value uniform 0.01; }
    outlet      { type zeroGradient; }
    walls       { type kqRWallFunction;   value uniform 0.01; }
    airfoil     { type kqRWallFunction;   value uniform 0.01; }
    frontAndBack { type empty; }
}
```

**Turbulent dissipation (epsilon)** — from `projects/03_naca_airfoil_analysis/0/epsilon`:

```c
dimensions      [0 2 -3 0 0 0 0];
internalField   uniform 0.01;

boundaryField
{
    inlet       { type fixedValue;              value uniform 0.01; }
    outlet      { type zeroGradient; }
    walls       { type epsilonWallFunction;     value uniform 0.01; }
    airfoil     { type epsilonWallFunction;     value uniform 0.01; }
    frontAndBack { type empty; }
}
```

**Turbulent viscosity (nut)** — from `projects/03_naca_airfoil_analysis/0/nut`:

```c
dimensions      [0 2 -1 0 0 0 0];
internalField   uniform 0;

boundaryField
{
    inlet       { type calculated;          value uniform 0; }
    outlet      { type calculated;          value uniform 0; }
    walls       { type nutkWallFunction;    value uniform 0; }
    airfoil     { type nutkWallFunction;    value uniform 0; }
    frontAndBack { type empty; }
}
```

### 10.3 — Common Boundary Condition Types

| BC Type | Category | Description |
|---------|----------|-------------|
| `fixedValue` | Dirichlet | Fixed value at the boundary |
| `zeroGradient` | Neumann | Zero normal gradient ($\partial\phi/\partial n = 0$) |
| `noSlip` | Dirichlet | Velocity = 0 at walls (shorthand for fixedValue (0 0 0)) |
| `empty` | Constraint | 2D cases — no solution on front/back faces |
| `symmetry` | Constraint | Symmetry plane — zero normal component |
| `totalPressure` | Mixed | $p_0 = p + \frac{1}{2}\rho|U|^{2}$ — adjusts with velocity |
| `fixedGradient` | Neumann | Specified normal gradient (e.g., heat flux) |
| `inletOutlet` | Mixed | Switches between fixedValue (inflow) and zeroGradient (outflow) |
| `calculated` | Derived | Value computed from other fields (e.g., nut from k, $\varepsilon$) |
| `kqRWallFunction` | Wall func. | Wall function for k, q, R fields |
| `epsilonWallFunction` | Wall func. | Wall function for $\varepsilon$ near walls |
| `nutkWallFunction` | Wall func. | Wall function for $\nu_t$ based on k |

> **Tip:** For a complete treatment of boundary conditions, including inlet turbulence
> estimation and wall function theory, see [05_boundary_conditions.md](05_boundary_conditions.md).

---

## 11 — `blockMeshDict` — Mesh Generation

The `blockMeshDict` (located in `constant/polyMesh/` or `system/`) defines simple structured
meshes from vertices, blocks, and boundary patches. It is the simplest meshing tool in
OpenFOAM.

```
 ┌───────────────────────────────────────────────────────────────┐
 │                  blockMeshDict Structure                      │
 │                                                               │
 │   vertices ( ... );       ◄── List of (x y z) coordinates    │
 │   blocks   ( ... );       ◄── Hex blocks with cell counts    │
 │   edges    ( ... );       ◄── Curved edges (optional)        │
 │   boundary ( ... );       ◄── Named patches with face lists  │
 │   mergePatchPairs ( ... );◄── Stitch patches (optional)      │
 └───────────────────────────────────────────────────────────────┘
```

| Keyword | Description |
|---------|-------------|
| `vertices` | Ordered list of 3D points defining block corners |
| `blocks` | `hex (v0 v1 ... v7) (nx ny nz) simpleGrading (gx gy gz)` |
| `edges` | `arc`, `spline`, `polyLine` for curved block edges |
| `boundary` | Named patches: `{ type patch; faces ( (f0 f1 f2 f3) ... ); }` |
| `mergePatchPairs` | Pairs of patches to merge for multi-block topologies |

> For detailed meshing workflows including `snappyHexMesh`, see [04_meshing.md](04_meshing.md).

---

## 12 — `decomposeParDict` — Parallel Decomposition

When running in parallel with `mpirun`, OpenFOAM splits the mesh and fields across processors.
The `decomposeParDict` (in `system/`) controls how.

```c
numberOfSubdomains  4;       // must match -np in mpirun

method              scotch;  // decomposition method

// Alternative: simple geometric decomposition
// method          simple;
// simpleCoeffs
// {
//     n           (2 2 1);
//     delta       0.001;
// }
```

| Method | Description | When to use |
|--------|-------------|-------------|
| `scotch` | Graph-based; minimises inter-processor communication | Default; best for most cases |
| `hierarchical` | Geometric splitting along axes in order | Simple geometries |
| `simple` | Geometric splitting with fixed divisions | When you know the geometry well |
| `manual` | User specifies processor per cell | Special cases |

> For a full treatment of parallel execution, see [07_parallelization.md](07_parallelization.md).

---

## 13 — Master Dictionary Reference Table

| Dictionary | Location | Class | Purpose | Used in cavity? | Used in NACA? |
|-----------|----------|-------|---------|:-:|:-:|
| `controlDict` | `system/` | `dictionary` | Solver selection, time control, I/O | ✓ | ✓ |
| `fvSchemes` | `system/` | `dictionary` | Discretization schemes for PDE terms | ✓ | ✓ |
| `fvSolution` | `system/` | `dictionary` | Linear solvers, algorithms, relaxation | ✓ | ✓ |
| `decomposeParDict` | `system/` | `dictionary` | Parallel domain decomposition | — | — |
| `snappyHexMeshDict` | `system/` | `dictionary` | Complex mesh generation | — | — |
| `blockMeshDict` | `constant/polyMesh/` | `dictionary` | Simple structured mesh | ✓ | — |
| `transportProperties` | `constant/` | `dictionary` | Fluid transport model and viscosity | ✓ | ✓ |
| `turbulenceProperties` | `constant/` | `dictionary` | Turbulence model selection | — | ✓ |
| `U` | `0/` | `volVectorField` | Velocity initial & boundary conditions | ✓ | ✓ |
| `p` | `0/` | `volScalarField` | Pressure initial & boundary conditions | ✓ | ✓ |
| `k` | `0/` | `volScalarField` | Turbulent kinetic energy IC/BC | — | ✓ |
| `epsilon` | `0/` | `volScalarField` | Turbulent dissipation rate IC/BC | — | ✓ |
| `nut` | `0/` | `volScalarField` | Turbulent viscosity IC/BC | — | ✓ |
| `omega` | `0/` | `volScalarField` | Specific dissipation rate IC/BC (k-ω models) | — | — |
| `nuTilda` | `0/` | `volScalarField` | Modified viscosity IC/BC (Spalart-Allmaras) | — | — |

---

## Quick-Reference: How All Dictionaries Connect

```
 ┌─────────────────────────────────────────────────────────────────────────────┐
 │                                                                            │
 │   controlDict ─────▶ selects SOLVER (icoFoam / simpleFoam / ...)          │
 │        │                    │                                              │
 │        │                    ▼                                              │
 │        │         Solver reads ALL other dictionaries:                      │
 │        │                                                                   │
 │        │    ┌──────────────────┐  ┌───────────────────┐  ┌────────────┐   │
 │        │    │  fvSchemes       │  │  fvSolution       │  │  0/ files  │   │
 │        │    │  HOW to          │  │  HOW ACCURATELY   │  │  WHAT the  │   │
 │        │    │  discretize      │  │  to solve the     │  │  initial & │   │
 │        │    │  the PDEs        │  │  linear systems   │  │  boundary  │   │
 │        │    └──────────────────┘  └───────────────────┘  │  values are│   │
 │        │                                                  └────────────┘   │
 │        │    ┌──────────────────┐  ┌───────────────────┐                   │
 │        │    │transportProperties│ │turbulenceProperties│                   │
 │        │    │  WHAT the fluid  │  │  WHICH turbulence │                   │
 │        │    │  properties are  │  │  model to use     │                   │
 │        │    └──────────────────┘  └───────────────────┘                   │
 │        │                                                                   │
 │        └─────▶ WHEN to write, HOW LONG to run, WHAT format to use         │
 │                                                                            │
 └─────────────────────────────────────────────────────────────────────────────┘
```

---

*All code examples in this document are taken directly from the project files in this
repository. Cavity = `projects/01_lid_driven_cavity`, NACA = `projects/03_naca_airfoil_analysis`.*
