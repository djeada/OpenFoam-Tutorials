# `icoFoam` Solver — Complete Analysis

> **Cross-references:**
> - [01 — Short Intro to CFD](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md)
> - [02 — OpenFOAM Cases](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md)
> - [03 — OpenFOAM Dictionaries](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md)
> - [04 — Meshing](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/04_meshing.md)
> - [05 — Boundary Conditions](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md)
> - [08 — CFL Number](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md)
> - [09 — Linear Solvers](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/09_linear_solvers.md)

---

## 1. Introduction — What icoFoam Solves

`icoFoam` is the simplest transient CFD solver in OpenFOAM. It solves the
**incompressible Navier-Stokes equations** for **laminar**, **Newtonian** fluids using
the **PISO** (Pressure-Implicit with Splitting of Operators) algorithm.

```
┌──────────────────────────────────────────────────────────────────────┐
│                      icoFoam AT A GLANCE                             │
├──────────────────────────────────────────────────────────────────────┤
│  Flow regime  :  Incompressible  (constant density ρ)               │
│  Viscosity    :  Newtonian       (constant viscosity μ)             │
│  Turbulence   :  Laminar ONLY    (no turbulence model)              │
│  Time         :  Transient       (time-dependent)                   │
│  Algorithm    :  PISO            (pressure-velocity coupling)       │
│  Typical Re   :  < ~2000         (before turbulence onset)          │
├──────────────────────────────────────────────────────────────────────┤
│  Source file  :  applications/solvers/incompressible/icoFoam/       │
│                  icoFoam.C                                           │
└──────────────────────────────────────────────────────────────────────┘
```

**What do these terms mean physically?**

- **Incompressible** — Density does not change with pressure. Valid for liquids and
  low-speed gas flows (Mach < 0.3). This eliminates the energy equation entirely.
- **Laminar** — The flow is smooth and orderly. No turbulent eddies. No turbulence
  model is needed or included. See [06 — Turbulence Models](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md)
  for what to use when turbulence matters.
- **Newtonian** — Shear stress is linearly proportional to strain rate. The viscosity
  $\nu$ is a constant. Water, air, and most common fluids are Newtonian. Blood, polymers,
  and ketchup are *not*.

**When to use icoFoam:** Teaching, validation cases, low-Re flows (pipe flow at
Re < 2000, Stokes flow, lid-driven cavity, creeping flows). For anything more complex,
see Section 8 (Limitations) and Section 9 (Related Solvers).

---

## 2. The Incompressible Navier-Stokes Equations

icoFoam solves two coupled equations — the **continuity** (mass conservation) and
**momentum** equations for an incompressible Newtonian fluid:

### Continuity Equation (mass conservation)

$$\nabla \cdot u = 0$$

This says: the velocity field is **divergence-free**. Fluid is neither created nor
destroyed. For an incompressible fluid, this replaces the full mass conservation
equation since $\rho$ = constant.

### Momentum Equation

$$\frac{\partial u}{\partial t} + \nabla \cdot (uu) - \nabla \cdot (\nu \nabla u) = -\frac{1}{\rho} \nabla p$$

In OpenFOAM's "kinematic" form (dividing through by $\rho$), pressure `p` has units
of $\text{m}^{2}/\text{s}^{2}$ (kinematic pressure = $p/\rho$), so:

$$\frac{\partial u}{\partial t} + \nabla \cdot (uu) - \nabla \cdot (\nu \nabla u) = -\nabla p$$

### Physical Meaning of Each Term

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│    ∂u/∂t          ∇·(uu)          ν ∇²u             −∇p                    │
│   ────────      ──────────      ──────────        ──────────               │
│   TIME           CONVECTION      DIFFUSION         PRESSURE                │
│   DERIVATIVE     (inertia)       (viscous)         GRADIENT                │
│                                                                             │
│   How fast       Momentum        Viscous           Pressure                │
│   velocity       carried by      friction          force that              │
│   changes        the flow        smoothing         drives/                 │
│   at a point     itself          the flow          decelerates             │
│                                                    the flow                │
│                                                                             │
│   fvm::ddt(U)    fvm::div(       fvm::laplacian(   fvc::grad(p)           │
│                   phi,U)          nu,U)                                    │
│                                                                             │
│   [Implicit]     [Implicit]      [Implicit]        [Explicit]              │
│   → matrix A     → matrix A      → matrix A        → source b             │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Why no time derivative for pressure?**
In incompressible flow, pressure is *not* a thermodynamic variable — it is a
**constraint field** that enforces $\nabla \cdot u = 0$ at every time step. There is no
pressure wave propagation (speed of sound $\to \infty$). Pressure adjusts
instantaneously to maintain divergence-free velocity. This is why we need a
special algorithm (PISO) to couple p and u.

---

## 3. The PISO Algorithm — Full Flowchart

PISO (Pressure-Implicit with Splitting of Operators) was developed by Issa (1986).
It is a **non-iterative** predictor-corrector approach for transient flows. Within
each time step, it performs ONE momentum prediction followed by multiple pressure
corrections (typically 2).

### Complete Algorithm Flowchart

```
┌═══════════════════════════════════════════════════════════════════════════┐
║                         icoFoam MAIN ALGORITHM                          ║
╠═══════════════════════════════════════════════════════════════════════════╣
║                                                                         ║
║   ┌─────────────────────────────────────────────────────────────────┐   ║
║   │  INITIALIZATION                                                 │   ║
║   │  • createTime    → set up time controls                        │   ║
║   │  • createMesh    → read polyMesh (see note 04)                 │   ║
║   │  • createFields  → read U, p, phi, nu from files               │   ║
║   │  • pisoControl   → read PISO settings from fvSolution          │   ║
║   └──────────────────────────┬──────────────────────────────────────┘   ║
║                              ▼                                          ║
║   ┌──────────────────────────────────────────────────────────────┐      ║
║   │ ┌────────────────────────────────────────────────────────┐   │      ║
║   │ │           TIME LOOP: while (runTime.loop())            │   │      ║
║   │ │                                                        │   │      ║
║   │ │  ┌──────────────────────────────────────────────────┐  │   │      ║
║   │ │  │  STEP 0: Courant Number Check                    │  │   │      ║
║   │ │  │  Co = |U| · Δt / Δx   (must be < 1 for PISO)    │  │   │      ║
║   │ │  │  See: note 08_cfl_number.md                      │  │   │      ║
║   │ │  └────────────────────┬─────────────────────────────┘  │   │      ║
║   │ │                       ▼                                │   │      ║
║   │ │  ┌──────────────────────────────────────────────────┐  │   │      ║
║   │ │  │  STEP 1: MOMENTUM PREDICTOR                      │  │   │      ║
║   │ │  │                                                  │  │   │      ║
║   │ │  │  Assemble:  A·U = H − ∇p_old                    │  │   │      ║
║   │ │  │                                                  │  │   │      ║
║   │ │  │  Where:                                          │  │   │      ║
║   │ │  │    A = diagonal coefficients of UEqn matrix      │  │   │      ║
║   │ │  │    H = off-diagonal + source terms               │  │   │      ║
║   │ │  │                                                  │  │   │      ║
║   │ │  │  Solve →  Intermediate velocity U*               │  │   │      ║
║   │ │  │  (U* does NOT satisfy continuity yet!)           │  │   │      ║
║   │ │  └────────────────────┬─────────────────────────────┘  │   │      ║
║   │ │                       ▼                                │   │      ║
║   │ │  ┌──────────────────────────────────────────────────┐  │   │      ║
║   │ │  │  STEP 2: PISO CORRECTOR LOOP                     │  │   │      ║
║   │ │  │  (repeats nCorrectors times, typically 2)        │  │   │      ║
║   │ │  │                                                  │  │   │      ║
║   │ │  │  ┌────────────────────────────────────────────┐  │  │   │      ║
║   │ │  │  │ 2a. Compute rAU = 1/A                      │  │  │   │      ║
║   │ │  │  │ 2b. Compute HbyA = H/A                     │  │  │   │      ║
║   │ │  │  │ 2c. Compute face flux phiHbyA              │  │  │   │      ║
║   │ │  │  │ 2d. Adjust phi for global conservation     │  │  │   │      ║
║   │ │  │  └──────────────────┬─────────────────────────┘  │  │   │      ║
║   │ │  │                     ▼                            │  │   │      ║
║   │ │  │  ┌────────────────────────────────────────────┐  │  │   │      ║
║   │ │  │  │ 2e. PRESSURE EQUATION                      │  │  │   │      ║
║   │ │  │  │     ∇·(rAU · ∇p) = ∇·(phiHbyA)            │  │  │   │      ║
║   │ │  │  │                                            │  │  │   │      ║
║   │ │  │  │   (+ non-orthogonal correction sub-loop)   │  │  │   │      ║
║   │ │  │  │     Solve → corrected pressure p           │  │  │   │      ║
║   │ │  │  └──────────────────┬─────────────────────────┘  │  │   │      ║
║   │ │  │                     ▼                            │  │   │      ║
║   │ │  │  ┌────────────────────────────────────────────┐  │  │   │      ║
║   │ │  │  │ 2f. VELOCITY CORRECTION                    │  │  │   │      ║
║   │ │  │  │     U = HbyA − rAU · ∇p                   │  │  │   │      ║
║   │ │  │  │     → divergence-free velocity U           │  │  │   │      ║
║   │ │  │  └──────────────────┬─────────────────────────┘  │  │   │      ║
║   │ │  │                     ▼                            │  │   │      ║
║   │ │  │  Check continuity error: ∇·U ≈ 0                │  │   │      ║
║   │ │  └──────────────────────────────────────────────────┘  │   │      ║
║   │ │                       ▼                                │   │      ║
║   │ │  ┌──────────────────────────────────────────────────┐  │   │      ║
║   │ │  │  STEP 3: Write results if write interval reached │  │   │      ║
║   │ │  │  Advance time: t = t + Δt                        │  │   │      ║
║   │ │  └──────────────────────────────────────────────────┘  │   │      ║
║   │ │                                                        │   │      ║
║   │ └──────────────────────────┬─────────────────────────────┘   │      ║
║   │                            │  Loop back if t < endTime       │      ║
║   │                            └─────────────────────────────────│      ║
║   └──────────────────────────────────────────────────────────────┘      ║
║                                                                         ║
╚═══════════════════════════════════════════════════════════════════════════╝
```

### Why PISO Works for Transient Flows

PISO is efficient for transient problems because it performs **no outer iterations**
within a time step. It relies on a small time step (CFL < 1, see
[08 — CFL Number](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md)) so that the change per step is small enough
that one momentum prediction + a few corrections suffice. This makes it faster
per-time-step than iterative methods like SIMPLE, but it demands small $\Delta t$.

---

## 4. Source Code Analysis — Deep Dive

The full source lives in `applications/solvers/incompressible/icoFoam/icoFoam.C`.
Below we analyze every section with line-by-line annotations.

### 4a. Headers and Includes

```cpp
#include "fvCFD.H"       // Master include: finite volume method classes,
                          //   field types (volScalarField, volVectorField, ...),
                          //   fvm:: and fvc:: namespaces, mesh, time, I/O
                          //   Basically: everything you need for a FV solver.

#include "pisoControl.H"  // PISO loop controller — reads nCorrectors and
                          //   nNonOrthogonalCorrectors from fvSolution,
                          //   provides .correct(), .correctNonOrthogonal(),
                          //   .momentumPredictor(), etc.
```

`fvCFD.H` is the umbrella header. It pulls in:

```
fvCFD.H
 ├── fvMesh.H          → the computational mesh
 ├── volFields.H       → cell-centered fields (volScalarField, volVectorField)
 ├── surfaceFields.H   → face-centered fields (surfaceScalarField)
 ├── fvm.H             → implicit discretization operators (fvm::ddt, fvm::div, ...)
 ├── fvc.H             → explicit discretization operators (fvc::grad, fvc::div, ...)
 ├── fvMatrices.H      → linear system classes (fvScalarMatrix, fvVectorMatrix)
 ├── Time.H            → time management
 └── argList.H         → command-line argument parsing
```

### 4b. Initialization

```cpp
int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"    // Parse command-line args, set up case path
    #include "createTime.H"          // Create the runTime object from controlDict
    #include "createMesh.H"          // Read polyMesh/ → create fvMesh object

    pisoControl piso(mesh);          // Read PISO settings from fvSolution dict

    #include "createFields.H"        // Read U, p, phi, nu from case files
    #include "initContinuityErrs.H"  // Initialize cumulative continuity error tracker
```

**What `createFields.H` does** (for icoFoam):

```cpp
// Read p (pressure) from 0/p
volScalarField p(...);

// Read U (velocity) from 0/U
volVectorField U(...);

// Create phi (face flux) = interpolate(U) · face area vectors
surfaceScalarField phi("phi", fvc::flux(U));

// Read nu (kinematic viscosity) from constant/transportProperties
dimensionedScalar nu(...);
```

**How fields map to the mesh:**

```
┌──────────────────────────────────────────────────────────────────┐
│                     MESH → FIELD MAPPING                         │
│                                                                  │
│           cell centers                face centers               │
│          ┌────┬────┬────┐          ┌──┬──┬──┬──┬──┐             │
│          │ p₁ │ p₂ │ p₃ │          │f₁│f₂│f₃│f₄│f₅│             │
│          │ U₁ │ U₂ │ U₃ │          │φ₁│φ₂│φ₃│φ₄│φ₅│             │
│          └────┴────┴────┘          └──┴──┴──┴──┴──┘             │
│                                                                  │
│   volScalarField  p  → one scalar per cell    (Nc values)        │
│   volVectorField  U  → one vector per cell    (Nc values)        │
│   surfaceScalarField phi → one scalar per face (Nf values)       │
│                                                                  │
│   phi = U · Sf  (volumetric flux through each face)              │
│   Sf = face area vector (outward normal × area)                  │
└──────────────────────────────────────────────────────────────────┘
```

### 4c. Time Loop

```cpp
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())                       // Advance t by Δt each iteration
    {                                            // Loop until t >= endTime
        Info<< "Time = " << runTime.timeName()   // Print current time
            << nl << endl;

        #include "CourantNo.H"                   // Compute and print Courant number
```

### 4d. Courant Number Check

`CourantNo.H` computes the CFL number (see [08 — CFL Number](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md)):

$$Co = \max_f \frac{|\phi_f|}{|d| \cdot |S_f|} \cdot \Delta t$$

where:
- $\phi_f$ = face flux (phi)
- $d$ = distance between neighboring cell centers
- $S_f$ = face area

For PISO to be stable, **Co must be < 1** (typically $\leq 0.5$). This is the main
restriction of the PISO algorithm — it requires small time steps.

### 4e. Momentum Predictor — Deep Dive

```cpp
        // Assemble the momentum equation matrix
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)              // ∂U/∂t   — implicit time derivative
          + fvm::div(phi, U)         // ∇·(φU)  — implicit convection
          - fvm::laplacian(nu, U)    // ν∇²U    — implicit diffusion
        );
```

**What `fvm::` vs `fvc::` means:**

```
┌──────────────────────────────────────────────────────────────────┐
│               fvm:: (IMPLICIT)  vs  fvc:: (EXPLICIT)             │
├──────────────────────────────────────────────────────────────────┤
│                                                                  │
│   The discretized equation is:   [A] · {x} = {b}                │
│                                                                  │
│   fvm:: operators   → contribute to MATRIX [A]                   │
│     fvm::ddt(U)        → time derivative term in A               │
│     fvm::div(phi,U)    → convection term in A                    │
│     fvm::laplacian()   → diffusion term in A                     │
│                                                                  │
│   fvc:: operators   → contribute to SOURCE VECTOR {b}            │
│     fvc::grad(p)       → pressure gradient (known from old p)    │
│     fvc::div(phi)      → explicit divergence                     │
│     fvc::flux(U)       → compute face flux from volume field     │
│                                                                  │
│   WHY?  Implicit terms (fvm) are solved simultaneously → stable  │
│         Explicit terms (fvc) use known values → may be unstable  │
│         Rule: put the "hard" terms (convection, diffusion) in A  │
│               and the "easy" terms (pressure grad) in b          │
└──────────────────────────────────────────────────────────────────┘
```

```cpp
        // Optionally solve for U* (the intermediate velocity)
        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
            //            ^^^^^^^^^^^^^^
            //   Pressure gradient uses OLD p (from previous time step)
            //   This is the "prediction" — U* won't satisfy ∇·U = 0
        }
```

**How the matrix maps to the mesh:**

```
         cell 0      cell 1      cell 2
        ┌────────┬──────────┬──────────┐
        │  a₀₀   │   a₀₁    │          │
        │        │          │          │  UEqn matrix [A]:
        ├────────┼──────────┼──────────┤
        │  a₁₀   │   a₁₁    │   a₁₂    │  a_ii = diagonal  → UEqn.A()
        │        │          │          │  a_ij = off-diag  → part of UEqn.H()
        ├────────┼──────────┼──────────┤
        │        │   a₂₁    │   a₂₂    │  Each row = one cell
        │        │          │          │  Neighbors = off-diagonal entries
        └────────┴──────────┴──────────┘
                                          Non-zero pattern follows mesh topology
```

### 4f. PISO Corrector Loop — Step by Step

This is the heart of icoFoam. Each sub-step is annotated:

```cpp
        // --- PISO loop (repeats nCorrectors times)
        while (piso.correct())
        {
```

#### Step 2a–2b: Compute rAU and HbyA

```cpp
            volScalarField rAU(1.0/UEqn.A());
            //  rAU = 1/A = reciprocal of diagonal coefficients
            //  This is a scalar field: one value per cell
            //  Used to scale the pressure equation

            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            //  H = "off-diagonal + source" part of the momentum equation
            //  HbyA = H/A = what velocity would be WITHOUT pressure gradient
            //  constrainHbyA: applies boundary constraints
```

**Decomposition of the momentum equation:**

$$A \cdot U = H - \nabla p$$

Therefore:

$$U = \frac{H}{A} - \frac{1}{A} \cdot \nabla p$$
$$U = HbyA - rAU \cdot \nabla p$$

This is the KEY IDENTITY that connects velocity and pressure!

#### Step 2c–2d: Flux calculation

```cpp
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)                              // Face flux from HbyA
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)   // Temporal correction
            );

            adjustPhi(phiHbyA, U, p);
            //  Adjust fluxes for global mass conservation
            //  (ensures total mass in = total mass out for closed domains)

            constrainPressure(p, U, phiHbyA, rAU);
            //  Update pressure BCs for flux consistency
```

#### Step 2e: Pressure Equation

```cpp
            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                    //  ∇·(rAU · ∇p) = ∇·(phiHbyA)
                    //
                    //  This is a POISSON EQUATION for pressure!
                    //  Derived by taking ∇· of the velocity correction:
                    //    ∇·U = ∇·HbyA − ∇·(rAU · ∇p)
                    //    Since ∇·U = 0 (continuity):
                    //    ∇·(rAU · ∇p) = ∇·HbyA
                );

                pEqn.setReference(pRefCell, pRefValue);
                //  Fix pressure at one point (since p equation is singular —
                //  only ∇p appears in momentum, so p is determined up to a constant)

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                //  Use the "p" solver for intermediate iterations,
                //  "pFinal" solver for the last PISO corrector step
                //  (see fvSolution — pFinal has relTol = 0 for tighter solve)

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                    //  Correct face fluxes so they satisfy continuity exactly
                    //  pEqn.flux() = face flux from the pressure gradient
                }
            }
```

**Why non-orthogonal corrections?**

```
┌────────────────────────────────────────────────────────────────────────┐
│          ORTHOGONAL MESH              NON-ORTHOGONAL MESH             │
│                                                                        │
│       ┌────┬────┬────┐             ┌────┬─────┬────┐                  │
│       │    │    │    │            /    /     /    /                    │
│       ├────┼────┼────┤           /────/─────/────/                     │
│       │    │    │    │          /    /     /    /                      │
│       └────┴────┴────┘         /────/─────/────/                       │
│                                                                        │
│   d ∥ Sf  → perfect            d ∦ Sf → error in gradient             │
│   No correction needed         Need extra iterations to correct        │
│                                                                        │
│   d = vector connecting        Sf = face area normal vector            │
│       cell centers                                                     │
│                                                                        │
│   nNonOrthogonalCorrectors = 0     nNonOrthogonalCorrectors = 2       │
│   (cavity case)                    (elbow case — non-orthogonal mesh)  │
└────────────────────────────────────────────────────────────────────────┘
```

#### Step 2f: Velocity Correction and Continuity Check

```cpp
            #include "continuityErrs.H"
            //  Compute and report: sum(|∇·phi|) — should be ~machine epsilon

            U = HbyA - rAU*fvc::grad(p);
            //  THE velocity correction equation:
            //  U = H/A − (1/A)·∇p_new
            //  Now U satisfies BOTH momentum and continuity (approximately)

            U.correctBoundaryConditions();
            //  Re-apply boundary conditions to the corrected velocity
        }
```

### 4g. Output Writing

```cpp
        runTime.write();
        //  Write U, p, phi to the current time directory
        //  (only if the write interval is reached — see controlDict)

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}
```

---

## 5. PISO vs SIMPLE vs PIMPLE — Algorithm Comparison

OpenFOAM has three main pressure-velocity coupling algorithms. Understanding when
to use each is critical.

### SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)

```
┌─────────────────────────────────────────────────────────────────┐
│                SIMPLE Algorithm (STEADY-STATE)                   │
│                                                                  │
│   ┌───────────────────────────────────────────────────────┐     │
│   │  1. Solve momentum with under-relaxation              │     │
│   │     A·U* = H − ∇p_old                                │     │
│   │     U = αᵤ·U* + (1−αᵤ)·U_old   (relax velocity)     │     │
│   └──────────────────────┬────────────────────────────────┘     │
│                          ▼                                       │
│   ┌───────────────────────────────────────────────────────┐     │
│   │  2. Solve pressure correction equation                 │     │
│   │     ∇·(1/A · ∇p') = ∇·U*                              │     │
│   │     p = p_old + αₚ·p'          (relax pressure)       │     │
│   └──────────────────────┬────────────────────────────────┘     │
│                          ▼                                       │
│   ┌───────────────────────────────────────────────────────┐     │
│   │  3. Correct velocity                                   │     │
│   │     U = U* − (1/A)·∇p'                                │     │
│   └──────────────────────┬────────────────────────────────┘     │
│                          ▼                                       │
│          Check convergence → if not converged, go to 1           │
│                                                                  │
│   Key: ITERATES to convergence within EACH "time step"           │
│   Uses under-relaxation factors (αᵤ, αₚ) for stability          │
└─────────────────────────────────────────────────────────────────┘
```

### PIMPLE (merged PISO-SIMPLE)

```
┌─────────────────────────────────────────────────────────────────┐
│              PIMPLE Algorithm (TRANSIENT, large Δt)              │
│                                                                  │
│   ┌───────────────────────────────────────────────────────┐     │
│   │  OUTER LOOP (nOuterCorrectors times):                  │     │
│   │                                                        │     │
│   │   ┌─────────────────────────────────────────────────┐  │     │
│   │   │  1. Solve momentum (with under-relaxation)      │  │     │
│   │   └────────────────────┬────────────────────────────┘  │     │
│   │                        ▼                               │     │
│   │   ┌─────────────────────────────────────────────────┐  │     │
│   │   │  INNER LOOP (nCorrectors times):                │  │     │
│   │   │    2. Solve pressure                            │  │     │
│   │   │    3. Correct velocity                          │  │     │
│   │   │    (same as PISO correctors)                    │  │     │
│   │   └─────────────────────────────────────────────────┘  │     │
│   │                                                        │     │
│   └───────────────────────────────────────────────────────┘     │
│                                                                  │
│   Key: Combines SIMPLE outer iterations + PISO inner corrections │
│   Allows CFL >> 1 (large time steps for transient flows)         │
└─────────────────────────────────────────────────────────────────┘
```

### Comparison Table

```
┌──────────────┬──────────────┬──────────────────┬──────────────────┐
│              │    PISO      │     SIMPLE       │    PIMPLE        │
├──────────────┼──────────────┼──────────────────┼──────────────────┤
│ Time         │ Transient    │ Steady-state     │ Transient        │
│ CFL limit    │ Co < 1       │ N/A (no real Δt) │ Co can be >> 1   │
│ Outer iters  │ 1 (none)     │ Many (converge)  │ nOuterCorrectors │
│ P-corrections│ nCorrectors  │ 1                │ nCorrectors      │
│ Relaxation   │ No           │ Yes (required)   │ Yes (outer loop) │
│ Cost/step    │ Low          │ High             │ Medium-High      │
│ Total cost   │ Many steps   │ Fewer steps      │ Fewer steps      │
│ OpenFOAM     │ icoFoam,     │ simpleFoam       │ pimpleFoam       │
│ solver       │ pisoFoam     │                  │                  │
│ Best for     │ Small Δt,    │ Only care about  │ Large Δt,        │
│              │ accuracy in  │ final converged  │ practical        │
│              │ time         │ solution         │ transient sims   │
└──────────────┴──────────────┴──────────────────┴──────────────────┘
```

---

## 6. icoFoam in Practice — Real Project Configuration

Below are the ACTUAL configuration files from our tutorial projects that use icoFoam.
See [02 — OpenFOAM Cases](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md) for the general case structure and
[03 — OpenFOAM Dictionaries](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md) for dictionary syntax.

### 6.1 Lid-Driven Cavity (`projects/01_lid_driven_cavity/`)

**The lid-driven cavity** is the classic icoFoam validation case. A square box with
the top wall sliding at constant velocity drives a recirculating flow.

```
    U = (1, 0, 0) ───────────────────→   movingWall (fixedValue)
    ┌─────────────────────────────────┐
    │  ↗  →   →   →   →   →   ↘     │
    │  ↑                        ↓     │
    │  ↑      PRIMARY           ↓     │   fixedWalls
    │  ↑       VORTEX           ↓     │   (noSlip)
    │  ↑                        ↓     │
    │  ↑                        ↓     │
    │  ↖  ←   ←   ←   ←   ←   ↙     │
    │                                 │
    │  ↘ secondary    secondary ↙     │
    │    vortex          vortex       │
    └─────────────────────────────────┘
              fixedWalls (noSlip)

    frontAndBack faces: type "empty" (2D simulation)
```

#### controlDict — Solver and Time Settings

```c
// From: projects/01_lid_driven_cavity/system/controlDict

application     icoFoam;          // ← This tells OpenFOAM which solver to run

startFrom       startTime;
startTime       0;                // Start at t = 0
stopAt          endTime;
endTime         0.5;              // Run for 0.5 seconds

deltaT          0.001;            // Time step = 1 ms
                                  // For a 1 m/s lid on a 0.1m cavity:
                                  // Co ≈ 1.0 × 0.001 / 0.005 = 0.2 ✓ (< 1)

writeControl    timeStep;
writeInterval   1;                // Write EVERY time step (500 output dirs!)
purgeWrite      0;                // Keep all output directories

writeFormat     ascii;            // Human-readable output
writePrecision  6;
writeCompression off;
runTimeModifiable true;           // Can change settings while running
```

#### fvSolution — Linear Solvers and PISO Settings

```c
// From: projects/01_lid_driven_cavity/system/fvSolution

solvers
{
    p                              // Pressure solver
    {
        solver          PCG;       // Preconditioned Conjugate Gradient
                                   // (symmetric matrix → use CG)
        preconditioner  DIC;       // Diagonal Incomplete Cholesky
        tolerance       1e-06;     // Absolute convergence tolerance
        relTol          0.05;      // Relative tolerance (5% reduction)
    }                              // See: note 09_linear_solvers.md

    pFinal                         // Pressure solver for LAST PISO corrector
    {
        $p;                        // Inherit all settings from 'p' above
        relTol          0;         // But solve to FULL tolerance (no relative shortcut)
    }                              // This ensures the final correction is accurate

    U                              // Velocity solver
    {
        solver          smoothSolver;      // Iterative smoother
        smoother        symGaussSeidel;    // Symmetric Gauss-Seidel
        tolerance       1e-05;
        relTol          0;
    }
}

PISO
{
    nCorrectors     2;                     // 2 pressure-velocity correction steps
                                           // (standard for PISO — rarely change this)
    nNonOrthogonalCorrectors 0;            // 0 because cavity mesh is orthogonal
    pRefCell        0;                     // Fix pressure in cell 0
    pRefValue       0;                     // Set reference pressure to 0
}
```

#### fvSchemes — Discretization Schemes

```c
// From: projects/01_lid_driven_cavity/system/fvSchemes

ddtSchemes
{
    default         Euler;             // First-order implicit time stepping
}                                      // (simple, stable, first-order accurate)

gradSchemes
{
    default         Gauss linear;      // Standard gradient: Gauss theorem + linear interp
    grad(p)         Gauss linear;      // Pressure gradient (same as default)
    grad(U)         cellLimited Gauss linear 1;
    //              ^^^^^^^^^^^^^^^^^^^^^^^^^
    //  Cell-limited gradient for velocity — prevents overshoots near
    //  boundaries. The "1" means full limiting. Important for stability
    //  with linearUpwind convection scheme.
}

divSchemes
{
    default         none;              // No default — force explicit specification
    div(phi,U)      Gauss linearUpwind grad(U);
    //              ^^^^^^^^^^^^^^^^^^^^^^^^
    //  Convection: Gauss integration + linearUpwind interpolation
    //  Second-order accurate, bounded (uses grad(U) for upwind direction)
    //  More accurate than pure upwind, more stable than pure linear
}

laplacianSchemes
{
    default         Gauss linear corrected;
    //              ^^^^^^^^^^^^^^^^^^^^^^^^
    //  Diffusion: Gauss integration + linear interpolation + non-orthogonal correction
    //  "corrected" adds explicit correction for non-orthogonal meshes
}

interpolationSchemes
{
    default         linear;            // Linear (central) interpolation cell→face
}

snGradSchemes
{
    default         corrected;         // Surface-normal gradient with correction
}
```

#### Boundary Conditions — U (velocity)

```c
// From: projects/01_lid_driven_cavity/0/U

dimensions      [0 1 -1 0 0 0 0];     // m/s  (length/time)

internalField   uniform (0 0 0);       // Start with zero velocity everywhere

boundaryField
{
    movingWall                         // TOP WALL — the "lid"
    {
        type            fixedValue;
        value           uniform (1 0 0);   // Slides at 1 m/s in x-direction
    }

    fixedWalls                         // BOTTOM + SIDES — stationary walls
    {
        type            noSlip;        // U = (0,0,0) — no-slip condition
    }                                  // See: note 05_boundary_conditions.md

    frontAndBack                       // FRONT + BACK faces
    {
        type            empty;         // 2D simulation — ignore z-direction
    }
}
```

#### Boundary Conditions — p (pressure)

```c
// From: projects/01_lid_driven_cavity/0/p

dimensions      [0 2 -2 0 0 0 0];     // m²/s²  (kinematic pressure = p/ρ)
                                       // NOT Pascals! OpenFOAM incompressible
                                       // solvers use p/ρ.

internalField   uniform 0;            // Initial pressure = 0 everywhere

boundaryField
{
    movingWall
    {
        type            zeroGradient;  // ∂p/∂n = 0 at walls
    }                                  // (standard for walls in incompressible flow)

    fixedWalls
    {
        type            zeroGradient;  // Same for fixed walls
    }

    frontAndBack
    {
        type            empty;         // 2D simulation
    }
}
```

**Boundary condition pairing rules for incompressible flow:**

```
┌──────────────────────────────────────────────────────────────────┐
│           BOUNDARY CONDITION PAIRING (U and p)                   │
├──────────────────────────────────────────────────────────────────┤
│                                                                  │
│   Boundary Type     │    U                │    p                 │
│   ──────────────────┼─────────────────────┼──────────────────────│
│   Wall              │  fixedValue/noSlip   │  zeroGradient       │
│   Inlet (velocity)  │  fixedValue          │  zeroGradient       │
│   Outlet (pressure) │  zeroGradient        │  fixedValue         │
│   Symmetry          │  symmetry            │  symmetry           │
│   2D (front/back)   │  empty               │  empty              │
│                                                                  │
│   RULE: You must fix either U or p on each boundary, not both!   │
│   (One Dirichlet + one Neumann per boundary)                     │
└──────────────────────────────────────────────────────────────────┘
```

### 6.2 Elbow Case (`projects/02_elbow/`)

The elbow case also uses icoFoam but with a **non-orthogonal mesh** (from Fluent),
two inlets, and one outlet.

```
                    velocity-inlet-6
                    U = (0, 3, 0)
                         │
                         ▼
         ┌───────────────┐
         │               │
         │               │
         │    ELBOW      │
         │   (mixing     │
         │    region)    │──────→  pressure-outlet-7
         │               │         p = 0
         │               │
         └───────────────┘
                 ▲
                 │
         velocity-inlet-5
         U = (1, 0, 0)

    Walls: wall-4, wall-8 (noSlip)
    2D:    frontAndBackPlanes (empty)
```

Key differences from the cavity case:

```c
// From: projects/02_elbow/system/controldict
application     icoFoam;
endTime         75;        // Much longer simulation (75s vs 0.5s)
deltaT          0.05;      // Larger time step (non-orthogonal mesh is coarser)
writeInterval   20;        // Write every 20 time steps

// From: projects/02_elbow/system/fvsolution — PISO settings
PISO
{
    nCorrectors     2;
    nNonOrthogonalCorrectors 2;   // ← NON-ZERO! Mesh is non-orthogonal
}                                 //   (Fluent mesh import → skewed cells)
```

```c
// From: projects/02_elbow/system/fvschemes — different convection scheme
divSchemes
{
    default         none;
    div(phi,U)      Gauss limitedLinearV 1;
    //              ^^^^^^^^^^^^^^^^^^^^^^
    //  limitedLinearV: vector-limited linear scheme
    //  More robust than linearUpwind for non-orthogonal meshes
    //  The "V" variant limits each vector component independently
}
```

---

## 7. Key OpenFOAM Programming Concepts

Understanding these concepts is essential for reading any OpenFOAM solver source.

### 7.1 fvm:: vs fvc:: — Implicit vs Explicit

```
┌────────────────────────────────────────────────────────────────┐
│                                                                │
│   The linear system:  [A] · {φ} = {b}                         │
│                                                                │
│   ┌──────────────────────────────────────────────────────┐     │
│   │  fvm::ddt(U)           → adds to [A] and {b}        │     │
│   │  fvm::div(phi, U)      → adds to [A] and {b}        │     │
│   │  fvm::laplacian(nu, U) → adds to [A] and {b}        │     │
│   │                                                      │     │
│   │  These create an fvMatrix (system to be solved)      │     │
│   └──────────────────────────────────────────────────────┘     │
│                                                                │
│   ┌──────────────────────────────────────────────────────┐     │
│   │  fvc::grad(p)          → returns a field (vector)    │     │
│   │  fvc::div(phi)         → returns a field (scalar)    │     │
│   │  fvc::laplacian(nu, U) → returns a field (vector)    │     │
│   │                                                      │     │
│   │  These evaluate immediately using CURRENT field vals  │     │
│   └──────────────────────────────────────────────────────┘     │
│                                                                │
│   Example in icoFoam:                                          │
│     fvm::ddt(U) + fvm::div(phi,U) - fvm::laplacian(nu,U)      │
│     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       │
│     All implicit → builds the UEqn matrix                      │
│                                                                │
│     solve(UEqn == -fvc::grad(p))                               │
│                    ^^^^^^^^^^^^                                 │
│     Explicit → computed from current p, goes into source {b}   │
│                                                                │
└────────────────────────────────────────────────────────────────┘
```

### 7.2 Field Types

```
┌────────────────────────────────────────────────────────────────────────┐
│                     OpenFOAM FIELD TYPES                               │
├────────────────────────────────────────────────────────────────────────┤
│                                                                        │
│  TYPE                      │ LOCATION  │ EXAMPLE     │ USED FOR        │
│  ─────────────────────────-┼───────────┼─────────────┼─────────────────│
│  volScalarField            │ cell      │ p (press.)  │ One scalar/cell │
│  volVectorField            │ cell      │ U (veloc.)  │ One vector/cell │
│  volTensorField            │ cell      │ stress      │ One tensor/cell │
│  surfaceScalarField        │ face      │ phi (flux)  │ One scalar/face │
│  surfaceVectorField        │ face      │ Sf (area)   │ One vector/face │
│                                                                        │
│  "vol"     → values stored at cell CENTERS                             │
│  "surface" → values stored at cell FACES                               │
│                                                                        │
│  In icoFoam:                                                           │
│    p   : volScalarField        → pressure at each cell center          │
│    U   : volVectorField        → velocity at each cell center          │
│    phi : surfaceScalarField    → volume flux through each face         │
│    nu  : dimensionedScalar     → single constant (kinematic viscosity) │
│                                                                        │
│  Dimensions format: [kg  m  s  K  mol  A  cd]                         │
│    p:   [0 2 -2 0 0 0 0]  →  m²/s²  (kinematic pressure)             │
│    U:   [0 1 -1 0 0 0 0]  →  m/s    (velocity)                       │
│    phi: [0 3 -1 0 0 0 0]  →  m³/s   (volumetric flow rate)           │
│    nu:  [0 2 -1 0 0 0 0]  →  m²/s   (kinematic viscosity)            │
└────────────────────────────────────────────────────────────────────────┘
```

### 7.3 Extending Solvers

icoFoam is the starting point for writing custom solvers. Common extensions:

1. **Add turbulence** → becomes pisoFoam (add `turbulence->correct()`)
2. **Add energy equation** → add `volScalarField T` + `fvm::ddt(T) + fvm::div(phi,T) - fvm::laplacian(DT,T)`
3. **Add scalar transport** → same pattern as energy, different field
4. **Non-Newtonian** → replace constant `nu` with a viscosity model

The structure is always the same:
```
1. Read fields
2. Time loop
3. Assemble matrix (fvm:: terms)
4. Add source (fvc:: terms)
5. Solve
6. Correct (PISO/SIMPLE/PIMPLE)
7. Write
```

---

## 8. Limitations of icoFoam

```
┌──────────────────────────────────────────────────────────────────────┐
│                    icoFoam LIMITATIONS                                │
├───────────────────┬──────────────────────────────────────────────────┤
│ LIMITATION        │ WHAT TO USE INSTEAD                              │
├───────────────────┼──────────────────────────────────────────────────┤
│ Laminar only      │ pisoFoam (adds RANS/LES turbulence models)      │
│ (no turbulence)   │ See: note 06_turbulence_models.md                │
├───────────────────┼──────────────────────────────────────────────────┤
│ Newtonian only    │ nonNewtonianIcoFoam (power-law, Bird-Carreau,   │
│ (constant ν)      │ Cross model viscosities)                         │
├───────────────────┼──────────────────────────────────────────────────┤
│ Incompressible    │ rhoPimpleFoam / sonicFoam (compressible)        │
│ (constant ρ)      │ buoyantPimpleFoam (buoyancy-driven)             │
├───────────────────┼──────────────────────────────────────────────────┤
│ CFL < 1 required  │ pimpleFoam (allows CFL >> 1 with outer iters)  │
│ (small Δt)        │ See: note 08_cfl_number.md                       │
├───────────────────┼──────────────────────────────────────────────────┤
│ No steady-state   │ simpleFoam (SIMPLE algorithm, no time stepping) │
│ capability        │                                                  │
├───────────────────┼──────────────────────────────────────────────────┤
│ Single phase      │ interFoam (VOF two-phase), multiphaseInterFoam  │
│ only              │                                                  │
└───────────────────┴──────────────────────────────────────────────────┘
```

---

## 9. Related Solvers — Quick Reference

```
┌───────────────────────┬────────────┬────────────┬──────────────┬───────────┐
│ Solver                │ Steady /   │ Turbulence │ Compressible │ Algorithm │
│                       │ Transient  │            │              │           │
├───────────────────────┼────────────┼────────────┼──────────────┼───────────┤
│ icoFoam               │ Transient  │ No         │ No           │ PISO      │
│ pisoFoam              │ Transient  │ Yes        │ No           │ PISO      │
│ pimpleFoam            │ Transient  │ Yes        │ No           │ PIMPLE    │
│ simpleFoam            │ Steady     │ Yes        │ No           │ SIMPLE    │
│ nonNewtonianIcoFoam   │ Transient  │ No         │ No           │ PISO      │
│ adjointShapeOptFoam   │ Steady     │ No         │ No           │ SIMPLE    │
│ rhoPimpleFoam         │ Transient  │ Yes        │ Yes          │ PIMPLE    │
│ rhoSimpleFoam         │ Steady     │ Yes        │ Yes          │ SIMPLE    │
│ sonicFoam             │ Transient  │ Yes        │ Yes          │ PISO      │
│ interFoam             │ Transient  │ Yes        │ No (2-phase) │ PIMPLE    │
└───────────────────────┴────────────┴────────────┴──────────────┴───────────┘
```

**Decision flowchart:**

```
                     Is the flow steady-state?
                    /                          \
                 YES                            NO (transient)
                  │                              │
            simpleFoam                   Is turbulence important?
            (or rhoSimpleFoam             /                    \
             if compressible)           NO                     YES
                                         │                      │
                                    Is CFL < 1?           Is CFL < 1?
                                    /         \           /         \
                                  YES          NO       YES          NO
                                   │            │        │            │
                               icoFoam     pimpleFoam  pisoFoam  pimpleFoam
                               (laminar,   (laminar,   (RANS/    (RANS/LES,
                                low Re)     large Δt)   LES)      large Δt)
```

---

## 10. Summary — The Complete icoFoam Picture

```
┌═══════════════════════════════════════════════════════════════════════════════┐
║                                                                             ║
║    icoFoam: Transient + Incompressible + Laminar + Newtonian + PISO         ║
║                                                                             ║
║    EQUATIONS:    ∇·u = 0   and   ∂u/∂t + ∇·(uu) − ν∇²u = −∇p              ║
║                                                                             ║
║    ALGORITHM:    1. Predict U*  (momentum with old p)                       ║
║                  2. Solve pressure Poisson equation                         ║
║                  3. Correct U from new p                                    ║
║                  4. Repeat 2-3 (nCorrectors times)                         ║
║                                                                             ║
║    REQUIRES:     Co < 1  (small time steps)                                 ║
║                  Laminar flow (Re < ~2000)                                  ║
║                  Constant density and viscosity                             ║
║                                                                             ║
║    FILES:        controlDict  → application icoFoam; deltaT; endTime        ║
║                  fvSolution   → PISO { nCorrectors 2; }                     ║
║                  fvSchemes    → ddtSchemes, divSchemes, laplacianSchemes     ║
║                  0/U, 0/p    → boundary conditions                          ║
║                  constant/transportProperties → nu                          ║
║                                                                             ║
║    UPGRADE TO:   pisoFoam (turbulence), pimpleFoam (large CFL),            ║
║                  simpleFoam (steady-state)                                  ║
║                                                                             ║
╚═══════════════════════════════════════════════════════════════════════════════╝
```
