# A Comprehensive Introduction to Computational Fluid Dynamics (CFD)

[All notes](../README.md#-notes-overview) · [Next: Note 02](02_openfoam_cases.md)

> **Version scope:** Classic OpenFOAM Foundation 6, matching the project launchers.
> Other releases/forks may need different dictionaries or executables; see
> [version guidance](../README.md#version-and-example-scope) and the
> [current reference cases](02_openfoam_cases.md#current-reference-cases).

<details>
<summary>Contents</summary>

- [What Is CFD and Why Does It Matter?](#what-is-cfd-and-why-does-it-matter)
  - [Where CFD Is Used](#where-cfd-is-used)
  - [A Brief History](#a-brief-history)
  - [CFD vs Other Approaches](#cfd-vs-other-approaches)
  - [CFD Tools — The Software Landscape](#cfd-tools--the-software-landscape)
- [The CFD Workflow](#the-cfd-workflow)
  - [Step-by-Step Breakdown](#step-by-step-breakdown)
- [Governing Equations](#governing-equations)
  - [The Control Volume Concept](#the-control-volume-concept)
  - [Conservation of Mass (Continuity Equation)](#conservation-of-mass-continuity-equation)
  - [Conservation of Momentum (Navier-Stokes Equations)](#conservation-of-momentum-navier-stokes-equations)
  - [Conservation of Energy](#conservation-of-energy)
  - [Turbulence and the Closure Problem](#turbulence-and-the-closure-problem)
- [Discretization — From Calculus to Algebra](#discretization--from-calculus-to-algebra)
  - [Mesh: Dividing Space Into Cells](#mesh-dividing-space-into-cells)
  - [Three Discretization Philosophies](#three-discretization-philosophies)
  - [Comparison Table](#comparison-table)
  - [From PDE to Linear System](#from-pde-to-linear-system)
- [5. The Finite Volume Method in Detail](#5-the-finite-volume-method-in-detail)
  - [Integral Form and Gauss's Divergence Theorem](#integral-form-and-gausss-divergence-theorem)
  - [Cell-Centered vs Vertex-Centered](#cell-centered-vs-vertex-centered)
  - [Discrete Operators in OpenFOAM](#discrete-operators-in-openfoam)
- [Numerical Schemes — Accuracy, Stability, and Trade-Offs](#numerical-schemes--accuracy-stability-and-trade-offs)
  - [Time Discretization](#time-discretization)
  - [Spatial Discretization (Convection Schemes)](#spatial-discretization-convection-schemes)
  - [Gradient Schemes](#gradient-schemes)
  - [Laplacian (Diffusion) Schemes](#laplacian-diffusion-schemes)
- [Pressure-Velocity Coupling](#pressure-velocity-coupling)
  - [The Problem](#the-problem)
  - [The Major Algorithms](#the-major-algorithms)
  - [SIMPLE vs PISO vs PIMPLE](#simple-vs-piso-vs-pimple)
- [Convergence and Accuracy](#convergence-and-accuracy)
  - [Residuals — How We Know the Solution Is Converging](#residuals--how-we-know-the-solution-is-converging)
  - [Grid Independence](#grid-independence)
  - [The CFL Number](#the-cfl-number)
  - [Sources of Error in CFD](#sources-of-error-in-cfd)
- [OpenFOAM in the CFD Landscape](#openfoam-in-the-cfd-landscape)
  - [Comparison with Commercial Solvers](#comparison-with-commercial-solvers)
  - [Why Open-Source Matters](#why-open-source-matters)
  - [OpenFOAM's Strengths and Limitations](#openfoams-strengths-and-limitations)
- [Quick Reference — Key Takeaways](#quick-reference--key-takeaways)
- [Where to Go Next](#where-to-go-next)

</details>

> **What this note covers:** The fundamental theory behind CFD — what it is, why it
> works, and how the mathematics connects to the OpenFOAM cases in this repository.
> This is the starting point before diving into OpenFOAM-specific notes.

## What Is CFD and Why Does It Matter?

Computational Fluid Dynamics (CFD) is a branch of fluid mechanics that uses numerical methods, algorithms, and computers to solve and analyze problems involving fluid flows. Instead of building expensive physical prototypes or running wind-tunnel experiments for every design iteration, engineers can simulate fluid behavior digitally — predicting pressures, velocities, temperatures, and forces with remarkable accuracy.

### Where CFD Is Used

| Industry       | Application Examples                                              |
|----------------|-------------------------------------------------------------------|
| **Aerospace**  | Wing design, turbine blade cooling, re-entry heating              |
| **Automotive** | External aerodynamics, engine cooling, HVAC cabin comfort         |
| **Energy**     | Wind turbine siting, nuclear reactor cooling, combustion modeling  |
| **Biomedical** | Blood flow in arteries, drug delivery, respiratory airflow        |
| **Marine**     | Hull resistance, propeller cavitation, wave loading               |
| **Civil Eng.** | Building wind loads, pollution dispersion, ventilation design      |
| **Electronics**| Chip thermal management, data center cooling                      |

### A Brief History

- **1920s–1940s:** Lewis Fry Richardson attempted the first numerical weather prediction by hand — dividing the atmosphere into cells and solving balance equations.
- **1960s:** Digital computers enabled practical CFD methods such as Harlow and Welch’s [Marker-and-Cell (MAC) method](https://doi.org/10.1063/1.1761178) at Los Alamos (1965).
- **1970s:** Patankar and Spalding introduced the [SIMPLE pressure-velocity coupling procedure](https://doi.org/10.1016/0017-9310(72)90054-3) (1972); it was not part of the 1960s MAC development.
- **1980s–1990s:** Commercial CFD software emerged (Fluent, CFX, STAR-CD). Turbulence modeling matured ($k$-$\varepsilon$, $k$-$\omega$ SST).
- **2004:** OpenFOAM was released as open-source, making CFD accessible to researchers and students worldwide.
- **Today:** GPU-accelerated solvers, machine-learning-augmented turbulence models, and cloud-based HPC make CFD faster and more accessible than ever.

> **💡 Tip:** In this repository, our projects span the range from simple laminar
> benchmarks (`projects/01_lid_driven_cavity/`) all the way to turbulent external
> aerodynamics (`projects/06_ahmed_body_aerodynamics/`). Each one follows the same
> fundamental CFD workflow described next.

### CFD vs Other Approaches

Engineers have three fundamental ways to study fluid flow. Each has a role, and the best
engineering practice combines them.

| Approach | Description | Strengths | Weaknesses |
|----------|-------------|-----------|------------|
| **Analytical** | Derive closed-form solutions from governing equations (e.g., Hagen-Poiseuille pipe flow, Blasius boundary layer). | Exact, instant, deep physical insight. | Only works for highly simplified geometries and assumptions (steady, laminar, simple shapes). Very few real problems have analytical solutions. |
| **Experimental** | Build physical models and measure flow quantities in wind tunnels, water channels, or test rigs. | Captures real physics including turbulence, transition, and complex interactions. Provides ground truth for validation. | Expensive, time-consuming, limited measurement resolution, scaling issues (Reynolds number matching), safety constraints (e.g., nuclear, biomedical). |
| **CFD (Numerical)** | Discretize governing equations on a computational mesh and solve them iteratively on a computer. | Any geometry, any flow regime, full-field data (pressure, velocity, temperature everywhere), cheap parameter sweeps, safe to study dangerous conditions. | Results are only as good as the model (turbulence closure, mesh resolution), requires validation, demands expertise. |

```
  ┌───────────────────────────────────────────────────────────────────────────┐
  │                  THE THREE PILLARS OF FLUID MECHANICS                     │
  │                                                                           │
  │      Analytical              Experimental                CFD              │
  │     ┌──────────┐           ┌──────────────┐        ┌──────────────┐       │
  │     │  Exact   │           │  Wind tunnel │        │  Computer    │       │
  │     │  math    │           │  Water tank  │        │  simulation  │       │
  │     │ solution │           │  Test rig    │        │  (OpenFOAM)  │       │
  │     └─────┬────┘           └──────┬───────┘        └──────┬───────┘       │
  │           │                       │                       │               │
  │           └───────────────────────┼───────────────────────┘               │
  │                                   ▼                                       │
  │                          VALIDATED DESIGN                                 │
  │                                                                           │
  │   Best practice: use all three where possible. CFD for exploration,       │
  │   experiments for validation, analytical solutions for sanity checks.     │
  └───────────────────────────────────────────────────────────────────────────┘
```

**When to prefer CFD over experiments:**

- Early design stages when no prototype exists yet
- Parameter sweeps (e.g., testing 50 angles of attack vs. 5 in a wind tunnel)
- Hazardous or inaccessible environments (nuclear reactors, human arteries)
- When you need full-field data (not just point measurements)

**When experiments are still essential:**

- Validating your CFD model (you cannot trust a simulation that has never been checked)
- Turbulence transition and complex multiphase phenomena that models approximate poorly
- Regulatory certification (e.g., aerospace, automotive crash testing)

### CFD Tools — The Software Landscape

A variety of CFD tools exist, ranging from general-purpose commercial packages to
specialized open-source solvers:

| Tool | Type | Strengths | Typical Use |
|------|------|-----------|-------------|
| **OpenFOAM** | Open-source (GPL), C++ | Free, fully customizable, unlimited parallelization, huge model library | Academic research, automotive, energy, marine — anyone who needs flexibility or runs on HPC clusters |
| **ANSYS Fluent** | Commercial | Polished GUI, extensive documentation, robust meshing (Fluent Meshing), industry standard | Aerospace, automotive, consulting firms with license budgets |
| **ANSYS CFX** | Commercial | Strong multiphase and turbomachinery capabilities, coupled solver | Turbomachinery, rotating equipment |
| **Siemens STAR-CCM+** | Commercial | Polyhedral meshing, integrated workflow, good automation API (Java) | Automotive OEMs (BMW, VW, Toyota), marine |
| **COMSOL Multiphysics** | Commercial | Excellent multi-physics coupling (structural, thermal, electromagnetic + CFD) | Biomedical, MEMS, electronics cooling, coupled physics |
| **SU2** | Open-source (LGPL) | Purpose-built for aerodynamic shape optimization, adjoint methods | Aerospace research, shape optimization |
| **Gerris / Basilisk** | Open-source | Adaptive octree meshing, excellent for free-surface and multiphase flows | Academic multiphase research |
| **Code_Saturne** | Open-source (GPL) | Developed by EDF, strong in nuclear and industrial flows | Nuclear thermal-hydraulics, industrial flows |

> **💡 Tip:** This repository focuses exclusively on **OpenFOAM**. The reasons are
> detailed in [OpenFOAM in the CFD Landscape](#openfoam-in-the-cfd-landscape)
> later in this note.

## The CFD Workflow

Every CFD project — whether academic or industrial — follows essentially the same
pipeline. Here it is at a glance:

```plaintext
┌──────────────────────────────────────────────────────────────────────────┐
│                        THE CFD WORKFLOW PIPELINE                         │
│                                                                          │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐               │
│  │ PROBLEM  │   │ GEOMETRY │   │ MESHING  │   │ PHYSICS  │               │
│  │ DEFINI-  │──▶│ CREATION │──▶│          │──▶│  SETUP   │               │
│  │  TION    │   │  (CAD)   │   │          │   │          │               │
│  └──────────┘   └──────────┘   └──────────┘   └──────────┘               │
│       │                                              │                   │
│       │         What are we solving?                 │                   │
│       │         What accuracy do we need?            ▼                   │
│       │                                        ┌──────────┐              │
│       │                                        │ SOLVING  │              │
│       │                                        │ (Run the │              │
│       │                                        │  solver) │              │
│       │                                        └────┬─────┘              │
│       │                                             │                    │
│       │         ┌──────────┐   ┌──────────────┐     │                    │
│       │         │ VALIDA-  │   │    POST-     │     │                    │
│       └────────▶│  TION    │◀──│  PROCESSING  │◀────┘                    │
│                 │          │   │              │                          │
│                 └──────────┘   └──────────────┘                          │
│                                                                          │
│  Compare with experiments / analytical solutions / grid studies          │
└──────────────────────────────────────────────────────────────────────────┘
```

### Step-by-Step Breakdown

| Step | What Happens | OpenFOAM Equivalent |
|------|-------------|---------------------|
| **Problem Definition** | Identify the physics, domain, and objectives | Choosing a solver (e.g., `icoFoam`, `simpleFoam`) |
| **Geometry Creation** | Build or import the physical shape | STL files, `blockMeshDict` geometry entries |
| **Meshing** | Divide geometry into small cells | `blockMesh`, `snappyHexMesh` (see note 04) |
| **Physics Setup** | Define equations, boundary conditions, fluid properties | `0/`, `constant/`, `system/` directories |
| **Solving** | Run the numerical solver | `icoFoam`, `simpleFoam`, `pimpleFoam`, etc. |
| **Post-Processing** | Visualize and extract results | ParaView, `postProcess` utilities |
| **Validation** | Verify results against known data | Comparison scripts, grid independence studies |

> **📂 Repo reference:** Look at any project directory (e.g., `projects/01_lid_driven_cavity/`)
> and you will see this workflow reflected in the folder structure: geometry in
> `constant/polyMesh`, physics in `0/` and `constant/`, solver settings in `system/`.

## Governing Equations

All of CFD rests on three conservation laws from classical physics. These are partial differential equations (PDEs) that describe how mass, momentum, and energy are transported through a fluid.

### The Control Volume Concept

Before writing equations, we need to understand what we are writing them *for*. In CFD, we consider a small fixed region of space — a **control volume** — and track what flows in and out of it:

```plaintext
                    ┌─── mass flux in (top) ────┐
                    │     ṁ_top = ρ u · n A     │
                    ▼                           │
         ┌──────────────────────────┐           │
         │                          │           │
  mass   │                          │   mass    │
  flux ──▶   CONTROL VOLUME (CV)    ├──▶ flux   │
  in     │                          │   out     │
 (left)  │   • density ρ            │  (right)  │
         │   • velocity u           │           │
         │   • pressure p           │           │
         │   • energy e             │           │
         │                          │           │
         └──────────────────────────┘           │
                    │                           │
                    │     ṁ_bot = ρ u · n A     │
                    ▼                           │
                    └─── mass flux out (bot) ───┘

  The fundamental idea: Rate of change inside CV =
      (What flows IN) − (What flows OUT) + (Sources)
```

This is the **Reynolds Transport Theorem** in action. Every governing equation in CFD has this structure: *accumulation = flux in − flux out + sources*.

### Conservation of Mass (Continuity Equation)

**Physical meaning:** Mass cannot be created or destroyed. If fluid enters a region, it must either accumulate there or leave somewhere else.

**General (compressible) form:**

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

where $\rho$ is density, $t$ is time, $\mathbf{u}$ is the velocity vector, and $\nabla \cdot$ is the divergence operator.

**Incompressible simplification:** When density is constant ($\rho = \text{const}$), the equation reduces to:

$$\nabla \cdot \mathbf{u} = 0$$

This says that the velocity field must be **divergence-free** — fluid cannot pile up or thin out anywhere.

> **📂 Repo reference:** In our lid-driven cavity project
> (`projects/01_lid_driven_cavity/`), `icoFoam` assumes incompressible flow, so it
> enforces $\nabla \cdot \mathbf{u} = 0$ at every time step.

### Conservation of Momentum (Navier-Stokes Equations)

**Physical meaning:** Newton's second law applied to a fluid element — the rate of change of momentum equals the sum of all forces (pressure, viscous, gravity).

**General (compressible) form:**

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g}$$

where $p$ is pressure, $\boldsymbol{\tau}$ is the viscous stress tensor, and $\mathbf{g}$ is gravitational acceleration.

**Incompressible Navier-Stokes equations** (constant $\rho$, Newtonian fluid):

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{g}$$

where $\nu = \mu / \rho$ is the kinematic viscosity.

> **💡 Key insight:** The incompressible Navier-Stokes equations are what `icoFoam`
> solves. The left side captures how fluid accelerates (locally in time and by being
> carried along — *convection*). The right side captures the forces: pressure gradient
> pushing the fluid, viscous diffusion smoothing out velocity gradients, and gravity.

This is the general constant-density equation. The stock OpenFOAM 6 `icoFoam` implementation has no explicit gravity source and uses kinematic pressure, `p = p_physical/rho`, in its momentum equation. See the [solver analysis](10_icofoam_solver_analysis.md) before translating physical equations into dictionary units.

### Conservation of Energy

**Physical meaning:** Energy is conserved — the rate of change of energy inside a control volume equals the net heat and work transfer.

$$\frac{\partial (\rho e)}{\partial t} + \nabla \cdot (\rho e \mathbf{u}) = -p \nabla \cdot \mathbf{u} + \boldsymbol{\tau}:\nabla\mathbf{u} + \nabla\cdot(\kappa\nabla T) + Q$$

Here $e$ is specific internal energy, $\kappa$ is thermal conductivity, $T$ is temperature, and $Q$ is a volumetric heat source. This form assumes Fourier heat conduction. Viscous dissipation is $\boldsymbol{\tau}:\nabla\mathbf{u}$; stress-work flux and gravitational work appear in the corresponding mechanical/total-energy balance. Mixing those terms into an internal-energy equation double-counts mechanical contributions. See [conservation of energy](https://doc.cfd.direct/notes/cfd-general-principles/conservation-of-energy).

> **Note:** For the incompressible, isothermal flows in most of our tutorial projects,
> the energy equation is not solved — temperature is assumed constant. Solvers like
> `buoyantSimpleFoam` add the energy equation when thermal effects matter.

### Turbulence and the Closure Problem

When flow becomes turbulent (high Reynolds number), the Navier-Stokes equations are still valid, but resolving every tiny eddy is computationally prohibitive. The **Reynolds-Averaged Navier-Stokes (RANS)** approach decomposes each variable into a mean and a fluctuating part, then time-averages. This introduces unknown **Reynolds stress** terms that require a *turbulence model* to close the system.

Common turbulence models (see note 06 for details):

| Model      | Type    | Best For                                    |
|------------|---------|---------------------------------------------|
| $k$-$\varepsilon$        | RANS    | General industrial flows, free-shear layers |
| $k$-$\omega$ SST    | RANS    | Boundary layers, adverse pressure gradients  |
| Spalart-Allmaras | RANS | Aerospace external flows                |
| LES        | Scale-resolving | Unsteady flows needing eddy resolution |

> **📂 Repo reference:** The NACA airfoil project (`projects/04_naca_airfoil_analysis/`)
> uses the Spalart–Allmaras turbulence model with `simpleFoam`, while the lid-driven cavity
> (`projects/01_lid_driven_cavity/`) is solved as laminar.

## Discretization — From Calculus to Algebra

The governing equations are continuous PDEs — they apply at every point in space and every instant in time. Computers cannot solve continuous equations directly. We must **discretize**: convert the continuous domain into a finite set of cells and convert the PDEs into algebraic equations that can be solved with linear algebra.

### Mesh: Dividing Space Into Cells

```plaintext
  CONTINUOUS DOMAIN                     DISCRETIZED DOMAIN (MESH)

  ┌─────────────────────┐              ┌────┬────┬────┬────┬────┐
  │                     │              │    │    │    │    │    │
  │   Fluid fills the   │              ├────┼────┼────┼────┼────┤
  │   entire region     │    ──▶       │    │    │    │    │    │
  │   continuously      │              ├────┼────┼────┼────┼────┤
  │                     │              │    │    │    │    │    │
  └─────────────────────┘              ├────┼────┼────┼────┼────┤
                                       │    │    │    │    │    │
  Infinite degrees of freedom          └────┴────┴────┴────┴────┘

                                       Finite number of cells, each
                                       with one value of p, U, etc.
```

**1D mesh terminology:**

```plaintext
  face  cell center  face  cell center  face
   │         │        │        │          │
   ▼         ▼        ▼        ▼          ▼

   |         ●        |        ●          |     ●      |     ●     |     ●     |
   |       cell 1     |       cell 2      |  cell 3    |  cell 4   |  cell 5   |
   |                  |                   |            |           |           |
   |  ◀──── Δx ────▶  |  ◀──── Δx ─────▶  | ◀── Δx ──▶ |◀── Δx ──▶ |◀── Δx ──▶ |

   ●  = cell centroid (where field values are stored in FVM)
   |  = cell face     (where fluxes are computed)
```

### Three Discretization Philosophies

```plaintext
  ┌──────────────────────────────────────────────────────────────────────┐
  │                  DISCRETIZATION APPROACHES                           │
  │                                                                      │
  │  FINITE DIFFERENCE (FDM)      FINITE VOLUME (FVM)                    │
  │                                                                      │
  │   ●────●────●────●            ┌────┬────┬────┐     Values stored     │
  │   │    │    │    │            │ ●  │ ●  │ ●  │     at cell CENTERS   │
  │   ●────●────●────●            ├────┼────┼────┤     Fluxes across     │
  │   │    │    │    │            │ ●  │ ●  │ ●  │     cell FACES        │
  │   ●────●────●────●            └────┴────┴────┘                       │
  │                                                                      │
  │   Values at grid NODES        Integral form of                       │
  │   Derivatives ≈ differences   equations → inherently                 │
  │   between neighbors           CONSERVATIVE                           │
  │                                                                      │
  │  FINITE ELEMENT (FEM)                                                │
  │                                                                      │
  │   ●─────●─────●              Values represented by                   │
  │   │\    │    /│              basis functions over                    │
  │   │  \  │  /  │              each element.                           │
  │   │    \│/    │              Minimizes a weighted                    │
  │   ●─────●─────●              residual (weak form).                   │
  │   │    /│\    │              Excellent for complex                   │
  │   │  /  │  \  │              geometry & solid mech.                  │
  │   │/    │    \│                                                      │
  │   ●─────●─────●                                                      │
  │                                                                      │
  └──────────────────────────────────────────────────────────────────────┘
```

### Comparison Table

| Feature                  | FDM                      | FVM                         | FEM                          |
|--------------------------|--------------------------|-----------------------------|------------------------------|
| **Formulation**          | Differential (pointwise) | Integral (over volumes)     | Weak / variational           |
| **Conservation**         | Not inherent             | Inherently conservative ✔  | Depends on formulation       |
| **Mesh flexibility**     | Structured grids only    | Unstructured OK             | Very flexible, unstructured  |
| **Complex geometry**     | Poor                     | Good                        | Excellent                    |
| **Implementation**       | Simple                   | Moderate                    | Complex                      |
| **Typical use**          | Research, simple domains | CFD (most codes)            | Structural, multiphysics     |
| **Example codes**        | Custom academic codes    | **OpenFOAM**, Fluent, CFX   | COMSOL, Abaqus, FEniCS       |

> **⚠️ Why OpenFOAM uses FVM:** The finite volume method naturally conserves mass,
> momentum, and energy because it works with the *integral* form of the equations.
> Fluxes leaving one cell enter the neighbor exactly — nothing is lost or created at
> interfaces. This is critical for physically meaningful CFD results.

### From PDE to Linear System

The discretization process converts a PDE like:

$$\frac{\partial \phi}{\partial t} + \nabla \cdot (\mathbf{u} \phi) = \Gamma \nabla^2 \phi$$

into a system of algebraic equations, one per cell:

$$a_P \phi_P = \sum_{N} a_N \phi_N + b_P$$

where $\phi_P$ is the value at cell $P$, $\phi_N$ are values at neighboring cells, and $a_P$, $a_N$, $b_P$ are coefficients that come from the discretization. Assembled over all cells, this becomes a matrix equation:

$$\mathbf{A} \boldsymbol{\phi} = \mathbf{b}$$

This sparse linear system is what the linear solvers (PCG, smoothSolver, GAMG — see note 09) actually solve at every iteration or time step.

## 5. The Finite Volume Method in Detail

Since OpenFOAM is built entirely on FVM, it is worth understanding how this method works
in more depth.

### Integral Form and Gauss's Divergence Theorem

Starting from the general transport equation in differential form:

$$\frac{\partial \phi}{\partial t} + \nabla \cdot (\mathbf{u} \phi) = \nabla \cdot (\Gamma \nabla \phi) + S$$

We integrate over a control volume $V$ bounded by surface $\partial V$:

$$\int_V \frac{\partial \phi}{\partial t} \, dV + \int_V \nabla \cdot (\mathbf{u} \phi) \, dV = \int_V \nabla \cdot (\Gamma \nabla \phi) \, dV + \int_V S \, dV$$

Applying **Gauss's Divergence Theorem** ($\int_V \nabla \cdot \mathbf{F} \, dV = \oint_{\partial V} \mathbf{F} \cdot \mathbf{n} \, dA$), the volume integrals of divergence terms become surface integrals:

$$\frac{d}{dt} \int_V \phi \, dV + \oint_{\partial V} (\mathbf{u} \phi) \cdot \mathbf{n} \, dA = \oint_{\partial V} (\Gamma \nabla \phi) \cdot \mathbf{n} \, dA + \int_V S \, dV$$

```plaintext
  ┌─────────────────────────────────────────────────────────────────┐
  │              FINITE VOLUME: FACE FLUX BALANCE                   │
  │                                                                 │
  │                        n_top ▲                                  │
  │                    ┌─────────┼─────────┐                        │
  │                    │ F_top   │         │                        │
  │                    │ ════════╪══════   │                        │
  │                    │         │         │                        │
  │      n_left        │    cell center    │       n_right          │
  │    ◀───────────────│─── ● (P) ────────│───────────────▶         │
  │       F_left       │         │         │      F_right           │
  │                    │         │         │                        │
  │                    │ ════════╪══════   │                        │
  │                    │ F_bot   │         │                        │
  │                    └─────────┼─────────┘                        │
  │                              ▼ n_bot                            │
  │                                                                 │
  │   F = flux through face = (ρ u φ) · n · A                       │
  │   Sum of all face fluxes = sources inside cell                  │
  │   ΣF_faces = S_cell · V_cell                                    │
  │                                                                 │
  │   KEY PROPERTY: Flux leaving cell P through a shared face       │
  │   is EXACTLY the flux entering the neighbor → CONSERVATION      │
  └─────────────────────────────────────────────────────────────────┘
```

### Cell-Centered vs Vertex-Centered

| Approach | Where values are stored | Used by |
|----------|------------------------|---------|
| **Cell-centered** | At the centroid of each cell | **OpenFOAM**, Fluent, STAR-CCM+ |
| Vertex-centered | At the vertices (corners) of each cell | SU2, some FEM-based codes |

OpenFOAM stores volume fields such as p, U, k, and $\varepsilon$ at cell centers. Surface fields such as the face flux `phi` are stored on faces. When a volume-field value is needed at a face, an interpolation scheme supplies it; the data location depends on the field type.

### Discrete Operators in OpenFOAM

OpenFOAM's `fvm::` and `fvc::` namespaces map directly to the FVM discretization:

| Operator             | Continuous          | OpenFOAM Code                  |
|----------------------|---------------------|--------------------------------|
| Time derivative      | $\partial/\partial t$ | `fvm::ddt(U)`                |
| Convection           | $\nabla \cdot (\mathbf{u}\phi)$ | `fvm::div(phi, U)` |
| Diffusion (Laplacian)| $\nabla^2 \phi$     | `fvm::laplacian(nu, U)`        |
| Gradient             | $\nabla p$          | `fvc::grad(p)`                 |

> **📂 Repo reference:** In note 10 (`notes/10_icofoam_solver_analysis.md`), the
> icoFoam source code is analyzed line by line showing exactly how these operators
> appear in the momentum equation assembly.

## Numerical Schemes — Accuracy, Stability, and Trade-Offs

Discretization schemes control *how* the integrals and derivatives are approximated.
The choice of scheme directly affects accuracy, stability, and convergence.

### Time Discretization

```plaintext
  TIME DISCRETIZATION: HOW WE MARCH THROUGH TIME

  t^n                    t^(n+1)                  t^(n+2)
   │                       │                       │
   ●━━━━━━━━━━━━━━━━━━━━━━━●━━━━━━━━━━━━━━━━━━━━━━━●━━━━━▶ time
   │       Δt              │        Δt             │
   │                       │                       │
   │  ┌─── Euler ────┐     │                       │
   │  │ Use value at │     │                       │
   │  │ t^n to get   │     │                       │
   │  │ t^(n+1)      │     │                       │
   │  └──────────────┘     │                       │
   │                       │                       │
   │  ┌─── Backward ──────────────┐                │
   │  │ Uses t^n AND t^(n-1)      │                │
   │  │ for 2nd-order accuracy    │                │
   │  └───────────────────────────┘                │
   │                       │                       │
   │  ┌─── Crank-Nicolson ─┐                       │
   │  │ Averages t^n and   │                       │
   │  │ t^(n+1) → 2nd order│                       │
   │  │ but can oscillate  │                       │
   │  └────────────────────┘                       │
```

The linear stability entries below refer to the standard scalar test problem;
they do not guarantee a stable coupled CFD calculation. Boundedness also
depends on the spatial discretization. In OpenFOAM, `CrankNicolson 1` is the
full Crank–Nicolson scheme, while `CrankNicolson 0` reduces to Euler; `0.5`
blends toward Euler and is not the full second-order setting. See the
[v6 Crank–Nicolson implementation notes](https://cpp.openfoam.org/v6/classFoam_1_1fv_1_1CrankNicolsonDdtScheme.html).

| Scheme | Order | Linear stability | Bounded? | OpenFOAM keyword |
|--------|-------|-----------|----------|------------------|
| **Euler** (implicit) | 1st | A-stable | With a monotone spatial operator | `Euler` |
| **backward** | 2nd at constant step | A-stable at constant step | No (can overshoot) | `backward` |
| **Crank-Nicolson** | 2nd when fully centered | A-stable | No (can oscillate) | `CrankNicolson 1` (full); `0.5` (damped blend) |
| **steadyState** | N/A | N/A | N/A | `steadyState` |

> **📂 Repo reference:** In our lid-driven cavity
> (`projects/01_lid_driven_cavity/system/fvSchemes`), Euler time stepping is used:
> ```
> ddtSchemes { default Euler; }
> ```
> This is first-order but very robust — a good default for getting a simulation running.
> The airfoil project (`projects/04_naca_airfoil_analysis/system/fvSchemes`) uses
> `steadyState` because we seek a time-independent solution.

### Spatial Discretization (Convection Schemes)

The convection term $\nabla \cdot (\mathbf{u} \phi)$ is the trickiest to discretize because it can introduce numerical oscillations (wiggles) or excessive smearing (artificial diffusion).

```plaintext
  ACCURACY vs. STABILITY TRADE-OFF

  More Stable                                          More Accurate
  (diffusive)                                          (can oscillate)
  ◀━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━▶

  ┌──────────┐    ┌──────────────┐    ┌────────┐    ┌──────────────┐
  │ 1st-order│    │ 2nd-order    │    │ linear │    │ Higher-order │
  │  upwind  │    │ linearUpwind │    │(central│    │ (QUICK, etc.)│
  │          │    │              │    │  diff) │    │              │
  └──────────┘    └──────────────┘    └────────┘    └──────────────┘
   Very stable     Good balance        Accurate      Best accuracy
   Very diffusive  Most common!        Can wiggle    Complex, costly
```

| Scheme | Order | Bounded? | Diffusion | OpenFOAM keyword |
|--------|-------|----------|-----------|------------------|
| **upwind** | 1st | Yes ✔ | High (smears gradients) | `Gauss upwind` |
| **linearUpwind** | 2nd | Mostly ✔ | Low | `Gauss linearUpwind grad(U)` |
| **linear** (central) | 2nd | No ✘ | None | `Gauss linear` |
| **QUICK** | 3rd | No ✘ | Very low | `Gauss QUICK` |
| **limitedLinear** | 2nd | Yes ✔ | Low | `Gauss limitedLinear 1` |

> **📂 Repo reference:** The cavity project uses `Gauss linearUpwind grad(U)` for
> velocity convection — a good compromise of accuracy and stability. The airfoil project
> uses `bounded Gauss upwind` for `U` and the SA transport variable `nuTilda`.
> The $k$ and $\varepsilon$ examples elsewhere illustrate two-equation models, rather than fields in this current airfoil case.

### Gradient Schemes

Gradients are needed for pressure gradient evaluation and for reconstructing face values in higher-order convection schemes. Common choices:

- `Gauss linear` — compute the gradient using Gauss's theorem with linear interpolation to faces. This is the most common default.
- `cellLimited Gauss linear 1` — same, but limit the gradient so that reconstructed face values do not exceed/undershoot cell neighbor values. Prevents oscillations.

### Laplacian (Diffusion) Schemes

The diffusion term $\nabla \cdot (\Gamma \nabla \phi)$ is discretized as:

$$\sum_f \Gamma_f (\nabla \phi)_f \cdot \mathbf{S}_f$$

In OpenFOAM: `Gauss linear corrected` — linear interpolation of $\Gamma$ to the face, with a non-orthogonal correction for meshes whose face normals are not aligned with the vector between cell centers.

> **⚠️ Warning:** On highly non-orthogonal meshes (angle > 70°), you may need
> `nNonOrthogonalCorrectors` > 0 in `fvSolution` to get accurate diffusion terms.

## Pressure-Velocity Coupling

### The Problem

In incompressible flow, there is no explicit equation for pressure — the continuity equation $\nabla \cdot \mathbf{u} = 0$ acts as a *constraint* on the velocity field. The momentum equation gives us velocity, but it contains the pressure gradient. We need a way to couple pressure and velocity so both the momentum equation and the continuity constraint are satisfied simultaneously.

### The Major Algorithms

```plaintext
  ┌─────────────────────────────────────────────────────────────────────┐
  │                     PISO ALGORITHM FLOWCHART                        │
  │                                                                     │
  │     ┌───────────────────────────┐                                   │
  │     │  Start new time step t^n  │                                   │
  │     └─────────────┬─────────────┘                                   │
  │                   ▼                                                 │
  │     ┌───────────────────────────┐                                   │
  │     │  Solve MOMENTUM equation  │◀─── Assemble H(U), A matrix       │
  │     │  for predicted U*         │                                   │
  │     └─────────────┬─────────────┘                                   │
  │                   ▼                                                 │
  │     ┌───────────────────────────┐                                   │
  │     │  Solve PRESSURE equation  │◀─── From continuity constraint    │
  │     │  (Poisson equation)       │     ∇·U = 0 → ∇²p = ∇·(H/A)       │
  │     └─────────────┬─────────────┘                                   │
  │                   ▼                                                 │
  │     ┌───────────────────────────┐                                   │
  │     │  CORRECT velocity:        │                                   │
  │     │  U = H/A − (1/A)∇p        │                                   │
  │     └─────────────┬─────────────┘                                   │
  │                   ▼                                                 │
  │            ┌──────────────┐                                         │
  │            │ Corrector    │ NO    (typically 2 correctors           │
  │            │ loops done?  ├──────▶ for standard PISO)               │
  │            └──────┬───────┘         Loop back to pressure solve     │
  │                   │ YES                                             │
  │                   ▼                                                 │
  │     ┌───────────────────────────┐                                   │
  │     │  Advance to t^(n+1)       │                                   │
  │     │  Write fields if needed   │                                   │
  │     └───────────────────────────┘                                   │
  └─────────────────────────────────────────────────────────────────────┘
```

### SIMPLE vs PISO vs PIMPLE

| Feature | SIMPLE | PISO | PIMPLE |
|---------|--------|------|--------|
| **Full name** | Semi-Implicit Method for Pressure-Linked Equations | Pressure-Implicit with Splitting of Operators | Merged PISO-SIMPLE |
| **Type** | Steady-state (iterative) | Transient (time-accurate) | Transient or pseudo-steady |
| **Outer loops** | Many iterations to convergence | 1 per time step | Configurable (1 = PISO, >1 = SIMPLE-like) |
| **Pressure correctors** | 1 | 2 (typically) | Configurable |
| **Under-relaxation** | Commonly needed; factors depend on case | Often omitted at small steps | Optional |
| **Time step** | N/A (no real time) | Small Co is a common coupling/accuracy target | Outer corrections can improve coupling at larger Co |
| **OpenFOAM solver** | `simpleFoam` | `icoFoam`, `pisoFoam` | `pimpleFoam` |
| **Best for** | Steady-state flows | Small $\Delta t$ transient | Large $\Delta t$ transient |

> **📂 Repo reference:** Our cavity project uses **PISO** (see
> `projects/01_lid_driven_cavity/system/fvSolution`) with `nCorrectors 2`, while the
> airfoil project uses **SIMPLE** (see `projects/04_naca_airfoil_analysis/system/fvSolution`)
> with under-relaxation factors of 0.2 for pressure, 0.3 for velocity, and 0.2 for `nuTilda`.
> The often-used 0.3/0.7 pressure/velocity pair is a starting example, not a theoretical requirement or a sum-to-one rule.

> **💡 Tip:** If you are unsure which algorithm to use, `pimpleFoam` with PIMPLE is the
> most flexible — it degenerates to PISO with `nOuterCorrectors 1` and behaves like
> SIMPLE with many outer correctors and under-relaxation.

## Convergence and Accuracy

### Residuals — How We Know the Solution Is Converging

At each iteration (or time step), the linear solver reduces the **residual** — the
imbalance in the discretized equation $\mathbf{A}\boldsymbol{\phi} = \mathbf{b}$.
The initial residual is:

$$r = ||\mathbf{b} - \mathbf{A}\boldsymbol{\phi}||$$

As the solver iterates, this should decrease. OpenFOAM prints residuals every iteration:

```
Time = 0.5
smoothSolver:  Solving for Ux, Initial residual = 3.2e-04, Final residual = 1.1e-06, No Iterations 4
smoothSolver:  Solving for Uy, Initial residual = 2.8e-04, Final residual = 9.3e-07, No Iterations 4
PCG:           Solving for p,  Initial residual = 5.1e-03, Final residual = 2.4e-06, No Iterations 12
```

**Rules of thumb:**

- Initial residuals dropping below ~$10^{-4}$ for velocity and ~$10^{-3}$ for pressure typically indicates decent convergence for SIMPLE-based steady solvers.
- For transient solvers, inspect both absolute and relative stopping criteria. A linear solve can stop on `tolerance`, `relTol`, or its iteration limit; `relTol 0` disables the relative criterion but does not guarantee convergence before `maxIter`.

OpenFOAM reports normalized linear residuals, not simply the unnormalized norm written above. Low residuals measure the algebraic solve; also check conservation, forces or other outputs, and mesh/time-step sensitivity before calling a CFD result accurate.

### Grid Independence

A CFD result is only meaningful if it does not change significantly when the mesh is refined. A **grid independence study** (or mesh convergence study) works as follows:

1. Run the simulation on a coarse mesh.
2. Refine the mesh (double the resolution in each direction).
3. Compare a key quantity (drag coefficient, pressure drop, etc.).
4. If the quantity changes by less than ~1–2%, the coarser mesh is sufficient.

> **💡 Tip:** Always report grid independence in academic work. A beautiful CFD picture
> on a coarse mesh can be completely wrong!

### The CFL Number

The Courant-Friedrichs-Lewy (CFL) number constrains the time step relative to the mesh and flow velocity:

$$\text{CFL} = \frac{u \, \Delta t}{\Delta x}$$

where $u$ is the local velocity, $\Delta t$ is the time step, and $\Delta x$ is the local cell size.

- **Explicit schemes:** the allowable CFL depends on the complete spatial/time discretization. Forward Euler with first-order upwind for 1D linear advection is stable for $0 \leq |u|\Delta t/\Delta x \leq 1$; CFL below one is not sufficient for every explicit scheme.
- **Implicit schemes:** some linear discretizations have no advective stability limit, but implicit time stepping does not make a nonlinear CFD calculation unconditionally stable or time accurate. PISO/PIMPLE corrections, explicit terms, and physical time scales still matter.

> **📖 See also:** Note 08 (`notes/08_cfl_number.md`) provides an in-depth treatment
> of the CFL number, including practical guidelines for choosing time steps in OpenFOAM.

### Sources of Error in CFD

| Error Type | Description | How to Minimize |
|------------|-------------|-----------------|
| **Modeling error** | Wrong physics (e.g., laminar model for turbulent flow) | Choose appropriate models |
| **Discretization error** | Truncation from finite mesh/time step | Refine mesh and time step |
| **Iterative error** | Solver not fully converged | Tighten tolerances, more iterations |
| **Round-off error** | Finite precision of floating-point arithmetic | Usually negligible (use double precision) |
| **Input error** | Wrong boundary conditions, geometry, or properties | Careful setup and validation |

## OpenFOAM in the CFD Landscape

### Comparison with Commercial Solvers

| Feature | OpenFOAM | ANSYS Fluent | Siemens STAR-CCM+ |
|---------|----------|-------------|-------------------|
| **Cost** | Free (GPL) | $$$$ (tens of thousands/year) | $$$$ |
| **Source code** | Fully open ✔ | Closed | Closed |
| **Customization** | Unlimited (C++) | Limited (UDFs) | Limited (Java macros) |
| **GUI** | Minimal (ParaView for post) | Integrated GUI | Integrated GUI |
| **Learning curve** | Steep | Moderate | Moderate |
| **Meshing** | blockMesh, snappyHexMesh | Fluent Meshing, ICEM | Built-in polyhedral mesher |
| **Parallelization** | Built-in (OpenMPI) ✔ | Built-in | Built-in |
| **Community** | Large, active ✔ | Vendor support | Vendor support |
| **Industry adoption** | Growing (automotive, energy) | Dominant | Strong (automotive) |

### Why Open-Source Matters

- **Transparency:** You can read every line of the solver. If it gives a wrong answer, you can find out *why* at the source-code level.
- **Reproducibility:** Anyone can run your exact simulation without purchasing a license. This is critical for academic publications.
- **Customization:** Need a custom boundary condition, a new turbulence model, or a coupled multi-physics solver? Write it in C++ and compile.
- **Cost:** Running 1000-core HPC jobs? No per-core license fees.

### OpenFOAM's Strengths and Limitations

**Strengths:**

- Comprehensive solver library (incompressible, compressible, multiphase, combustion, electromagnetics, solid mechanics).
- Flexible meshing with `snappyHexMesh` for complex geometries.
- Mature parallelization with domain decomposition.
- Extensive turbulence model library.

**Limitations:**

- No integrated GUI — everything is done through text dictionaries and the command line.
- Steep learning curve for beginners (which is why this repo exists!).
- Documentation can be sparse for advanced features.
- Two main forks (openfoam.org and openfoam.com) can cause confusion.

##  Quick Reference — Key Takeaways

| # | Concept                      | Description                                                                                                                    |
| - | ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| 1 | CFD                          | Solves fluid flow problems numerically by discretizing governing PDEs into algebraic equations on a mesh.                      |
| 2 | Conservation Laws            | Three laws govern fluid flow: mass (continuity), momentum (Navier–Stokes), and energy.                                         |
| 3 | Incompressible Navier–Stokes | Defined by $\nabla \cdot \mathbf{u} = 0$ and the momentum equation with $\nu \nabla^2 \mathbf{u}$; solved by `icoFoam`.        |
| 4 | Finite Volume Method         | Converts PDEs into integral form over control volumes, ensuring conservation of transported quantities.                        |
| 5 | Numerical Schemes            | Trade-off between accuracy and stability: upwind (stable but diffusive) vs. central differencing (accurate but may oscillate). |
| 6 | Pressure-Velocity Coupling   | Algorithms like SIMPLE, PISO, and PIMPLE handle coupling since no explicit pressure equation exists in incompressible flow.    |
| 7 | Validation                   | Essential step: check residuals, perform grid independence studies, and compare with analytical or experimental data.          |
| 8 | OpenFOAM                     | Open-source FVM-based CFD toolbox; free, customizable, and widely used in academia and industry.                               |

## Where to Go Next

| Topic | Note |
|-------|------|
| OpenFOAM case structure | `notes/02_openfoam_cases.md` |
| Dictionary files explained | `notes/03_openfoam_dictionaries.md` |
| Meshing in depth | `notes/04_meshing.md` |
| Boundary conditions | `notes/05_boundary_conditions.md` |
| Turbulence models | `notes/06_turbulence_models.md` |
| Parallel computing | `notes/07_parallelization.md` |
| CFL number deep dive | `notes/08_cfl_number.md` |
| Linear solvers | `notes/09_linear_solvers.md` |
| icoFoam source analysis | `notes/10_icofoam_solver_analysis.md` |
| Viscosity models | `notes/11_viscosity_models.md` |
| Multiphase flows | `notes/12_multiphase_flows.md` |
| **Hands-on:** Lid-driven cavity | `projects/01_lid_driven_cavity/` |
| **Hands-on:** NACA airfoil | `projects/04_naca_airfoil_analysis/` |


*Original revision: 2025; corrections and clarifications: 2026-10-16. Part of the [OpenFOAM-Tutorials](https://github.com/djeada/OpenFoam-Tutorials/blob/main/README.md) repository.*
