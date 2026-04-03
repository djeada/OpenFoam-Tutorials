```
 ╔═══════════════════════════════════════════════════════════════╗
 ║              OpenFOAM Tutorials & Learning Hub                ║
 ║         From Zero to CFD — A Hands-On Journey                 ║
 ╚═══════════════════════════════════════════════════════════════╝
```

# OpenFOAM Tutorials & Learning Hub

> **A complete, self-study curriculum for learning Computational Fluid Dynamics (CFD) with OpenFOAM — from first principles to production-grade simulations.**

Whether you're a student encountering CFD for the first time, an engineer switching from commercial solvers, or a researcher looking for a structured OpenFOAM reference — this repository has you covered. It combines **12 in-depth theory notes** with **12 progressively challenging simulation projects** so you can learn by *reading* and *doing*.

---

## Table of Contents

- [Why This Repository?](#why-this-repository)
- [Repository Structure](#repository-structure)
- [Learning Path](#learning-path)
- [Notes Overview](#-notes-overview)
- [Projects Overview](#-projects-overview)
- [Quick Start](#-quick-start)
- [Installation](#installation)
- [Key Concepts at a Glance](#-key-concepts-at-a-glance)
- [How to Contribute](#how-to-contribute)
- [References & Resources](#references--resources)
- [License](#license)

---

## Why This Repository?

Most OpenFOAM resources fall into one of two traps: they're either *too theoretical* (pages of equations with no hands-on work) or *too cookbook-style* (copy-paste commands with no understanding). This repository bridges that gap:

| What You Get | Why It Matters |
|---|---|
| 📚 **12 Theory Notes** | Understand the *why* behind every setting and parameter |
| 🔬 **12 Simulation Projects** | Build real muscle memory by running actual CFD cases |
| 📈 **Progressive Difficulty** | Start with laminar cavity flow, finish with multiphase & rotating machinery |
| 🧩 **Cross-Referenced** | Notes link to projects; projects reference the theory they rely on |
| 🐳 **Docker-Ready** | Get running in minutes without dependency headaches |

---

## Repository Structure

```
OpenFoam-Tutorials/
│
├── notes/                              ← 📚 Theory & Reference (12 comprehensive guides)
│   ├── 01_short_intro_to_cfd.md           CFD fundamentals & governing equations
│   ├── 02_openfoam_cases.md               Case structure & setup workflow
│   ├── 03_openfoam_dictionaries.md        Dictionary format & key dictionaries
│   ├── 04_meshing.md                      blockMesh, snappyHexMesh, mesh quality
│   ├── 05_boundary_conditions.md          BC types, initial conditions & physical meaning
│   ├── 06_turbulence_models.md            RANS, LES, DNS & model selection
│   ├── 07_parallelization.md              Domain decomposition, MPI & scaling
│   ├── 08_cfl_number.md                   CFL concept, stability & time stepping
│   ├── 09_linear_solvers.md               Solvers, preconditioners & convergence
│   ├── 10_icofoam_solver_analysis.md      PISO algorithm & source code walkthrough
│   ├── 11_viscosity_models.md             Newtonian & non-Newtonian viscosity models
│   └── 12_multiphase_flows.md             VOF method, interFoam & multiphase setup
│
├── projects/                           ← 🔬 Hands-On Simulations (12 projects)
│   ├── 01_lid_driven_cavity/              Classic benchmark — icoFoam, blockMesh
│   ├── 02_elbow/                          Curved pipe flow — Fluent mesh import
│   ├── 03_flow_around_cylinder/           External flow — snappyHexMesh
│   ├── 04_naca_airfoil_analysis/          Airfoil aerodynamics — k-epsilon
│   ├── 05_heated_pipe_flow/               Conjugate heat transfer — buoyantSimpleFoam
│   ├── 06_ahmed_body_aerodynamics/        Automotive aero — k-omega SST
│   ├── 07_backward_facing_step/           Separated flow validation — k-omega SST
│   ├── 08_dam_break/                      Multiphase VOF — interFoam
│   ├── 09_wind_turbine_blade/             Rotating machinery — MRF + simpleFoam
│   ├── 10_mixing_t_junction/              Thermal mixing — buoyantSimpleFoam
│   ├── 11_boat_hull_fixed/                Fixed hull in water — interFoam + snappyHexMesh
│   └── 12_boat_hull_floating/             Floating body dynamics — interDyMFoam + 6-DoF
│
├── document.pdf                        ← 📄 Supplementary reference document
├── LICENSE                             ← ⚖️  MIT License
└── README.md                           ← 📖 You are here
```

---

## Learning Path

Follow this roadmap from left to right. Each stage builds on the previous one.

```
 🟢 BEGINNER                      🟡 INTERMEDIATE                 🔴 ADVANCED
 ┌──────────────────────┐         ┌──────────────────────┐        ┌──────────────────────┐
 │                      │         │                      │        │                      │
 │  📚 Notes 01–03      │         │  📚 Notes 04–07,     │        │  📚 Notes 08–10      │
 │  ─────────────────   │         │      11–12           │        │  ─────────────────   │
 │  • CFD Fundamentals  │         │  ─────────────────   │        │  • CFL & Stability   │
 │  • Case Structure    │────────▶│  • Meshing           │───────▶│  • Linear Solvers    │
 │  • Dictionaries      │         │  • Boundary Conds.   │        │  • icoFoam Deep Dive │
 │                      │         │  • Turbulence Models │        │                      │
 │                      │         │  • Parallelization   │        │                      │
 │                      │         │  • Viscosity Models  │        │                      │
 │                      │         │  • Multiphase Flows  │        │                      │
 ├──────────────────────┤         ├──────────────────────┤        ├──────────────────────┤
 │                      │         │                      │        │                      │
 │  🔬 Projects         │         │  🔬 Projects         │        │  🔬 Projects         │
 │  ─────────────────   │         │  ─────────────────   │        │  ─────────────────   │
 │  01 Lid-Driven       │         │  03 Cylinder Flow    │        │  05 Ahmed Body       │
 │     Cavity           │         │  04 NACA Airfoil     │        │  06 Backward Step    │
 │  02 Elbow Flow       │         │  04 Heated Pipe      │        │  07 Dam Break (VOF)  │
 │                      │         │                      │        │  08 Wind Turbine     │
 │                      │         │                      │        │  09 T-Junction Mix   │
 │                      │         │                      │        │  10 Boat Hull Fixed  │
 │                      │         │                      │        │  11 Boat Hull Float  │
 │                      │         │                      │        │                      │
 └──────────────────────┘         └──────────────────────┘        └──────────────────────┘

 Est. 1–2 weeks                    Est. 2–3 weeks                  Est. 3–4 weeks
```

> **💡 Tip:** Don't skip the theory notes! Understanding *why* a setting exists is the difference between debugging a simulation in 10 minutes vs. 10 hours.

---

## 📚 Notes Overview

Twelve self-contained guides covering every major concept you'll need. Read them in order for the best experience, or jump to a specific topic as a reference.

| # | Note | Key Topics | Level |
|---|------|------------|-------|
| 01 | [**A Short Introduction to CFD**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md) | Navier-Stokes equations, CFD vs other approaches, CFD tools, finite volume method, discretization | 🟢 Beginner |
| 02 | [**OpenFOAM Cases**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md) | Case directory structure (`0/`, `constant/`, `system/`), setup workflow, running a case | 🟢 Beginner |
| 03 | [**OpenFOAM Dictionaries**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md) | Dictionary syntax, `controlDict`, `fvSchemes`, `fvSolution`, `transportProperties` | 🟢 Beginner |
| 04 | [**Meshing**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/04_meshing.md) | `blockMesh` fundamentals, `snappyHexMesh` workflow, mesh quality metrics, refinement strategies | 🟡 Intermediate |
| 05 | [**Boundary & Initial Conditions**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/05_boundary_conditions.md) | Initial conditions, Dirichlet/Neumann types, `fixedValue`, `zeroGradient`, `inlet/outlet`, wall functions | 🟡 Intermediate |
| 06 | [**Turbulence Models**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md) | RANS vs LES vs DNS, k-epsilon, k-omega SST, model selection guidelines, y+ requirements | 🟡 Intermediate |
| 07 | [**Parallelization**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/07_parallelization.md) | Domain decomposition methods, `decomposePar`, MPI execution, scaling, `reconstructPar` | 🟡 Intermediate |
| 08 | [**CFL Number**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/08_cfl_number.md) | Courant-Friedrichs-Lewy condition, stability criteria, adaptive time stepping, `maxCo` | 🔴 Advanced |
| 09 | [**Linear Solvers**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/09_linear_solvers.md) | PCG, PBiCGStab, GAMG, preconditioners, tolerance settings, convergence troubleshooting | 🔴 Advanced |
| 10 | [**icoFoam Solver Analysis**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/10_icofoam_solver_analysis.md) | PISO algorithm step-by-step, C++ source code walkthrough, pressure-velocity coupling | 🔴 Advanced |
| 11 | [**Viscosity Models**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/11_viscosity_models.md) | Newtonian & non-Newtonian fluids, power-law, Bird-Carreau, Herschel-Bulkley, `transportProperties` | 🟡 Intermediate |
| 12 | [**Multiphase Flows**](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/12_multiphase_flows.md) | Volume of Fluid (VOF), `interFoam`, surface tension, `alpha` field, `setFields`, multiphase setup | 🟡 Intermediate |

---

## 🔬 Projects Overview

Twelve simulation projects arranged by increasing complexity. Each project directory contains the complete OpenFOAM case structure ready to run.

| # | Project | Description | Solver | Mesh | Turbulence | Level |
|---|---------|-------------|--------|------|------------|-------|
| 01 | [**Lid-Driven Cavity**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/01_lid_driven_cavity) | The "Hello World" of CFD — a square cavity with a moving top wall. Validates against Ghia et al. benchmark data. | `icoFoam` | `blockMesh` | Laminar | 🟢 Beginner |
| 02 | [**Elbow**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/02_elbow) | Flow through a curved pipe elbow. Demonstrates importing a Fluent mesh into OpenFOAM and post-processing in ParaView. | `icoFoam` | Fluent import | Laminar | 🟢 Beginner |
| 03 | [**Flow Around a Cylinder**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/03_flow_around_cylinder) | External flow past a circular cylinder. Introduces `snappyHexMesh` for geometry-based meshing and drag/lift computation. | `simpleFoam` | `snappyHexMesh` | Laminar/Low-Re | 🟡 Intermediate |
| 04 | [**NACA Airfoil Analysis**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis) | Aerodynamic analysis of a NACA airfoil profile. Covers angle-of-attack sweeps, lift/drag coefficients, and pressure distributions. | `simpleFoam` | Structured | k-epsilon | 🟡 Intermediate |
| 05 | [**Heated Pipe Flow**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/05_heated_pipe_flow) | Conjugate heat transfer in a heated pipe. Introduces the energy equation, buoyancy effects, and thermal boundary conditions. | `buoyantSimpleFoam` | `blockMesh` | k-epsilon | 🟡 Intermediate |
| 06 | [**Ahmed Body Aerodynamics**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/06_ahmed_body_aerodynamics) | Full 3D automotive aerodynamics simulation of the Ahmed body — the standard reference geometry for vehicle aero studies. | `simpleFoam` | `snappyHexMesh` | k-omega SST | 🔴 Advanced |
| 07 | [**Backward-Facing Step**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/07_backward_facing_step) | Classic turbulent separated flow benchmark. Validates reattachment length against NASA experimental data (Xr/H ≈ 6.0). | `simpleFoam` | `blockMesh` | k-omega SST | 🔴 Advanced |
| 08 | [**Dam Break**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/08_dam_break) | Multiphase free-surface flow using Volume of Fluid (VOF). Water column collapse with gravity — transient interFoam simulation. | `interFoam` | `blockMesh` | Laminar | 🔴 Advanced |
| 09 | [**Wind Turbine Blade**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/09_wind_turbine_blade) | Rotating machinery simulation using the Multiple Reference Frame (MRF) approach. Introduces `topoSet` and `MRFProperties`. | `simpleFoam` | `blockMesh` | k-epsilon | 🔴 Advanced |
| 10 | [**Mixing T-Junction**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/10_mixing_t_junction) | Thermal mixing of hot (350 K) and cold (300 K) streams in a T-junction. Multi-block structured mesh with buoyancy coupling. | `buoyantSimpleFoam` | `blockMesh` | k-epsilon | 🔴 Advanced |
| 11 | [**Boat Hull — Fixed**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/11_boat_hull_fixed) | Fixed toy boat hull in calm water with free surface. Uses your own STL geometry with snappyHexMesh for drag and wave pattern analysis. | `interFoam` | `snappyHexMesh` | Laminar | 🔴 Advanced |
| 12 | [**Boat Hull — Floating**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/12_boat_hull_floating) | Floating rigid-body boat with 6-DoF motion. Dynamic mesh morphing with heave, pitch, and roll — the capstone project. | `interDyMFoam` | `snappyHexMesh` | Laminar | 🔴 Advanced |

### Project Progression Map

```
  FUNDAMENTALS           EXTERNAL AERO          HEAT TRANSFER          ADVANCED PHYSICS
  ┌──────────┐          ┌──────────┐           ┌──────────┐          ┌──────────┐
  │ 01 Cavity│          │03 Cylinder│           │05 Heated│           │08 Dam    │
  │ icoFoam  │─────────▶│ simpleFoam│──────────▶│   Pipe   │     ┌──▶│  Break   │
  │ blockMesh│          │ snappyHex │           │ buoyant  │     │   │ interFoam│
  └──────────┘          └──────────┘           └──────────┘     │   │ VOF      │
       │                     │                      │            │   └──────────┘
       ▼                     ▼                      ▼            │        │
  ┌──────────┐          ┌──────────┐           ┌──────────┐     │        ▼
  │ 02 Elbow │          │04 Airfoil │           │10 T-Junct│     │   ┌──────────┐
  │ icoFoam  │          │ simpleFoam│           │ buoyant  │     └──▶│11 Boat   │
  │ Fluent   │          │ k-epsilon │           │ thermal  │         │ Fixed    │
  └──────────┘          └──────────┘           └──────────┘         │interFoam │
                             │                                       └──────────┘
                             ▼                                            │
                        ┌──────────┐           ┌──────────┐              ▼
                        │06 Ahmed  │           │07 BFS    │         ┌──────────┐
                        │ Body 3D  │           │ k-ω SST  │         │12 Boat   │
                        │ k-ω SST  │           │validation│         │ Floating │
                        └──────────┘           └──────────┘         │interDyM  │
                                                                    │ 6-DoF    │
                        ┌──────────┐                                └──────────┘
                        │09 Wind   │
                        │ Turbine  │
                        │ MRF      │
                        └──────────┘

  Learn:                Learn:                 Learn:                Learn:
  Basic setup           Complex meshes         Energy equation       Multiphase VOF
  Post-processing       Force coefficients     Buoyancy coupling     snappyHexMesh + STL
  Mesh import           Turbulence models      Thermal BCs           Floating body / 6-DoF
```

---

## 🚀 Quick Start

Want to run your first CFD simulation right now? Here's the fastest path.

### Prerequisites

- **OpenFOAM** installed (v6+ or ESI v2006+) — see [Installation](#installation) below
- **ParaView** for visualization (recommended)
- Basic terminal/command-line familiarity

### Your First Simulation in 5 Commands

Run the classic lid-driven cavity case:

```bash
# 1. Clone this repository
git clone https://github.com/your-username/OpenFoam-Tutorials.git
cd OpenFoam-Tutorials/projects/01_lid_driven_cavity

# 2. Generate the mesh
blockMesh

# 3. Run the solver
icoFoam

# 4. Visualize the results
paraFoam
```

That's it! You should see velocity contours of the lid-driven cavity flow.

> **📌 Next Steps:** Read through [Note 01 — Intro to CFD](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md) to understand what just happened, then explore [Note 02 — Case Structure](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md) to understand *why* the directories are organized the way they are.

---

## Installation

There are two recommended ways to get OpenFOAM up and running.

### Option A: Direct Installation (Ubuntu/Debian)

The simplest approach if you're on Ubuntu:

```bash
# Add the OpenFOAM repository
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
sudo add-apt-repository http://dl.openfoam.org/ubuntu

# Update and install (example: OpenFOAM 11)
sudo apt-get update
sudo apt-get install -y openfoam11

# Source the environment (add to ~/.bashrc for persistence)
source /opt/openfoam11/etc/bashrc

# Verify installation
simpleFoam -help
```

> **💡 Tip:** Add the `source` line to your `~/.bashrc` so OpenFOAM is available every time you open a terminal.

### Option B: Docker (Any Platform — Recommended for Portability)

Docker containers provide an isolated, consistent environment without modifying your host system. This is the best option for macOS/Windows users or if you need multiple OpenFOAM versions.

**1. Pull and run the OpenFOAM container:**

```bash
docker container run -ti cfdengine/openfoam
```

**2. Persist your data on the host (recommended):**

Mount a local directory so your simulation data survives container restarts:

```bash
# Create a directory for your OpenFOAM data
mkdir -p ~/openfoam_data

# Run with volume mount
docker container run -ti --rm \
  -v ~/openfoam_data:/home/foam \
  cfdengine/openfoam
```

**3. Handle SELinux restrictions (Fedora/RHEL/CentOS):**

If you encounter permission errors on SELinux-enabled systems:

```bash
docker container run --privileged -ti --rm \
  -v ~/openfoam_data:/home/foam \
  cfdengine/openfoam
```

### Setting Up the Container Environment

Once inside the Docker container:

```bash
# 1. Source the OpenFOAM environment
#    Check your version first:
ls /opt

#    Then source it (adjust version number as needed):
source /opt/openfoam6/etc/bashrc

# 2. Verify the FOAM_RUN variable
echo $FOAM_RUN

#    If empty, set it manually:
export FOAM_RUN=/home/foam

# 3. Create and enter the working directory
mkdir -p $FOAM_RUN && cd $FOAM_RUN

# 4. Set the OpenFOAM installation directory
export OPENFOAM_DIR=/opt/openfoam6/
```

### Verify Your Installation

Regardless of installation method, run these checks:

```bash
# Check that OpenFOAM is accessible
blockMesh -help

# Check solver availability
icoFoam -help
simpleFoam -help

# Run the built-in cavity tutorial as a smoke test
cp -r $FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity $FOAM_RUN/
cd $FOAM_RUN/cavity
blockMesh && icoFoam
```

If the solver runs without errors, you're good to go! 🎉

---

## 📋 Key Concepts at a Glance

A quick reference for the concepts you'll encounter most often. For deeper coverage, see the corresponding notes.

### The Three Essential Directories

Every OpenFOAM case follows the same structure:

```
myCase/
├── 0/                 ← Initial & boundary conditions
│   ├── U                 Velocity field
│   ├── p                 Pressure field
│   └── ...               Other fields (k, epsilon, omega, nut, etc.)
│
├── constant/          ← Physical properties & mesh
│   ├── transportProperties    Fluid properties (viscosity, density)
│   ├── turbulenceProperties   Turbulence model selection
│   └── polyMesh/              The computational mesh
│       ├── points
│       ├── faces
│       ├── owner
│       ├── neighbour
│       └── boundary
│
└── system/            ← Simulation controls
    ├── controlDict           Run control (time step, end time, output)
    ├── fvSchemes              Discretization schemes
    ├── fvSolution             Solver settings & tolerances
    ├── decomposeParDict       Parallel decomposition settings
    └── ...
```

> **📖 Full details:** [Note 02 — OpenFOAM Cases](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md)

### Common Solvers

| Solver | Flow Type | Time | Algorithm | Typical Use Case |
|--------|-----------|------|-----------|------------------|
| `icoFoam` | Incompressible, laminar | Transient | PISO | Learning, simple laminar flows |
| `simpleFoam` | Incompressible, turbulent | Steady-state | SIMPLE | Most external/internal aero problems |
| `pisoFoam` | Incompressible, turbulent | Transient | PISO | Unsteady turbulent flows |
| `pimpleFoam` | Incompressible, turbulent | Transient | PIMPLE | Large time steps, general purpose |
| `interFoam` | Multiphase (VOF) | Transient | PIMPLE | Free surface, wave, sloshing |
| `rhoSimpleFoam` | Compressible, turbulent | Steady-state | SIMPLE | High-speed aero, heat transfer |
| `buoyantSimpleFoam` | Buoyant, turbulent | Steady-state | SIMPLE | Natural convection, HVAC |

> **📖 Full details:** [Note 10 — icoFoam Solver Analysis](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/10_icofoam_solver_analysis.md)

### Common Commands

| Command | Purpose | When to Use |
|---------|---------|-------------|
| `blockMesh` | Generate structured hex mesh | Simple geometries, learning |
| `snappyHexMesh` | Generate complex unstructured mesh | Real-world geometries with STL files |
| `checkMesh` | Validate mesh quality | **Always** run after meshing |
| `decomposePar` | Split case for parallel runs | Before running on multiple processors |
| `mpirun -np N solver` | Run solver in parallel | Large cases needing multiple CPUs |
| `reconstructPar` | Recombine parallel results | After a parallel run, before post-processing |
| `paraFoam` | Launch ParaView for visualization | Viewing results |
| `foamToVTK` | Convert results to VTK format | Alternative post-processing path |
| `postProcess -func 'forceCoeffs'` | Compute force coefficients | Drag, lift, moment calculations |
| `foamLog` | Extract residuals from log | Monitoring convergence |

> **📖 Full details:** [Note 03 — Dictionaries](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/03_openfoam_dictionaries.md) and [Note 07 — Parallelization](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/07_parallelization.md)

### Turbulence Model Quick Selector

Not sure which turbulence model to use? Start here:

```
Is the flow laminar?
│
├── YES → No turbulence model needed (set to "laminar")
│
└── NO → Is it a simple internal flow or external aero?
    │
    ├── Internal (pipes, ducts) → k-epsilon (standard, robust)
    │
    └── External (vehicles, airfoils) → k-omega SST (better near-wall)
        │
        └── Need to capture vortex shedding? → LES (expensive!)
```

> **📖 Full details:** [Note 06 — Turbulence Models](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/06_turbulence_models.md)

---

## How to Contribute

Contributions that enhance this repository are very welcome! Here's how you can help:

### Ways to Contribute

- 🐛 **Bug fixes** — Typos, broken links, incorrect settings in cases
- 📝 **New notes** — Additional theory topics or deeper dives into existing ones
- 🔬 **New projects** — More simulation tutorials at any difficulty level
- 📊 **Validation data** — Experimental/benchmark data to compare against simulation results
- 🎨 **Visualizations** — ParaView screenshots, animations, or result plots

### Contribution Workflow

```bash
# 1. Fork the repository on GitHub

# 2. Clone your fork
git clone https://github.com/YOUR-USERNAME/OpenFoam-Tutorials.git
cd OpenFoam-Tutorials

# 3. Create a feature branch
git checkout -b feature/your-amazing-feature

# 4. Make your changes and commit
git add .
git commit -m "Add: description of your changes"

# 5. Push and open a Pull Request
git push origin feature/your-amazing-feature
```

Then open a Pull Request on GitHub with a clear description of what you've added or changed.

### Style Guidelines

- **Notes:** Use Markdown with clear headings, code blocks, and equations where relevant
- **Projects:** Include a local `README.md` in each project directory explaining the case setup, expected results, and how to run it
- **Code blocks:** Always specify the language (e.g., ` ```bash `, ` ```cpp `)
- **File references:** Use absolute links (e.g., `[Note 01](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md)`)

---

## References & Resources

### Official Documentation

| Resource | Link |
|----------|------|
| OpenFOAM Foundation User Guide | [openfoam.org/docs](https://openfoam.org/documentation/) |
| OpenFOAM ESI User Guide | [openfoam.com/documentation](https://www.openfoam.com/documentation/user-guide) |
| OpenFOAM C++ API Reference | [cpp.openfoam.org](https://cpp.openfoam.org/) |
| OpenFOAM Wiki | [wiki.openfoam.com](https://wiki.openfoam.com/) |
| OpenFOAM User Guide (local) | `$OPENFOAM_DIR/doc/Guides/OpenFOAMUserGuide-A4.pdf` |

### Video Tutorials

| Resource | Description | Link |
|----------|-------------|------|
| József Nagy's Series | Comprehensive YouTube course covering basics through advanced topics | [YouTube Playlist](https://youtube.com/playlist?list=PLcOe4WUSsMkH6DLHpsYyveaqjKxnEnQqB&si=Ddju9M-vhwHQWlnn) |
| CFD Direct Videos | Practical tutorials on specific OpenFOAM cases and features | [cfd.direct](https://cfd.direct/) |

### Courses & Written Tutorials

| Resource | Description | Link |
|----------|-------------|------|
| Wolf Dynamics Training | Extensive intro course — basics to custom boundary conditions | [figshare.com](https://figshare.com/) |
| Jibran Haider's Course | Beginner-friendly introduction to OpenFOAM fundamentals | [jibranhaider.com](https://jibranhaider.com/) |
| Artur Lidtke's Tutorials | C++ concepts in OpenFOAM, from basic to advanced | [GitHub](https://github.com/sayin/OpenFOAM-Tutorials) |

### Meshing Resources

| Resource | Description | Link |
|----------|-------------|------|
| ENGYS' snappyHexMesh Tour | In-depth guide to `snappyHexMesh` for complex mesh generation | [openfoamwiki.net](https://openfoamwiki.net/) |

### Recommended Books

- **"The OpenFOAM Technology Primer"** by T. Marić, J. Höpken, K. Mooney — Best overall OpenFOAM reference
- **"An Introduction to Computational Fluid Dynamics"** by H.K. Versteeg & W. Malalasekera — Excellent FVM fundamentals
- **"Computational Methods for Fluid Dynamics"** by J.H. Ferziger & M. Perić — Graduate-level CFD theory

---

## License

This project is licensed under the [MIT License](https://github.com/djeada/OpenFoam-Tutorials/blob/main/LICENSE) — see the LICENSE file for details.

You are free to use, modify, and distribute this material for any purpose. Attribution is appreciated but not required.

---

<div align="center">

**Happy simulating!** 🌊

*If this repository helped you, consider giving it a ⭐ on GitHub.*

</div>
