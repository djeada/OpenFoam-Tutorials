```
 ╔═══════════════════════════════════════════════════════════════╗
 ║              OpenFOAM Tutorials & Learning Hub                ║
 ║         From Zero to CFD — A Hands-On Journey                 ║
 ╚═══════════════════════════════════════════════════════════════╝
```

# OpenFOAM Tutorials & Learning Hub

> **A complete, self-study curriculum for learning Computational Fluid Dynamics (CFD) with OpenFOAM — from first principles to production-grade simulations.**

Whether you're a student encountering CFD for the first time, an engineer switching from commercial solvers, or a researcher looking for a structured OpenFOAM reference — this repository has you covered. It combines **12 in-depth theory notes** with **13 progressively challenging simulation projects** so you can learn by *reading* and *doing*.

---

## Version and Example Scope

The project launchers source **OpenFOAM Foundation 6** (`/opt/openfoam6/etc/bashrc`).
The notes explain that classic solver/dictionary interface; later Foundation
releases and OpenCFD releases can use different solver names, model dictionaries,
and boundary-condition syntax. A larger version number is not a guarantee that
a case runs unchanged. Use the project launcher’s matching environment or port
and verify the case explicitly.

The notes preserve illustrative k-epsilon and snappyHexMesh alternatives. The
current NACA project uses **Spalart–Allmaras with Gmsh**, as documented in the
[checked reference-case table](notes/02_openfoam_cases.md#current-reference-cases).
Each guide has a section index and previous/next links. Numerical examples are
teaching estimates; use the project files and validation scope when running them.

## Table of Contents

- [Version and Example Scope](#version-and-example-scope)
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
| 🔬 **13 Simulation Projects** | Build real muscle memory by running actual CFD cases |
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
├── projects/                           ← 🔬 Hands-On Simulations (13 projects)
│   ├── 01_lid_driven_cavity/              Classic benchmark — icoFoam, blockMesh
│   ├── 02_elbow/                          Curved pipe flow — Fluent mesh import
│   ├── 03_flow_around_cylinder/           External flow — Gmsh mesh, pimpleFoam
│   ├── 04_naca_airfoil_analysis/          Airfoil aerodynamics — Spalart-Allmaras, Gmsh mesh
│   ├── 05_heated_pipe_flow/               Heated turbulent pipe flow — buoyantSimpleFoam
│   ├── 06_ahmed_body_aerodynamics/        Automotive aero — k-omega SST
│   ├── 07_backward_facing_step/           Separated flow validation — k-omega SST
│   ├── 08_dam_break/                      Multiphase VOF — interFoam
│   ├── 09_wind_turbine_blade/             Rotating machinery — MRF + simpleFoam
│   ├── 10_mixing_t_junction/              Thermal mixing — buoyantSimpleFoam
│   ├── 11_boat_hull_fixed/                Fixed hull in water — interFoam + snappyHexMesh
│   ├── 12_boat_hull_floating/             Floating body dynamics — interFoam + 6-DoF
│   └── 13_window_ventilation/             Open versus tilted window — buoyancy + air replacement
│
├── tools/
│   └── make_reel.py                    ← 🎬 Vertical reels/shorts from any finished case
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
 │  01 Lid-Driven       │         │  03 Cylinder Flow    │        │  06 Ahmed Body       │
 │     Cavity           │         │  04 NACA Airfoil     │        │  07 Backward Step    │
 │  02 Elbow Flow       │         │  05 Heated Pipe      │        │  08 Dam Break (VOF)  │
 │                      │         │                      │        │  09 Wind Turbine     │
 │                      │         │                      │        │  10 T-Junction Mix   │
 │                      │         │                      │        │  11 Boat Hull Fixed  │
 │                      │         │                      │        │  12 Boat Hull Float  │
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

Thirteen simulation projects arranged by increasing complexity. Each project directory contains the complete OpenFOAM case structure ready to run.

| # | Project | Description | Solver | Mesh | Turbulence | Level |
|---|---------|-------------|--------|------|------------|-------|
| 01 | [**Lid-Driven Cavity**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/01_lid_driven_cavity) | The "Hello World" of CFD — a square cavity with a moving top wall. Validates against Ghia et al. benchmark data. | `icoFoam` | `blockMesh` | Laminar | 🟢 Beginner |
| 02 | [**Elbow**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/02_elbow) | Flow through a curved pipe elbow. Demonstrates importing a Fluent mesh into OpenFOAM and post-processing in ParaView. | `icoFoam` | Fluent import | Laminar | 🟢 Beginner |
| 03 | [**Flow Around a Cylinder**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/03_flow_around_cylinder) | External flow past a circular cylinder at Re = 200. Introduces importing an externally generated Gmsh mesh (`gmshToFoam`), vortex shedding and drag/lift/Strouhal evaluation. | `pimpleFoam` | Gmsh import | Laminar | 🟡 Intermediate |
| 04 | [**NACA Airfoil Analysis**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis) | Steady RANS of a NACA 0012 airfoil at a fixed angle of attack on a Gmsh-generated mesh. Covers lift/drag/moment coefficients and pressure distribution. | `simpleFoam` | Gmsh import | Spalart-Allmaras | 🟡 Intermediate |
| 05 | [**Heated Pipe Flow**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/05_heated_pipe_flow) | Turbulent heated pipe flow (single region, not conjugate). Introduces the energy equation, buoyancy effects, thermal boundary conditions and an energy-balance check. | `buoyantSimpleFoam` | `blockMesh` | k-epsilon | 🟡 Intermediate |
| 06 | [**Ahmed Body Aerodynamics**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/06_ahmed_body_aerodynamics) | Full 3D automotive aerodynamics simulation of the Ahmed body — the standard reference geometry for vehicle aero studies. | `simpleFoam` | `snappyHexMesh` | k-omega SST | 🔴 Advanced |
| 07 | [**Backward-Facing Step**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/07_backward_facing_step) | Classic turbulent separated flow benchmark. Reattachment length is measured from the wall shear stress (Xr/H $\approx 9.0$ for this setup, kept as a regression baseline). | `simpleFoam` | `blockMesh` | k-omega SST | 🔴 Advanced |
| 08 | [**Dam Break**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/08_dam_break) | Multiphase free-surface flow using Volume of Fluid (VOF). Water column collapse with gravity — transient interFoam simulation. | `interFoam` | `blockMesh` | Laminar | 🔴 Advanced |
| 09 | [**Wind Turbine Blade**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/09_wind_turbine_blade) | Demonstrates the Multiple Reference Frame (MRF) approach with a rotating cell zone (no blades). Introduces `topoSet` and `MRFProperties`. | `simpleFoam` | `blockMesh` | k-epsilon | 🔴 Advanced |
| 10 | [**Mixing T-Junction**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/10_mixing_t_junction) | Thermal mixing of hot (350 K) and cold (300 K) streams in a T-junction. Multi-block structured mesh with buoyancy coupling. | `buoyantSimpleFoam` | `blockMesh` | k-epsilon | 🔴 Advanced |
| 11 | [**Boat Hull — Fixed**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/11_boat_hull_fixed) | Fixed toy boat hull in calm water with free surface. Uses your own STL geometry with snappyHexMesh for drag and wave pattern analysis. | `interFoam` | `snappyHexMesh` | Laminar | 🔴 Advanced |
| 12 | [**Boat Hull — Floating**](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/12_boat_hull_floating) | Floating rigid-body boat with 6-DoF motion. Dynamic mesh morphing with constrained heave/yaw motion — the capstone project. | `interFoam` | `snappyHexMesh` | Laminar | 🔴 Advanced |
| 13 | [**Window Ventilation**](projects/13_window_ventilation) | Transient comparison of a wide-open window and a tilted sash, with a radiator, an outdoor region, and passive-tracer air-replacement measurements. | `buoyantBoussinesqPimpleFoam` | Gmsh import | Laminar (2D teaching model) | 🔴 Advanced |

### Project Progression Map

```
  FUNDAMENTALS           EXTERNAL AERO          HEAT TRANSFER          ADVANCED PHYSICS
  ┌──────────┐          ┌──────────┐           ┌──────────┐          ┌──────────┐
  │ 01 Cavity│          │03 Cylinder│           │05 Heated│           │08 Dam    │
  │ icoFoam  │─────────▶│ pimpleFoam│──────────▶│   Pipe   │     ┌──▶│  Break   │
  │ blockMesh│          │ Gmsh mesh │           │ buoyant  │     │   │ interFoam│
  └──────────┘          └──────────┘           └──────────┘     │   │ VOF      │
       │                     │                      │            │   └──────────┘
       ▼                     ▼                      ▼            │        │
  ┌──────────┐          ┌──────────┐           ┌──────────┐     │        ▼
  │ 02 Elbow │          │04 Airfoil │           │10 T-Junct│     │   ┌──────────┐
  │ icoFoam  │          │ simpleFoam│           │ buoyant  │     └──▶│11 Boat   │
  │ Fluent   │          │ S-A model │           │ thermal  │         │ Fixed    │
  └──────────┘          └──────────┘           └──────────┘         │interFoam │
                             │                                       └──────────┘
                             ▼                                            │
                        ┌──────────┐           ┌──────────┐              ▼
                        │06 Ahmed  │           │07 BFS    │         ┌──────────┐
                        │ Body 3D  │           │ k-ω SST  │         │12 Boat   │
                        │ k-ω SST  │           │validation│         │ Floating │
                        └──────────┘           └──────────┘         │interFoam │
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

Project [13 — Window Ventilation](projects/13_window_ventilation) extends the heat-transfer
branch with transient natural convection and passive-scalar measurements.

---

## 🚀 Quick Start

Want to run your first CFD simulation right now? Here's the fastest path.

### Prerequisites

- **A matching OpenFOAM environment** (the project launchers target Foundation 6) — see [Installation](#installation) below
- **ParaView** for visualization (recommended)
- Basic terminal/command-line familiarity

### Your First Simulation in 5 Commands

Run the classic lid-driven cavity case:

```bash
# 1. Clone this repository
git clone https://github.com/djeada/OpenFoam-Tutorials.git
cd OpenFoam-Tutorials/projects/01_lid_driven_cavity

# 2. Generate the mesh
blockMesh

# Inspect mesh quality before solving
checkMesh

# 3. Run the solver
icoFoam

# 4. Visualize the results
paraFoam
```

That's it! You should see velocity contours of the lid-driven cavity flow.

### Running every project with Docker

You do not need a native OpenFOAM install. Each project has a `run_docker.sh` that runs the case in
the `cfdengine/openfoam` image (OpenFOAM 6) and finishes with `python3 verify.py`, which prints
PASS/FAIL checks against the targets in that project's README:

```bash
cd projects/01_lid_driven_cavity
./run_docker.sh        # clean, run the case in Docker, then verify
./open_paraview.sh     # optional: open the results in ParaView
```

Requirements on the host: `docker` and `python3`. Projects **03**, **04**, and **13** also need `gmsh`
(the image has no Gmsh, so their meshes are generated on the host). Gmsh 4.12.x, which is the
`apt` package on Ubuntu 24.04, crashes on the project 03 geometry. Use 4.11.1, 4.13.1 or 4.15.2,
for example `python3 -m venv ~/gmsh-venv && ~/gmsh-venv/bin/pip install gmsh==4.15.2`, and put
`~/gmsh-venv/bin` first in `PATH`. Rough run times: seconds to a few minutes for 01, 02, 04, 05,
07, 08, 09 and 10; 20 to 30 minutes for 03 and 11; hours for the full runs of 06 and 12. Each
project README lists its measured time where available. Project **13** takes about 30 minutes for
its default 300-second run (both cases in parallel) and also provides a 2-second smoke-run command.

### Turning results into reels / shorts

`tools/make_reel.py` renders a vertical 1080×1920 MP4 (Instagram Reels, YouTube Shorts,
TikTok) from any project that has finished results. Each reel has an intro card, about 35
seconds of the simulation played back smoothly (fields are interpolated between written
times), timed captions, and a closing card. It uses a calm ambient soundtrack generated by
`tools/ambient_music.py`, so there are no licensing issues, and frames are rendered at 2× and
downscaled. It runs under ParaView's `pvbatch` and needs `ffmpeg` with libx264:

```bash
pvbatch tools/make_reel.py projects/08_dam_break                 # one project
pvbatch tools/make_reel.py projects/03_flow_around_cylinder --field vorticity --range -2 2
pvbatch tools/make_reel.py --all                                  # every project with results
```

Videos are written to `output/artifacts/<project>/<reel>.mp4`, which git ignores. Without
configuration the script picks a field (`alpha.water`, `T`, `U`, `p_rgh`, `p`), shows 2D cases
face-on and slices 3D cases through the centre, and fixes one colour range for the whole clip.
The intro text defaults to the first sentence of the project README.
`--field`, `--range`, `--colormap`, `--normal`/`--up`, `--zoom`/`--center`, `--duration`,
`--skip-before`/`--skip-after` and `--streamlines N` override the defaults. `--field vorticity`
colours by how fast the fluid spins about the viewing axis (computed from `U`), which makes
vortex streets and shear layers stand out. `--streamlines N` draws thin white streamlines from
an N-wide grid of fixed seeds; in `reel.json`, `"streamlines": {"lines": [[start, end, seeds]]}`
seeds them along lines instead, for example across an inlet or from a vortex centre to a wall.
`--lic` adds a flow-texture overlay, which suits steady or slowly changing flow; in fast
transients it shimmers and greys out the colours. `--music track.mp3` uses your own audio
and `--no-music` makes the video silent. `--supersample 1` gives a fast draft, and `--still`
writes a PNG of the final frame. A project can also ship a `reel.json` with stacked comparison
panels, labels, live numbers from `postProcessing`, captions timed to simulation time, and
intro/outro text that can quote final values. See
[project 13](projects/13_window_ventilation/reel.json) for an example, and projects
[01](projects/01_lid_driven_cavity/reel.json), [02](projects/02_elbow/reel.json) and
[03](projects/03_flow_around_cylinder/reel.json) for streamlines and vorticity. Interpolation
cross-fades between written times, so fast events such as splashes look better when results
are written more often. Use `--no-interpolate` if the mesh topology changes between time
steps. Add `--force-offscreen-rendering` to `pvbatch` on headless machines.

> **📌 Next Steps:** Read through [Note 01 — Intro to CFD](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/01_short_intro_to_cfd.md) to understand what just happened, then explore [Note 02 — Case Structure](https://github.com/djeada/OpenFoam-Tutorials/blob/main/notes/02_openfoam_cases.md) to understand *why* the directories are organized the way they are.

---

## Installation

There are two recommended ways to get OpenFOAM up and running.

<a id="option-a-direct-installation-ubuntudebian"></a>

### Option A: Direct Installation (Release-Specific)

Follow the [Foundation’s official installation instructions](https://openfoam.org/download/)
for your release and supported operating system. The former example installed
OpenFOAM 11; that is a different interface from these Foundation 6 cases (modern
Foundation releases use [modular solver workflows](https://doc.cfd.direct/openfoam/user-guide-v11/solvers-modules)). Do not assume `simpleFoam`
exists in every release, or combine one release’s package instructions with
another release’s dictionaries.

For an **existing Foundation 6 installation** at the path used by the launchers:

```bash
source /opt/openfoam6/etc/bashrc
echo "$WM_PROJECT_VERSION"
blockMesh -help
icoFoam -help
simpleFoam -help
```

Source the matching environment in each shell that runs a case. You may add that
line to your shell startup file once you have confirmed the installed path.
If your host cannot run the required release natively, use the project’s Docker
workflow. Package availability depends on the OS; there is no universal apt
command that installs a compatible release on every Ubuntu/Debian version.

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

For SELinux label denials, relabel a dedicated simulation directory for the
container using `:Z` (`:z` when deliberately shared by multiple containers):

```bash
docker container run -ti --rm \
  -v "$HOME/openfoam_data:/home/foam:Z" \
  cfdengine/openfoam
```

`--privileged` grants broad host/device access and is not the usual fix for a
bind-mount label problem. Limit relabeling to your simulation directory; also
check normal ownership/permissions. See [Docker’s bind-mount documentation](https://docs.docker.com/engine/storage/bind-mounts/).
The image name in the launchers is unpinned; record the actual image digest and
`WM_PROJECT_VERSION` when reproducibility matters.

### Setting Up the Container Environment

Once inside the Docker container:

```bash
# 1. Source the OpenFOAM environment
#    Check your version first:
ls /opt

#    Then source it (adjust version number as needed):
source /opt/openfoam6/etc/bashrc

# 2. Verify the FOAM_RUN variable
echo "$FOAM_RUN"

#    If empty, set it manually:
export FOAM_RUN=/home/foam

# 3. Create and enter the working directory
mkdir -p "$FOAM_RUN" && cd "$FOAM_RUN"

# 4. Set the OpenFOAM installation directory
export OPENFOAM_DIR=/opt/openfoam6/
```

### Verify Your Installation

In the matching Foundation 6 shell, run these checks. Use a fresh destination
for the tutorial copy to avoid overwriting an existing cavity case:

```bash
# Check that OpenFOAM is accessible
blockMesh -help

# Check solver availability
icoFoam -help
simpleFoam -help

# Run the built-in cavity tutorial as a smoke test
mkdir -p "$FOAM_RUN"
cp -r "$FOAM_TUTORIALS/incompressible/icoFoam/cavity/cavity" "$FOAM_RUN/"
cd "$FOAM_RUN/cavity"
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
| `mpirun -np N solver -parallel` | Run solver in parallel | Large cases needing multiple CPUs |
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
        └── Need resolved turbulent eddies? → LES (costly)
```

Laminar transient solvers can capture shedding at suitable Reynolds numbers
(e.g. project 03); URANS can capture modeled large-scale unsteadiness. LES is
a choice for resolving turbulent structures and needs adequate 3D mesh/time
resolution, not a prerequisite for every shedding calculation.

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
