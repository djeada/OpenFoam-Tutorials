# 02 — Elbow Flow (icoFoam)

![Difficulty: Beginner–Intermediate](https://img.shields.io/badge/difficulty-beginner--intermediate-yellowgreen)
![Solver: icoFoam](https://img.shields.io/badge/solver-icoFoam-blue)
![Mesh: Fluent Import](https://img.shields.io/badge/mesh-Fluent%20import-orange)
![Flow: Laminar Transient](https://img.shields.io/badge/flow-laminar%20transient-lightgrey)
![Algorithm: PISO](https://img.shields.io/badge/algorithm-PISO-yellow)

This tutorial simulates **incompressible, laminar, transient** flow through a
two-inlet elbow geometry using the **icoFoam** solver and the **PISO** algorithm.
It builds on the cavity tutorial by introducing an imported unstructured mesh,
multiple inlet boundaries, and Fluent interoperability tools, making it the natural
second step in learning OpenFOAM.

---

## Table of Contents

1. [Problem Description](#problem-description)
2. [Physics and Theory](#physics-and-theory)
3. [Case Structure](#case-structure)
4. [Mesh Details](#mesh-details)
5. [Boundary Conditions](#boundary-conditions)
6. [Fluid Properties](#fluid-properties)
7. [Numerical Schemes and Solver Settings](#numerical-schemes-and-solver-settings)
8. [Automation Scripts](#automation-scripts)
9. [How to Run](#how-to-run)
10. [Expected Results](#expected-results)
11. [Key Differences from the Cavity Tutorial](#key-differences-from-the-cavity-tutorial)
12. [Exercises](#exercises)
13. [References](#references)

---

## Problem Description

The **elbow flow** problem models two fluid streams entering a curved channel from
perpendicular directions, mixing in a shared passage, and exiting through a common
outlet. This configuration is ubiquitous in piping systems, HVAC ductwork, and
chemical-process manifolds, making it a practical step up from the simple lid-driven
cavity.

### Why Study Elbow Flow?

- **Curved-geometry effects** — pressure gradients across the bend drive secondary
  flow patterns (Dean vortices) that do not appear in straight ducts.
- **Multi-inlet mixing** — two streams with different velocity vectors collide,
  creating shear layers, recirculation zones, and complex transient behaviour.
- **Pressure drop estimation** — engineers must predict losses through bends to size
  pumps and blowers correctly.
- **External mesh import** — unlike the cavity case (which uses `blockMesh`), this
  case imports a Fluent `.msh` file, introducing a real-world meshing workflow.
- **Fluent round-trip** — the case demonstrates exporting results back to Fluent
  format, a common requirement in mixed-tool engineering teams.

### Geometry Overview

```
                     pressure-outlet-7
                      ┌───────────┐
                      │  ↑  ↑  ↑  │
                      │  ↑  ↑  ↑  │
                      │           │
              ┌───────┘           └───────┐
              │                           │
              │        Elbow Region       │
              │     (mixing + turning)    │
              │                           │
    ──────────┘    ←── secondary ──→      │
                       flow               │
    velocity-inlet-5                      │  wall-8
    ─→ ─→ ─→ ─→ ─→ (1, 0, 0)            │
    ──────────┐                           │
      wall-4  │                           │
              │           ↑  ↑  ↑         │
              └───────┐   ↑  ↑  ↑   ┌────┘
                      │   ↑  ↑  ↑   │
                      │             │
                      └──────┬──────┘
                             │
                      velocity-inlet-6
                        (0, 3, 0)

    frontAndBackPlanes: empty (2-D constraint)
```

**What is happening physically:**

- **velocity-inlet-5** drives fluid horizontally at 1 m/s into the lower-left
  branch of the elbow.
- **velocity-inlet-6** drives fluid vertically at 3 m/s into the bottom branch.
- The two streams collide in the elbow region, exchange momentum, and exit upward
  through **pressure-outlet-7**.
- **wall-4** and **wall-8** are solid no-slip walls bounding the channel.
- **frontAndBackPlanes** enforces a 2-D (single-cell-thick) simulation by applying
  the `empty` boundary condition.

The velocity ratio of 3:1 means the vertical stream dominates the combined flow,
deflecting the horizontal stream strongly upward. Where the two streams meet, a
shear layer forms. Momentum diffusion (controlled by the kinematic viscosity $\nu$)
slowly smooths the velocity differences over time, but the transient development of
this mixing region is one of the key features to observe.

This setup is representative of many practical engineering situations:

| Application area      | Example                                              |
|-----------------------|------------------------------------------------------|
| Piping systems        | Mixing tees, manifolds, header pipes                 |
| HVAC                  | Duct junctions where supply streams merge            |
| Chemical engineering  | Inline mixing of reagent streams                     |
| Hydraulics            | River confluences, irrigation channel merges         |

---

## Physics and Theory

### Governing Equations

`icoFoam` solves the incompressible Navier–Stokes equations for a Newtonian fluid
with constant kinematic viscosity:

$$\nabla \cdot U = 0$$ (continuity)

$$\partial U / \partial t + \nabla \cdot (UU) = -\nabla p + \nu \nabla^2 U$$ (momentum)

where **U** is the velocity vector, **p** is the kinematic pressure (pressure
divided by density, p/$\rho$), and **$\nu$** is the kinematic viscosity. Because the fluid
is incompressible and the viscosity is constant, no energy equation or equation of
state is needed.

### Secondary Flow and the Dean Number

When fluid negotiates a bend, centrifugal effects push the faster-moving fluid
toward the outer wall while the slower fluid near the walls recirculates inward.
This establishes a pair of counter-rotating vortices (Dean vortices) superimposed
on the primary axial flow. The strength of this secondary flow is characterised by
the **Dean number**:

$$De = Re \times \sqrt{D / 2R}$$

where `Re` is the Reynolds number, `D` is the hydraulic diameter of the duct, and
`R` is the bend radius of curvature. Higher Dean numbers produce stronger secondary
flows and larger pressure losses through the bend.

### Mixing of Two Inlet Streams

The two inlet velocities differ in both magnitude and direction. Where the streams
meet, a shear layer forms. Viscous diffusion and convective transport gradually
blend the two streams. The time scale for complete mixing depends on the velocity
ratio, the channel geometry, and the Reynolds number.

### Reynolds Number Estimate

With $\nu = 0.01$ m$^2$/s, a characteristic velocity of ~3 m/s (the faster inlet), and a
channel width on the order of 1 m, the Reynolds number is roughly:

$$Re = U \times L / \nu \approx 3 \times 1 / 0.01 = 300$$

This is well within the laminar regime (Re ≪ 2300 for internal flow), which
justifies the use of `icoFoam` without a turbulence model.

> **See also:**
> [notes/04_meshing.md](../../notes/04_meshing.md) for general meshing theory,
> [notes/05_boundary_conditions.md](../../notes/05_boundary_conditions.md) for a
> deeper dive into boundary condition types, and
> [notes/10_icofoam_solver_analysis.md](../../notes/10_icofoam_solver_analysis.md)
> for a detailed walkthrough of the icoFoam solver algorithm.

---

## Case Structure

All filenames in this project are **lowercase** — this is important when editing
files or writing scripts. Note in particular: `controldict`, `fvschemes`,
`fvsolution`, `transportproperties`, `foamdatatofluentdict`, `allrun`, `allclean`.

```
02_elbow/
├── 0/
│   ├── p                        # Pressure initial & boundary conditions
│   └── u                        # Velocity initial & boundary conditions
├── constant/
│   ├── polymesh/
│   │   ├── boundary             # Patch definitions (6 patches)
│   │   ├── cellZones            # Cell zone information
│   │   ├── faceZones            # Face zone information
│   │   ├── faces                # Face-to-vertex connectivity
│   │   ├── neighbour            # Neighbour cell for each internal face
│   │   ├── owner                # Owner cell for each face
│   │   ├── pointZones           # Point zone information
│   │   └── points               # Vertex coordinates
│   └── transportproperties      # Kinematic viscosity (nu = 0.01)
├── system/
│   ├── controldict              # Time control, write interval, I/O
│   ├── foamdatatofluentdict     # Field-to-Fluent-ID export mapping
│   ├── fvschemes                # Discretisation schemes
│   └── fvsolution               # Linear solvers, PISO settings
├── allclean                     # Script: remove all generated files
├── allrun                       # Script: full automated run
├── elbow.msh                    # Fluent-format mesh (~32 KB)
└── readme.md                    # This file
```

---

## Mesh Details

### External Mesh — Fluent Format

Unlike the cavity tutorial (which generates its mesh internally with `blockMesh`),
the elbow geometry is too complex for simple block decomposition. The mesh was
originally created in ANSYS Fluent's meshing environment, which excels at generating
unstructured meshes for curved and irregular domains.

The source file `elbow.msh` (~32 KB) is included in the case directory and must be
converted to OpenFOAM's native `polyMesh` format before running:

```bash
fluentMeshToFoam elbow.msh
```

This command:

1. Reads the Fluent mesh (nodes, faces, cells, boundary zones).
2. Reorders cells for OpenFOAM's face-based addressing scheme.
3. Writes `constant/polymesh/` with `points`, `faces`, `owner`, `neighbour`, and
   `boundary` files.
4. Maps Fluent zone names (e.g. `velocity-inlet-5`) to OpenFOAM patch names.

### Mesh Statistics

The converted mesh has the following characteristics (from `constant/polymesh/boundary`):

| Patch                | Type    | Faces |
|----------------------|---------|------:|
| wall-4               | wall    |   100 |
| velocity-inlet-5     | patch   |     8 |
| velocity-inlet-6     | patch   |     4 |
| pressure-outlet-7    | patch   |     8 |
| wall-8               | wall    |    34 |
| frontAndBackPlanes   | empty   |  1836 |

- **Internal faces:** 1300
- **Total cells:** 918 (2-D; each cell has one front and one back `empty` face,
  so 1836 $\div$ 2 = 918)
- **Mesh type:** Unstructured (triangular/quadrilateral mix typical of Fluent meshes)

### Unstructured Mesh Characteristics

Because the mesh is unstructured, the solver must handle:

- **Non-orthogonality** — corrected by the 2 non-orthogonal corrector loops in PISO
  (configured in `fvsolution`).
- **Variable cell size** — finer near walls and in the mixing zone, coarser away
  from regions of interest.
- **Arbitrary connectivity** — each cell may have a different number of neighbours,
  unlike the regular structure produced by `blockMesh`.

### Mesh Quality Check

After converting, always verify mesh quality:

```bash
checkMesh
```

Look for the summary at the end: no failed checks should be reported. Key metrics
to inspect include non-orthogonality (the 2 non-orthogonal correctors help manage
moderate values), skewness, and aspect ratio.

> **Further reading:**
> [notes/04_meshing.md](../../notes/04_meshing.md) covers meshing strategies,
> quality metrics, `blockMesh` vs external tools, and troubleshooting.

---

## Boundary Conditions

The boundary conditions are defined in `0/u` (velocity) and `0/p` (pressure).

### Velocity Field — `0/u`

| Patch                | BC Type        | Value / Setting          | Physical Meaning                           |
|----------------------|----------------|--------------------------|---------------------------------------------|
| wall-4               | `noSlip`       | —                        | Zero velocity at the wall (no-slip)         |
| velocity-inlet-5     | `fixedValue`   | `uniform (1 0 0)` m/s   | Horizontal inflow at 1 m/s                  |
| velocity-inlet-6     | `fixedValue`   | `uniform (0 3 0)` m/s   | Vertical inflow at 3 m/s                    |
| pressure-outlet-7    | `zeroGradient` | —                        | Velocity extrapolated to outlet             |
| wall-8               | `noSlip`       | —                        | Zero velocity at the wall (no-slip)         |
| frontAndBackPlanes   | `empty`        | —                        | 2-D constraint (no solution in z-direction) |

### Pressure Field — `0/p`

| Patch                | BC Type        | Value / Setting     | Physical Meaning                                  |
|----------------------|----------------|---------------------|---------------------------------------------------|
| wall-4               | `zeroGradient` | —                   | No pressure flux through the wall                 |
| velocity-inlet-5     | `zeroGradient` | —                   | Pressure adjusts to satisfy continuity at inlet   |
| velocity-inlet-6     | `zeroGradient` | —                   | Pressure adjusts to satisfy continuity at inlet   |
| pressure-outlet-7    | `fixedValue`   | `uniform 0`         | Reference pressure fixed at the outlet            |
| wall-8               | `zeroGradient` | —                   | No pressure flux through the wall                 |
| frontAndBackPlanes   | `empty`        | —                   | 2-D constraint                                    |

### Internal Field Initialisation

Both fields are initialised to uniform zero:

- **p:** `internalField uniform 0;`
- **U:** `internalField uniform (0 0 0);`

This means the simulation starts from quiescence and develops flow as the inlet
conditions take effect from the first time step.

### Key Points

- Fixing pressure at the outlet and velocity at the inlets is a **well-posed**
  combination for incompressible flow. Never fix both p and U on the same boundary.
- The `noSlip` boundary for walls is a shorthand that sets velocity to
  `uniform (0 0 0)` without requiring an explicit `value` entry.
- The `empty` type tells OpenFOAM to ignore those faces entirely — essential for 2-D
  cases modelled as a single-cell-thick 3-D slice.

---

## Fluid Properties

Defined in `constant/transportproperties`:

```
nu    [0 2 -1 0 0 0 0]    0.01;
```

| Property            | Symbol | Value   | Dimensions | Dimensional formula        |
|---------------------|--------|---------|------------|----------------------------|
| Kinematic viscosity | $\nu$      | 0.01    | m$^2$/s       | [0 2 -1 0 0 0 0]          |

This is a relatively high viscosity (for context, water at 20 °C has
$\nu \approx 1 \times 10^{-6}$ m$^2$/s). The large value keeps the Reynolds number low (~300),
ensuring the flow remains laminar and the simulation converges without a
turbulence model.

---

## Numerical Schemes and Solver Settings

### Time Control — `system/controldict`

| Parameter          | Value         | Notes                                       |
|--------------------|---------------|---------------------------------------------|
| `application`      | `icoFoam`     | Incompressible laminar transient solver      |
| `startFrom`        | `latestTime`  | Resume from the last written time directory  |
| `startTime`        | `0`           | Fallback if no time directories exist        |
| `endTime`          | `75`          | Simulation runs for 75 seconds               |
| `deltaT`           | `0.05`        | Fixed time step of 50 ms                     |
| `writeControl`     | `timeStep`    | Write output every N time steps              |
| `writeInterval`    | `20`          | Write every 20 steps → every 1.0 s           |
| `purgeWrite`       | `0`           | Keep all time directories                    |
| `writeFormat`      | `ascii`       | Human-readable output                        |
| `writePrecision`   | `6`           | 6 significant digits                         |

Total time steps: 75 / 0.05 = **1500**. With `writeInterval 20`, this produces
**75 output directories** (plus the initial `0` directory).

### Discretisation Schemes — `system/fvschemes`

| Category               | Scheme                      | Notes                                    |
|------------------------|-----------------------------|------------------------------------------|
| `ddtSchemes`           | `Euler`                     | First-order implicit, unconditionally stable |
| `gradSchemes`          | `Gauss linear`              | Second-order, cell-centred gradient      |
| `divSchemes`           | `Gauss limitedLinearV 1`    | TVD-limited linear for `div(phi,U)`      |
| `laplacianSchemes`     | `Gauss linear corrected`    | Second-order with non-orthogonal correction |
| `interpolationSchemes` | `linear`                    | Linear (central) face interpolation      |
| `snGradSchemes`        | `corrected`                 | Explicit non-orthogonal correction       |

The `limitedLinearV 1` divergence scheme applies a TVD slope limiter to the
convective term, preventing spurious oscillations that can arise with pure central
differencing on the unstructured mesh — particularly near the inlet junction where
velocity gradients are steep.

### Linear Solvers and PISO Algorithm — `system/fvsolution`

**Pressure solver (`p`):**

| Setting          | Value       |
|------------------|-------------|
| Solver           | `PCG`       |
| Preconditioner   | `DIC`       |
| Tolerance        | 1e-06       |
| Relative tol.    | 0.05        |

**Final pressure corrector (`pFinal`):** Same as `p` but with `relTol 0` — solved
to the absolute tolerance on the last PISO corrector step.

**Velocity solver (`U`):**

| Setting          | Value             |
|------------------|-------------------|
| Solver           | `smoothSolver`    |
| Smoother         | `symGaussSeidel`  |
| Tolerance        | 1e-05             |
| Relative tol.    | 0                 |

**PISO algorithm settings:**

| Parameter                       | Value | Purpose                                      |
|---------------------------------|------:|----------------------------------------------|
| `nCorrectors`                   |     2 | Pressure–velocity coupling correction steps  |
| `nNonOrthogonalCorrectors`      |     2 | Extra corrections for mesh non-orthogonality |

Two PISO correctors enforce the pressure–velocity coupling within each time step.
The two non-orthogonal correctors improve accuracy on the unstructured mesh, where
face normals do not perfectly align with the line connecting neighbouring cell
centres.

### Fluent Export Dictionary — `system/foamdatatofluentdict`

This file maps OpenFOAM field names to Fluent variable IDs so that results can be
exported back to Fluent format using `foamDataToFluent`:

| OpenFOAM Field | Fluent ID |
|----------------|----------:|
| `p`            |         1 |
| `U`            |         2 |
| `T`            |         3 |
| `h`            |         4 |
| `k`            |         5 |
| `epsilon`      |         6 |
| `alpha1`       |       150 |

Only `p` and `U` are active in this case; the other entries are placeholders that
allow the dictionary to be reused in more complex simulations.

---

## Automation Scripts

### `allrun` — Run the Full Case

```bash
#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Get application directory
application=$(getApplication)

runApplication fluentMeshToFoam elbow.msh   # 1. Convert Fluent mesh → polyMesh
runApplication "$application"               # 2. Run icoFoam
runApplication foamMeshToFluent             # 3. Convert mesh back to Fluent
runApplication foamDataToFluent             # 4. Export field data to Fluent
```

**Step-by-step breakdown:**

1. **`fluentMeshToFoam elbow.msh`** — reads the Fluent `.msh` file and writes
   OpenFOAM's `constant/polymesh/` directory.
2. **`icoFoam`** — runs the transient solver for 1500 time steps (75 s at $\Delta t$ = 0.05).
   The `getApplication` function reads the solver name from `controldict`.
3. **`foamMeshToFluent`** — converts the OpenFOAM mesh back to Fluent format and
   writes it to `fluentInterface/`.
4. **`foamDataToFluent`** — exports field data (p, U, etc.) using the IDs from
   `foamdatatofluentdict`.

The `runApplication` wrapper logs each command's output to a `log.<command>` file
and **skips** re-running a step if its log file already exists. To force a re-run,
delete the relevant log file or run `./allclean` first.

### `allclean` — Reset the Case

```bash
#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial clean functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions

rm -f constant/polyMesh/boundary > /dev/null 2>&1
rm -rf fluentInterface
cleanCase
```

This removes:

- All generated time-step directories (1/, 2/, … 75/)
- All `log.*` files from each utility
- The `constant/polyMesh/boundary` file (regenerated by `fluentMeshToFoam`)
- The `fluentInterface/` export directory
- Processor directories (from parallel runs, if any)
- Any VTK output

After running `allclean`, the case is back to its original state and ready for a
fresh run.

---

## How to Run

### Prerequisites

- OpenFOAM 6+ installed and sourced (`. $WM_PROJECT_DIR/etc/bashrc`).
- `icoFoam`, `fluentMeshToFoam`, `checkMesh`, and `paraFoam` (or `foamToVTK` +
  ParaView) available on your `$PATH`.

### Quick Start (Automated)

```bash
cd projects/02_elbow
./allclean            # ensure a clean state
./allrun              # convert mesh → solve → export to Fluent
```

Check progress by tailing the solver log:

```bash
tail -f log.icoFoam
```

### Manual Step-by-Step

If you prefer to run each stage individually for learning purposes:

```bash
# 1. Navigate to the case directory
cd projects/02_elbow

# 2. Clean any previous results
./allclean

# 3. Convert the Fluent mesh to OpenFOAM format
fluentMeshToFoam elbow.msh

# 4. Inspect the mesh (optional but recommended)
checkMesh

# 5. Run the solver
icoFoam

# 6. View results in ParaView
#    Option A — native OpenFOAM reader:
paraFoam
#    Option B — convert to VTK first:
foamToVTK
paraview VTK/

# 7. (Optional) Export results back to Fluent format
foamMeshToFluent
foamDataToFluent
```

### Monitoring Convergence

While `icoFoam` runs, watch the residuals printed to the terminal. For a healthy
simulation the initial residuals for `p` and `Ux`/`Uy` should decrease over the
early time steps and then stabilise. If residuals diverge, reduce `deltaT` or
check boundary conditions.

You can also extract residual data from the log file for plotting:

```bash
foamLog log.icoFoam
# Creates a logs/ directory with per-field residual data
# Plot with gnuplot, matplotlib, or your preferred tool
```

### Post-Processing with ParaView

In ParaView:

1. Open the case directory with the built-in OpenFOAM reader, or load VTK files
   from the `VTK/` directory after running `foamToVTK`.
2. Apply a **Surface** representation coloured by velocity magnitude (`U`) to
   visualise the flow field.
3. Use **Glyph** filter on `U` to see velocity vectors and flow direction.
4. Apply a **Contour** or change colouring to `p` to visualise pressure distribution.
5. Create a **Plot Over Line** from the inner to outer wall at the bend midpoint
   to see the velocity profile and pressure gradient across the duct.
6. Use **Stream Tracer** to follow fluid paths through the elbow.

### Cleaning Up

To reset the case to its original state:

```bash
./allclean
```

---

## Expected Results

### Flow Development

As the simulation evolves from t = 0 to t = 75 s, expect the following phases:

1. **Early transient (0–5 s):** Both inlet streams begin to fill the domain. The
   faster vertical stream (3 m/s from inlet-6) reaches the outlet region first.
   A clear shear layer forms where the horizontal and vertical streams meet.

2. **Development phase (5–30 s):** The flow negotiates the bend. Secondary flow
   structures develop. Recirculation bubbles appear near the inner wall of the
   elbow where the horizontal stream curves upward. The faster vertical stream
   dominates momentum in the mixing zone and deflects the slower horizontal
   stream.

3. **Quasi-steady state (30–75 s):** The velocity field reaches an approximately
   steady pattern. The pressure field shows a clear gradient from the inlets to
   the outlet, with the highest pressure at the outer wall of the bend. The outlet
   velocity profile is skewed toward the wall-8 side.

### Features to Look For

| Feature                | Where to Observe                                  |
|------------------------|---------------------------------------------------|
| Recirculation zone     | Inner wall of the elbow, near the stream junction |
| High-velocity core     | Centre of the outlet passage                      |
| Shear layer            | Interface between the two merging streams         |
| Outer-wall pressure    | Higher pressure on the outer wall of the bend     |
| Velocity overshoot     | Inner wall, just downstream of the bend           |
| Dean vortices          | Cross-sectional views through the curved passage  |
| Gradual mixing         | Downstream of the junction toward the outlet      |

### Quantitative Checks

- **Maximum velocity** should be roughly 3–4 m/s (near the vertical inlet region).
- **Outlet pressure** should be approximately 0 (fixed reference).
- **Mass conservation:** total mass outflow should match the sum of inflows. Verify
  with `postProcess -func 'patchFlowRate(name=velocity-inlet-5)'` (if available in
  your OpenFOAM version).

---

## Key Differences from the Cavity Tutorial

| Aspect                   | Cavity (01)                          | Elbow (02)                                       |
|--------------------------|--------------------------------------|--------------------------------------------------|
| **Geometry**             | Simple square box                    | Curved channel with two inlets                   |
| **Mesh generation**      | `blockMesh` (built-in, structured)   | `fluentMeshToFoam` (external import)             |
| **Mesh type**            | Structured hexahedral                | Unstructured (mixed tri/quad)                    |
| **Number of inlets**     | 0 (lid-driven)                       | 2 (velocity-inlet-5, velocity-inlet-6)           |
| **Number of outlets**    | 0 (enclosed)                         | 1 (pressure-outlet-7)                            |
| **Driving mechanism**    | Moving wall (lid)                    | Prescribed inlet velocities                      |
| **Non-orth. correctors** | 0 (perfect hex grid)                 | 2 (needed for unstructured mesh)                 |
| **Divergence scheme**    | `Gauss linear`                       | `Gauss limitedLinearV 1` (TVD limiter)           |
| **Fluent interop**       | None                                 | Import `.msh` + export via `foamDataToFluent`    |
| **Simulation duration**  | ~0.5 s                               | 75 s (to reach quasi-steady state)               |
| **Key physics**          | Single central vortex                | Dean vortices, multi-stream mixing, recirculation|
| **New skills learned**   | Basic OpenFOAM workflow              | External meshing, Fluent conversion, PISO tuning |

---

## Exercises

These exercises progressively increase in difficulty. Try them in order.

### Exercise 1 — Change Inlet Velocities

Modify `0/u` to set both inlets to the same speed, e.g. `(1 0 0)` and `(0 1 0)`.

- How does the mixing pattern change when the velocity ratio is 1:1 instead of 1:3?
- Does the recirculation zone at the inner wall grow or shrink?
- Is the outlet velocity profile more symmetric?

### Exercise 2 — Reduce Viscosity

In `constant/transportproperties`, change `nu` from `0.01` to `0.001`.

- What happens to the flow? Does it become unsteady or show oscillations?
- You may need to decrease `deltaT` to maintain stability (check CFL number).
- Compare velocity profiles at the outlet for both viscosities.

### Exercise 3 — Increase PISO Correctors

Change `nCorrectors` from 2 to 4 in `system/fvsolution`. Re-run and compare:

- Does the residual history improve (faster convergence per time step)?
- Is the solution noticeably different?
- What is the trade-off in computation time?

### Exercise 4 — Refine the Mesh

Open `elbow.msh` in a meshing tool (e.g. Gmsh or ANSYS Meshing) and generate a
finer mesh with roughly double the cell count. Export as `.msh` and re-run:

```bash
fluentMeshToFoam elbow_fine.msh
icoFoam
```

- Is the solution noticeably different? This is a basic **mesh independence** check.
- Pay attention to the recirculation zone size and the outlet velocity profile.

### Exercise 5 — Try a Different Divergence Scheme

Replace `Gauss limitedLinearV 1` with `Gauss upwind` in `system/fvschemes`.

- How does the additional numerical diffusion from upwind affect the sharpness of
  the shear layer between the two streams?
- Is the recirculation zone still visible, or has it been smeared out?

### Exercise 6 — Add a Passive Scalar

Add a scalar transport equation to track how the fluid from each inlet mixes:

1. Create a file `0/T` (copy structure from `0/u` but use `volScalarField`).
2. Set inlet-5 to `T = 0` and inlet-6 to `T = 1`.
3. Add appropriate `divSchemes`, `laplacianSchemes`, and solver entries for `T`.
4. Use `scalarTransportFoam` or couple a transport equation into the solver.

This exercise goes beyond `icoFoam` and is a stepping-stone to multi-species
transport simulations.

### Exercise 7 — Export and Compare in Fluent

Run the full `allrun` script to generate the Fluent export files in
`fluentInterface/`. If you have access to ANSYS Fluent:

- Import the mesh and field data.
- Compare visualisations and residual values between OpenFOAM and Fluent.
- Are there differences due to the different solver implementations?

---

## References

1. **OpenFOAM User Guide — icoFoam**
   [https://openfoam.org/](https://openfoam.org/) — Official documentation for the
   incompressible laminar transient solver.

2. **OpenFOAM Tutorial: Elbow**
   The original tutorial shipped with OpenFOAM under
   `$FOAM_TUTORIALS/incompressible/icoFoam/elbow/`.

3. **Dean, W. R. (1928)**
   *The stream-line motion of fluid in a curved pipe.* Phil. Mag. Series 7, 5(30),
   pp. 673–695. The foundational paper on secondary flow in curved channels.

4. **Versteeg, H. K. & Malalasekera, W.**
   *An Introduction to Computational Fluid Dynamics: The Finite Volume Method.*
   Pearson, 2nd edition. Excellent reference for the PISO algorithm and finite
   volume discretisation.

5. **Repository notes** — deeper background for concepts used in this tutorial:

   | Note file                        | Relevance                                            |
   |----------------------------------|------------------------------------------------------|
   | [01_short_intro_to_cfd.md](../../notes/01_short_intro_to_cfd.md) | Fundamental CFD concepts, governing equations |
   | [02_openfoam_cases.md](../../notes/02_openfoam_cases.md)         | OpenFOAM case directory structure (0/, constant/, system/) |
   | [03_openfoam_dictionaries.md](../../notes/03_openfoam_dictionaries.md) | Dictionary syntax used in controldict, fvschemes, etc. |
   | [04_meshing.md](../../notes/04_meshing.md)                       | Mesh generation; blockMesh vs imported meshes |
   | [05_boundary_conditions.md](../../notes/05_boundary_conditions.md) | fixedValue, zeroGradient, noSlip, empty types |
   | [08_cfl_number.md](../../notes/08_cfl_number.md)                 | CFL condition and time-step selection |
   | [09_linear_solvers.md](../../notes/09_linear_solvers.md)         | PCG, smoothSolver, DIC, symGaussSeidel |
   | [10_icofoam_solver_analysis.md](../../notes/10_icofoam_solver_analysis.md) | Detailed analysis of the icoFoam algorithm |
