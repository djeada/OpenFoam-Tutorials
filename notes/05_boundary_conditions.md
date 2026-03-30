# Boundary Conditions in OpenFOAM

## 1. Introduction — Why Boundary Conditions Matter

Boundary conditions (BCs) are the **single most important input** to any CFD simulation. They
define how the computational domain interacts with the outside world: where fluid enters, where it
leaves, what happens at walls, and how the simulation "knows" about physical reality.

```
  ┌─────────────────────────────────────────────────────────────────────┐
  │                WHY BOUNDARY CONDITIONS MATTER                      │
  │                                                                     │
  │   ┌──────────────┐     BCs define     ┌──────────────────────┐     │
  │   │  GOVERNING   │ ──────────────────► │  UNIQUE SOLUTION     │     │
  │   │  EQUATIONS   │  what makes the     │  to the PDEs         │     │
  │   │  (N-S eqns)  │  problem solvable   │  (your CFD result)   │     │
  │   └──────────────┘                     └──────────────────────┘     │
  │                                                                     │
  │   Without BCs: infinite possible solutions (ill-posed problem)      │
  │   Wrong BCs:   converged but WRONG solution (garbage in/out)        │
  │   Right BCs:   physically meaningful, accurate results              │
  └─────────────────────────────────────────────────────────────────────┘
```

**From a mathematical perspective**, the Navier-Stokes equations are partial differential equations
(PDEs). A unique solution requires boundary conditions — without them, there are infinitely many
solutions. The type of BC (Dirichlet, Neumann, mixed) determines which mathematical problem you
are actually solving.

**From a physical perspective**, BCs represent real-world constraints: the speed of incoming air,
the pressure at an outlet, the no-penetration condition at a solid wall. Getting them right means
your simulation reflects reality.

> **⚠️ Key Insight**: A CFD simulation with perfect mesh, perfect numerics, but wrong boundary
> conditions will give you a perfectly converged **wrong answer**. Always validate your BCs first.

---

## 1b. Initial Conditions — Setting the Starting Point

While boundary conditions define what happens at the **edges** of your domain, **initial
conditions** (ICs) define the field values at **every cell inside the domain** at the very
start of the simulation (time = 0). Together, BCs and ICs form a complete **initial-boundary
value problem** — the mathematical formulation that makes the Navier-Stokes equations solvable
for a unique solution over time.

```
  ┌───────────────────────────────────────────────────────────────────────┐
  │              INITIAL vs BOUNDARY CONDITIONS                          │
  │                                                                       │
  │   BOUNDARY CONDITIONS (BCs)        INITIAL CONDITIONS (ICs)           │
  │   ─────────────────────────        ────────────────────────           │
  │   WHERE: domain boundaries         WHERE: every internal cell         │
  │   WHEN:  all time steps            WHEN:  t = 0 only                  │
  │   HOW:   patch-by-patch in 0/      HOW:   internalField in 0/        │
  │   WHY:   define interaction        WHY:   give solver a starting      │
  │          with outside world               point for iteration         │
  └───────────────────────────────────────────────────────────────────────┘
```

### How Initial Conditions Work in OpenFOAM

In OpenFOAM, initial conditions are set via the `internalField` entry in each field file
inside the `0/` directory. The same file that defines boundary conditions also defines the
initial conditions:

```c
// File: 0/U
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    object      U;
}

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0);      // ← INITIAL CONDITION: all cells start at rest

boundaryField                          // ← BOUNDARY CONDITIONS below
{
    inlet
    {
        type            fixedValue;
        value           uniform (1 0 0);
    }
    // ...
}
```

The `internalField` keyword accepts two forms:

| Form | Syntax | When to Use |
|------|--------|-------------|
| **uniform** | `internalField uniform (0 0 0);` | Same value everywhere — the most common choice |
| **nonuniform** | `internalField nonuniform List<vector> N ( ... );` | Different values per cell — used when restarting from a previous simulation or when `setFields` or `mapFields` has been applied |

### Why Initial Conditions Matter

**For transient simulations** (e.g., `icoFoam`, `pimpleFoam`, `interFoam`):
- ICs represent the **physical starting state** of the flow at t = 0.
- A poor initial field can cause divergence in the first few time steps.
- Physical consistency matters: if you set `U = (100 0 0)` everywhere but the inlet is `(1 0 0)`,
  the solver will struggle with the massive initial discontinuity.

**For steady-state simulations** (e.g., `simpleFoam`):
- There is no physical "time zero" — the solver iterates toward a converged state.
- ICs serve as the **initial guess** for the iterative solver.
- A good initial guess can dramatically speed up convergence.
- A bad guess will not change the final answer (if the solver converges), but may cause
  divergence or require many more iterations.

### Practical Tips for Initial Conditions

1. **Start simple:** `uniform (0 0 0)` for velocity and `uniform 0` for pressure is fine for
   most problems. The solver will iterate toward the correct solution.

2. **Use `potentialFoam` for better initialization:** Running `potentialFoam` before your main
   solver computes an irrotational velocity field that satisfies continuity. This gives a much
   better starting point, especially for external aerodynamics:
   ```bash
   potentialFoam    # Initialise velocity field
   simpleFoam       # Run the actual solver
   ```

3. **Use `setFields` for non-uniform ICs:** For multiphase simulations (e.g., dam break),
   you need to initialize the volume fraction field (`alpha.water`) with different values in
   different regions:
   ```bash
   setFields    # Sets alpha.water = 1 in the water region, 0 in the air region
   interFoam    # Run the multiphase solver
   ```

4. **Restart from a previous solution:** Copy a time directory from a coarser mesh or a
   similar case. Use `mapFields` to interpolate fields between different meshes:
   ```bash
   mapFields ../coarse_case -sourceTime latestTime
   ```

5. **Turbulent fields need sensible ICs:** If using RANS turbulence models, initialize `k`,
   `epsilon` (or `omega`), and `nut` with physically reasonable values — not zero. See
   [notes/06_turbulence_models.md](06_turbulence_models.md) for how to calculate appropriate
   initial values.

> **⚠️ Common Pitfall:** Setting `internalField uniform 0` for turbulent kinetic energy `k`
> will cause most turbulence models to crash at the very first iteration, because the turbulent
> viscosity `nut ∝ k²/ε` becomes undefined. Always use a small positive value.

---

## 2. How Boundary Conditions Work in OpenFOAM

In OpenFOAM, every boundary face in the mesh is assigned to a **patch**. Each patch then gets a
BC type for **every field variable** (U, p, k, epsilon, nut, etc.) in the `0/` directory.

```
  Physical Domain                    OpenFOAM Mapping
  ┌──────────────────────┐
  │                      │          constant/polyMesh/boundary
  │  INLET ═══►  ►  ►   │            → defines patch names & types
  │                      │
  │      FLUID DOMAIN    │          0/U  (velocity BCs per patch)
  │      (cells)         │          0/p  (pressure BCs per patch)
  │                      │          0/k  (turbulent KE BCs per patch)
  │  WALL ░░░░░░░░░░░░░  │          0/epsilon  (dissipation BCs)
  │                      │          0/nut (turbulent viscosity BCs)
  │  ◄  ◄  ◄  OUTLET    │
  └──────────────────────┘
```

Here is how the mapping works for a typical incompressible flow:

```
  ┌───────────────────────────────────────────────────────────────────────────┐
  │  Physical Boundary        0/U file                0/p file               │
  │  ─────────────────        ──────────              ──────────             │
  │                                                                          │
  │  INLET  ═══►              fixedValue              zeroGradient           │
  │  (fluid enters)           value (Ux Uy Uz)        (pressure floats)      │
  │                                                                          │
  │  OUTLET ◄═══              zeroGradient            fixedValue             │
  │  (fluid leaves)           (velocity floats)       value 0                │
  │                                                                          │
  │  WALL   ░░░░              noSlip                  zeroGradient           │
  │  (solid surface)          (U = 0 at wall)         (no pressure flux)     │
  │                                                                          │
  │  FRONT/BACK               empty                   empty                  │
  │  (2D simulation)          (no variation in z)     (no variation in z)    │
  └───────────────────────────────────────────────────────────────────────────┘
```

### Dimensional Consistency

Every field file starts with a `dimensions` specification using the standard 7-element array:

```
dimensions  [ kg  m  s  K  mol  A  cd ];
```

| Field      | Dimensions              | Meaning                        |
|------------|-------------------------|--------------------------------|
| U          | `[0 1 -1 0 0 0 0]`     | m/s (velocity)                 |
| p (incomp) | `[0 2 -2 0 0 0 0]`    | $m^{2}/s^{2}$ (kinematic pressure)     |
| p (comp)   | `[1 -1 -2 0 0 0 0]`   | $kg/(m \cdot s^{2})$ = Pa      |
| k          | `[0 2 -2 0 0 0 0]`     | $m^{2}/s^{2}$ (turbulent KE)           |
| epsilon    | `[0 2 -3 0 0 0 0]`     | $m^{2}/s^{3}$ (dissipation rate)       |
| nut        | `[0 2 -1 0 0 0 0]`     | $m^{2}/s$ (turbulent viscosity)     |
| omega      | `[0 0 -1 0 0 0 0]`     | 1/s (specific dissipation)     |

> **💡 Tip**: OpenFOAM uses **kinematic** pressure ($p/\rho$) for incompressible solvers. This is why
> the dimensions are `[0 2 -2 ...]` rather than `[1 -1 -2 ...]` (Pascals).

---

## 3. BC Classification — Dirichlet, Neumann, and Mixed

All boundary conditions in any PDE solver fall into three mathematical categories:

```
  ┌─────────────────────────────────────────────────────────────────────┐
  │            THE THREE TYPES OF BOUNDARY CONDITIONS                  │
  │                                                                     │
  │  1. DIRICHLET (fixedValue)         φ = specified value              │
  │     ─────────────────────                                           │
  │     value ●━━━━━━━━━━━━━━━━         "I know the VALUE at the       │
  │           ┃                          boundary"                      │
  │           ┃  φ(x)                                                   │
  │           ┗━━━━━━━━━━► x                                            │
  │                                                                     │
  │  2. NEUMANN (zeroGradient/fixedGradient)   ∂φ/∂n = specified       │
  │     ──────────────────────────────────                              │
  │              ━━━━━━━━━━● slope = 0   "I know the GRADIENT at the   │
  │           ┃ ╱                         boundary"                     │
  │           ┃╱  φ(x)                                                  │
  │           ┗━━━━━━━━━━► x                                            │
  │                                                                     │
  │  3. MIXED / ROBIN                   a·φ + b·∂φ/∂n = c              │
  │     ────────────                                                    │
  │           ━━━━━━━━━━●~~ blend       "A combination of value and    │
  │           ┃ ╱                         gradient"                     │
  │           ┃╱  φ(x)                                                  │
  │           ┗━━━━━━━━━━► x                                            │
  └─────────────────────────────────────────────────────────────────────┘
```

| Mathematical Type | OpenFOAM BC            | What It Sets        | Example Use               |
|-------------------|------------------------|---------------------|---------------------------|
| Dirichlet         | `fixedValue`           | Value at boundary   | Inlet velocity            |
| Neumann (zero)    | `zeroGradient`         | Zero gradient       | Outlet velocity           |
| Neumann (nonzero) | `fixedGradient`        | Specified gradient  | Heat flux at wall         |
| Mixed             | `mixed`                | Value + gradient    | inletOutlet               |
| Special           | `noSlip`, `slip`, etc. | Physical constraint | Wall conditions           |

> **⚠️ Warning**: You must specify **one BC per patch per field**. Specifying both value and
> gradient (over-determined) or neither (under-determined) leads to solver failure.

---

## 4. Complete BC Reference

### 4.1 fixedValue (Dirichlet)

Sets the field to a **specific constant value** at the boundary.

```
  Boundary
     │
     ●━━━━━ φ = specified value (e.g., U = 1 m/s)
     │
     │  ╱
     │ ╱  interior values develop freely
     │╱
     ┗━━━━━━━━━━━━━━━━► distance from boundary
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Dirichlet                                       |
| **Sets**        | Value of the field at the boundary face          |
| **Used for**    | Inlet velocity, outlet pressure, fixed temp      |
| **Variables**   | U, p, T, k, epsilon, omega                      |

```cpp
inlet
{
    type            fixedValue;
    value           uniform (1 0 0);    // vector for U
}
// or for scalar:
outlet
{
    type            fixedValue;
    value           uniform 0;          // scalar for p
}
```

### 4.2 zeroGradient (Neumann, zero flux)

Sets the **normal gradient to zero** at the boundary — the field value at the boundary equals the
value in the adjacent cell. The field "floats" to whatever value the interior solution produces.

```
  Boundary
     │
     │━━━━━━━━━━● value "floats" (∂φ/∂n = 0)
     │         ╱
     │        ╱  slope approaches zero at boundary
     │       ╱
     │      ╱
     ┗━━━━━━━━━━━━━━━━► distance from boundary
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Neumann (zero)                                  |
| **Sets**        | $\partial\phi/\partial n = 0$ at the boundary    |
| **Used for**    | Outlet velocity, inlet pressure, wall pressure   |
| **Variables**   | U, p, T, k, epsilon                             |

```cpp
outlet
{
    type            zeroGradient;
}
```

> **💡 Tip**: `zeroGradient` for pressure at walls is almost always correct for incompressible
> flows. The pressure adjusts itself; you don't need to specify it.

### 4.3 fixedGradient (Neumann, specified flux)

Sets a **specific non-zero gradient** at the boundary. Useful when you know the flux (e.g., a
known heat flux through a wall) but not the value itself.

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Neumann (non-zero)                              |
| **Sets**        | $\partial\phi/\partial n$ = specified gradient    |
| **Used for**    | Heat flux walls, specified mass flux              |
| **Variables**   | T, scalar fields                                 |

```cpp
heatedWall
{
    type            fixedGradient;
    gradient        uniform 500;    // e.g., heat flux
}
```

### 4.4 noSlip (Wall velocity)

The **no-slip condition** is the fundamental wall BC for viscous flows. It enforces **zero
velocity** at a solid wall — fluid "sticks" to the surface.

```
     Wall (noSlip)
  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░
  U=0 ─┤
       ─┤──►
       ─┤────────►
       ─┤──────────────►
       ─┤────────────────────►   Free-stream velocity
       ─┤──────────────────────────►
         │
         └──────────────────────────► y (distance from wall)

  Velocity profile develops from zero at the wall
  to the free-stream value far from the wall
  (this is the BOUNDARY LAYER)
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Dirichlet (U = 0)                               |
| **Sets**        | All velocity components to zero                  |
| **Used for**    | All solid walls in viscous flow                  |
| **Variables**   | U only                                           |

```cpp
walls
{
    type            noSlip;
}
// Equivalent to:
//   type  fixedValue;
//   value uniform (0 0 0);
```

### 4.5 slip (Frictionless wall)

Allows **tangential flow** but prevents **normal flow** through the boundary. Models a
frictionless wall or symmetry plane for velocity.

```
     Slip wall
  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░
       ──────────────────────►   velocity tangent to wall is FREE
       ──────────────────────►   (no friction)
       ──────────────────────►
     ↕ U_normal = 0              but NO flow through the wall
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Mixed (tangential free, normal zero)             |
| **Sets**        | Normal velocity = 0, tangential = free           |
| **Used for**    | Far-field walls, inviscid walls, free-slip       |
| **Variables**   | U only                                           |

```cpp
farWall
{
    type            slip;
}
```

### 4.6 empty (2D Simulations)

Tells OpenFOAM that a patch has **no physical extent** in one direction. Required for 2D
simulations, which in OpenFOAM are actually thin 3D slabs with one cell in the z-direction.

```
  ┌─────────────────────────────────────────────────────────────┐
  │  3D MESH (one cell thick)         2D SIMULATION             │
  │                                                             │
  │       ┌──────────┐                ┌──────────┐              │
  │      ╱          ╱│               │          │              │
  │     ╱  front   ╱ │    empty      │  solved  │              │
  │    ╱  (empty) ╱  │   ══════►     │   2D     │              │
  │   ┌──────────┐   │               │  domain  │              │
  │   │          │   │               │          │              │
  │   │          │  ╱                └──────────┘              │
  │   │   back   │ ╱   z-direction                             │
  │   │  (empty) │╱    has NO effect                           │
  │   └──────────┘     on solution                             │
  └─────────────────────────────────────────────────────────────┘
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Special (geometric)                              |
| **Sets**        | No solution in that direction                    |
| **Used for**    | Front/back patches in 2D simulations             |
| **Variables**   | ALL fields (U, p, k, epsilon, nut, etc.)         |

```cpp
frontAndBack
{
    type            empty;
}
```

> **⚠️ Warning**: The mesh **must** be exactly one cell thick in the empty direction, and the
> patch type in `constant/polyMesh/boundary` must also be set to `empty`.

### 4.7 symmetryPlane / symmetry

Used to model a **plane of symmetry**, halving the computational domain. The solution is mirrored
across this plane, saving half the computational cost.

```
  ┌────────────────────────────────────────────────────────────────┐
  │  FULL DOMAIN                    HALF DOMAIN + SYMMETRY         │
  │                                                                │
  │  ┌──────────────┐              ┌──────────────┐               │
  │  │   ●          │              │   ●          │               │
  │  │  ╱ ╲         │              │  ╱ ╲         │               │
  │  │ ╱   ╲  flow  │   ═════►    │ ╱   ╲  flow  │               │
  │  ├──────────────┤  cut here   ╞══════════════╡ symmetryPlane  │
  │  │ ╲   ╱  flow  │              (mirror image                  │
  │  │  ╲ ╱         │               is implied)                   │
  │  │   ●          │                                              │
  │  └──────────────┘              Saves 50% of cells!            │
  └────────────────────────────────────────────────────────────────┘
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Special (geometric)                              |
| **Sets**        | Zero normal gradient, zero normal velocity       |
| **Used for**    | Symmetric geometries (half-planes, quarter)      |
| **Variables**   | ALL fields                                       |

```cpp
midPlane
{
    type            symmetryPlane;
}
// or the more general version:
midPlane
{
    type            symmetry;
}
```

### 4.8 wedge (Axisymmetric Simulations)

For **axisymmetric problems** (pipes, nozzles, rotating bodies), OpenFOAM uses a thin wedge-shaped
slice instead of a full 3D domain.

```
  ┌──────────────────────────────────────────────────────────────────┐
  │  FULL 3D (expensive)              WEDGE SLICE (cheap)            │
  │                                                                  │
  │       ╭─────────╮                    ╱│                          │
  │      ╱    ●     ╲                   ╱ │ ← wedge patch (back)     │
  │     │   axis     │    ═════►       ╱  │                          │
  │      ╲           ╱               ●────┤ axis                     │
  │       ╰─────────╯                 ╲  │                          │
  │    (millions of cells)             ╲ │ ← wedge patch (front)    │
  │                                      ╲│                          │
  │                                   ~5° wedge angle               │
  │                                   (typically 2.5° each side)     │
  └──────────────────────────────────────────────────────────────────┘
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Special (geometric)                              |
| **Sets**        | Axisymmetric constraint on thin wedge            |
| **Used for**    | Pipes, jets, nozzles, axisymmetric bodies        |
| **Variables**   | ALL fields                                       |

```cpp
front
{
    type            wedge;
}
back
{
    type            wedge;
}
```

### 4.9 Inlet/Outlet Types

These **mixed** BCs switch behavior depending on flow direction — crucial for cases where
backflow might occur at boundaries.

**inletOutlet**: Acts as `fixedValue` when flow enters and `zeroGradient` when flow exits.

```cpp
outlet
{
    type            inletOutlet;
    inletValue      uniform (0 0 0);    // value if flow reverses
    value           uniform (0 0 0);    // initial guess
}
```

**pressureInletOutletVelocity**: For pressure boundaries where flow direction may reverse.

```cpp
outlet
{
    type            pressureInletOutletVelocity;
    value           uniform (0 0 0);
}
```

| BC Type                           | Outflow Behavior | Inflow Behavior    | Used For       |
|-----------------------------------|------------------|--------------------|----------------|
| `inletOutlet`                     | zeroGradient     | fixedValue         | U at outlets   |
| `pressureInletOutletVelocity`     | zeroGradient     | from pressure      | U at openings  |
| `totalPressure`                   | adjusts p        | adjusts p          | p at inlets    |

### 4.10 freestream (Far-field)

Used for **external aerodynamics** where boundaries are far from the body. Applies a freestream
value when flow enters the domain and zeroGradient when flow leaves.

```cpp
farfield
{
    type            freestream;
    freestreamValue uniform (25.75 0 0);    // freestream velocity
}
```

### 4.11 Turbulent Wall Functions

When using RANS turbulence models, the near-wall region requires special treatment. Wall functions
**bridge** the viscous sublayer, allowing coarser meshes near walls.

```
  ┌──────────────────────────────────────────────────────────────────┐
  │  NEAR-WALL VELOCITY PROFILE AND WALL FUNCTIONS                  │
  │                                                                  │
  │  u+                                                              │
  │  │                                    ╱ log law region           │
  │  │                              ╱╱╱╱╱                            │
  │  │                         ╱╱╱╱╱                                 │
  │  │                    ╱╱╱╱   ← wall function bridges this gap    │
  │  │              ╱╱╱╱╱                                            │
  │  │         ╱╱╱╱                                                  │
  │  │    ╱╱╱╱                                                       │
  │  │  ╱╱  buffer layer                                             │
  │  │╱╱                                                             │
  │  ╱ viscous sublayer (u+ = y+)                                    │
  │  ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ y+         │
  │  0    5        30         100        1000                        │
  └──────────────────────────────────────────────────────────────────┘
```

| Wall Function            | Field   | Description                                     |
|--------------------------|---------|-------------------------------------------------|
| `kqRWallFunction`        | k, q, R | Zero-gradient for turbulent kinetic energy       |
| `epsilonWallFunction`    | epsilon | Calculates epsilon in near-wall cell             |
| `omegaWallFunction`      | omega   | Calculates omega in near-wall cell               |
| `nutkWallFunction`       | nut     | Calculates nut from k and log-law                |
| `nutUWallFunction`       | nut     | Calculates nut from velocity (no k needed)       |

```cpp
// From projects/04_naca_airfoil_analysis/0/k:
airfoil
{
    type            kqRWallFunction;
    value           uniform 0.01;
}

// From projects/04_naca_airfoil_analysis/0/epsilon:
airfoil
{
    type            epsilonWallFunction;
    value           uniform 0.01;
}

// From projects/04_naca_airfoil_analysis/0/nut:
airfoil
{
    type            nutkWallFunction;
    value           uniform 0;
}
```

> **💡 Tip**: The `value` entry in wall functions is an initial guess, not a fixed value. The
> wall function overrides it during the simulation.

### 4.12 calculated (Derived fields)

The `calculated` type means the field is **computed from other fields**, not directly solved. The
solver calculates the value rather than applying a BC.

```cpp
// From projects/04_naca_airfoil_analysis/0/nut:
inlet
{
    type            calculated;
    value           uniform 0;
}
```

| Property        | Detail                                         |
|-----------------|-------------------------------------------------|
| **Type**        | Derived                                          |
| **Sets**        | Nothing — value comes from other fields          |
| **Used for**    | nut at non-wall patches, derived quantities      |
| **Variables**   | nut, alphat, mut                                 |

### 4.13 cyclic (Periodic Boundaries)

Models **periodically repeating** geometries. The solution on one patch is mapped directly to its
partner patch.

```
  ┌────────────────────────────────────────────────────────────────┐
  │  PERIODIC / CYCLIC BOUNDARIES                                  │
  │                                                                │
  │   cyclic_left              cyclic_right                        │
  │      │  ┌──────────────────────┐  │                            │
  │      │  │ → → → → → → → → → → │  │                            │
  │      │  │ → → → → → → → → → → │  │                            │
  │      │  │ → → → → → → → → → → │  │                            │
  │      │  └──────────────────────┘  │                            │
  │      │◄─────── matched ──────────►│                            │
  │                                                                │
  │   What exits the right enters the left (and vice versa)        │
  │   Used for: fully-developed channel flow, repeating geometry   │
  └────────────────────────────────────────────────────────────────┘
```

```cpp
cyclic_left
{
    type            cyclic;
    neighbourPatch  cyclic_right;
}
cyclic_right
{
    type            cyclic;
    neighbourPatch  cyclic_left;
}
```

---

## 5. Boundary Conditions by Variable — Quick Reference

This table shows the **typical** BC choice for each variable at each boundary type:

| Variable   | Inlet              | Outlet            | Wall (laminar)   | Wall (turbulent)         | Symmetry        |
|------------|--------------------|--------------------|------------------|--------------------------|-----------------|
| **U**      | `fixedValue`       | `zeroGradient`     | `noSlip`         | `noSlip`                 | `symmetry`      |
| **p**      | `zeroGradient`     | `fixedValue (0)`   | `zeroGradient`   | `zeroGradient`           | `symmetry`      |
| **k**      | `fixedValue`       | `zeroGradient`     | —                | `kqRWallFunction`        | `symmetry`      |
| **epsilon**| `fixedValue`       | `zeroGradient`     | —                | `epsilonWallFunction`    | `symmetry`      |
| **omega**  | `fixedValue`       | `zeroGradient`     | —                | `omegaWallFunction`      | `symmetry`      |
| **nut**    | `calculated (0)`   | `calculated (0)`   | —                | `nutkWallFunction`       | `symmetry`      |
| **T**      | `fixedValue`       | `zeroGradient`     | `fixedValue` or `fixedGradient` | `fixedValue` or `fixedGradient` | `symmetry` |

> **💡 Tip**: For incompressible solvers, pressure is typically `zeroGradient` at the inlet and
> `fixedValue uniform 0` at the outlet. For compressible solvers, this may differ.

---

## 6. Laminar vs Turbulent Boundary Conditions

The biggest difference in BC setup is between **laminar** and **turbulent** simulations. Turbulent
cases require additional field variables (k, epsilon or omega, nut) with their own BCs.

### 6.1 Laminar Case — Lid-Driven Cavity

From `projects/01_lid_driven_cavity/0/U`:

```cpp
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
}
```

From `projects/01_lid_driven_cavity/0/p`:

```cpp
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
}
```

**Only two files needed**: `U` and `p`. No turbulence variables.

### 6.2 Turbulent Case — NACA Airfoil

The turbulent case requires **five** field files. Here are the additional turbulence files:

From `projects/04_naca_airfoil_analysis/0/k` (turbulent kinetic energy):

```cpp
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

From `projects/04_naca_airfoil_analysis/0/epsilon` (dissipation rate):

```cpp
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

From `projects/04_naca_airfoil_analysis/0/nut` (turbulent viscosity):

```cpp
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

### 6.3 Side-by-Side Comparison

```
  ┌──────────────────────────────┬────────────────────────────────────┐
  │  LAMINAR (cavity)            │  TURBULENT (airfoil, k-epsilon)    │
  │                              │                                    │
  │  Files needed:               │  Files needed:                     │
  │    0/U  ✓                    │    0/U       ✓                     │
  │    0/p  ✓                    │    0/p       ✓                     │
  │                              │    0/k       ✓  ← NEW             │
  │                              │    0/epsilon ✓  ← NEW             │
  │                              │    0/nut     ✓  ← NEW             │
  │                              │                                    │
  │  Wall BCs:                   │  Wall BCs:                         │
  │    U: noSlip                 │    U:       noSlip                 │
  │    p: zeroGradient           │    p:       zeroGradient           │
  │                              │    k:       kqRWallFunction        │
  │                              │    epsilon: epsilonWallFunction    │
  │                              │    nut:     nutkWallFunction       │
  └──────────────────────────────┴────────────────────────────────────┘
```

---

## 7. Calculating Turbulent BC Values

Getting the right inlet values for k, epsilon, and omega is critical. Here are the standard
formulas based on **turbulence intensity** (I) and **turbulent length scale** (l).

### Formulas

$$k = 1.5 \times (U \times I)^{2}$$

$$\varepsilon = \frac{C_{\mu}^{0.75} \times k^{1.5}}{l} \quad \text{where } C_{\mu} = 0.09$$

$$\omega = \frac{k^{0.5}}{C_{\mu}^{0.25} \times l}$$

$$\nu_t = \frac{C_{\mu} \times k^{2}}{\varepsilon} \quad \text{(or computed automatically)}$$

Where:
- $U$ = mean flow velocity (m/s)
- $I$ = turbulence intensity (fraction, e.g., 0.05 for 5%)
- $l$ = turbulent length scale (m), often estimated as $0.07 \times D$ ($D$ = characteristic length)

### Typical Values

| Flow Type                  | Turbulence Intensity ($I$) | Length Scale ($l$)     |
|----------------------------|--------------------------|------------------------|
| Low-turbulence wind tunnel | 0.5% – 1%               | $0.01 \times D$        |
| Medium turbulence          | 1% – 5%                 | $0.07 \times D$        |
| High turbulence (urban)    | 5% – 20%                | $0.1 \times D$         |
| Pipe flow (fully dev.)     | ~5%                      | $0.07 \times D_{pipe}$ |
| Free-stream (external)     | 0.1% – 1%               | $0.01 \times chord$    |

### Example Calculation

For the NACA airfoil case with $U = 1$ m/s, $I = 10\%$, $l = 0.01$ m:


$$k = 1.5 \times (1.0 \times 0.1)^{2} = 1.5 \times 0.01 = 0.015 \; m^{2}/s^{2}$$

$$\varepsilon = \frac{0.09^{0.75} \times 0.015^{1.5}}{0.01} \approx 0.028 \; m^{2}/s^{3}$$

> **💡 Tip**: If you are unsure about turbulence intensity, 5% is a reasonable default for
> most industrial applications. The solution is often not very sensitive to inlet turbulence
> values if the domain is long enough for turbulence to develop naturally.

---

## 8. Real Case Studies from This Repository

### 8.1 Lid-Driven Cavity — The Moving Wall Drives Everything

In the lid-driven cavity, the **movingWall** BC is the sole driver of the entire flow. There is
no inlet or outlet — the top wall moves and drags fluid via viscous forces.

```
  ┌─────────────────────────────────────────────────────┐
  │                                                     │
  │  movingWall: fixedValue (1 0 0)                     │
  │  ══════════════════════════════►  U = 1 m/s         │
  │  ┌─────────────────────────────┐                    │
  │  │  ╭──────────────────────╮   │                    │
  │  │  │     primary vortex   │   │  fixedWalls:       │
  │  │  │    ╭──╮              │   │  noSlip            │
  │  │  │    │  │  secondary   │   │  (U = 0)           │
  │  │  │    ╰──╯   vortex    │   │                    │
  │  │  ╰──────────────────────╯   │                    │
  │  └─────────────────────────────┘                    │
  │  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░                     │
  │  fixedWalls: noSlip (U = 0)                         │
  │                                                     │
  │  frontAndBack: empty (2D simulation)                │
  └─────────────────────────────────────────────────────┘
```

### 8.2 Elbow — Two Inlets with Different Velocity Directions

The elbow case demonstrates how to handle **multiple inlets** with **different flow directions**.

From `projects/02_elbow/0/u`:

```cpp
boundaryField
{
    wall-4
    {
        type            noSlip;
    }
    velocity-inlet-5
    {
        type            fixedValue;
        value           uniform (1 0 0);    // horizontal flow, 1 m/s
    }
    velocity-inlet-6
    {
        type            fixedValue;
        value           uniform (0 3 0);    // vertical flow, 3 m/s
    }
    pressure-outlet-7
    {
        type            zeroGradient;
    }
    wall-8
    {
        type            noSlip;
    }
    frontAndBackPlanes
    {
        type            empty;
    }
}
```

```
  ┌──────────────────────────────────────────────────────┐
  │                                                      │
  │  velocity-inlet-6                                    │
  │  fixedValue (0 3 0)                                  │
  │         │ │ │                                        │
  │         ▼ ▼ ▼   3 m/s downward                      │
  │    ┌────────────┐                                    │
  │    │            │                                    │
  │    │   ELBOW    ├──────► pressure-outlet-7           │
  │    │   mixing   │        zeroGradient                │
  │    │   region   │                                    │
  │    └────────────┘                                    │
  │         ▲                                            │
  │    ═════╪═══════► 1 m/s rightward                    │
  │  velocity-inlet-5                                    │
  │  fixedValue (1 0 0)                                  │
  │                                                      │
  │  Two streams meet and mix in the elbow               │
  └──────────────────────────────────────────────────────┘
```

### 8.3 NACA Airfoil — Wall Functions on the Airfoil Surface

The airfoil case shows turbulent wall functions applied to a curved aerodynamic surface.

From `projects/04_naca_airfoil_analysis/0/U`:

```cpp
boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (1 0 0);
    }
    outlet
    {
        type            zeroGradient;
    }
    walls
    {
        type            noSlip;
    }
    airfoil
    {
        type            noSlip;
    }
    frontAndBack
    {
        type            empty;
    }
}
```

```
  ┌─────────────────────────────────────────────────────────────┐
  │  NACA AIRFOIL — BOUNDARY CONDITION MAP                      │
  │                                                             │
  │  inlet               walls (noSlip)              outlet     │
  │  fixedValue    ┌──────────────────────────┐    zeroGradient │
  │  (1 0 0)  ►   │         airfoil           │   ►            │
  │  ══════════►   │    ╱‾‾‾‾‾‾‾‾‾‾‾‾╲       │   ══════════►  │
  │  ══════════►   │  ╱                ╲      │   ══════════►  │
  │  ══════════►   │ ●──────────────────●     │   ══════════►  │
  │  ══════════►   │  ╲________________╱      │   ══════════►  │
  │  ══════════►   │                          │   ══════════►  │
  │  ══════════►   │    noSlip + wall funcs   │   ══════════►  │
  │                └──────────────────────────┘                 │
  │                                                             │
  │  airfoil patch gets:                                        │
  │    U:       noSlip                                          │
  │    p:       zeroGradient                                    │
  │    k:       kqRWallFunction                                 │
  │    epsilon: epsilonWallFunction                              │
  │    nut:     nutkWallFunction                                │
  │                                                             │
  │  outlet patch (p):  fixedValue uniform 0  (reference p)     │
  │  frontAndBack:      empty (2D simulation)                   │
  └─────────────────────────────────────────────────────────────┘
```

---

## 9. Common BC Mistakes and Troubleshooting

### 9.1 Ill-Posed BC Combinations

| Mistake                                    | Symptom                              | Fix                          |
|--------------------------------------------|--------------------------------------|------------------------------|
| `fixedValue` on U at both inlet AND outlet | Solver diverges or wrong mass flow   | Use `zeroGradient` at outlet |
| `zeroGradient` on p at ALL boundaries      | Floating pressure (no reference)     | Fix p at one boundary        |
| `fixedValue` on p at both inlet and outlet | Over-constrained pressure            | Use `zeroGradient` at inlet  |
| Missing `empty` on front/back (2D)         | Solver crashes immediately           | Add `empty` to both patches  |

### 9.2 Missing Turbulence BCs

```
  ┌──────────────────────────────────────────────────────────────────┐
  │  ERROR: Missing turbulence BCs                                   │
  │                                                                  │
  │  If you switch from laminar to turbulent (e.g., add kEpsilon)    │
  │  you MUST create files for ALL turbulence variables:             │
  │                                                                  │
  │    0/k        ← turbulent kinetic energy                         │
  │    0/epsilon  ← dissipation rate (k-epsilon model)               │
  │    0/omega    ← specific dissipation (k-omega model, INSTEAD)    │
  │    0/nut      ← turbulent viscosity                              │
  │                                                                  │
  │  Missing any of these → immediate crash with error like:         │
  │  "Cannot find patchField entry for patch <name>"                 │
  └──────────────────────────────────────────────────────────────────┘
```

### 9.3 Empty BC Misuse

> **⚠️ Warning**: The `empty` BC must match the mesh topology. Common mistakes:
> - Using `empty` on a patch that is not exactly one cell thick → **crash**
> - Using `empty` for one field but not others on the same patch → **crash**
> - Forgetting to set `type empty;` in `constant/polyMesh/boundary` → **crash**
> - Using `empty` on an internal patch → **crash**

### 9.4 Dimensional Mismatch

If you see errors like `"dimensions do not match"`, check that the `value` entry in your BC
matches the field dimensions:

```
// WRONG — pressure is scalar, not vector:
inlet
{
    type    fixedValue;
    value   uniform (1 0 0);    // ← vector in a scalar field!
}

// CORRECT:
inlet
{
    type    fixedValue;
    value   uniform 0;          // ← scalar for pressure
}
```

---

## 10. Choosing the Right BC — Decision Flowchart

```
  START: What kind of boundary is this?
    │
    ├── INLET (fluid enters)
    │     │
    │     ├── Know the velocity?
    │     │     YES → U: fixedValue, p: zeroGradient
    │     │     NO  → Know the pressure?
    │     │             YES → p: fixedValue/totalPressure, U: pressureInletVelocity
    │     │
    │     └── Turbulent?
    │           YES → Also set k, epsilon/omega: fixedValue
    │                  nut: calculated
    │           NO  → Done (only U and p needed)
    │
    ├── OUTLET (fluid leaves)
    │     │
    │     ├── Backflow possible?
    │     │     YES → U: inletOutlet, p: fixedValue 0
    │     │     NO  → U: zeroGradient, p: fixedValue 0
    │     │
    │     └── Turbulent?
    │           YES → k, epsilon/omega: zeroGradient, nut: calculated
    │           NO  → Done
    │
    ├── WALL (solid surface)
    │     │
    │     ├── Viscous (most cases)?
    │     │     YES → U: noSlip
    │     │     NO  → U: slip (inviscid/far-field)
    │     │
    │     ├── p: zeroGradient (almost always)
    │     │
    │     └── Turbulent?
    │           YES → k: kqRWallFunction
    │                  epsilon: epsilonWallFunction  (or omega: omegaWallFunction)
    │                  nut: nutkWallFunction
    │           NO  → Done
    │
    ├── SYMMETRY PLANE
    │     └── ALL fields: symmetryPlane (or symmetry)
    │
    ├── 2D FRONT/BACK
    │     └── ALL fields: empty
    │
    ├── AXISYMMETRIC WEDGE
    │     └── ALL fields: wedge
    │
    └── PERIODIC (repeating geometry)
          └── ALL fields: cyclic (with neighbourPatch)
```

> **💡 Final Tip**: When in doubt, start with the simplest valid BC combination (fixedValue inlet,
> zeroGradient outlet for U; zeroGradient inlet, fixedValue 0 outlet for p). Get the simulation
> running first, then refine the BCs for physical accuracy.
