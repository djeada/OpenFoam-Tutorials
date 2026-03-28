# The CFL Number: Stability, Time Stepping, and OpenFOAM Practice

## Introduction — What the CFL Number Really Means

The Courant-Friedrichs-Lewy (CFL) number answers one critical question in every
CFD simulation:

> **"How far does information travel through the grid relative to the cell size
> during a single time step?"**

Named after Richard Courant, Kurt Friedrichs, and Hans Lewy, who published the
condition in 1928, the CFL number is a dimensionless ratio that governs whether
a numerical simulation will produce meaningful results or blow up spectacularly.

Think of it this way: your simulation marches forward in discrete time steps. At
each step, the solver looks at neighboring cells to compute what happens next. If
the flow moves information *faster* than the solver can track — skipping over
cells entirely — the solver is blind to what actually happened. The result?
Instability, divergence, and garbage output.

---

## Physical Intuition — Information Travel on a Grid

The CFL number connects three quantities: the flow velocity, the time step, and
the grid spacing. Here is what that looks like on a 1D grid:

```
                        u = velocity
                        ─────────────→

    Time step Δt
    ┌─────────┐
    │         │       During one Δt, flow travels a distance u·Δt
    │  u · Δt │
    │         │       The grid has cells of width Δx
    └─────────┘
                      CFL = u·Δt / Δx

    ←──Δx──→←──Δx──→←──Δx──→←──Δx──→←──Δx──→
    ┌───────┬───────┬───────┬───────┬───────┐
    │       │       │       │       │       │
    │Cell 1 │Cell 2 │Cell 3 │Cell 4 │Cell 5 │
    │       │       │       │       │       │
    └───────┴───────┴───────┴───────┴───────┘

    CFL < 1 :  Information travels LESS than one cell       ✓ SAFE
               per time step. The solver "sees" everything.

    CFL = 1 :  Information travels EXACTLY one cell         ⚠ Boundary
               per time step. Just barely OK.

    CFL > 1 :  Information SKIPS cells — the solver         ✗ UNSTABLE
               misses what happened in between!
```

> **Intuition**: Imagine you are filming a bouncing ball with a camera. If the
> frame rate is high enough (small Δt), you see every bounce. If the frame rate
> is too low (large Δt), the ball teleports between frames and you cannot
> reconstruct its path. CFL < 1 means your "frame rate" is fast enough to
> capture the physics.

---

## The CFL Condition — Derived

### One-Dimensional Formula

For a 1D problem with uniform grid spacing Δx, velocity u, and time step Δt:

```
              u · Δt
    CFL  =  ─────────
                Δx

    Stability requires:  CFL  ≤  CFL_max

    For most explicit schemes:  CFL_max = 1
```

Where:
- **u** — the local flow velocity (m/s)
- **Δt** — the simulation time step (s)
- **Δx** — the grid cell size (m)

### Two-Dimensional Extension

In 2D, information can travel in both the x and y directions simultaneously.
The CFL condition becomes:

```
              u · Δt     v · Δt
    CFL  =  ───────── + ─────────  ≤  CFL_max
                Δx         Δy
```

Where u and v are the velocity components and Δx, Δy are the cell dimensions
in each direction.

### General Three-Dimensional Form

For 3D flows (what most OpenFOAM cases use):

```
              u · Δt     v · Δt     w · Δt
    CFL  =  ───────── + ───────── + ─────────  ≤  CFL_max
                Δx         Δy         Δz
```

In OpenFOAM's finite volume framework, this is computed per-cell using the
face fluxes and cell volumes, which naturally handles non-uniform and
unstructured meshes:

```
                  Σ_f |φ_f|
    Co_cell  =  ────────────  · Δt
                  2 · V_cell

    Where:
      φ_f    = face flux (volumetric flow rate through face f)
      V_cell = cell volume
      Σ_f    = sum over all faces of the cell
```

### Connection to Von Neumann Stability Analysis

The CFL condition can be derived rigorously via von Neumann stability analysis.
The idea: decompose the numerical error into Fourier modes and check whether
each mode grows or decays over time. For the first-order upwind scheme applied
to the linear advection equation:

```
    Amplification factor:  |g| = |1 - CFL·(1 - e^(-ikΔx))|

    For stability:  |g| ≤ 1  for all wavenumbers k
    This yields:    0 ≤ CFL ≤ 1
```

The key takeaway: the CFL limit is not a rule of thumb — it is a **mathematical
requirement** for stability of explicit time integration.

---

## Visual: What Happens When CFL > 1

```
    CFL ≤ 1  (STABLE)                     CFL > 1  (UNSTABLE)

    The domain of dependence of           The numerical domain of dependence
    the PDE is INSIDE the numerical       is SMALLER than the physical one.
    stencil.                              Information is LOST.

    t+Δt  ┌─────┬─────┬─────┬─────┐      t+Δt  ┌─────┬─────┬─────┬─────┐
          │  ?  │  ?  │  ?  │  ?  │            │  ?  │  ?  │  ?  │  ?  │
          └──┬──┴──┬──┴──┬──┴─────┘            └──┬──┴─────┴──┬──┴─────┘
             │╲    │    ╱│                        │╲          ╱│
             │ ╲   │   ╱ │                        │ ╲  GAP!  ╱ │
             │  ╲  │  ╱  │                        │  ╲      ╱  │
    t     ┌──┴──┬──┴──┬──┴──┬─────┐      t     ┌──┴──┬─────┬──┴──┬─────┐
          │ old │ old │ old │ old │            │ old │ old │ old │ old │
          └─────┴─────┴─────┴─────┘            └─────┴─────┴─────┴─────┘

    Each new value depends on             Solver only looks at immediate
    neighbors that contain all            neighbors, but information has
    the needed information.  ✓            traveled further → MISSED!  ✗


    Result: Smooth, physical solution     Result: Oscillations → divergence

    ──────────────────────────────────    ──────────────────────────────────
    Value                                 Value
    │     ╱──────                         │    ╱╲  ╱╲
    │    ╱                                │   ╱  ╲╱  ╲  ╱╲
    │   ╱                                 │  ╱        ╲╱  ╲ ╱╲
    │──╱                                  │─╱              ╲╱
    └──────────────→ x                    └──────────────────→ x
    Stable propagation                    Oscillating, divergent
```

> **Warning**: Even CFL slightly above 1 can cause instability in explicit
> schemes. The simulation may run for a while before oscillations grow large
> enough to crash it. Always monitor your Courant number!

---

## Explicit vs Implicit Time Integration Schemes

This is one of the most important practical distinctions in CFD. The CFL
condition applies *strictly* only to explicit schemes, but understanding why
helps you choose the right approach.

### How the Schemes Differ

```
    EXPLICIT SCHEME                        IMPLICIT SCHEME
    (Forward Euler, RK4, etc.)             (Backward Euler, Crank-Nicolson)

    Known values at time t                 Known values at time t
    are used DIRECTLY to                   AND unknown values at t+Δt
    compute values at t+Δt.               are used TOGETHER.

    t+Δt:   [  ?  ]                        t+Δt:   [  ?  ]──[  ?  ]──[  ?  ]
              ↑                                      ↑    ╲   ↑   ╱    ↑
    t:     [old]──[old]──[old]             t:     [old]──[old]──[old]

    Stencil at t only                      Stencil spans t AND t+Δt
    → Limited information reach            → Implicit coupling captures
    → CFL ≤ 1 required                      distant information
                                           → No strict CFL limit

    Cheap per step (no system solve)       Expensive per step (matrix solve)
    Many small steps needed                Fewer large steps possible
```

### Comparison Table

| Property              | Explicit Schemes          | Implicit Schemes             |
|-----------------------|---------------------------|------------------------------|
| **CFL limit**         | Strict: CFL ≤ 1           | No strict stability limit    |
| **Cost per time step**| Low (direct evaluation)   | High (linear system solve)   |
| **Time accuracy**     | Depends on scheme order   | Depends on scheme order      |
| **Memory usage**      | Low                       | Higher (stores matrix)       |
| **Best for**          | Waves, acoustics, LES     | Steady/slow transients       |
| **OpenFOAM solvers**  | icoFoam, interFoam        | simpleFoam, pimpleFoam       |
| **Parallelization**   | Excellent                 | Requires parallel solvers    |
| **Robustness**        | Fragile if CFL violated   | More forgiving               |

### Why Implicit Schemes Tolerate Higher CFL

Implicit schemes include unknown future values in the discretization. This means
the solution at t+Δt depends on *all* cells simultaneously (through the matrix
system), not just immediate neighbors. The domain of dependence of the numerical
scheme is effectively the entire domain — it always contains the physical domain
of dependence, regardless of CFL.

> **Trade-off**: Implicit schemes let you take bigger time steps, but each step
> costs more (you must assemble and solve a large sparse linear system). The
> total wall-clock time may or may not be less — it depends on the problem!

### The PIMPLE Algorithm — Best of Both Worlds?

OpenFOAM's `pimpleFoam` uses the PIMPLE algorithm (PISO + SIMPLE), which
performs outer correction loops within each time step. This allows it to run
at CFL > 1 while maintaining reasonable accuracy:

```
    Standard PISO (CFL ≤ 1):              PIMPLE (CFL can exceed 1):

    ┌─ Time step ──────────────┐           ┌─ Time step ──────────────────┐
    │  1× momentum predictor   │           │  Outer loop (2-5 iterations) │
    │  2× pressure corrector   │           │  ┌─ Each iteration ────────┐ │
    │  (no outer iterations)   │           │  │ momentum predictor      │ │
    └──────────────────────────┘           │  │ pressure corrector(s)   │ │
                                           │  │ check convergence       │ │
                                           │  └─────────────────────────┘ │
                                           └──────────────────────────────┘
    Fast per step, needs small Δt          Slower per step, allows large Δt
```

---

## CFL in OpenFOAM — Practical Details

### How OpenFOAM Reports the Courant Number

Every transient OpenFOAM solver computes and prints the Courant number at each
time step. A typical log output looks like:

```
    Time = 0.005

    Courant Number mean: 0.0823416 max: 0.852394
    smoothSolver: Solving for Ux, Initial residual = 0.0367, Final residual = 2.18e-06
    smoothSolver: Solving for Uy, Initial residual = 0.142,  Final residual = 8.93e-06
    ...
```

What these mean:
- **mean**: The volume-weighted average CFL across all cells
- **max**: The highest CFL in any single cell (this is the critical one!)

> **Tip**: Always watch the **max** Courant number. A mean of 0.1 is useless if
> the max is 5.0 — one unstable cell can crash the entire simulation.

### The CourantNo.H Include

Inside OpenFOAM's transient solvers, the Courant number calculation is handled
by an included header file `CourantNo.H`. The logic is essentially:

```cpp
// Simplified from OpenFOAM's CourantNo.H
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

scalarField sumPhi
(
    fvc::surfaceSum(mag(phi))().primitiveField()
);

CoNum = 0.5 * gMax(sumPhi / mesh.V().field()) * runTime.deltaTValue();
meanCoNum = 0.5 * (gSum(sumPhi) / gSum(mesh.V().field())) * runTime.deltaTValue();

Info<< "Courant Number mean: " << meanCoNum
    << " max: " << CoNum << endl;
```

This computes per-cell Courant numbers using face fluxes (φ) and cell volumes
(V), then reports the global maximum and mean.

### Real Code: Lid-Driven Cavity (Transient, Fixed Δt)

From our project `projects/01_lid_driven_cavity/system/controlDict`:

```c
application     icoFoam;

// Start time of the simulation
startFrom       startTime;
startTime       0;

// Stop time of the simulation
stopAt          endTime;
endTime         0.5;

// Time step size
deltaT          0.001;

// Control output based on time step
writeControl    timeStep;

// Write interval for output data
writeInterval   1;

// Purge older time directories
purgeWrite      0;

// Output format
writeFormat     ascii;

// Precision of the output data
writePrecision  6;

// Compression for output data
writeCompression off;

// Format of time in output files
timeFormat      general;
timePrecision   6;

// Allow runtime modifications
runTimeModifiable true;
```

Key observations for CFL:
- **`application icoFoam`** — transient, laminar, incompressible (explicit-like PISO)
- **`deltaT 0.001`** — fixed time step of 1 ms
- **No `adjustTimeStep`** — the user pre-calculated Δt to keep CFL safe
- For a 1 m/s lid velocity and typical cavity mesh (~0.01 m cells):
  CFL ≈ 1.0 × 0.001 / 0.01 = 0.1 ✓

### Real Code: NACA Airfoil (Steady-State, No CFL Concern)

From our project `projects/03_naca_airfoil_analysis/system/controlDict`:

```c
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

Key observations:
- **`application simpleFoam`** — steady-state RANS solver (SIMPLE algorithm)
- **`endTime 1000`** — these are *iteration* counts, not physical time
- **`deltaT 0.1`** — acts as under-relaxation pseudo-time step
- **No CFL concern** — SIMPLE iterates toward steady state; there is no physical
  time marching so the Courant number has no stability meaning here

> **Key Insight**: Steady-state solvers like simpleFoam do not have a CFL
> stability constraint. The `deltaT` in simpleFoam controls iteration stepping,
> not physical time advancement. Stability is governed by under-relaxation
> factors in `fvSolution` instead.

---

## Automatic Time Step Adjustment

For transient cases where flow conditions change significantly (startup,
vortex shedding, sudden inflows), a fixed time step is often either too
conservative (slow) or too aggressive (unstable). OpenFOAM can adjust Δt
automatically.

### Configuration

Add these entries to `controlDict`:

```c
adjustTimeStep  yes;       // Enable automatic time step adjustment
maxCo           0.5;       // Target maximum Courant number
maxDeltaT       0.01;      // Upper bound on time step (seconds)
```

OpenFOAM will compute the Courant number at each step and scale Δt so that
the maximum Courant number stays at or below `maxCo`.

### How It Works — Adaptive Δt Over Time

```
    Δt
    │
    │         ╱╲
    │        ╱  ╲        ╱╲
    │       ╱    ╲      ╱  ╲
    │      ╱      ╲    ╱    ╲       ╱╲
    │     ╱        ╲  ╱      ╲     ╱  ╲
    │    ╱          ╲╱        ╲   ╱    ╲
    │   ╱                      ╲─╱      ╲────────
    │──╱
    └─────────────────────────────────────────────→ Time
         ↑              ↑                  ↑
      Startup:       Vortex shedding:   Flow settles:
      Small Δt       Δt oscillates      Δt stabilizes
      (high CFL      with flow          at maxDeltaT
       risk)          dynamics           or maxCo limit

    At every step:  Δt_new = maxCo × Δt_old / Co_max_current
                    Δt_new = min(Δt_new, maxDeltaT)
```

### When to Use Fixed vs Adjustable Time Step

| Scenario                        | Fixed Δt          | Adjustable Δt      |
|---------------------------------|--------------------|--------------------|
| Flow velocity is well-known     | ✓ Good choice      | Works but overkill |
| Complex startup transients      | Risk of instability| ✓ Recommended      |
| Uniform mesh                    | ✓ Easy to estimate | Either works       |
| Highly non-uniform mesh         | Hard to pick right | ✓ Recommended      |
| Consistent output intervals     | ✓ Simpler          | Need writeControl  |
| Vortex shedding / oscillations  | May waste time     | ✓ Adapts well      |
| Validation / benchmarking       | ✓ Reproducible     | Harder to compare  |

> **Tip**: When using `adjustTimeStep yes`, set `writeControl adjustableRunTime`
> (not `timeStep`) so that output is written at regular physical time intervals
> rather than at every variable-sized step.

---

## CFL Guidelines for Different OpenFOAM Solvers

Not all solvers have the same CFL requirements. The algorithm used determines
how strictly the CFL condition must be enforced:

| Solver        | Algorithm | Type                 | Typical CFL | Notes                                       |
|---------------|-----------|----------------------|-------------|---------------------------------------------|
| icoFoam       | PISO      | Transient, laminar   | ≤ 1         | Strict CFL required for stability           |
| pisoFoam      | PISO      | Transient, turbulent | ≤ 1         | Same as icoFoam, with turbulence models     |
| pimpleFoam    | PIMPLE    | Transient, turbulent | 1–10+       | Outer corrections allow higher CFL          |
| simpleFoam    | SIMPLE    | Steady-state         | N/A         | No physical time stepping; no CFL concern   |
| interFoam     | PISO/VOF  | Multiphase transient | ≤ 0.5       | Interface CFL (maxAlphaCo) also applies     |
| interIsoFoam  | PISO/isoAdvector | Multiphase    | ≤ 0.5       | Geometric interface reconstruction          |
| rhoPimpleFoam | PIMPLE    | Compressible trans.  | ≤ 1         | Acoustic CFL may dominate                   |
| rhoSimpleFoam | SIMPLE    | Compressible steady  | N/A         | Steady-state, no CFL                        |
| sonicFoam     | PISO      | Compressible trans.  | ≤ 0.5       | Acoustic CFL is very restrictive            |
| buoyantPimpleFoam | PIMPLE | Buoyant transient  | 1–5         | Check both flow and buoyancy scales         |

```
    CFL Tolerance Spectrum:

    STRICT                                                 RELAXED
    CFL ≤ 0.5          CFL ≤ 1           CFL ~ 1-10       No CFL
    ├──────────────────┼─────────────────┼────────────────┤
    │                  │                 │                │
    interFoam       icoFoam          pimpleFoam      simpleFoam
    sonicFoam       pisoFoam                         rhoSimpleFoam
                    rhoPimpleFoam
```

---

## Acoustic CFL — Compressible Flows

For incompressible flows, the CFL condition involves only the flow velocity u.
But for compressible flows, pressure waves (sound) travel at the speed of sound
c, which is often *much faster* than the flow itself.

### The Acoustic CFL Condition

```
                  (|u| + c) · Δt
    CFL_acoustic = ────────────────
                        Δx

    Where:
      u = flow velocity
      c = speed of sound (343 m/s in air at 20°C!)
      Δt = time step
      Δx = cell size
```

### Why This Matters

```
    Incompressible (water in a pipe):     Compressible (air around a wing):

    u = 2 m/s                             u = 100 m/s
    Δx = 0.01 m                           Δx = 0.01 m
    CFL ≤ 1                               c = 343 m/s
                                          CFL uses (u + c) = 443 m/s!

    Δt ≤ Δx/u                             Δt ≤ Δx/(u+c)
    Δt ≤ 0.01/2                            Δt ≤ 0.01/443
    Δt ≤ 0.005 s                           Δt ≤ 0.0000226 s  ← 220× smaller!

    This is why compressible explicit solvers need TINY time steps
    compared to incompressible solvers on the same mesh.
```

> **Warning**: If you switch from an incompressible solver (icoFoam) to a
> compressible solver (sonicFoam) on a similar mesh, expect to need a time step
> that is **orders of magnitude smaller**. Plan your compute budget accordingly.

---

## Interface CFL for Multiphase Flows

Multiphase solvers like `interFoam` track the interface between fluids (e.g.,
water and air) using the Volume of Fluid (VOF) method. The interface is
represented by the volume fraction field α (alpha), and it requires its own
CFL condition.

### Why a Separate Interface CFL?

The α field is advected by a special equation with a compression term to keep
the interface sharp. This equation is more sensitive to CFL violations than the
momentum equation:

```
    ┌─────────────────────────────────────────────────────┐
    │                                                     │
    │   CFL too high for interface:                       │
    │                                                     │
    │   α = 1 (water)      α = 0 (air)                   │
    │   ████████████░░░░░░░░░░░░░░░░░░    BEFORE          │
    │                                                     │
    │   ████████░░░░████░░░░░░░░░░░░░░    AFTER           │
    │          ↑         ↑                                │
    │        Smeared    Fragments                         │
    │        interface  of water                          │
    │                   in air!                           │
    │                                                     │
    │   The interface becomes diffuse and unphysical.     │
    │   Droplets and bubbles appear as numerical          │
    │   artifacts, not real physics.                      │
    └─────────────────────────────────────────────────────┘
```

### OpenFOAM Configuration for Multiphase CFL

In `interFoam`, you control both the flow CFL and the interface CFL:

```c
adjustTimeStep  yes;
maxCo           0.5;        // Maximum flow Courant number
maxAlphaCo      0.3;        // Maximum interface Courant number (stricter!)
maxDeltaT       1e-3;
```

The solver takes the *more restrictive* of the two:
`Δt = min(Δt from maxCo, Δt from maxAlphaCo)`

> **Tip**: For free-surface flows (waves, sloshing, dam breaks), keep
> `maxAlphaCo` at 0.2–0.4. A sharp interface is critical for accurate results.
> The extra cost of smaller time steps is usually worth it.

---

## Practical Guidelines Summary

Here are recommended CFL values for common simulation scenarios:

| Scenario                             | Recommended CFL | Justification                                     |
|--------------------------------------|:---------------:|---------------------------------------------------|
| Laminar pipe flow (icoFoam)          | 0.5 – 1.0      | Well-behaved flow, standard PISO limit             |
| Turbulent external flow (pisoFoam)   | 0.5 – 0.8      | Leave margin for velocity fluctuations             |
| LES / DES simulations                | 0.3 – 0.5      | Accuracy of turbulent structures requires low CFL  |
| PIMPLE with outer corrections        | 1.0 – 5.0      | Outer loops recover stability; start conservative  |
| Steady-state RANS (simpleFoam)       | N/A             | No CFL constraint; use under-relaxation instead    |
| Free-surface / VOF (interFoam)       | 0.2 – 0.5      | Interface sharpness demands tight CFL              |
| Compressible subsonic (rhoPimpleFoam)| 0.5 – 1.0      | Account for acoustic CFL                           |
| Compressible transonic/supersonic    | 0.3 – 0.5      | Shock capturing needs conservative CFL             |
| Combustion / reacting flows          | 0.2 – 0.5      | Chemical time scales may be very short             |
| Startup / initialization phase       | 0.1 – 0.3      | Flow field is far from physical; reduce CFL early  |

> **Rule of Thumb**: When in doubt, start with CFL = 0.5 and `adjustTimeStep
> yes`. Monitor the simulation for a few hundred steps. If it's stable and the
> Courant number is well below your limit, you can cautiously increase maxCo.

---

## Troubleshooting CFL-Related Problems

### Problem: Simulation Diverges

```
    Symptom:  "FOAM FATAL ERROR" or residuals suddenly jump to 1e+20

    Diagnostic steps:
    ┌──────────────────────────────────────────────────────┐
    │  1. Check the last Courant number before crash       │
    │     → grep "Courant Number" log.solver | tail -20    │
    │                                                      │
    │  2. Is max CFL >> 1?                                 │
    │     YES → Reduce deltaT or enable adjustTimeStep     │
    │     NO  → Problem is elsewhere (BCs, mesh, etc.)     │
    │                                                      │
    │  3. Did CFL suddenly spike?                          │
    │     YES → Check for mesh quality issues near the     │
    │           spike location (very small cells)          │
    │     NO  → Gradual growth suggests scheme instability │
    └──────────────────────────────────────────────────────┘
```

### Problem: CFL is Too Low (Wasting Compute Time)

If your max Courant number is consistently 0.001, you are taking time steps
that are 1000× smaller than necessary:

```
    Solution options:
    ┌──────────────────────────────────────────────────────┐
    │  • Increase deltaT (if using fixed time step)        │
    │  • Enable adjustTimeStep with appropriate maxCo      │
    │  • Check if your mesh has a few extremely small       │
    │    cells that force tiny Δt for the whole domain     │
    │    → Coarsen those cells or use local time stepping  │
    └──────────────────────────────────────────────────────┘
```

### Problem: CFL Varies Wildly Across the Domain

This usually means highly non-uniform mesh or large velocity gradients:

```
    ┌─────────────────────────────────────────────────┐
    │                                                 │
    │  Near wall (tiny cells):    CFL = 8.5  ← !!    │
    │  In freestream (big cells): CFL = 0.01          │
    │                                                 │
    │  The solver reports max CFL = 8.5               │
    │  but mean CFL = 0.05                            │
    │                                                 │
    │  Fix: Use adjustTimeStep (lets Δt adapt)        │
    │       or switch to pimpleFoam (tolerates high   │
    │       CFL with outer corrections)               │
    │       or improve mesh grading (smoother          │
    │       cell size transitions)                    │
    └─────────────────────────────────────────────────┘
```

> **Tip**: Use `foamMonitor` or plot the Courant number field in ParaView
> (`Co` field if written) to visualize *where* CFL is highest. This often
> reveals mesh problems you did not expect.

---

## Calculating CFL Before Running — Estimate Your Time Step

Do not guess your time step. Estimate it before running:

### The Formula

```
                CFL_target × Δx_min
    Δt_est  =  ────────────────────
                    u_max

    Where:
      CFL_target = your desired maximum CFL (e.g., 0.5)
      Δx_min     = your smallest cell size (check with checkMesh)
      u_max      = expected maximum velocity in the domain
```

### Step-by-Step Estimation

```
    ┌─────────────────────────────────────────────────────────────┐
    │                                                             │
    │  Step 1: Find smallest cell size                            │
    │  ─────────────────────────────                              │
    │  Run:  checkMesh | grep "Min volume"                        │
    │  Δx_min ≈ (Min volume)^(1/3)    for hex cells              │
    │                                                             │
    │  Example: Min volume = 1e-9 m³                              │
    │           Δx_min ≈ (1e-9)^(1/3) = 0.001 m = 1 mm          │
    │                                                             │
    │  Step 2: Estimate maximum velocity                          │
    │  ────────────────────────────────                           │
    │  From boundary conditions, inlet velocity, or               │
    │  analytical estimates.                                      │
    │                                                             │
    │  Example: Lid-driven cavity with U_lid = 1 m/s             │
    │           u_max ≈ 1 m/s (may be slightly higher locally)   │
    │                                                             │
    │  Step 3: Compute Δt                                         │
    │  ──────────────────                                         │
    │  Δt = CFL_target × Δx_min / u_max                          │
    │  Δt = 0.5 × 0.001 / 1.0                                    │
    │  Δt = 0.0005 s                                              │
    │                                                             │
    │  Step 4: Round down for safety                              │
    │  ────────────────────────────                               │
    │  Use Δt = 0.0005 s or Δt = 0.001 s with CFL monitoring     │
    │                                                             │
    └─────────────────────────────────────────────────────────────┘
```

### Quick Reference Calculation

```
    Example: Flow around a cylinder
    ─────────────────────────────────
    Inlet velocity:     u = 10 m/s
    Smallest cell:      Δx = 0.0005 m  (near cylinder wall)
    Target CFL:         0.5

    Δt = 0.5 × 0.0005 / 10 = 2.5e-5 s = 25 μs

    At endTime = 1 s, that is 40,000 time steps.
    Consider: Is this feasible? If not, coarsen the mesh
    or use pimpleFoam with higher CFL tolerance.
```

---

## Key Takeaways

```
    ╔══════════════════════════════════════════════════════════════╗
    ║                    CFL NUMBER CHEAT SHEET                   ║
    ╠══════════════════════════════════════════════════════════════╣
    ║                                                              ║
    ║  CFL = u·Δt/Δx     Keep ≤ 1 for explicit solvers            ║
    ║                                                              ║
    ║  ALWAYS check max Courant number, not just the mean          ║
    ║                                                              ║
    ║  Use adjustTimeStep for complex transient flows              ║
    ║                                                              ║
    ║  Steady-state solvers (SIMPLE) have NO CFL constraint        ║
    ║                                                              ║
    ║  PIMPLE can exceed CFL = 1 thanks to outer corrections       ║
    ║                                                              ║
    ║  Compressible flows: CFL uses (u + c), not just u            ║
    ║                                                              ║
    ║  Multiphase flows: check BOTH maxCo AND maxAlphaCo           ║
    ║                                                              ║
    ║  Estimate Δt BEFORE running: Δt = CFL × Δx_min / u_max      ║
    ║                                                              ║
    ║  When in doubt: CFL = 0.5, adjustTimeStep yes                ║
    ║                                                              ║
    ╚══════════════════════════════════════════════════════════════╝
```
