# Linear Solvers, Preconditioners, and Solution Control in OpenFOAM

## Table of Contents

1. [Introduction: Why Linear Solvers Matter](#1-introduction-why-linear-solvers-matter)
2. [From PDEs to Linear Systems](#2-from-pdes-to-linear-systems)
3. [Solver Categories Overview](#3-solver-categories-overview)
4. [OpenFOAM Solvers in Detail](#4-openfoam-solvers-in-detail)
5. [Preconditioners](#5-preconditioners)
6. [Smoothers](#6-smoothers)
7. [Convergence Control](#7-convergence-control)
8. [Solver Decision Tree](#8-solver-decision-tree)
9. [Real Project Examples](#9-real-project-examples)
10. [Relaxation Factors](#10-relaxation-factors)
11. [Pressure-Velocity Coupling Algorithms](#11-pressure-velocity-coupling-algorithms)
12. [Monitoring Convergence](#12-monitoring-convergence)
13. [Troubleshooting](#13-troubleshooting)
14. [Quick Reference Tables](#14-quick-reference-tables)

---

## 1. Introduction: Why Linear Solvers Matter

Every CFD simulation ultimately boils down to solving enormous systems of linear equations.
The Navier-Stokes equations are **nonlinear partial differential equations**, but after
discretization (finite volume method) and linearization, each time step or outer iteration
requires solving one or more **linear systems** of the form **Ax = b**.

The linear solver is often the **single most expensive** part of a CFD simulation. Choosing
the wrong solver can mean the difference between a simulation that converges in minutes
versus one that takes hours — or never converges at all.

```
┌─────────────────────────────────────────────────────────────────────┐
│                     THE CFD SIMULATION PIPELINE                     │
│                                                                     │
│  Physical       Governing        Discretized       Linear           │
│  Problem   ───► Equations   ───► Equations    ───► System  ───► x   │
│                 (PDEs)           (algebraic)       Ax = b           │
│                                                                     │
│  Fluid flow     Navier-Stokes    Finite Volume     PCG, GAMG,      │
│  around a       + turbulence     Method on mesh    smoothSolver     │
│  wing           models                             solve this!      │
│                                                                     │
│  ⏱ Time spent:  ~0%              ~5-10%            ~85-95%          │
│                                                                     │
│  The linear solver dominates total computation time.                │
└─────────────────────────────────────────────────────────────────────┘
```

This document covers everything you need to know about configuring linear solvers in
OpenFOAM: which solver to pick, what preconditioner to pair it with, how to set
convergence tolerances, and how to troubleshoot when things go wrong.

---

## 2. From PDEs to Linear Systems

### The Discretization Process

Consider a simple 1D heat conduction (Laplace) equation. The finite volume method turns
the continuous PDE into a sparse matrix system:

```
  CONTINUOUS PDE                    DISCRETIZED LINEAR SYSTEM

  ∂²T                              ┌                    ┐ ┌    ┐   ┌        ┐
  ─── = 0                          │  2  -1   0   0   0 │ │ T₁ │   │ T_left │
  ∂x²                              │ -1   2  -1   0   0 │ │ T₂ │   │   0    │
                                   │  0  -1   2  -1   0 │ │ T₃ │ = │   0    │
  on a 1D domain                   │  0   0  -1   2  -1 │ │ T₄ │   │   0    │
  with boundary                    │  0   0   0  -1   2 │ │ T₅ │   │ T_right│
  conditions                       └                    ┘ └    ┘   └        ┘

                                              A        ·  x    =     b

  Continuous world                  Discrete world: sparse linear system
  (infinite DOF)                    (finite DOF, one per cell)
```

### Matrix Properties in CFD

The properties of matrix **A** determine which solvers we can use:

| Property | Meaning | Which equations? |
|----------|---------|------------------|
| **Symmetric** | A = Aᵀ (aᵢⱼ = aⱼᵢ) | Pressure (Laplacian) |
| **Asymmetric** | A ≠ Aᵀ | Velocity, turbulence (convection terms) |
| **Positive-definite** | All eigenvalues > 0 | Pressure (with proper BCs) |
| **Sparse** | Mostly zeros | ALL CFD matrices (only neighbors are non-zero) |
| **Banded** | Non-zeros near diagonal | Structured meshes |
| **Diagonally dominant** | |aᵢᵢ| ≥ Σⱼ≠ᵢ |aᵢⱼ| | Well-posed CFD problems |

### Why Not Direct Solvers?

Direct solvers (Gaussian elimination, LU decomposition) give the **exact** answer, but:

```
  Problem size (cells)    Memory for direct solve    Time complexity
  ─────────────────────   ───────────────────────    ──────────────
       1,000              ~8 MB                      ~1 second
      10,000              ~800 MB                    ~minutes
     100,000              ~80 GB                     ~hours
   1,000,000              ~8 TB  ← impossible!       ~days
  10,000,000              forget about it            ~heat death

  Iterative solvers: memory scales as O(N), time as O(N) to O(N log N)
  Direct solvers:    memory scales as O(N²), time as O(N³)
```

For any realistic CFD problem (100k+ cells), **iterative solvers are the only option**.

---

## 3. Solver Categories Overview

```
  Linear Solvers in OpenFOAM
  │
  ├── Krylov Subspace Methods (build solution in expanding subspace)
  │   │
  │   ├── PCG ─────────── Preconditioned Conjugate Gradient
  │   │                    For: symmetric positive-definite (pressure)
  │   │
  │   ├── PBiCGStab ───── Preconditioned Bi-Conjugate Gradient Stabilized
  │   │                    For: asymmetric matrices (velocity, turbulence)
  │   │
  │   └── (PBiCG) ─────── Legacy, replaced by PBiCGStab
  │                        Avoid in new setups
  │
  ├── Smoothing-Based Solvers
  │   │
  │   └── smoothSolver ── Uses iterative smoothers (GS, symGS, DIC, DILU)
  │                        Versatile, works for symmetric & asymmetric
  │
  ├── Multigrid Methods
  │   │
  │   └── GAMG ────────── Geometric-Algebraic Multi-Grid
  │                        Best for large pressure systems
  │                        Uses coarsening + smoothing cycles
  │
  └── Diagonal Solver
      │
      └── diagonal ────── Direct diagonal solve (trivial systems only)
                           Used for density in compressible flows
```

### Solver Comparison Summary

| Solver | Matrix Type | Best For | Speed | Memory | Robustness |
|--------|------------|----------|-------|--------|------------|
| **PCG** | Symmetric PD | Pressure (small-medium) | ★★★☆ | ★★★★ | ★★★★ |
| **PBiCGStab** | Asymmetric | Velocity, turbulence | ★★★☆ | ★★★☆ | ★★★☆ |
| **smoothSolver** | Any | Velocity, turbulence, scalars | ★★★☆ | ★★★★ | ★★★★ |
| **GAMG** | Symmetric PD | Pressure (large cases) | ★★★★ | ★★★☆ | ★★★★ |
| **diagonal** | Diagonal | Trivial systems | ★★★★★ | ★★★★★ | ★★★★★ |

---

## 4. OpenFOAM Solvers in Detail

### 4.1 PCG — Preconditioned Conjugate Gradient

**What it does**: The Conjugate Gradient method builds an orthogonal search direction at
each iteration, guaranteeing convergence in at most N iterations (where N = matrix size).
With preconditioning, convergence is dramatically accelerated.

**Requirements**: Matrix must be **symmetric positive-definite** (SPD).

**When to use**: Pressure equation on small to medium meshes (< ~500k cells).

**How it works** (simplified):

```
  PCG Algorithm (conceptual)
  ═══════════════════════════════════════════════════════════

  Initialize: x₀ = initial guess, r₀ = b - A·x₀
               z₀ = M⁻¹·r₀       (apply preconditioner)
               p₀ = z₀             (first search direction)

  For k = 0, 1, 2, ...
  ┌─────────────────────────────────────────────────────────┐
  │  αₖ = (rₖᵀ·zₖ) / (pₖᵀ·A·pₖ)    ← step size           │
  │  xₖ₊₁ = xₖ + αₖ·pₖ              ← update solution     │
  │  rₖ₊₁ = rₖ - αₖ·A·pₖ            ← update residual     │
  │                                                         │
  │  Check: ‖rₖ₊₁‖ < tolerance?  → STOP (converged!)       │
  │                                                         │
  │  zₖ₊₁ = M⁻¹·rₖ₊₁                ← precondition        │
  │  βₖ = (rₖ₊₁ᵀ·zₖ₊₁)/(rₖᵀ·zₖ)    ← conjugacy coeff     │
  │  pₖ₊₁ = zₖ₊₁ + βₖ·pₖ            ← new search dir      │
  └─────────────────────────────────────────────────────────┘

  Convergence behavior (symmetric matrix):

  Residual ‖r‖
  │
  │╲
  │ ╲
  │  ╲
  │   ╲
  │    ╲
  │     ╲╲
  │       ╲╲╲
  │          ╲╲╲╲╲___________  ← monotonic decrease (guaranteed!)
  │
  └────────────────────────────── Iterations

  Key property: residual decreases MONOTONICALLY for SPD matrices
```

**Typical configuration** (from our lid-driven cavity project):

```c
p
{
    solver          PCG;        // Conjugate Gradient
    preconditioner  DIC;        // Diagonal Incomplete Cholesky
    tolerance       1e-06;      // Absolute residual threshold
    relTol          0.05;       // Relative residual drop (5%)
}
```

### 4.2 PBiCGStab — Preconditioned Bi-Conjugate Gradient Stabilized

**What it does**: Extends CG to **non-symmetric** matrices using a bi-orthogonal approach.
The "Stabilized" variant (van der Vorst, 1992) fixes the erratic convergence of plain BiCG.

**Requirements**: Works on **any** matrix (symmetric or asymmetric).

**When to use**: Velocity, turbulence quantities, or any equation with strong convection.

```
  PBiCGStab Convergence (asymmetric matrix):

  Residual ‖r‖
  │
  │╲
  │ ╲ ╱╲
  │  ╳   ╲
  │ ╱     ╲╱╲
  │╱         ╲
  │           ╲╱╲
  │              ╲╲
  │                ╲╲╲╲________  ← may oscillate but trends down
  │
  └────────────────────────────── Iterations

  Note: convergence is NOT monotonic — residual can temporarily increase!
  This is normal for asymmetric systems.
```

**Typical configuration**:

```c
U
{
    solver          PBiCGStab;    // Bi-Conjugate Gradient Stabilized
    preconditioner  DILU;         // Diagonal Incomplete LU
    tolerance       1e-05;
    relTol          0.1;
}
```

> **Note**: `PBiCG` (without "Stab") is the older variant. OpenFOAM still accepts it but
> internally uses PBiCGStab. Always prefer `PBiCGStab` in new configurations.

### 4.3 smoothSolver

**What it does**: Rather than building Krylov subspaces, smoothSolver simply applies a
**smoother** (iterative method like Gauss-Seidel) repeatedly until convergence. This is
simpler but surprisingly effective, especially for well-conditioned systems.

**When to use**: Velocity, turbulence, scalar transport. Very popular default choice.

**Available smoothers**:

| Smoother | Full Name | For Matrix Type | Notes |
|----------|-----------|----------------|-------|
| `GaussSeidel` | Gauss-Seidel | Any | Basic, robust |
| `symGaussSeidel` | Symmetric Gauss-Seidel | Any | Forward + backward sweep, better for symmetric |
| `DIC` | Diagonal Incomplete Cholesky | Symmetric | Also works as smoother |
| `DILU` | Diagonal Incomplete LU | Asymmetric | Also works as smoother |
| `DICGaussSeidel` | Combined DIC + GS | Symmetric | Good for pressure in GAMG |

**Typical configuration** (from our lid-driven cavity project):

```c
U
{
    solver          smoothSolver;
    smoother        symGaussSeidel;   // symmetric Gauss-Seidel smoother
    tolerance       1e-05;
    relTol          0;                // solve to absolute tolerance (transient)
}
```

**From our NACA airfoil project** (turbulence quantities):

```c
k
{
    solver          smoothSolver;
    smoother        symGaussSeidel;
    tolerance       1e-05;
    relTol          0.1;              // relative tolerance OK (steady-state)
}

epsilon
{
    solver          smoothSolver;
    smoother        symGaussSeidel;
    tolerance       1e-05;
    relTol          0.1;
}
```

### 4.4 GAMG — Geometric-Algebraic Multi-Grid

**What it does**: GAMG is the **most powerful solver** for pressure equations in large cases.
It exploits the fact that iterative smoothers (Gauss-Seidel) are great at removing
**high-frequency** errors but terrible at **low-frequency** errors. By coarsening the grid,
low-frequency errors on the fine grid become high-frequency errors on the coarse grid.

```
  THE MULTIGRID IDEA
  ══════════════════

  Problem: Gauss-Seidel quickly kills high-freq error but stalls on low-freq

  High-frequency error          Low-frequency error
  (oscillates cell-to-cell)     (smooth, varies slowly)

  ┊╱╲╱╲╱╲╱╲╱╲╱╲╱╲╱╲┊          ┊                    ┊
  ┊                  ┊          ┊         ╱╲         ┊
  ┊                  ┊          ┊       ╱    ╲       ┊
  ┊                  ┊          ┊     ╱        ╲     ┊
  ┊                  ┊          ┊   ╱            ╲   ┊
  ┊                  ┊          ┊ ╱                ╲ ┊
  ┊                  ┊          ┊╱                  ╲┊

  GS kills this fast ✓          GS barely touches this ✗

  Solution: Transfer to coarser grid where low-freq → high-freq!

  Fine grid (original)       Coarse grid             Coarsest grid
  ┌─┬─┬─┬─┬─┬─┬─┬─┐        ┌──┬──┬──┬──┐           ┌────┬────┐
  │ │ │ │ │ │ │ │ │        │  │  │  │  │           │    │    │
  ├─┼─┼─┼─┼─┼─┼─┼─┤        ├──┼──┼──┼──┤           ├────┼────┤
  │ │ │ │ │ │ │ │ │        │  │  │  │  │           │    │    │
  ├─┼─┼─┼─┼─┼─┼─┼─┤        ├──┼──┼──┼──┤           └────┴────┘
  │ │ │ │ │ │ │ │ │        │  │  │  │  │           4 cells: direct
  ├─┼─┼─┼─┼─┼─┼─┼─┤        ├──┼──┼──┼──┤           solve is cheap!
  │ │ │ │ │ │ │ │ │        │  │  │  │  │
  └─┴─┴─┴─┴─┴─┴─┴─┘        └──┴──┴──┴──┘
  64 cells                   16 cells
```

**The V-Cycle**:

```
  Level 0 (fine)     Smooth ──────────────────────────── Smooth
                          \                             /
  Level 1 (coarse)         Smooth ────────────── Smooth
                                \               /
  Level 2 (coarser)              Smooth ── Smooth
                                      \   /
  Level 3 (coarsest)                  SOLVE
                                  (direct or few
                                   iterations)

       RESTRICT ──────────────────►   ◄─────────────── PROLONGATE
    (fine → coarse)                                  (coarse → fine)
    transfer residual                                interpolate correction
```

**Why GAMG excels for pressure**: The pressure equation (Laplacian) is the hardest to solve
because it is **elliptic** — information propagates everywhere. A change at one boundary
affects cells far away. GAMG's coarsening captures these long-range interactions efficiently.

**Typical GAMG configuration**:

```c
p
{
    solver          GAMG;
    smoother        GaussSeidel;       // smoother at each level
    tolerance       1e-06;
    relTol          0.01;

    // GAMG-specific settings
    nCellsInCoarsestLevel  10;         // stop coarsening at ~10 cells
    agglomerator    faceAreaPair;      // how to group cells
    mergeLevels     1;                 // aggressiveness of coarsening
    cacheAgglomeration true;           // reuse coarsening (faster)
    nPreSweeps      0;                 // smoothing before restriction
    nPostSweeps     2;                 // smoothing after prolongation
    nFinestSweeps   2;                 // smoothing on finest level
}
```

**GAMG-specific settings explained**:

| Setting | Default | Purpose |
|---------|---------|---------|
| `nCellsInCoarsestLevel` | 10 | Min cells before direct solve on coarsest level |
| `agglomerator` | `faceAreaPair` | Cell grouping strategy for coarsening |
| `mergeLevels` | 1 | How aggressively to coarsen (higher = more aggressive) |
| `cacheAgglomeration` | `true` | Cache the coarsening hierarchy (saves time) |
| `nPreSweeps` | 0 | Smoother sweeps before restricting to coarser level |
| `nPostSweeps` | 2 | Smoother sweeps after prolongating from coarser level |
| `nFinestSweeps` | 2 | Smoother sweeps on the finest (original) level |
| `maxIter` | 1000 | Maximum GAMG cycles |
| `minIter` | 1 | Minimum GAMG cycles (ensures at least 1 solve) |

---

## 5. Preconditioners

### What Preconditioners Do

Preconditioning transforms the linear system into one that is easier to solve. Instead of
solving **Ax = b** directly, we solve **M⁻¹Ax = M⁻¹b** where **M** approximates **A**.

```
  Without preconditioning:         With preconditioning:

  Solve: Ax = b                    Solve: M⁻¹Ax = M⁻¹b
  Condition number: κ(A) = 10⁶     Condition number: κ(M⁻¹A) = 10²

  Iterations needed: ~1000          Iterations needed: ~10
  ─────────────────────────         ──────────────────────────
  Slow!                             Fast!

  The ideal preconditioner: M = A  →  M⁻¹A = I  →  1 iteration!
  But computing M⁻¹ = A⁻¹ is as hard as solving the original problem.

  So we use CHEAP APPROXIMATIONS of A:
  ┌──────────────────────────────────────────────────────────┐
  │  M ≈ A,  but M⁻¹ is cheap to compute                    │
  │                                                          │
  │  Trade-off: better approximation = fewer iterations      │
  │             but more expensive per iteration              │
  └──────────────────────────────────────────────────────────┘
```

### OpenFOAM Preconditioners

| Preconditioner | Full Name | Matrix Type | Pairs With | Typical Use |
|---------------|-----------|-------------|-----------|-------------|
| **DIC** | Diagonal Incomplete Cholesky | Symmetric | PCG | Pressure |
| **FDIC** | Faster DIC | Symmetric | PCG | Pressure (faster, less robust) |
| **DILU** | Diagonal Incomplete LU | Asymmetric | PBiCGStab | Velocity, turbulence |
| **diagonal** | Diagonal (Jacobi) | Any | Any Krylov | Very cheap, less effective |
| **GAMG** | GAMG as preconditioner | Symmetric | PCG | Pressure (large cases) |
| **none** | No preconditioning | Any | Any Krylov | Almost never useful |

### DIC — Diagonal Incomplete Cholesky

For **symmetric** matrices. Approximates the Cholesky decomposition (A = LLᵀ) but only
fills in entries where A already has non-zeros (hence "incomplete"). The "diagonal" variant
modifies only the diagonal to maintain the correct row sums.

```c
p
{
    solver          PCG;
    preconditioner  DIC;      // ← standard choice for pressure with PCG
    tolerance       1e-06;
    relTol          0.05;
}
```

### DILU — Diagonal Incomplete LU

For **asymmetric** matrices. Approximates the LU decomposition (A = LU) but only on the
existing sparsity pattern. The "diagonal" variant is a simplified form that is cheaper
to compute.

```c
U
{
    solver          PBiCGStab;
    preconditioner  DILU;     // ← standard choice for velocity with PBiCGStab
    tolerance       1e-05;
    relTol          0.1;
}
```

### FDIC — Faster Diagonal Incomplete Cholesky

A performance-optimized variant of DIC. Uses a simplified computation that is faster per
iteration but may require slightly more iterations. Good trade-off for large meshes.

### diagonal — Diagonal (Jacobi) Preconditioner

The simplest possible preconditioner: M = diag(A). Very cheap to apply, but provides
minimal improvement in condition number. Rarely the best choice.

### GAMG as Preconditioner

GAMG can be used as a preconditioner for PCG, combining the optimal convergence properties
of CG with the multi-level error reduction of GAMG:

```c
p
{
    solver          PCG;
    preconditioner
    {
        preconditioner  GAMG;
        smoother        DICGaussSeidel;
        tolerance       1e-06;
        relTol          0.01;
        nCellsInCoarsestLevel 10;
    }
    tolerance       1e-06;
    relTol          0.01;
}
```

---

## 6. Smoothers

Smoothers are iterative methods used within smoothSolver and GAMG. They "smooth" the error,
removing high-frequency components:

| Smoother | Full Name | Description |
|----------|-----------|-------------|
| `GaussSeidel` | Gauss-Seidel | Forward sweep, uses newest values immediately |
| `symGaussSeidel` | Symmetric Gauss-Seidel | Forward + backward sweep, better smoothing |
| `DIC` | Diagonal Incomplete Cholesky | Used as smoother (symmetric matrices) |
| `DILU` | Diagonal Incomplete LU | Used as smoother (asymmetric matrices) |
| `DICGaussSeidel` | Combined DIC + Gauss-Seidel | Popular for GAMG pressure |

### How Gauss-Seidel Smoothing Works

```
  Original error (before smoothing):
  ┊                                         ┊
  ┊   ╱╲   ╱╲      ╱╲     ╱╲╱╲╱╲           ┊  ← mix of high-freq
  ┊  ╱  ╲ ╱  ╲    ╱  ╲   ╱      ╲          ┊    and low-freq error
  ┊ ╱    ╳    ╲  ╱    ╲ ╱        ╲         ┊
  ┊╱          ╲╱      ╲╱          ╲────────┊

  After a few GS sweeps:
  ┊                                         ┊
  ┊                                         ┊
  ┊       ╱──╲         ╱──╲                 ┊  ← high-freq gone!
  ┊     ╱     ╲      ╱     ╲               ┊    only smooth (low-freq)
  ┊   ╱        ╲   ╱        ╲──────────── ┊    error remains
  ┊─╱           ╲╱                         ┊

  This is exactly what GAMG needs: smooth remaining error transfers
  well to the coarse grid for efficient removal there.
```

---

## 7. Convergence Control

### tolerance vs relTol

These two settings work together to determine when the solver stops iterating:

```
  Residual
  │
  │ r₀ (initial residual)
  │╲
  │ ╲
  │  ╲         relTol threshold = r₀ × relTol
  │   ╲         (relative drop from initial)
  │    ╲        ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─
  │     ╲╲
  │       ╲╲
  │         ╲╲╲╲╲
  │              ╲╲╲╲╲╲___
  │                       ╲╲╲╲╲╲______________ tolerance (absolute floor)
  │─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─
  │
  └──────────────────────────────────────────── Iterations

  Solver STOPS when EITHER condition is met:
    ‖residual‖ < tolerance         (absolute convergence)
    ‖residual‖ < r₀ × relTol      (relative convergence)

  For FINAL correctors (pFinal), set relTol = 0
  to force solving to the absolute tolerance.
```

**Key parameters**:

| Parameter | Type | Meaning | Typical Values |
|-----------|------|---------|---------------|
| `tolerance` | Absolute | Stop if residual falls below this | 1e-06 (pressure), 1e-05 (velocity) |
| `relTol` | Relative | Stop if residual drops by this factor | 0.01–0.1 (steady), 0 (final corrector) |
| `minIter` | Integer | Minimum iterations (always do at least this many) | 1 |
| `maxIter` | Integer | Maximum iterations (give up after this) | 1000 (default) |

### Steady-State vs Transient Convergence Strategy

```
  STEADY-STATE (SIMPLE):                  TRANSIENT (PISO/PIMPLE):

  ┌────────────────────────┐              ┌─────────────────────────────┐
  │ Each outer iteration   │              │ Each time step MUST be      │
  │ improves solution      │              │ fully converged for         │
  │ incrementally.         │              │ time-accuracy.              │
  │                        │              │                             │
  │ relTol = 0.01–0.1     │              │ relTol = 0 for final solve  │
  │ (no need to fully     │              │ (use tolerance as target)   │
  │  converge each step)  │              │                             │
  │                        │              │ Use pFinal with relTol = 0  │
  │ Convergence comes from │              │ to ensure tight solve on    │
  │ outer SIMPLE loop      │              │ last PISO corrector         │
  └────────────────────────┘              └─────────────────────────────┘
```

### Reading Solver Output

When OpenFOAM runs, each solver reports its progress:

```
  Typical solver output line:

  smoothSolver:  Solving for Ux, Initial residual = 0.0123, Final residual = 4.56e-07, No Iterations 3
  │              │                │                         │                          │
  │              │                │                         │                          └─ iterations used
  │              │                │                         └─ residual at end
  │              │                └─ residual at start of this solve
  │              └─ which field variable
  └─ which solver

  GAMG:  Solving for p, Initial residual = 0.456, Final residual = 0.00345, No Iterations 7
```

**What the numbers mean**:
- **Initial residual**: How far off the current solution is before solving
- **Final residual**: How far off after solving (should be much smaller)
- **No Iterations**: How many iterations the solver needed

---

## 8. Solver Decision Tree

Use this flowchart to choose the right solver for each equation:

```
  ╔═══════════════════════════════════════════════════════════════════╗
  ║              WHICH SOLVER SHOULD I USE?                          ║
  ╠═══════════════════════════════════════════════════════════════════╣
  ║                                                                   ║
  ║   What equation are you solving?                                  ║
  ║   │                                                               ║
  ║   ├── PRESSURE (p, p_rgh) ── symmetric positive-definite          ║
  ║   │   │                                                           ║
  ║   │   ├── Small mesh (<500k cells)?                               ║
  ║   │   │   │                                                       ║
  ║   │   │   ├── YES ──► PCG + DIC                                   ║
  ║   │   │   │           Simple, reliable, low memory                ║
  ║   │   │   │                                                       ║
  ║   │   │   └── NO ───► GAMG + GaussSeidel                          ║
  ║   │   │               Scales well, fastest for large meshes       ║
  ║   │   │                                                           ║
  ║   │   └── Struggling with convergence?                            ║
  ║   │       └── Try: PCG + GAMG (as preconditioner)                 ║
  ║   │                                                               ║
  ║   ├── VELOCITY (U) ── asymmetric                                  ║
  ║   │   │                                                           ║
  ║   │   ├── Default ──────► smoothSolver + symGaussSeidel           ║
  ║   │   │                   Most popular, good all-rounder          ║
  ║   │   │                                                           ║
  ║   │   └── Alternative ──► PBiCGStab + DILU                        ║
  ║   │                       Sometimes faster, less robust           ║
  ║   │                                                               ║
  ║   ├── TURBULENCE (k, ε, ω, nuTilda) ── asymmetric                ║
  ║   │   │                                                           ║
  ║   │   ├── Default ──────► smoothSolver + symGaussSeidel           ║
  ║   │   │                                                           ║
  ║   │   └── Alternative ──► PBiCGStab + DILU                        ║
  ║   │                                                               ║
  ║   └── SCALAR TRANSPORT (T, species) ── usually asymmetric         ║
  ║       │                                                           ║
  ║       └── Default ──────► smoothSolver + symGaussSeidel           ║
  ║                           OR PBiCGStab + DILU                      ║
  ║                                                                   ║
  ╚═══════════════════════════════════════════════════════════════════╝
```

---

## 9. Real Project Examples

### 9.1 Lid-Driven Cavity — Transient (PISO)

**Source**: `projects/01_lid_driven_cavity/system/fvSolution`

This is a transient, incompressible flow solved with PISO. The mesh is small, so simple
solvers suffice.

```c
// ─────────────────────────────────────────────────────────────
// fvSolution — Lid-Driven Cavity (icoFoam with PISO)
// ─────────────────────────────────────────────────────────────

solvers
{
    // ── PRESSURE ──────────────────────────────────────────────
    // Symmetric positive-definite → use PCG + DIC
    // Small mesh → no need for GAMG
    p
    {
        solver          PCG;            // Conjugate Gradient
        preconditioner  DIC;            // Incomplete Cholesky
        tolerance       1e-06;          // Absolute convergence target
        relTol          0.05;           // 5% relative drop is enough
    }                                   // (not the final corrector)

    // ── PRESSURE (final PISO corrector) ──────────────────────
    // Inherits from p, but with relTol=0 → solve to tolerance
    pFinal
    {
        $p;                             // Copy all settings from p
        relTol          0;              // Must reach absolute tolerance!
    }

    // ── VELOCITY ──────────────────────────────────────────────
    // Asymmetric matrix → smoothSolver with symGaussSeidel
    // relTol=0 because transient: each time step must converge
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0;
    }
}

// ── PRESSURE-VELOCITY COUPLING ────────────────────────────────
// PISO: predictor-corrector for transient flows
PISO
{
    nCorrectors     2;                  // 2 pressure corrector steps
    nNonOrthogonalCorrectors 0;         // 0 = orthogonal mesh assumed
    pRefCell        0;                  // Reference pressure cell
    pRefValue       0;                  // Reference pressure value
}
```

**Key observations**:
- Small transient case → PCG+DIC is sufficient for pressure
- `pFinal` uses `$p` to inherit settings, then overrides `relTol` to 0
- Velocity uses `relTol 0` because every time step must be properly converged
- PISO with 2 correctors is standard for transient incompressible flow
- `pRefCell` and `pRefValue` set the pressure reference (needed when no pressure BC)

### 9.2 NACA Airfoil — Steady-State (SIMPLE)

**Source**: `projects/04_naca_airfoil_analysis/system/fvSolution`

This is a steady-state, turbulent (RANS) simulation using SIMPLE. The approach is different:
relaxation factors are critical, and solvers don't need to fully converge each outer iteration.

```c
// ─────────────────────────────────────────────────────────────
// fvSolution — NACA Airfoil (simpleFoam with SIMPLE)
// ─────────────────────────────────────────────────────────────

solvers
{
    // ── PRESSURE ──────────────────────────────────────────────
    // Even for steady-state, pressure is the hardest equation.
    // PCG + DIC works well for moderate mesh sizes.
    p
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-06;
        relTol          0.05;           // Don't fully converge each
    }                                   // SIMPLE iteration

    // ── VELOCITY ──────────────────────────────────────────────
    // smoothSolver is the workhorse for velocity
    U
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;            // 10% relative drop per iteration
    }

    // ── TURBULENT KINETIC ENERGY ──────────────────────────────
    k
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }

    // ── TURBULENT DISSIPATION RATE ────────────────────────────
    epsilon
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-05;
        relTol          0.1;
    }
}

// ── SIMPLE ALGORITHM ──────────────────────────────────────────
// Steady-state pressure-velocity coupling with under-relaxation
SIMPLE
{
    nNonOrthogonalCorrectors 0;         // Increase for non-orthogonal meshes
    consistent yes;                     // SIMPLEC: modified pressure correction
                                        // allows higher relaxation factors

    // Stop the simulation when all residuals drop below these:
    residualControl
    {
        p       1e-3;                   // Pressure converged at 1e-3
        U       1e-4;                   // Velocity converged at 1e-4
        "(k|epsilon)" 1e-4;             // Turbulence converged at 1e-4
    }
}

// ── PISO (not used by simpleFoam, but included for reference) ──
PISO
{
    nCorrectors     2;
    nNonOrthogonalCorrectors 0;
}

// ── RELAXATION FACTORS ────────────────────────────────────────
// Critical for SIMPLE stability! Without these, SIMPLE diverges.
relaxationFactors
{
    fields
    {
        p       0.3;                    // Pressure: conservative (0.3)
    }
    equations
    {
        U       0.7;                    // Velocity: moderate (0.7)
        k       0.7;                    // Turb. KE: moderate (0.7)
        epsilon 0.7;                    // Turb. dissipation: moderate (0.7)
    }
}
```

**Key observations**:
- Steady-state → `relTol > 0` for all solvers (partial convergence each iteration is fine)
- `consistent yes` enables SIMPLEC, which is more stable and allows larger relaxation factors
- `residualControl` tells the SIMPLE loop when to stop the entire simulation
- All turbulence quantities use the same solver configuration (common practice)
- Relaxation factors are **essential** — $p$=0.3, $U$=$k$=$\varepsilon$=0.7 are standard starting values

---

## 10. Relaxation Factors

### What Under-Relaxation Does

In steady-state SIMPLE, the solution is updated iteratively. Without relaxation, the updates
can overshoot and cause divergence. Under-relaxation **dampens** the update:

```
  x_new = x_old + α · (x_computed - x_old)

  where α = relaxation factor (0 < α ≤ 1)

  α = 1.0:  Full update (no relaxation) — fast but may diverge
  α = 0.5:  Half update — more stable, slower convergence
  α = 0.3:  30% update — very stable, much slower

  ┌────────────────────────────────────────────────────────────┐
  │  Convergence speed vs stability trade-off:                 │
  │                                                            │
  │  Diverges ←──────────────────────────────→ Never converges │
  │             α=1.0    α=0.7    α=0.3    α=0.01              │
  │             (fast    (good    (safe    (too slow,           │
  │              but     balance)  but      wastes time)        │
  │              risky)           slow)                         │
  │                                                            │
  │  Sweet spot: fastest α that doesn't diverge                │
  └────────────────────────────────────────────────────────────┘
```

### Recommended Starting Values

| Variable | Relaxation Factor | Notes |
|----------|------------------|-------|
| **p** (pressure) | 0.3 | Most sensitive — always low |
| **U** (velocity) | 0.7 | Standard starting point |
| **k** (turb. KE) | 0.7 | Match velocity |
| **ε** (turb. dissipation) | 0.7 | Match velocity |
| **ω** (specific dissipation) | 0.7 | Match velocity |
| **nuTilda** (SA model) | 0.7 | Match velocity |
| **T** (temperature) | 0.7–1.0 | Often less sensitive |

### Fields vs Equations

OpenFOAM distinguishes between relaxing **fields** and **equations**:

```c
relaxationFactors
{
    fields      // Applied to the FIELD VALUES after solving
    {
        p   0.3;    // Dampens pressure field update
    }
    equations   // Applied to the EQUATION (matrix) before solving
    {
        U   0.7;    // Modifies the diagonal of the U matrix
        k   0.7;    // This implicitly provides under-relaxation
    }
}
```

- **Field relaxation** (`fields`): `p_new = p_old + α·(p_solved - p_old)`
  Used for pressure because the pressure equation itself is not relaxed.
- **Equation relaxation** (`equations`): Modifies the diagonal dominance of the matrix.
  More robust than field relaxation for transport equations.

### Important: p + U relaxation factors should sum to ~1.0

A common rule of thumb for SIMPLE: `α_p + α_U ≈ 1.0`. So if p=0.3, then U=0.7.
This is not a strict requirement but is a good guideline. For SIMPLEC (`consistent yes`),
higher factors may work (e.g., p=0.5, U=0.9).

---

## 11. Pressure-Velocity Coupling Algorithms

The coupling between pressure and velocity is the central challenge in incompressible CFD.
OpenFOAM implements several algorithms, each configured in `fvSolution`:

### SIMPLE — Semi-Implicit Method for Pressure-Linked Equations

```
  ╔═══════════════════════════════════════════════════════╗
  ║  SIMPLE Algorithm (Steady-State)                     ║
  ║                                                       ║
  ║  ┌─────────────────────────────────────────────────┐  ║
  ║  │  1. Solve momentum (U) with guessed pressure     │  ║
  ║  │  2. Solve pressure correction equation (p)       │  ║
  ║  │  3. Correct velocity with pressure gradient      │  ║
  ║  │  4. Under-relax all fields                       │  ║
  ║  │  5. Check convergence → if not, go to step 1     │  ║
  ║  └─────────────────────────────────────────────────┘  ║
  ║                                                       ║
  ║  Key: Under-relaxation is ESSENTIAL for stability     ║
  ║  Config: relaxationFactors { ... }                    ║
  ╚═══════════════════════════════════════════════════════╝
```

```c
SIMPLE
{
    nNonOrthogonalCorrectors 0;   // Extra p solves for non-ortho meshes
    consistent      yes;           // SIMPLEC variant (recommended)

    residualControl
    {
        p       1e-3;
        U       1e-4;
        "(k|epsilon)" 1e-4;
    }
}
```

### PISO — Pressure-Implicit with Splitting of Operators

```
  ╔═══════════════════════════════════════════════════════╗
  ║  PISO Algorithm (Transient)                          ║
  ║                                                       ║
  ║  For each time step Δt:                               ║
  ║  ┌─────────────────────────────────────────────────┐  ║
  ║  │  1. Solve momentum predictor (U*)                │  ║
  ║  │  2. PISO corrector loop (nCorrectors times):     │  ║
  ║  │     a. Solve pressure equation (p)               │  ║
  ║  │     b. Correct velocity (U)                      │  ║
  ║  │  3. Advance to next time step                    │  ║
  ║  └─────────────────────────────────────────────────┘  ║
  ║                                                       ║
  ║  Key: NO under-relaxation needed                      ║
  ║  Config: nCorrectors (typically 2-4)                  ║
  ╚═══════════════════════════════════════════════════════╝
```

```c
PISO
{
    nCorrectors     2;            // Pressure corrector steps per time step
    nNonOrthogonalCorrectors 0;   // Extra p solves for non-ortho meshes
    pRefCell        0;            // Reference pressure (closed domains)
    pRefValue       0;
}
```

### PIMPLE — Merged PISO-SIMPLE

```
  ╔═══════════════════════════════════════════════════════╗
  ║  PIMPLE Algorithm (Transient, allows large Δt)       ║
  ║                                                       ║
  ║  For each time step Δt:                               ║
  ║  ┌─────────────────────────────────────────────────┐  ║
  ║  │  Outer loop (nOuterCorrectors times):            │  ║
  ║  │  │  1. Solve momentum (U) with under-relaxation  │  ║
  ║  │  │  2. Inner PISO loop (nCorrectors times):      │  ║
  ║  │  │     a. Solve pressure (p)                     │  ║
  ║  │  │     b. Correct velocity (U)                   │  ║
  ║  │  │  3. Solve turbulence equations                │  ║
  ║  │  └───────────────────────────────────────────────│  ║
  ║  │  4. Advance to next time step                    │  ║
  ║  └─────────────────────────────────────────────────┘  ║
  ║                                                       ║
  ║  Key: Combines SIMPLE outer loop + PISO inner loop    ║
  ║  Allows larger Courant numbers (Co >> 1)              ║
  ╚═══════════════════════════════════════════════════════╝
```

```c
PIMPLE
{
    nOuterCorrectors    2;        // SIMPLE-like outer iterations
    nCorrectors         1;        // PISO correctors per outer
    nNonOrthogonalCorrectors 0;

    // Under-relaxation (only if nOuterCorrectors > 1)
    // Often disabled for final outer iteration
}
```

### Algorithm Comparison

| Feature | SIMPLE | PISO | PIMPLE |
|---------|--------|------|--------|
| **Type** | Steady-state | Transient | Transient (flexible) |
| **Solver** | simpleFoam | icoFoam, pisoFoam | pimpleFoam |
| **Under-relaxation** | Required | Not needed | Optional |
| **Time step size** | N/A (no time) | Small (Co < 1) | Large OK (Co >> 1) |
| **Outer iterations** | Until convergence | 1 per time step | nOuterCorrectors |
| **Inner p corrections** | 1 | nCorrectors | nCorrectors |
| **Best for** | Steady RANS | Small Δt transient | LES, large Δt |

---

## 12. Monitoring Convergence

### What Good Convergence Looks Like

```
  GOOD: Steady monotonic decrease         BAD: Oscillating / stalled
                                          
  log(residual)                           log(residual)
  │                                       │
  │╲                                      │    ╱╲    ╱╲    ╱╲
  │ ╲                                     │╲  ╱  ╲  ╱  ╲  ╱
  │  ╲                                    │ ╲╱    ╲╱    ╲╱
  │   ╲                                   │
  │    ╲                                  │
  │     ╲                                 │
  │      ╲╲                               │
  │        ╲╲╲                             │
  │           ╲╲╲╲╲_____                   │
  └──────────────────────                 └──────────────────────
  Iterations                              Iterations

  GOOD: All residuals dropping            BAD: Diverging
  together (p, U, k, ε)                  
                                          log(residual)
                                          │           ╱
                                          │         ╱
                                          │       ╱
                                          │     ╱
                                          │   ╱
                                          │ ╱
                                          │╱
                                          └──────────────────────
                                          Iterations
```

### Reading Log Files

Extract residuals from OpenFOAM log files with `grep`:

```bash
# Extract pressure residuals
grep "Solving for p" log.simpleFoam | awk '{print $8}' | tr -d ','

# Extract velocity residuals
grep "Solving for Ux" log.simpleFoam | awk '{print $8}' | tr -d ','

# Quick convergence check: show final residual of last 10 iterations
grep "Solving for p" log.simpleFoam | tail -10
```

### Plotting Residuals

Use `gnuplot` or the built-in OpenFOAM `foamMonitor` utility:

```bash
# Using foamMonitor (if available)
foamMonitor -l postProcessing/residuals/0/residuals.dat

# Using gnuplot
gnuplot -e "
  set logscale y;
  set xlabel 'Iteration';
  set ylabel 'Residual';
  plot 'postProcessing/residuals/0/residuals.dat' using 1:2 title 'p' with lines,
       '' using 1:3 title 'Ux' with lines,
       '' using 1:4 title 'Uy' with lines;
  pause -1
"
```

### Convergence Criteria for Different Problem Types

| Problem Type | p residual | U residual | Turbulence | Notes |
|-------------|-----------|-----------|-----------|-------|
| **Engineering estimate** | 1e-3 | 1e-4 | 1e-4 | Quick, approximate |
| **Standard accuracy** | 1e-4 | 1e-5 | 1e-5 | Most applications |
| **High accuracy** | 1e-5 | 1e-6 | 1e-6 | Validation, research |
| **Machine precision** | 1e-8 | 1e-8 | 1e-8 | Usually unnecessary |

---

## 13. Troubleshooting

### Problem: Solution Not Converging

```
  Symptom: Residuals stall at a plateau, never reaching target

  log(residual)
  │╲
  │ ╲
  │  ╲╲╲╲
  │       ──────────────────────────── stuck!
  │
  └──────────────────────────────────── Iterations
```

**Possible causes and fixes**:

| Cause | Fix |
|-------|-----|
| Mesh too coarse | Refine mesh, especially near walls |
| Poor mesh quality | Check with `checkMesh`, fix non-orthogonality |
| Wrong boundary conditions | Verify BCs match physics |
| Turbulence model mismatch | Check y+ values, wall functions |
| relTol too large | Reduce relTol (e.g., 0.1 → 0.01) |
| Relaxation too aggressive | Reduce relaxation factors |
| Non-orthogonal mesh | Increase `nNonOrthogonalCorrectors` to 1 or 2 |

### Problem: Oscillating Residuals

```
  Symptom: Residuals bounce up and down without converging

  log(residual)
  │
  │  ╱╲    ╱╲    ╱╲    ╱╲
  │ ╱  ╲  ╱  ╲  ╱  ╲  ╱  ╲
  │╱    ╲╱    ╲╱    ╲╱    ╲╱
  │
  └──────────────────────────── Iterations
```

**Possible causes and fixes**:

| Cause | Fix |
|-------|-----|
| Relaxation too high | Reduce p to 0.2, U to 0.5 |
| Numerical scheme too aggressive | Switch from `linearUpwind` to `upwind` temporarily |
| Unsteady flow solved as steady | Switch to transient solver (pimpleFoam) |
| Mesh quality issues | Check aspect ratio, skewness |

### Problem: Divergence (Blow-Up)

```
  Symptom: Residuals shoot up, solver crashes with floating-point exception

  #0  Foam::error::printStack at ...
  --> FOAM FATAL ERROR: Floating point exception
```

**Possible causes and fixes**:

| Cause | Fix |
|-------|-----|
| CFL number too high | Reduce time step (Co < 1 for PISO) |
| Bad initial conditions | Start from a simpler solution or potentialFoam |
| Mesh problems | Run `checkMesh`, look for negative volumes |
| Incompatible BC combination | Check that BCs are physically consistent |
| Relaxation too high | Start with p=0.1, U=0.3 then increase gradually |
| Wrong solver type | Verify symmetric solvers for pressure, asymmetric for velocity |

### Problem: Solver Taking Too Many Iterations

```
  Symptom: "No Iterations 1000" — solver hits maxIter without converging
```

**Possible causes and fixes**:

| Cause | Fix |
|-------|-----|
| tolerance too tight | Relax tolerance (1e-8 → 1e-6) |
| Poor preconditioner | Try stronger preconditioner (DIC → GAMG) |
| Wrong solver for equation | Switch e.g., PCG → GAMG for large pressure |
| Bad mesh quality | Improve mesh, check non-orthogonality |
| Physics not well-posed | Review boundary conditions and models |

---

## 14. Quick Reference Tables

### Complete Solver Configuration Reference

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `solver` | word | — | Solver name (PCG, PBiCGStab, smoothSolver, GAMG) |
| `preconditioner` | word | — | Preconditioner (DIC, DILU, FDIC, diagonal, GAMG, none) |
| `smoother` | word | — | Smoother for smoothSolver/GAMG |
| `tolerance` | scalar | 0 | Absolute convergence criterion |
| `relTol` | scalar | 0 | Relative convergence criterion |
| `minIter` | int | 0 | Minimum number of iterations |
| `maxIter` | int | 1000 | Maximum number of iterations |
| `nCellsInCoarsestLevel` | int | 10 | GAMG: cells on coarsest level |
| `agglomerator` | word | `faceAreaPair` | GAMG: coarsening method |
| `cacheAgglomeration` | bool | `true` | GAMG: cache coarse levels |
| `nPreSweeps` | int | 0 | GAMG: pre-smoothing sweeps |
| `nPostSweeps` | int | 2 | GAMG: post-smoothing sweeps |
| `nFinestSweeps` | int | 2 | GAMG: finest-level sweeps |
| `mergeLevels` | int | 1 | GAMG: coarsening aggressiveness |

### Common Solver Recipes

**Recipe 1: Small transient case (icoFoam/pisoFoam)**

```c
solvers
{
    p     { solver PCG;          preconditioner DIC;            tolerance 1e-06; relTol 0.05; }
    pFinal{ $p;                                                                  relTol 0;    }
    U     { solver smoothSolver; smoother       symGaussSeidel; tolerance 1e-05; relTol 0;    }
}
PISO { nCorrectors 2; nNonOrthogonalCorrectors 0; }
```

**Recipe 2: Large transient case (pimpleFoam)**

```c
solvers
{
    p     { solver GAMG;         smoother GaussSeidel; tolerance 1e-06; relTol 0.01; }
    pFinal{ $p;                                                         relTol 0;    }
    U     { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-05; relTol 0;  }
}
PIMPLE { nOuterCorrectors 2; nCorrectors 1; nNonOrthogonalCorrectors 1; }
```

**Recipe 3: Steady-state RANS (simpleFoam)**

```c
solvers
{
    p       { solver GAMG;         smoother GaussSeidel;    tolerance 1e-06; relTol 0.1;  }
    U       { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-05; relTol 0.1;  }
    k       { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-05; relTol 0.1;  }
    epsilon { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-05; relTol 0.1;  }
}
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; }
relaxationFactors
{
    fields    { p 0.3; }
    equations { U 0.7; k 0.7; epsilon 0.7; }
}
```

### Solver-Preconditioner-Smoother Pairing Rules

```
  ┌─────────────┬───────────────────┬──────────────────────┐
  │  Solver     │  Preconditioner   │  Smoother            │
  ├─────────────┼───────────────────┼──────────────────────┤
  │  PCG        │  DIC, FDIC,       │  (N/A — uses         │
  │             │  diagonal, GAMG   │   preconditioner)    │
  ├─────────────┼───────────────────┼──────────────────────┤
  │  PBiCGStab  │  DILU, diagonal   │  (N/A — uses         │
  │             │                   │   preconditioner)    │
  ├─────────────┼───────────────────┼──────────────────────┤
  │  smoothSolver│ (N/A — uses      │  GaussSeidel,        │
  │             │  smoother)        │  symGaussSeidel,     │
  │             │                   │  DIC, DILU           │
  ├─────────────┼───────────────────┼──────────────────────┤
  │  GAMG       │  (N/A — uses      │  GaussSeidel,        │
  │             │  smoother)        │  symGaussSeidel,     │
  │             │                   │  DICGaussSeidel      │
  └─────────────┴───────────────────┴──────────────────────┘

  PCG / PBiCGStab  →  need "preconditioner" keyword
  smoothSolver / GAMG  →  need "smoother" keyword
  Don't mix them up!
```

---

*Last updated based on project configurations from `01_lid_driven_cavity` and
`04_naca_airfoil_analysis`. For more on mesh generation, see `03_meshGeneration.md`.
For numerical schemes, see `08_numerical_schemes.md`.*
