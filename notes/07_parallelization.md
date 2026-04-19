# Parallelization in OpenFOAM

## Introduction: Why Parallelize?

Computational fluid dynamics simulations can require enormous amounts of processing time. A mesh with millions of cells solving transient turbulent flow may take days or even weeks on a single processor core. **Parallelization** splits the work across multiple processor cores, dramatically reducing wall-clock time.

OpenFOAM uses **MPI (Message Passing Interface)** for parallel communication. Each processor works on its own portion of the mesh and exchanges data with neighboring processors at shared boundaries. This is the **domain decomposition** approach — the single most important concept in OpenFOAM parallelization.

### Amdahl's Law — The Fundamental Limit

Not all parts of a simulation can be parallelized. Amdahl's law states that the maximum speedup is limited by the serial (non-parallelizable) fraction of the code:

$$
\text{Speedup} = \frac{1}{(1 - P) + \frac{P}{N}}
$$

Where:

- $P$ = fraction of code that is parallelizable  
- $N$ = number of processors

Even with infinite processors, if 5% of your code is serial, maximum speedup is only 20×. This is why efficient decomposition and minimizing communication overhead are critical.

## Domain Decomposition Concept

Domain decomposition divides the computational mesh into **N subdomains**, one per processor. Each processor solves the governing equations only on its subdomain, then communicates boundary data with its neighbors.

```
  Original Domain (1 processor)              Decomposed Domain (4 processors)
  ┌─────────────────────────────┐            ┌──────────────┬──────────────┐
  │                             │            │              │              │
  │                             │            │   Proc 0     │   Proc 1     │
  │                             │            │  (cells:     │  (cells:     │
  │       Full Mesh             │   ──────→  │   0..2499)   │   2500..4999 │
  │       10,000 cells          │            ├──────────────┼──────────────┤
  │                             │            │              │              │
  │                             │            │   Proc 2     │   Proc 3     │
  │                             │            │  (cells:     │  (cells:     │
  │                             │            │  5000..7499) │  7500..9999  │
  └─────────────────────────────┘            └──────────────┴──────────────┘
```

### How Data Is Split

When `decomposePar` runs, it creates **N directories** (`processor0/`, `processor1/`, ...) each containing a complete OpenFOAM case with:

- Its portion of the mesh (in `constant/polyMesh/`)
- Its portion of the field data (in `0/` and subsequent time directories)
- Processor boundary patch definitions for inter-processor communication

### Processor Boundaries — Ghost Cells and Halo Exchange

At the interfaces between subdomains, OpenFOAM creates special **processor patches**. These patches have **ghost cells** (also called halo cells) that store copies of values from the neighboring processor's cells. Before each solver iteration, processors exchange these values via MPI.

```
         Processor 0                    Processor 1
  ┌────────────────────────┐     ┌────────────────────────┐
  │                        │     │                        │
  │  ┌────┬────┬────┬─ ─ ─ │─ ─ ─│─ ─  ┬────┬────┬────┐   │
  │  │    │    │    │ghost │     │ghost│    │    │    │   │
  │  │ c1 │ c2 │ c3 │cell  │◄───►│cell │ c4 │ c5 │ c6 │   │
  │  │    │    │    │(c4') │     │(c3')│    │    │    │   │
  │  └────┴────┴────┴─ ─ ─ │─ ─ ─│─ ─  ┴────┴────┴────┘   │
  │                        │     │                        │
  └────────────────────────┘     └────────────────────────┘
              ▲          MPI Exchange           ▲
              │    (boundary values synced      │
              │     every iteration)            │
              └─────────────────────────────────┘
```

The ghost cells ensure that gradient calculations and flux computations at processor boundaries remain accurate, as if the mesh were still a single domain.

## Decomposition Methods

OpenFOAM supports several decomposition methods, each with different trade-offs. The method is specified in the `system/decomposeParDict` file.

### Simple — Geometric Slicing Along Axes

The **simple** method divides the domain by splitting it geometrically along the x, y, and z axes according to specified counts. It is fast but produces poor load balancing for non-uniform meshes.

```
  Simple decomposition: n = (4, 1, 1) — split into 4 along X-axis

  ┌────────┬────────┬────────┬────────┐
  │        │        │        │        │
  │  P0    │  P1    │  P2    │  P3    │
  │        │        │        │        │
  └────────┴────────┴────────┴────────┘
       X ──────────────────────→

  Simple decomposition: n = (2, 2, 1) — 2×2 grid

  ┌──────────────┬──────────────┐
  │              │              │
  │     P0       │     P1       │
  │              │              │
  ├──────────────┼──────────────┤
  │              │              │
  │     P2       │     P3       │
  │              │              │
  └──────────────┴──────────────┘
```

### Hierarchical — Multi-Level Splitting

The **hierarchical** method is similar to simple but applies the splits in a specified order (e.g., first along X, then Y, then Z). This gives slightly more control over the decomposition shape and is useful when the domain has a preferred direction.

```
  Hierarchical decomposition: n = (2, 2, 1), order = xyz

  Step 1: Split along X         Step 2: Split each half along Y
  ┌────────────┬────────────┐   ┌────────────┬────────────┐
  │            │            │   │    P0      │    P1      │
  │   Left     │   Right    │   ├────────────┼────────────┤
  │            │            │   │    P2      │    P3      │
  └────────────┴────────────┘   └────────────┴────────────┘
```

### Scotch — Graph-Based Optimal Partitioning (Recommended)

The **scotch** method uses the PT-Scotch graph partitioning library to minimize the number of inter-processor cell faces (communication volume) while keeping the cell count balanced across processors. It works well for complex geometries and non-uniform meshes, and is generally the **recommended method** for production runs.

```
  Scotch decomposition of an irregular mesh:

  ┌────────────────────────────────────────┐
  │  ┌─────────┐                           │
  │  │  P0     │    ┌──────────────┐       │
  │  │  (dense │    │              │       │
  │  │  region)│    │    P1        │       │
  │  └─────────┘    │              │       │
  │                 └──────────────┘       │
  │  ┌──────────────────┐  ┌────────────┐  │
  │  │                  │  │            │  │
  │  │      P2          │  │    P3      │  │
  │  │                  │  │            │  │
  │  └──────────────────┘  └────────────┘  │
  └────────────────────────────────────────┘

  Note: Scotch adapts partition shapes to mesh density,
  keeping cell counts balanced even with local refinement.
```

### Metis — Similar to Scotch

The **metis** method uses the METIS graph partitioning library. It produces similar results to scotch and can sometimes be faster or produce slightly different partitions. Requires METIS to be installed.

### Comparison Table

| Method       | Partition Quality | Speed     | Load Balancing | Best Used For                              |
|--------------|-------------------|-----------|----------------|--------------------------------------------|
| simple       | Low               | Very fast | Poor           | Quick tests, uniform hex meshes            |
| hierarchical | Low–Medium        | Fast      | Fair           | Structured meshes with a preferred axis    |
| scotch       | High              | Moderate  | Excellent      | **Production runs**, complex geometries    |
| metis        | High              | Moderate  | Excellent      | Alternative to scotch, large-scale runs    |

## The `decomposeParDict` — Deep Dive

The `system/decomposeParDict` controls how `decomposePar` splits the mesh. Below is a complete, annotated example.

### Complete Example: Scotch Method (Recommended)

```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Website:  https://openfoam.org                  |
|   \\  /    A nd           | Version:  6                                     |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains  4;

method              scotch;

// Optional: scotch-specific coefficients
scotchCoeffs
{
    // processorWeights — assign relative weights to processors
    // (useful if processors have different speeds)
    // processorWeights ( 1 1 1 1 );

    // strategy — scotch partitioning strategy string (advanced)
    // strategy "b{job=2,map=T,poli=S,sep=(/((|h{pass=10})x(|h{pass=10}))S)}";
}

// Distributed data settings (for multi-node clusters)
distributed     no;
roots           ();

// ************************************************************************* //
```

### Simple Method Configuration

```
numberOfSubdomains  4;

method              simple;

simpleCoeffs
{
    // Number of splits along each axis: (nX nY nZ)
    // Product must equal numberOfSubdomains: 4 × 1 × 1 = 4
    n               (4 1 1);

    // Cell skew factor — perturbation to help with load balancing
    // 0.001 is the default; larger values add more randomness
    delta           0.001;
}
```

### Hierarchical Method Configuration

```
numberOfSubdomains  8;

method              hierarchical;

hierarchicalCoeffs
{
    // Number of splits along each axis: (nX nY nZ)
    // Product must equal numberOfSubdomains: 2 × 2 × 2 = 8
    n               (2 2 2);

    // Cell skew factor
    delta           0.001;

    // Order in which axes are split: xyz, xzy, yxz, yzx, zxy, zyx
    order           xyz;
}
```

### Keyword Reference

| Keyword              | Required | Description                                                      |
|----------------------|----------|------------------------------------------------------------------|
| `numberOfSubdomains` | Yes      | Total number of processor subdomains (must match `-np` in mpirun)|
| `method`             | Yes      | Decomposition method: `simple`, `hierarchical`, `scotch`, `metis`|
| `simpleCoeffs`       | If simple| Sub-dictionary with `n` (split counts) and `delta`              |
| `hierarchicalCoeffs` | If hier. | Sub-dictionary with `n`, `delta`, and `order`                    |
| `scotchCoeffs`       | Optional | Scotch-specific settings (processor weights, strategy)           |
| `distributed`        | No       | `yes` if case data is distributed across network nodes           |
| `roots`              | No       | Root paths for distributed data on each node                     |

## The Complete Parallel Workflow

The parallel workflow in OpenFOAM follows a strict sequence of steps:

```
  ┌──────────────┐     ┌───────────────┐     ┌──────────────────┐     ┌─────────────────┐
  │   1. Mesh    │     │ 2. decompose  │     │  3. mpirun       │     │ 4. reconstruct  │
  │   Generation │────→│    Par        │────→│     solver       │────→│    Par          │
  │              │     │               │     │     -parallel    │     │                 │
  │  blockMesh   │     │ Splits mesh   │     │ Solves on each   │     │ Merges results  │
  │  snappyHex   │     │ into N parts  │     │ processor        │     │ into one        │
  └──────────────┘     └───────────────┘     └──────────────────┘     └─────────────────┘
         │                    │                      │                        │
         ▼                    ▼                      ▼                        ▼
    constant/            processor0/            processor0/             Unified time
    polyMesh/            processor1/            1/ 2/ 3/ ...           directories:
    (single mesh)        processor2/            (time dirs in          1/ 2/ 3/ ...
                         processor3/            each procDir)
```

### Step 1: Mesh Generation

Generate the mesh as usual with `blockMesh`, `snappyHexMesh`, or by importing from external tools (e.g., `fluentMeshToFoam` as used in the elbow tutorial). The mesh exists as a single domain at this point.

### Step 2: Domain Decomposition

Run `decomposePar` to split the mesh:

```bash
decomposePar
```

This creates the processor directories:

```
  case/
  ├── 0/                          ← Original boundary/initial conditions
  ├── constant/
  │   └── polyMesh/               ← Original full mesh
  ├── system/
  │   ├── controlDict
  │   ├── fvSchemes
  │   ├── fvSolution
  │   └── decomposeParDict        ← Decomposition configuration
  ├── processor0/                 ← Subdomain for processor 0
  │   ├── 0/                      ← Field data (portion)
  │   └── constant/
  │       └── polyMesh/           ← Mesh (portion + processor patches)
  ├── processor1/
  │   ├── 0/
  │   └── constant/
  │       └── polyMesh/
  ├── processor2/
  │   └── ...
  └── processor3/
      └── ...
```

### Step 3: Run Solver in Parallel

```bash
mpirun -np 4 icoFoam -parallel
```

Each processor writes its own time directories inside its `processorN/` folder. During the run, all processors communicate via MPI to exchange boundary data at processor interfaces.

### Step 4: Reconstruct Results

```bash
reconstructPar
```

This merges the field data from all `processorN/` directories back into the top-level case directory, producing unified time directories for post-processing.

## Running in Parallel — Detailed Guide

### Basic `mpirun` Syntax

```bash
mpirun -np <N> <solver> -parallel
```

- `-np <N>` — number of processors (must match `numberOfSubdomains`)
- `<solver>` — the OpenFOAM solver executable (e.g., `icoFoam`, `simpleFoam`, `pimpleFoam`)
- `-parallel` — tells the solver to run in decomposed mode

### Common `mpirun` Flags

| Flag               | Description                                                        |
|--------------------|--------------------------------------------------------------------|
| `-np 4`            | Run on 4 processors                                               |
| `--hostfile hosts` | Specify a file listing available nodes for cluster computing       |
| `--oversubscribe`  | Allow more processes than physical cores (useful for testing)      |
| `--bind-to core`   | Bind each process to a specific CPU core (better performance)      |
| `--map-by socket`  | Distribute processes across CPU sockets before filling one         |

### Running on Multiple Nodes (Cluster Computing)

For cluster jobs, create a **hostfile** listing available machines:

```
# hostfile
node01 slots=8
node02 slots=8
node03 slots=8
```

Then run:

```bash
mpirun -np 24 --hostfile hostfile simpleFoam -parallel
```

Ensure the case directory is accessible on all nodes (shared filesystem like NFS or Lustre), and set `distributed no;` in `decomposeParDict` when using a shared filesystem.

### Example: Parallel Allrun Script

Based on the elbow tutorial's `allrun` script structure, here is how you would adapt a case for parallel execution:

```bash
#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Get the solver name from system/controlDict
application=$(getApplication)

# Step 1: Generate or import the mesh
runApplication fluentMeshToFoam elbow.msh

# Step 2: Decompose the domain
runApplication decomposePar

# Step 3: Run the solver in parallel
runParallel $application

# Step 4: Reconstruct the results
runApplication reconstructPar

#------------------------------------------------------------------------------
```

The key difference from a serial script is the addition of `decomposePar`, the use of `runParallel` instead of `runApplication` for the solver, and `reconstructPar` at the end.

### Monitoring Parallel Runs

During execution, each processor writes its own log. You can monitor progress with:

```bash
# Watch the log file in real time
tail -f log.icoFoam

# Check residuals across all processors
grep "Solving for Ux" log.icoFoam | tail -5

# Monitor system resources
htop    # or top — check CPU usage across cores
```

## Reconstruction — Options and Strategies

### Full Reconstruction

Reconstruct all time steps from all processors:

```bash
reconstructPar
```

### Reconstruct Specific Time Steps

```bash
# Reconstruct only time 0.5
reconstructPar -time 0.5

# Reconstruct a range of times
reconstructPar -time 0.1:0.5

# Reconstruct only the latest time step
reconstructPar -latestTime
```

### When NOT to Reconstruct

For very large cases, reconstruction can take significant time and disk space. Alternatives:

**`paraFoam -builtin`**

ParaView can read decomposed data directly without reconstruction.

Launch with:

```bash
paraFoam -builtin -case .
```

Then select the `processor0/..processorN/` directories as a decomposed case.

**Post-process in parallel**

Many OpenFOAM post-processing utilities support `-parallel`:

```bash
mpirun -np 4 postProcess -func wallShearStress -parallel
```

### Reconstruction Reference

| Command                              | Description                              |
|--------------------------------------|------------------------------------------|
| `reconstructPar`                     | Reconstruct all time steps               |
| `reconstructPar -time 0.5`           | Reconstruct only time 0.5               |
| `reconstructPar -time 0.1:0.5`      | Reconstruct time range 0.1 to 0.5       |
| `reconstructPar -latestTime`         | Reconstruct the latest time step only    |
| `reconstructPar -newTimes`           | Reconstruct only times not yet done      |
| `reconstructPar -fields '(U p)'`    | Reconstruct only velocity and pressure   |

## Performance Optimization

### Optimal Number of Processors

The ideal processor count depends on the mesh size. A common rule of thumb:

```
  Recommended: 50,000 – 200,000 cells per processor

  Example:
    1,000,000 cell mesh  →  5 to 20 processors
    5,000,000 cell mesh  →  25 to 100 processors
   50,000,000 cell mesh  →  250 to 1000 processors
```

Going below 50k cells/processor usually degrades performance because communication overhead dominates. Going above 200k cells/processor leaves performance on the table.

### Speedup and Scaling

```text
  Speedup
    │
    │                                   ╱  Ideal (linear speedup: S = N)
    │                                 ╱
    │                               ╱
    │                             ╱
    │                           ╱
    │                        ──╯      Actual (good decomposition)
    │                     ──╯
    │                  ──╯
    │               ──╯
    │            ──╯              Actual (poor decomposition / too few cells)
    │         ──╯        ─────────╯
    │      ──╯    ───────╯
    │   ──╯  ─────╯
    │──╯──────────────────────────────────── Diminishing returns
    └─────────────────────────────────────── Number of Processors
      1   2   4   8   16   32   64
```

In practice, you rarely achieve linear speedup. The gap between ideal and actual comes from:

1. **Communication overhead** — processors exchange data at every iteration
2. **Load imbalance** — some processors finish earlier and wait for others
3. **Serial fractions** — file I/O, initialization, and reconstruction are serial
4. **Memory bandwidth** — multiple cores sharing the same memory bus

### Amdahl's Law vs Gustafson's Law

**Amdahl's law** says: for a fixed problem size, speedup is limited by the serial fraction. If 10% of the work is serial, max speedup is 10× no matter how many  processors you add.

**Gustafson's law** says: as you add processors, you can increase the problem size proportionally, and the serial fraction becomes negligible. This is the more practical view for CFD — you typically run bigger meshes when you have more processors available.

### Load Balancing Considerations

Poor load balancing means some processors do more work than others. The slowest processor
determines the overall speed.

```
  Good balance (scotch):          Poor balance (simple on non-uniform mesh):
  ┌──────┬──────┬──────┬──────┐   ┌──────┬──────┬──────┬──────┐
  │ 2501 │ 2500 │ 2499 │ 2500 │   │ 1200 │ 1800 │ 3500 │ 3500 │
  │cells │cells │cells │cells │   │cells │cells │cells │cells │
  └──────┴──────┴──────┴──────┘   └──────┴──────┴──────┴──────┘
  Max work: 2501 → balanced       Max work: 3500 → P0 idle 65% of time!
```

Use `scotch` or `metis` for best load balancing on non-uniform meshes.

### Network vs Computation

For **shared-memory** systems (single workstation, multi-core), communication is fast through shared memory. For **distributed-memory** clusters (multiple nodes connected by a network), communication speed is limited by network bandwidth and latency.

| System Type              | Communication | Ideal For                                |
|--------------------------|---------------|------------------------------------------|
| Multi-core workstation   | Shared memory | Small to medium cases (< 10M cells)     |
| Cluster (Ethernet)       | ~1 Gbps       | Medium cases, limited scaling            |
| Cluster (InfiniBand)     | ~100 Gbps     | Large-scale production, good scaling     |

## Practical Tips

### When NOT to Parallelize

- **Small meshes** (< 50,000 cells): Overhead from decomposition, communication, and reconstruction often exceeds the time saved. Run serial instead.
- **Very short runs**: If your simulation finishes in seconds, the startup and teardown time of MPI will dominate.
- **Debugging**: Always debug your case setup on a serial run first. Parallel error messages are harder to read and often less informative.

### Debugging Parallel Runs

If a parallel run crashes:

I. **Check all processor logs** — the error may only appear in one processor's output:

```bash
grep -r "FOAM FATAL" processor*/
grep -i "error" log.simpleFoam
```
II. **Run on 1 processor in parallel mode** to isolate MPI issues from solver issues:

```bash
mpirun -np 1 icoFoam -parallel
```

III. **Check decomposition quality** — visualize the decomposition in ParaView to ensure all processors have reasonable cell counts and shapes.

### Common Errors and Solutions

| Error                                          | Cause                                    | Solution                                    |
|------------------------------------------------|------------------------------------------|---------------------------------------------|
| `numberOfSubdomains != nProcs`                | Mismatch between dict and `-np`          | Ensure `-np N` matches `numberOfSubdomains` |
| `Could not find processor* directories`        | Forgot to run `decomposePar`             | Run `decomposePar` before `mpirun`          |
| `Floating point exception` on some processors | Load imbalance or bad mesh partition     | Switch to `scotch` method                   |
| `MPI_ABORT was invoked`                        | One processor crashed                    | Check individual processor logs             |
| `bus error` or segfault                        | Out of memory on a node                  | Use fewer procs per node or more memory     |
| Reconstruction fails or produces empty fields  | Run was interrupted before writing       | Use `reconstructPar -latestTime`            |

### MPI Environment Setup

Ensure your MPI environment is correctly configured:

```bash
# Check which MPI is available
which mpirun
mpirun --version

# For OpenFOAM, source the environment first
source /opt/openfoam6/etc/bashrc    # adjust path for your installation

# Verify OpenFOAM parallel tools are available
which decomposePar
which reconstructPar
```

## Reference: Common Parallel Commands

| Command                                          | Description                                        |
|--------------------------------------------------|----------------------------------------------------|
| `decomposePar`                                   | Decompose mesh and fields into processor dirs      |
| `decomposePar -force`                            | Overwrite existing decomposition                   |
| `decomposePar -copyZero`                         | Copy `0/` directory instead of decomposing fields  |
| `mpirun -np 4 icoFoam -parallel`                | Run icoFoam on 4 processors                       |
| `mpirun -np 8 simpleFoam -parallel`             | Run simpleFoam on 8 processors                    |
| `mpirun -np 4 postProcess -func yPlus -parallel`| Run post-processing in parallel                    |
| `reconstructPar`                                 | Reconstruct all time directories                   |
| `reconstructPar -latestTime`                     | Reconstruct only the latest time step              |
| `reconstructPar -time 0.1:0.5`                  | Reconstruct a specific time range                  |
| `reconstructPar -fields '(U p)'`                | Reconstruct only specific fields                   |
| `reconstructPar -newTimes`                       | Reconstruct only new (not yet reconstructed) times |
| `paraFoam -builtin`                              | View decomposed results without reconstruction     |

## Putting It All Together — Example for the Lid-Driven Cavity

Applying parallelization to the `01_lid_driven_cavity` project (which uses `icoFoam`):

**1. Create `system/decomposeParDict`:**

```
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}

numberOfSubdomains  4;
method              scotch;
```

**2. Run the parallel workflow:**

```bash
# Generate mesh
blockMesh

# Decompose
decomposePar

# Solve in parallel
mpirun -np 4 icoFoam -parallel

# Reconstruct
reconstructPar

# Visualize
paraFoam
```

**3. Clean up processor directories (optional):**

```bash
# Remove processor directories after reconstruction to save disk space
rm -rf processor*
```

This workflow applies identically to any OpenFOAM case — just change the solver name and adjust `numberOfSubdomains` based on your mesh size and available cores.
