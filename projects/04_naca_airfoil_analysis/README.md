# NACA 0012 Airfoil Analysis

This project contains two OpenFOAM airfoil cases built around a NACA 0012 profile:

- a default steady RANS case in the project root
- a separate transient animation case in [`transient_animation/`](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis/transient_animation)

The current default case is intended to be runnable and visually useful. It is not yet a
validation-grade reproduction of the NASA NACA 0012 benchmark.

## What Is In This Project

### Default case

- Solver: `simpleFoam`
- Flow model: incompressible, steady, RANS
- Turbulence model: `SpalartAllmaras`
- Mesh path: Gmsh `.geo` -> `gmshToFoam` -> OpenFOAM
- Geometry setup: airfoil rotated by about `4°` while freestream remains horizontal
- Intended output: attached steady field for pressure / velocity visualization

### Transient case

- Directory: [`transient_animation/`](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis/transient_animation)
- Solver: `pimpleFoam`
- Intended output: many saved time steps for animation and unsteady visualization
- This is the case to use if you want actual time evolution instead of SIMPLE iterations

## Current Status

The root case is stable and runs end to end through Docker.

From the latest completed steady run:

- `Cl ≈ 0.426`
- `Cd ≈ 0.044`

Results file:

- [`postProcessing/forceCoeffs/0/forceCoeffs.dat`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/postProcessing/forceCoeffs/0/forceCoeffs.dat)

Important limitation:

- the current wall treatment is still coarse for an airfoil case
- latest average `y+` is about `476`
- that is acceptable for a runnable tutorial-style wall-function case, but still too high
  for validation-quality drag prediction

`yPlus` output:

- [`postProcessing/yPlus/0/yPlus.dat`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/postProcessing/yPlus/0/yPlus.dat)

## Case Layout

```text
04_naca_airfoil_analysis/
├── 0/                              # active initial / boundary fields
├── 0.orig/                         # clean fields restored by Allrun
├── constant/
│   ├── transportProperties
│   ├── turbulenceProperties
│   ├── airfoil_benchmark.geo       # generated Gmsh geometry
│   ├── airfoil_benchmark.msh       # generated Gmsh mesh
│   └── triSurface/
├── mesh_generation_scripts/
│   ├── generate_naca0012_benchmark_geo.py
│   ├── generate_naca_0012_airfoil.py
│   └── process_naca0012.sh
├── system/
│   ├── changeDictionaryDict
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── transient_animation/            # separate transient case
├── Allclean
├── Allrun
├── run_docker.sh
├── open_paraview_streamlines.sh
├── paraview_airfoil.py
└── paraview_guide.md
```

## Physics And Setup

### Default steady case

- Freestream velocity: `U = 50 m/s`
- Kinematic viscosity: `nu = 8.333e-6 m^2/s`
- Chord length: `1 m`
- Reference span used for coefficients: `0.02 m`
- Approximate Reynolds number: `6e6`

Relevant files:

- [`constant/transportProperties`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/constant/transportProperties)
- [`constant/turbulenceProperties`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/constant/turbulenceProperties)
- [`system/controlDict`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/system/controlDict)

### Boundary conditions

Root case fields:

- [`0/U`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/U)
- [`0/p`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/p)
- [`0/nuTilda`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/nuTilda)
- [`0/nut`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/nut)

Patch types:

- `inlet`: fixed freestream velocity
- `outlet`: fixed pressure
- `top`, `bottom`: `symmetryPlane`
- `airfoil`: `wall`
- `frontAndBack`: `empty`

## Mesh Workflow

The root case no longer uses the older `blockMesh + snappyHexMesh` path as its main run
path.

Current root workflow:

1. [`Allrun`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/Allrun) restores clean fields from [`0.orig/`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0.orig)
2. [`mesh_generation_scripts/generate_naca0012_benchmark_geo.py`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/mesh_generation_scripts/generate_naca0012_benchmark_geo.py)
   writes a Gmsh geometry file
3. Gmsh generates `constant/airfoil_benchmark.msh`
4. `gmshToFoam` imports the mesh
5. `changeDictionary` fixes patch types
6. `renumberMesh`, `checkMesh`, `simpleFoam` run in sequence

The currently stable mesh is a practical compromise:

- it runs reliably
- it gives a useful attached-flow field
- it is still not fine enough near the wall for benchmark drag accuracy

## How To Run

### Recommended: Docker

Run the steady case from the project directory:

```bash
cd projects/04_naca_airfoil_analysis
./run_docker.sh
```

This script:

1. cleans the case through the OpenFOAM container
2. generates the Gmsh geometry on the host
3. generates the `.msh` mesh with host `gmsh`
4. runs the OpenFOAM workflow in Docker

### Native OpenFOAM shell

If OpenFOAM is installed locally and sourced:

```bash
cd projects/04_naca_airfoil_analysis
./Allclean
./Allrun
```

### Transient animation case

To run the transient case instead:

```bash
cd projects/04_naca_airfoil_analysis/transient_animation
./Allrun
```

That case uses `pimpleFoam` and writes many time steps for animation.

## ParaView

### One-command launcher

From the project directory:

```bash
./open_paraview_streamlines.sh
```

This opens ParaView with [`paraview_airfoil.py`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/paraview_airfoil.py) and builds the
default scene automatically.

Current default ParaView view:

- mid-plane slice
- `Surface LIC`
- colored by `U Magnitude`
- white airfoil surface
- tight 2D camera

### Useful launch options

The shell wrapper forwards extra arguments to the ParaView script:

```bash
./open_paraview_streamlines.sh --u-min 25 --u-max 65
./open_paraview_streamlines.sh --scalar p --pressure-limit 300
```

Use:

- `--scalar U` for velocity-magnitude style plots
- `--scalar p` for pressure plots

Manual notes are in [`paraview_guide.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/paraview_guide.md).

## Important Interpretation Notes

### The root case is steady

The root solver is `simpleFoam`, so the directories `10/`, `20/`, `100/`, `600/` are
SIMPLE iteration outputs, not physical seconds.

That means:

- the field becoming stationary is expected
- you will not get a true animation from the root case
- if you want real time evolution, use the transient case in
  [`transient_animation/`](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis/transient_animation)

### Why the internet images look stronger

Most reference-style airfoil images differ from this case in one or more of these ways:

- they are colored by velocity magnitude rather than pressure
- they use a transient or more strongly separated flow regime
- they use much finer boundary-layer resolution
- they use structured or higher-quality external-aero meshes

So ParaView can improve presentation, but it cannot create stronger physical effects than
the CFD solution actually contains.

## Known Limitations

- The steady root case is stable, but not validation-grade.
- Drag is still higher than expected for a high-quality NACA 0012 benchmark case.
- Near-wall resolution is still too coarse for trustworthy airfoil drag prediction.
- More aggressive boundary-layer meshes were tested, but in this repo’s current Gmsh
  workflow they either stalled during recombination or produced poor-quality tetrahedral
  meshes that destabilized the solver.

## Recommended Next Steps

If you want to improve physics rather than just visualization, the best next steps are:

1. move the steady case to a better airfoil mesh topology, ideally a structured or true
   C-grid style external-aero mesh
2. reduce wall-function dependence by improving first-layer control and `y+`
3. compare `Cl`, `Cd`, and surface pressure against NASA / TMR reference data
4. use the transient case for animation-oriented work rather than trying to animate the
   steady SIMPLE solution

## References

- NASA Turbulence Modeling Resource, NACA 0012 validation:
  <https://turbmodels.larc.nasa.gov/naca0012_val.html>
- NASA TMR, Spalart-Allmaras reference results:
  <https://turbmodels.larc.nasa.gov/naca0012_val_sa.html>
- NASA TMR, benchmark grids:
  <https://turbmodels.larc.nasa.gov/naca0012_grids.html>
- OpenFOAM `simpleFoam` documentation:
  <https://doc.openfoam.com/2306/tools/processing/solvers/rtm/incompressible/simpleFoam/>
- OpenFOAM `pimpleFoam` documentation:
  <https://doc.openfoam.com/2306/tools/processing/solvers/rtm/incompressible/pimpleFoam/>
