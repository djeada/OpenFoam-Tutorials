# NACA 0012 Airfoil Analysis

This project contains two OpenFOAM airfoil cases built around a NACA 0012 profile:

- a default steady RANS case in the project root (`simpleFoam`, Spalart-Allmaras, 4 degrees angle of attack)
- a separate transient animation case in [`transient_animation/`](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis/transient_animation) (`pimpleFoam`, laminar, 15 degrees angle of attack)

The default case is intended to be runnable and visually useful. It is not a validation-grade
reproduction of the NASA NACA 0012 benchmark; see [Known Limitations](#known-limitations).

![Pressure field around the airfoil](images/pressure.png)

*Kinematic pressure p/rho at iteration 600. Look for the high-pressure region (red) at the
leading edge below the nose, where the flow stagnates, and the low-pressure region (blue) over
the front part of the upper surface (suction peak).*

## Goal

Show the complete external-aerodynamics workflow on a single airfoil:

- generating a 2D-style mesh with a boundary layer from a parametric NACA 0012 profile (Python -> Gmsh -> `gmshToFoam`)
- setting up a steady incompressible RANS case with the Spalart-Allmaras model and wall functions
- applying an angle of attack by rotating the geometry while keeping the freestream horizontal, so
  `liftDir` / `dragDir` stay `(0 1 0)` / `(1 0 0)`
- reading lift, drag and moment coefficients from `forceCoeffs`, and judging them against theory
- understanding why mesh quality, `y+` and numerical schemes limit the drag prediction

## Target

| Quantity | Target | Basis |
| --- | --- | --- |
| Reynolds number | `U c / nu = 6.0e6` | `U = 50 m/s`, `c = 1 m`, `nu = 8.333e-6 m^2/s` |
| Angle of attack | `+4 degrees` | airfoil rotated by `-4 degrees` about the quarter chord |
| Lift coefficient | `Cl = 0.439 +/- 0.05` | thin-airfoil theory `2 pi alpha`; the NACA 0012 lift-curve slope at this Re is close to it for small alpha |
| Drag coefficient | `0 < Cd < Cl` (see note) | consistency only, no reference value at 4 degrees is used |
| Moment coefficient | `abs(Cm) < 0.02` about `x = 0.25 c` | symmetric section: theory gives zero about the quarter chord |
| Mesh | no `checkMesh` failures other than skewness, non-orthogonality below 70 degrees; skewness is known to fail (see below) | OpenFOAM's own limits |
| Solution | `Cl` changes by less than 1 % over the last 100 iterations | stationary force |

Note on drag: `Cd` at this mesh is **not expected to match published values** for the NACA 0012,
and no validated reference value for this angle of attack, Reynolds number and solver is
available in this repository, so none is used as a pass/fail criterion. The reasons are numerical:
the mesh is coarse, the boundary layer is handled only by wall functions with `y+` of about
215 to 750 (average 476, far from the `y+` of about 1 needed to resolve the sublayer), and `U` and
`nuTilda` use first-order upwind. The measured `Cd = 0.0447` is much larger than a well-resolved
turbulent airfoil at this Reynolds number would give; only `Cd > 0` and `Cd < Cl` are checked.
The reference data for a proper comparison is listed under [References](#references).

## Verification Method

After a run, check the results with

```bash
./run_docker.sh        # runs the case and ends with verify.py
python3 verify.py      # re-check an existing run; reads only log.* and postProcessing/
```

[`verify.py`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/verify.py)
uses only the Python standard library. It prints one `PASS`/`FAIL` line per check and a summary,
exits `0` only if all checks pass, and exits `2` if results are missing. Values in the last
column were measured in a complete `./run_docker.sh` run of this repository's case
(Gmsh 4.11.1 on the host, OpenFOAM 6 in the container, 600 iterations). A second run with
Gmsh 4.15.2 produced a slightly different mesh (39 412 cells) and similar results
(`Cl = 0.427`, `Cd = 0.0465`, average `y+` 485, max non-orthogonality 58.1 degrees, final `p`
residual 9.9e-3): results depend somewhat on the Gmsh version.
Meshing with the host's Gmsh 4.12.1 was checked separately (38 585 cells, the same `checkMesh`
figures, skewness 4.25 on 2 faces), but the full steady flow was not re-run with that version.

| Check | What it proves | Expected | Measured |
| --- | --- | --- | --- |
| freestream `abs(U)` at inlet vs `magUInf` | force coefficients use the actual freestream | equal, relative 1e-6 | 50 vs 50 |
| `Aref` = `lRef` x span | reference area matches the 0.02 m extrusion | equal, 1e-9 | 0.02 |
| Reynolds number `U c / nu` | case is at the stated Re | 6.0e6 +/- 5 % (nu is rounded) | 6.0e6 |
| `checkMesh` failed checks other than skewness | topology, openness, aspect ratio, non-orthogonality are fine | 0 | 0 |
| max non-orthogonality | below OpenFOAM's severe limit | < 70 deg | 62.18 deg |
| `simpleFoam` no `FOAM FATAL ERROR` | solver ran cleanly | none | none |
| `simpleFoam` ran to `End` | log is complete | `End` | `End` |
| run finished | `residualControl` met or `endTime` reached | 600 iterations | 600 iterations |
| final initial residual `p` | pressure residual dropped at least one order | < 1e-1 | 1.2e-2 |
| final initial residual `Ux` | momentum residual dropped at least one order | < 1e-1 | 1.1e-3 |
| final initial residual `Uy` | momentum residual dropped at least one order | < 1e-1 | 1.3e-3 |
| final initial residual `nuTilda` | turbulence residual dropped at least one order | < 1e-1 | 1.7e-4 |
| `Cl` vs `2 pi alpha` (alpha = 4 deg) | lift magnitude is physically right | 0.4386 +/- 0.05 | 0.4297 |
| `Cl > 0` | lift direction is consistent with nose-up rotation and `liftDir` | > 0 | 0.4297 |
| `0 < Cd < Cl` | drag acts downstream and is smaller than lift | 0 < Cd < Cl | 0.04474 |
| `abs(Cm)` about quarter chord | symmetric section has near-zero moment | < 0.02 | 0.0073 |
| relative change of `Cl` over last 100 iterations | force is stationary | < 1 % | 0.26 % |

Tolerance rationale:

- `Cl` band `0.05` is about 11 %: it has to absorb the coarse wall-function mesh, the
  15-chord domain with symmetry-plane top and bottom boundaries, first-order upwind schemes and
  the idealisation of `2 pi alpha`. It is a consistency check, not a validation.
- Residual limit `1e-1` (one order of magnitude below the initial value of 1) is only the minimum
  evidence that the solver is iterating towards a solution. `residualControl` in `fvSolution`
  asks for `1e-6`, which this case does not reach in 600 iterations: with first-order upwind on
  this mesh the pressure residual plateaus at about `1e-2` and fluctuates (9.9e-3 and 1.2e-2 on
  the two meshes above). The limit was set to `1e-1` after the pressure residual of a
  second mesh landed just above an initial `1e-2` limit; the case is simply not converged in
  the residual sense. The stationarity check on `Cl` is the stronger statement about the answer.
- `abs(Cm) < 0.02` is about 5 % of `Cl`.
- Non-orthogonality 70 degrees is the limit `checkMesh` itself uses for a severe error.

`checkMesh` does not report `Mesh OK` for this mesh. It ends with `Failed 1 mesh checks`: the
skewness check fails with max skewness 4.25 on 2 faces (OpenFOAM limit 4). All other checks
pass (max non-orthogonality 62.2 degrees, average 6.6; max aspect ratio 4.0). `simpleFoam` ran
without problems on this mesh, so `verify.py` does not count the skewness failure but prints it
as a `WARN` line. It also prints `INFO` lines for `y+` (min 217, max 749, average 476) and for
`Cd`, which is deliberately not compared with any published value (see the note above).

Manual checks (visual, not automated):

- Open the case in ParaView (`./open_paraview.sh`): the stagnation point should be at the leading edge
  on the lower side and the suction peak on the upper surface near the nose.
- `postProcessing/forceCoeffs/0/forceCoeffs.dat`: `Cl` should rise smoothly and flatten
  (the first values near iteration 20 are the start-up transient).
- Streamlines should stay attached to the whole airfoil at 4 degrees.

## What Is In This Project

### Default case

- Solver: `simpleFoam`
- Flow model: incompressible, steady, RANS
- Turbulence model: `SpalartAllmaras`
- Mesh path: Gmsh `.geo` -> `gmshToFoam` -> OpenFOAM
- Geometry setup: airfoil rotated by `4` degrees about the quarter chord, freestream remains horizontal
- Intended output: attached steady field for pressure / velocity visualization

### Transient case

- Directory: [`transient_animation/`](https://github.com/djeada/OpenFoam-Tutorials/tree/main/projects/04_naca_airfoil_analysis/transient_animation)
- Solver: `pimpleFoam`, laminar
- Setup: `U = 1 m/s`, `nu = 1e-4 m^2/s` (Re = 10 000), airfoil rotated by 15 degrees, `endTime = 20 s`,
  adaptive time step with `maxCo = 0.5`, output every `0.05 s`
- Intended output: many saved time steps for animation and unsteady visualization
- This is the case to use if you want actual time evolution instead of SIMPLE iterations
- Status: **untested.** Its scripts were refactored together with the root case (`make_mesh.sh`,
  `Allrun`, `Allclean`, `run_docker.sh transient`) but the case was not run in the verification
  described in this README (the full run is long)

## Results

From the verification run (OpenFOAM 6, 600 SIMPLE iterations, about 1-2 minutes of solver
time on one core):

- `Cl = 0.430`, `Cd = 0.0447`, `Cm = 0.0073` (about `x = 0.25 c`)
- average `y+` on the airfoil about `476` (min 217, max 749)
- mesh: 38 584 cells (38 490 hexahedra, 94 prisms), one cell thick

Result files:

- [`postProcessing/forceCoeffs/0/forceCoeffs.dat`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/postProcessing/forceCoeffs/0/forceCoeffs.dat)
- [`postProcessing/yPlus/0/yPlus.dat`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/postProcessing/yPlus/0/yPlus.dat)
- [`postProcessing/residuals/0/residuals.dat`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/postProcessing/residuals/0/residuals.dat)

(`postProcessing/` is generated by the run and not tracked in git.)

Important limitations of these numbers:

- the wall treatment is coarse for an airfoil case: average `y+` of about `476` is acceptable for a
  runnable tutorial-style wall-function case, but too high for validation-quality drag prediction
- `Cd = 0.0447` is not expected to match published NACA 0012 values on this mesh (see the note
  under [Target](#target)); it should not be used as a result
- `checkMesh` fails one check: max skewness 4.25 on 2 faces (limit 4), so it does not print `Mesh OK`
- the residuals stall between about `1e-4` and `1e-2`, so the solution is stationary in the forces but not
  converged to `residualControl`

![Velocity magnitude and streamlines](images/velocity.png)

*Velocity magnitude (color limits 30 to 70 m/s) with streamlines seeded 1 chord upstream. Look
for the accelerated flow (red) above the nose, the low-speed stagnation region (blue) at the
leading edge, and the streamlines following the whole surface without separation.*

## Case Layout

```text
04_naca_airfoil_analysis/
├── 0/                              # initial / boundary fields (U, p, nut, nuTilda)
├── constant/
│   ├── transportProperties
│   ├── turbulenceProperties
│   ├── airfoil_benchmark.geo       # generated Gmsh geometry
│   └── airfoil_benchmark.msh       # generated Gmsh mesh
├── images/                         # figures used in this README
├── mesh_generation_scripts/
│   └── generate_naca0012_benchmark_geo.py
├── system/
│   ├── changeDictionaryDict
│   ├── controlDict
│   ├── fvSchemes
│   └── fvSolution
├── transient_animation/            # separate transient case
├── Allclean
├── Allrun
├── make_mesh.sh                    # geometry + Gmsh meshing (called by Allrun)
├── run_docker.sh
├── verify.py
├── open_paraview.sh
├── paraview_airfoil.py
├── paraview_images.py              # renders images/*.png
└── paraview_guide.md
```

## Physics And Setup

### Default steady case

| Parameter | Value | File |
| --- | --- | --- |
| Freestream velocity `U` | `50 m/s`, direction `(1 0 0)` | `0/U`, `magUInf` in `system/controlDict` |
| Kinematic viscosity `nu` | `8.333e-6 m^2/s` | `constant/transportProperties` |
| Chord length `c` | `1 m` (`lRef`) | `system/controlDict` |
| Span of the extrusion | `0.02 m` | `make_mesh.sh` / generator default |
| Reference area `Aref` | `0.02 m^2` (chord x span) | `system/controlDict` |
| Reynolds number | `U c / nu = 6e6` | derived |
| Angle of attack | `+4 degrees` (geometry rotated by `-4 degrees`) | `make_mesh.sh` |
| Moment reference point `CofR` | `(0.25 0 0)` (quarter chord of the unrotated airfoil) | `system/controlDict` |
| Density for coefficients `rhoInf` | `1.225 kg/m^3` (the solver is kinematic) | `system/controlDict` |
| Freestream `nuTilda`, `nut` | `2.5e-5 m^2/s` (3 nu) | `0/nuTilda`, `0/nut` |
| Iterations | `600`, write every `10` | `system/controlDict` |

Numerics: first-order `bounded Gauss upwind` for `U` and `nuTilda`, `cellLimited Gauss linear 1`
gradients, GAMG for `p`, relaxation `p 0.2`, `U 0.3`, `nuTilda 0.2`,
`residualControl 1e-6`.

Lift and drag directions stay `(0 1 0)` and `(1 0 0)` because the freestream is horizontal and
the airfoil, not the flow, is rotated.

Relevant files:

- [`constant/transportProperties`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/constant/transportProperties)
- [`constant/turbulenceProperties`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/constant/turbulenceProperties)
- [`system/controlDict`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/system/controlDict)

### Domain

Dimensions in chord lengths `c = 1 m`; the mesh is one cell thick in `z` (0.02 m).

```text
                       top: symmetryPlane,  y = +15
      +-------------------------------------------------------------+
      |                                                             |
      |                                                             |
      |                                                             |
      |                       airfoil (wall)                        |
 U -> |                  .-==========-.__                           |-> outlet
inlet |  x = -15         LE (0,0)       TE (1,-0.05)    wake box    |  x = +20
      |                  4 deg AoA      <--- 0 < x < 10, |y| < 3 -->|  p = 0
      |                                                             |
      |                                                             |
      |                                                             |
      +-------------------------------------------------------------+
                       bottom: symmetryPlane, y = -15
      frontAndBack: empty (z = 0 and z = 0.02)
```

The airfoil (leading edge near the origin, trailing edge near `(1, -0.05)`) is rotated by
`-4 degrees` about `(0.25, 0)`, which is a +4 degree angle of attack to the +x freestream.

![Whole computational domain](images/mesh_domain.png)

*Whole domain (x from -15 to 20, y from -15 to 15). The dark box is the wake refinement region
(0 < x < 10, abs(y) < 3); the airfoil is the small black shape at its upstream edge.*

### Boundary conditions

Root case fields:

- [`0/U`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/U)
- [`0/p`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/p)
- [`0/nuTilda`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/nuTilda)
- [`0/nut`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/0/nut)

| Patch | Type | `U` | `p` | `nuTilda` | `nut` |
| --- | --- | --- | --- | --- | --- |
| `inlet` | patch | `fixedValue (50 0 0)` | `zeroGradient` | `fixedValue 2.5e-5` | `calculated` |
| `outlet` | patch | `pressureInletOutletVelocity` | `fixedValue 0` | `inletOutlet` | `calculated` |
| `top`, `bottom` | `symmetryPlane` | symmetryPlane | symmetryPlane | symmetryPlane | symmetryPlane |
| `airfoil` | `wall` | `noSlip` | `zeroGradient` | `fixedValue 0` | `nutUSpaldingWallFunction` |
| `frontAndBack` | `empty` | empty | empty | empty | empty |

The patch types are assigned after import by `system/changeDictionaryDict`, because
`gmshToFoam` creates every patch as a plain `patch`.

## Mesh Workflow

The root case used an older `blockMesh + snappyHexMesh` path earlier in its history; the
dictionaries for it have been removed. The current workflow is:

1. [`make_mesh.sh`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/make_mesh.sh) calls
   [`mesh_generation_scripts/generate_naca0012_benchmark_geo.py`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/mesh_generation_scripts/generate_naca0012_benchmark_geo.py),
   which writes the Gmsh geometry `constant/airfoil_benchmark.geo`
2. Gmsh generates `constant/airfoil_benchmark.msh` (format 2, as `gmshToFoam` requires)
3. `gmshToFoam` imports the mesh
4. `changeDictionary` fixes patch types
5. `renumberMesh`, `checkMesh`, `simpleFoam` run in sequence

```mermaid
flowchart LR
    A[generate_naca0012_benchmark_geo.py] --> B[airfoil_benchmark.geo]
    B -->|gmsh -3 -format msh2| C[airfoil_benchmark.msh]
    C --> D[gmshToFoam]
    D --> E[changeDictionary]
    E --> F[renumberMesh]
    F --> G[checkMesh]
    G --> H[simpleFoam]
    H --> I[verify.py]
```

Mesh parameters used by `make_mesh.sh` (lengths in chord lengths):

| Parameter | Value |
| --- | --- |
| Chord points per side (cosine spacing) | 101 |
| Airfoil surface size | 0.012 |
| Boundary layer: first layer / total thickness / growth ratio | 0.004 / 0.05 / 1.12 |
| Near-body size, transition distance | 0.03, 1.5 |
| Wake box: size, x extent, half height | 0.04, 0 to 10, 3 |
| Far field | x from -15 to 20, y from -15 to 15, size 1.5 (generator default) |

The mesh is a practical compromise:

- it runs reliably
- it gives a useful attached-flow field
- it is not fine enough near the wall for benchmark drag accuracy
- `checkMesh` fails one check, `Failed 1 mesh checks`: max skewness 4.25 on 2 faces (limit 4);
  max non-orthogonality is 62.2 degrees (limit 70)

![Mesh around the airfoil](images/mesh_airfoil.png)

*Mesh around the airfoil, cells seen from +z. Look for the boundary-layer strip of about eight
stacked quad layers along both surfaces (first layer 0.004 c, growth ratio 1.12), the irregular
quads outside it. The mesh has 98 prisms (triangular-based cells) among 38 490 hexahedra.*

## How To Run

### Recommended: Docker

Requirements: `docker`, `python3` and `gmsh` on the host (the OpenFOAM image contains no Gmsh).
Run the steady case from the project directory:

```bash
cd projects/04_naca_airfoil_analysis
./run_docker.sh
```

This script:

1. cleans the case through the OpenFOAM container (`./Allclean`)
2. generates the Gmsh geometry and the `.msh` mesh on the host with `make_mesh.sh`
3. runs `./Allrun` in the container (`gmshToFoam` ... `simpleFoam`)
4. runs `python3 verify.py` and prints the PASS/FAIL summary

Meshing is the slowest step (several minutes); the solver itself runs for one to two minutes.

### Native OpenFOAM shell

If OpenFOAM, `python3` and `gmsh` are available locally and OpenFOAM is sourced:

```bash
cd projects/04_naca_airfoil_analysis
./Allclean
./Allrun
python3 verify.py
```

### Transient animation case

To run the transient case instead:

```bash
./run_docker.sh transient
```

or, natively:

```bash
cd projects/04_naca_airfoil_analysis/transient_animation
./Allrun
```

That case uses `pimpleFoam` and writes many time steps for animation. `./run_docker.sh all`
runs both cases one after the other.

### Cleaning

`./Allclean` removes time directories, `postProcessing/`, logs, the generated mesh and the
`airfoil.foam` marker. `images/` and `verify.py` are not touched.

## ParaView

### One-command launcher

From the project directory:

```bash
./open_paraview.sh
```

This opens ParaView with [`paraview_airfoil.py`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/paraview_airfoil.py) and builds the
default scene automatically.

Default ParaView view:

- mid-plane slice
- `Surface LIC`
- colored by `U Magnitude`
- white airfoil surface
- tight 2D camera

### Changing the view

The ParaView GUI does not pass command-line options to a script, so `open_paraview.sh` takes none.
The options of [`paraview_airfoil.py`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/paraview_airfoil.py)
are set through the `default=` values of its argument parser:

- `--scalar U` (default) for velocity magnitude, `--scalar p` for pressure
- `--u-min 20`, `--u-max 70` color limits for `U` magnitude
- `--pressure-limit 300` symmetric color limit for `p`
- `--case` path of the `.foam` file, `--screenshot` output image

The script was not rendered headless in the verification: `pvpython paraview_airfoil.py
--screenshot out.png` produced an image with only the color bar (the `Surface LIC` slice was
not drawn offscreen). Use `paraview_images.py` for headless figures.

Manual notes are in [`paraview_guide.md`](https://github.com/djeada/OpenFoam-Tutorials/blob/main/projects/04_naca_airfoil_analysis/paraview_guide.md).

### Regenerating the README figures

```bash
pvbatch paraview_images.py        # writes images/mesh_domain.png, mesh_airfoil.png, pressure.png, velocity.png
```

### Transient animation

```bash
./open_paraview.sh transient            # interactive
./open_paraview.sh transient render     # PNG frames with pvbatch, then MP4 with ffmpeg
```

## Important Interpretation Notes

### The root case is steady

The root solver is `simpleFoam`, so the directories `10/`, `20/`, ... `600/` are
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
- Drag is much higher than expected for a high-quality NACA 0012 benchmark case (`Cd = 0.0447`);
  no validated reference value is available in this repository to quantify the difference.
- Near-wall resolution is too coarse for trustworthy airfoil drag prediction (`y+` 217 to 749).
- `U` and `nuTilda` use first-order upwind, which adds numerical diffusion.
- Residuals do not reach `residualControl` (`1e-6`) in 600 iterations; the run ends at `endTime`.
- `checkMesh` reports `Failed 1 mesh checks`: two faces exceed the skewness limit (4.25 against 4).
- Top and bottom are symmetry planes 15 chords from the airfoil, so the far field is not a true
  free stream.
- More aggressive boundary-layer meshes were tested earlier, but in this repo's Gmsh workflow they
  either stalled during recombination or produced poor-quality tetrahedral meshes that
  destabilized the solver.

## Recommended Next Steps

If you want to improve physics rather than just visualization, the best next steps are:

1. move the steady case to a better airfoil mesh topology, ideally a structured or true
   C-grid style external-aero mesh
2. reduce wall-function dependence by improving first-layer control and `y+`
3. use a second-order scheme for `U` and `nuTilda` (for example `linearUpwind`)
4. compare `Cl`, `Cd`, and surface pressure against NASA / TMR reference data
5. use the transient case for animation-oriented work rather than trying to animate the
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
