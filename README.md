<img width="720" height="258" alt="femeko_logo_v2" src="https://github.com/user-attachments/assets/2bb8ad52-ff73-4182-86b5-ddcad3cd55a4" />

# Femeko.jl

Creating a Finite Element simulation from scratch is unnecessarily complicated. From mesh generation, mesh quality, geometry handling, data structures, when does your FEM implementation start?

Femeko is my answer to "I want to make my own FEM, have full control over my implementation, but I want everything else related to FEM to be as easy as possible" (mesh generation, model handling and data structures).

Femeko streamlines Gmsh functionality, making geometry handling and mesh generation more straightforward

```julia
    # Import cad file 
    importCAD("STEP_Models/Fennec_Fox.step")

    # Generate a tetrahedral mesh
    mesh = Mesh(cells, meshSize, localSize, saveMesh) # typeof(mesh) -> struct
```
The `mesh` struct holds the mesh information in a digested format, such as node connectivity and node coordinates, element volume and surface element area, etc. For example, you can access the 4 nodes of the first tetrahedron with `mesh.t[:, 1]`, which would output a vector of 4 integers; access the 3 nodes of the 2nd surface triangle with `mesh.surfaceT[1:3, 2]` and its corresponding boundary ID with `mesh.surfaceT[4, 2]` (the 4th row stores the IDs of each triangle).

## Femeko.jl currently has full fledged implementations for
- Magnetostatics (non-linear magnetic materials under applied fields, permanent magnets, etc)
- Heat equation (with implicit time stepping)
- Viscous fluid (incompressible, static)
- Heat transfer with convection to a passing fluid (2D)
- Micromagnetics (both in time and steepest descent energy minimization)

Femeko has implementations for both 3D and 2D in most physics packages.

## Femeko.jl aims to have the following, in the future
- Full 3D Heat equation (with convection to a passing fluid)
- Elastostatics (stress)

## Table of Contents
- [Fully featured examples](#examples)
- [Basic model creation and mesh gen](#functionality)
- [Installation](#installation)
- [C++ available implementations](#current-c-alternatives-covered)

## Examples

Each physics package has its own folder: `Magnetostatics`, `Micromagnetics`, `Heat` and `Fluids`. 

### Magnetostatics
Example from `Magnetostatics`, the internal magnetic field of a plate aligned with the applied field.

<img width="551" height="443" alt="H_plate" src="https://github.com/user-attachments/assets/0b03a7a4-1872-4402-a10f-1654ce149a1f" />

### Heat
Most implementations in Femeko have both a 2D and 3D version. Here is a snapshot of the 2D heat simulation

<!--
<img width="551" height="443" alt="heat_2d" src="https://github.com/user-attachments/assets/232fba09-f23b-4201-9c46-4a996075fa89" />
-->

<img width="550" height="550" alt="Heat2D" src="https://github.com/user-attachments/assets/2c7c0250-78b3-484a-9f52-2062bc3cbef1" />


### Fluids
Femeko.jl can incorporate quadratic order elements, solving the stokes equation of a viscus fluid for both pressure and velocity.

<img width="1596" height="474" alt="Fluid_Simulation" src="https://github.com/user-attachments/assets/addc0521-1092-4262-8a01-dfd17651fdc9" />

### Micromagnetics
The Micromagnetics package has two distinct functionalities, based on the Landau-Lifshitz equation: the magnetization over time of the system; and an energy minimization by the steepest descent algorithm. The solver incorporates an external magnetic field, the demagnetizing field, the exchange field and the anisotropy field. This solver was validated against OOMMF, replicating Fig 2. of this article https://doi.org/10.1109/TMAG.2008.2001666 .

<img width="551" height="443" alt="M_time_permalloy" src="https://github.com/user-attachments/assets/5434942c-a6dd-4444-aadf-c945c17e593b" />


## Functionality
Make a high quality 3D mesh of your model and get all the properties you need, easily accessible in a simple MESH() object, powered by Gmsh. Make a local mesh refinement, based on the volume ID, or just set a target mesh size.

Automatically create a bounding shell for your 3D model, simplifying magnetic field simulations. You can define the scale of your bounding shell directly in import phase of your .step file, or keep the default "5x larger". The local refinement is automatically set for every cell that isn't the container volume.

You can import your geometry (and automatically create a bounding shell for open boundary problems) with
```
importCAD(fileName)
```

Or make your own geometry with cuboids
```
box = addCuboid([X, Y, Z], [W, D, H], cells)
```
And/or spheres as
```
addSphere([X, Y, Z], radius, cells)
```
where `cells` is an array of all current volume ID's. Each cell ID you add is tracked for you.
You can generate a mesh for your model simply by
```
mesh = Mesh(cells, meshSize, localSize, saveMesh)
```

`meshSize` sets the overall target element size, and `localSize` sets a target local element size, automatically defined for every volume ID in `cells`.

![twoBalls](https://github.com/user-attachments/assets/3b9549ba-3968-40f1-94a4-5c21ce37ca9e)

Both internal and bounding shell surfaces are preserved. You can neatly access the nodes of any surface element of your mesh directly: `mesh.surfaceT[:, i]` outputs a vector of length 3 with the node labels that make the `i` surface element. You also have direct access to mesh element volumes, surface triangle normals and the area of each surface triangle.

The output mesh object is optimized for Finite-Element simulations, see `meshExample.jl` for two simple examples of a) importing a cad file and b) making your own model with simple shapes.

### Installation
Main install:
- Open the Repl
- Press ']' key to switch to the package manager "pkg >"
- `add Gmsh LinearAlgebra SparseArrays IterativeSolvers Dierckx DelimitedFiles`

Most examples plot the results using Makie.jl. It is recommended that you do
- `add GLMakie` or/and `add CairoMakie`

Compiling C++ alternative implementations:
- First update your clone of the repository to include eigen by going to the terminal and run `git submodule update --init --recursive` while in the Femeko.jl folder
- Move to the `cFemeko/Magnetostatics/` folder and compile `magnetostatics.cpp` with `g++ -O3 -fPIC -shared -o magnetostatics.so magnetostatics.cpp`
- Adding this flag is recommended: `-fopenmp`

### Current C++ alternatives covered

- Magnetostatics has a complete C++ alternative available.
- Micromagnetics has a full implementation based on FEM and an older implementation based on FEM-BEM.

### Source code organization
Femeko is organized as:
| Directory         | Contents                                                         |
| -                 | -                                                                |
| `src/`            | source code to handle the geometry, mesh, and FEM functionality  |
| `Fluids/`         | Fluid simulation implementation                                  |
| `Heat/`           | Heat simulation implementation                                   |
| `Magnetostatics/` | Magnetostatics simulation implementation                         |
| `Micromagnetics/` | Micromagnetics simulation implementation                         |
| `STEP_Models/`    | Example .STEP files to import 3D models                          |
| `cFemeko/`        | C++ variants of Femeko for increased performance                 |
| `extern/`         | external libraries needed for C++ implementations                |

### License

Distributed under the MIT License. See `LICENSE` for more information.
