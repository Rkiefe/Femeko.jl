# Update:
Development migrated to [Femeko.jl](https://codeberg.org/rkiefe/Femeko.jl)

All future updates will only be available at codeberg. 


<img width="720" height="258" alt="femeko_logo_v2" src="https://github.com/user-attachments/assets/2bb8ad52-ff73-4182-86b5-ddcad3cd55a4" />

# Femeko.jl

Using Gmsh, Femeko simplifies the model creation and mesh generation, and provides the relevant data in a digested and simple structure.

Here is an example of how easy it is to create a model and mesh to simulate a magnet in free space with Femeko:
```julia
    cells = [] # Holds each entity (dim, tag) property     
    addRectangle([0.0,  0.0], [2.0, 1.0], cells) # 2 by 1 magnet
    ID = addDisk([0.0, 0.0], 10.0) # Disk of radius 10 with tag 'ID'
    unifyModel(cells, ID) # Combine the geometries to create a conforming mesh
    mesh = Mesh2D(cells, 2.0, 0.1) # Create a mesh with maximum element size 2.0
                                   # and with local 0.1 size mesh refinement on the magnet's boundary

    # or import your own model with
    importCAD("STEP_Models/Fennec_Fox.step")
```
Then, you can access the 3 nodes of the first element by `nodes = mesh.t[:, 1]`, and the x,y coordinates of each node of that element with `xy = mesh.p[:, nodes]`. You also have access to the boundary elements and their respective boundary ID with `mesh.surfaceT[:, elementTag]`, this would output the 2 nodes + boundary ID (in 2D) or 3 nodes + boundary ID (in 3D).

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

Femeko has two types of C++ implementations. A full 100% rewrite to C++ and a mix between Julia and C++ with Julia's `ccall()`.
The mixed language use has:
- Magnetostatics has a complete C++ alternative available.
- Micromagnetics has a full implementation with FEM and an older implementation with FEM-BEM.

The 100% rewrite has:
- Magnetostatics


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
| `cFemeko/`        | Combined use of C++ and Julia                 |
| `Femeko.cpp/`     | 100% rewrite to C++               |
| `extern/`         | external libraries needed for C++ implementations                |

### License

Distributed under the MIT License. See `LICENSE` for more information.
