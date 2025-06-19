# Femeko.jl
This repository aims to be as user-friendly and readable as possible so that people can adapt Femeko.jl to their needs or learn how to make their own implementations from a friendly starting point.

The mesh generation is handled by Gmsh. Femeko.jl has a wrapper to Gmsh, to facilitate model creation and mesh generation (with straight forward local refinement). The plots are powered by Makie.jl. That's it, everything else is Femeko.

## Femeko.jl currently has full fledged implementations for
- Magnetostatic fields from permanent magnets
- Calculate the magnetic field of paramagnets under an external field
- Full fledged micromagnetic simulations (both in time and steepest descent energy minimization)

## Femeko.jl aims to have the following, in the future
- Full Heat equation (with convection to a passing fluid)
- Viscous fluid dynamics (in-compressible)
- Elastostatics (stress)


## Table of Contents
- [Fully featured examples](#examples)
- [Basic model creation and mesh gen](#functionality)
- [Installation](#installation)
- [C++ available implementations](#current-c-alternatives-covered)

<!-- - [License](#license) -->

## Examples

This framework already has a few implementations such as `Magnetostatics` and `Micromagnetics`. The Micromagnetics package has two distinct functionalities, based on the Landau-Lifshitz equation: the magnetization over time of the system; and an energy minimization by the steepest descent algorithm. The solver incorporates an external magnetic field, the demagnetizing field, the exchange field and the anisotropy field. This solver was validated against OOMMF, replicating Fig 2. of this article https://doi.org/10.1109/TMAG.2008.2001666 .

![M_time_permalloy](https://github.com/user-attachments/assets/5434942c-a6dd-4444-aadf-c945c17e593b)


Here is a direct output of magneticField.jl example `Magnetostatics`, the internal magnetic field of a plate aligned with the aplied field.
![H_plate](https://github.com/user-attachments/assets/0b03a7a4-1872-4402-a10f-1654ce149a1f)

And now for a sphere. Both of these geometries were created with simple commands.
![magneticfield_example](https://github.com/user-attachments/assets/86fc8c7c-7e8a-4f6b-a807-0df6720f1a1b)

## Functionality
Make a high quality 3D mesh of your 3D model and get all the properties you need, easily accessible in a simple MESH() object, powered by Gmsh. Make a local refinement of your model, based on the volume ID, or just set a target mesh size.

Automatically create a bounding shell for your 3D model, simplifying magnetic field simulations. You can define the scale of your bounding shell directly in import phase of your .step file, or keep the default "5x larger". The local refinement is automatically set for every cell that isn't the container volume.

You can import your geometry (and automatically create a bounding shell for open boundary problems) with
```
importCAD(file)
```

Or make your own geometry with cuboids
```
box = addCuboid(position_center,[W,D,H])
```
And/or spheres as
```
addSphere(position_center,sphere_radius,cells)
```
Where `cells` is an array of volume ID's inside the bounding shell (considering you have solids inside a defined space by a bounding shell, such as with open boundary problems). Each cell ID you add is tracked for you.
You can generate a mesh for your volume simply by
```
mesh = Mesh(cells,meshSize,localSize,saveMesh)
```

![twoBalls](https://github.com/user-attachments/assets/3b9549ba-3968-40f1-94a4-5c21ce37ca9e)

Both internal and bounding shell surfaces are preserved. You can access the surface ID of each surface triangle of your mesh directly.

Automatically get the mesh element volumes, surface triangle normals and the area of each surface triangle.

The output mesh object is optimized for Finite-Element simulations, see ´meshExample.jl´ for two simple examples of a) importing a cad file and b) making your own model with simple shapes.

### Installation
Main install:
- Open the Repl
- Press ']' key to switch to the package manager "pkg >"
- `add Gmsh LinearAlgebra SparseArrays`

If you want the built in plots:
- `add GLMakie` or/and `add CairoMakie`

Compiling C++ alternative implementation:
- First update your clone of the repository to include eigen by going to the terminal and run `git submodule update --init --recursive` while in the Femeko.jl folder
- Move to the `cFemeko/Magnetostatics/` folder and compile `julia_wrapper.cpp` with `g++ -O3 -fPIC -shared -o julia_wrapper.so julia_wrapper.cpp`
- Note that additional flags are available such as `-fopenmp`

### Current C++ alternatives covered
This C++ alternatives are not for speed increases, suprisingly. They are to avoid the garbage collector and honestly some messy crashes that happen on Julia's side. This way, the mesh and plots are handled by Julia, while everything else (including the linear solver) is handled by C++.

- Magnetostatics has a complete C++ alternative available.
- Micromagnetics has a full implementation of the magnetostatic scalar potential with FEM-BEM in C++

Missing in Micromagnetics
- Exchange field
- Applied field (should be trivial)
- Anisotropy field
- Time step
- Steepest descent energy minimization


### License

Distributed under the MIT License. See `LICENSE` for more information.
