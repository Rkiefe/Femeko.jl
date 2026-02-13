# Femeko in C++
This is currently just a trial to check if I can rewrite the Femeko.jl logic in C++. The advantage would be creating a GUI and sharable binaries.

Currently it has
- basic 2D geometry handling (rectangles and circles)
- conforming mesh generation
- converts the gmsh mesh data to femeko format (a mesh struct)

## How-to
Download the gmsh sdk
```
wget http://gmsh.info/bin/Linux/gmsh-4.15.0-Linux64-sdk.tgz
```
unzip
```
tar -xzf gmsh-4.15.0-Linux64-sdk.tgz
```
Compile the file `mesh2Dexample.cpp` with
```
g++ mesh2Dexample.cpp -o mesh2Dexample.out -I gmsh-4.15.0-Linux64-sdk/include -L gmsh-4.15.0-Linux64-sdk/lib -l gmsh -Wl,-rpath,gmsh-4.15.0-Linux64-sdk/lib
```

