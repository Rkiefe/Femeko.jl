# Kick-starts the simulation environment

using Gmsh, LinearAlgebra

# Model and mesh generation
include("geometry.jl")
include("mesh.jl")

# Finite Element Method functions
include("FEM.jl")



