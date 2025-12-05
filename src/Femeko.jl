# Kick-starts the simulation environment

using LinearAlgebra, SparseArrays

# Model and mesh generation
include("geometry.jl")
include("mesh.jl")

# Finite Element Method functions
include("FEM.jl") 
