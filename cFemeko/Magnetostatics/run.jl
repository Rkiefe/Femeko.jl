#=
    Magnetostatic field simulation of a magnetic susceptible material
    under a uniform applied field

    All the FEM calculations are done in C++
=#

include("../../src/Femeko.jl")

# For plots | Uncomment the plot section of "main()"
using GLMakie

function main(meshSize=0, localSize=0, showGmsh=true, saveMesh=false)

    gmsh.initialize()

    mu0 = pi*4e-7                       # vacuum magnetic permeability
    Hext::Vector{Float64} = [1,0,0]     # T

    # Cuboid dimensions
    L::Vector{Float64} = [1.65, 1.65, 0.04]

    # List of volume cells
    cells = []

    # Add a cuboid
    addCuboid([0,0,0], L, cells) # position, dimensions, cell list, update cell list

    # Create a bounding shell
    box = addSphere([0,0,0], 5*maximum(L)) # Don't update the cell list and get the volume id

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id = unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize, saveMesh)

    println("\nOuter shell ID: ", shell_id)
    println("Number of elements ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of nodes ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ne)
    println("")

    if showGmsh # Show GUI
        gmsh.fltk.run()
    end
    gmsh.finalize()
    
    # Adjust to 0 indexing
    t::Matrix{Int32} = mesh.t .- 1
    surfaceT::Matrix{Int32} = mesh.surfaceT .- 1
    shell::Int32 = shell_id[1] - 1

    # Permeability
    mu::Vector{Float64} = ones(mesh.nt)
    mu[mesh.InsideElements] .= 1.0 + 2.0

    # Prepare the output
    u::Vector{Float64} = zeros(mesh.nv)

    # Run C++ version of FEM magnetostatic solver
    # Updates the input scalar potential 'u'
    @ccall "julia_wrapper.so".cMagnetoStatics(
        u::Ptr{Float64},
        mesh.p::Ptr{Float64},
        t::Ptr{Int32},
        surfaceT::Ptr{Int32},
        mesh.normal::Ptr{Float64},
        mesh.nv::Int32,
        mesh.nt::Int32,
        mesh.ne::Int32,
        mesh.VE::Ptr{Float64},
        mu::Ptr{Float64},
        Hext::Ptr{Float64},
        shell::Int32
    )::Cvoid

end # end of main

meshSize = 4
localSize = 0.1
showGmsh = false

main(meshSize, localSize, showGmsh)
