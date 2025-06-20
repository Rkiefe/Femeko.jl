#=
    Calculates the magnetostatic field (demagnetizing field)
    for a given magnetization, using FEM-BEM.

    All the FEM-BEM calculations are done in C++
=#

include("../../src/gmsh_wrapper.jl")
include("../../src/FEM.jl")

# For plots | Uncomment the plot section of "main()"
using GLMakie

function main(meshSize=0,showGmsh=true,saveMesh=false)
    #=
        This creates the mesh. C++ then handles the simulation
    =#
    
    mu0 = pi*4e-7                       # vacuum magnetic permeability
    Hext::Vector{Float64} = [1,0,0]     # T

    # Create a geometry
    gmsh.initialize()

    # >> Model
    L::Vector{Float64} = [1.65, 1.65, 0.04] # Dimensions of the object

    # Add the body
    addCuboid([0,0,0],L)

    # Generate Mesh
    mesh = Mesh([],meshSize,0.0,saveMesh)
    
    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))

    # Volume of elements of each mesh node | Needed for the demagnetizing field
    Vn::Vector{Float64} = zeros(mesh.nv)

    # Integral of basis function over the domain | Needed for the exchange field
    nodeVolume::Vector{Float64} = zeros(mesh.nv)
    
    for k in 1:mesh.nt
        Vn[mesh.t[:,k]]         .+= mesh.VE[k]
        nodeVolume[mesh.t[:,k]] .+= mesh.VE[k]/4
    end
    

    # Adjust to 0 indexing
    t::Matrix{Int32} = mesh.t .- 1
    surfaceT::Matrix{Int32} = mesh.surfaceT .- 1

    # Magnetization
    m::Matrix{Float64} = zeros(3,mesh.nv)
    m[1,:] .= 1

    # Prepare the output
    Hd::Matrix{Float64} = zeros(3,mesh.nv)

    # Calculate the magnetostatic field
    @ccall "julia_wrapper.so".demag(
        Hd::Ptr{Float64},
        mesh.p::Ptr{Float64},
        t::Ptr{Int32},
        surfaceT::Ptr{Int32},
        mesh.normal::Ptr{Float64},
        mesh.AE::Ptr{Float64},
        mesh.nv::Int32,
        mesh.nt::Int32,
        mesh.ne::Int32,
        mesh.VE::Ptr{Float64},
        Vn::Ptr{Float64},
        m::Ptr{Float64}
    )::Cvoid

    # Magnetic field intensity
    H::Vector{Float64} = zeros(mesh.nv)
    for i in 1:mesh.nv
        H[i] = norm(Hd[:,i])
    end

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Magnetic field H")
    
    scatterPlot = scatter!(ax, 
        mesh.p[1,:],
        mesh.p[2,:],
        mesh.p[3,:], 
        color = H, 
        colormap=:rainbow
        ) # markersize=20 .* mesh.VE./maximum(mesh.VE)

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    
    # save("H.png",fig)

end # end of main

meshSize = 4
showGmsh = false
saveMesh = false

main(meshSize,showGmsh,saveMesh)

