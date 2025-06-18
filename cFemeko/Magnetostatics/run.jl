#=
    Magnetostatic field simulation of a magnetic susceptible material
    under a uniform applied field
=#

include("../../src/gmsh_wrapper.jl")
include("../../src/FEM.jl")

# For plots | Uncomment the plot section of "main()"
using GLMakie

function main(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)
    #=
        This creates the mesh. C++ then handles the simulation
    =#
    
    mu0 = pi*4e-7                       # vacuum magnetic permeability
    Hext::Vector{Float64} = [1,0,0]     # T

    # Create a geometry
    gmsh.initialize()

    # >> Model
    L::Vector{Float64} = [1.65, 1.65, 0.04] # Dimensions of the object

    # Create a bounding shell
    box = addSphere([0,0,0],5*maximum(L))

    # Get how many surfaces compose the bounding shell
    temp = gmsh.model.getEntities(2)            # Get all surfaces of current model
    bounding_shell_n_surfaces = 1:length(temp)    # Get the number of surfaces in the bounding shell

    # List of cells inside the container
    cells = []

    # Add the magnetic body
    addCuboid([0,0,0],L,cells,true)

    # Fragment to make a unified geometry
    _, fragments = gmsh.model.occ.fragment([(3, box)], cells)
    gmsh.model.occ.synchronize()

    # Update container volume ID
    box = fragments[1][1][2]

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)
    
    # Get bounding shell surface id
    mesh.shell_id = gmsh.model.getAdjacencies(3, box)[2]

    # Must remove the surface Id of the interior surfaces
    mesh.shell_id = mesh.shell_id[bounding_shell_n_surfaces] # All other, are interior surfaces

    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

    # Adjust to 0 indexing
    t::Matrix{Int32} = mesh.t .- 1
    surfaceT::Matrix{Int32} = mesh.surfaceT .- 1
    shell_id::Int32 = mesh.shell_id[1] - 1

    # Permeability
    mu::Vector{Float64} = ones(mesh.nt)
    mu[mesh.InsideElements] .= 1.0 + 2.0

    # Prepare the output
    u::Vector{Float64} = zeros(mesh.nv)

    @time @ccall "../src/julia_wrapper.so".cMagnetoStatics(
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
        shell_id::Int32
    )::Cvoid

    # Magnetic field
    H_vectorField::Matrix{Float64} = zeros(mesh.nt,3)
    for k in 1:mesh.nt
        nds = mesh.t[:,k];

        # Sum the contributions
        for nd in nds
            # obtain the element parameters
            _,b,c,d = abcd(mesh.p,nds,nd)

            H_vectorField[k,1] -= u[nd]*b;
            H_vectorField[k,2] -= u[nd]*c;
            H_vectorField[k,3] -= u[nd]*d;
        end
    end

    # Magnetic field intensity
    H::Vector{Float64} = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = norm(H_vectorField[k,:])
    end

    # Element centroids
    centroids::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[1,k] = sum(mesh.p[1,nds])/4
        centroids[2,k] = sum(mesh.p[2,nds])/4
        centroids[3,k] = sum(mesh.p[3,nds])/4
    end

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Magnetic field H")
    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        centroids[3,mesh.InsideElements], 
        color = H[mesh.InsideElements], 
        colormap=:rainbow, 
        markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements]))

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    
    # # save("H.png",fig)

end # end of main

meshSize = 4
localSize = 0.01
showGmsh = false
saveMesh = false

main(meshSize,localSize,showGmsh,saveMesh)

