#=
    Benchmark comparison between C++ and Julia of my implementation

    Calculates the magnetostatic field (demagnetizing field)
    for a given magnetization, using FEM-BEM.
=#

include("../../src/gmsh_wrapper.jl")
include("../../src/BEM.jl")

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
    return

    # Element centroids
    centroids::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[1,k] = sum(mesh.p[1,nds])/4
        centroids[2,k] = sum(mesh.p[2,nds])/4
        centroids[3,k] = sum(mesh.p[3,nds])/4
    end

    # Adjust to 0 indexing
    t::Matrix{Int32} = mesh.t .- 1
    surfaceT::Matrix{Int32} = mesh.surfaceT .- 1

    # Magnetization
    m::Matrix{Float64} = zeros(3,mesh.nv)
    m[1,:] .= 1

    # Prepare the output
    u::Vector{Float64} = zeros(mesh.nv)


    # ----- Using Julia -----
    @time begin
    # FEM/BEM matrices
    A = denseStiffnessMatrix(mesh)  # ij
    B = Bmatrix(mesh)        # in
    C = Cmatrix(mesh)        # mj
    D = Dmatrix(mesh)        # mn

    LHS::Matrix{Float64} = [-A B; C D]; # Final BEM matrix

    RHS::Vector{Float64} = zeros(mesh.nv + mesh.ne)
    for s in 1:mesh.ne
        nds = @view mesh.surfaceT[1:3,s]
        RHS[nds] .+= dot(mesh.normal[:,s],mean(m[:,nds],2))*mesh.AE[s]/3
    end

    # Magnetic scalar potential
    u = LHS\-RHS
    end

    # Demag field | Elements
    Hdk::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k] # all nodes of that element

        # Sum the contributions
        for j in 1:4
            _, b::Float64, c::Float64, d::Float64 = abcd(mesh.p,nds,nds[j])

            Hdk[1,k] = Hdk[1,k] - u[nds[j]]*b;
            Hdk[2,k] = Hdk[2,k] - u[nds[j]]*c;
            Hdk[3,k] = Hdk[3,k] - u[nds[j]]*d;
        end
    end

    # Magnetic field intensity
    H::Vector{Float64} = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = norm(Hdk[:,k])
    end

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Julia")
    
    scatterPlot = scatter!(ax, 
        centroids[1,:],
        centroids[2,:],
        centroids[3,:], 
        color = H, 
        colormap=:rainbow, 
        markersize=20 .* mesh.VE./maximum(mesh.VE))

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    save("Julia.png",fig)
    
    @time begin
    # ---- Using C++ ----
    @ccall "julia_wrapper.so".demag(
        u::Ptr{Float64},
        mesh.p::Ptr{Float64},
        t::Ptr{Int32},
        surfaceT::Ptr{Int32},
        mesh.normal::Ptr{Float64},
        mesh.AE::Ptr{Float64},
        mesh.nv::Int32,
        mesh.nt::Int32,
        mesh.ne::Int32,
        mesh.VE::Ptr{Float64},
        m::Ptr{Float64}
    )::Cvoid
    end

    # Magnetic field
    H_vectorField::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k];

        # Sum the contributions
        for nd in nds
            # obtain the element parameters
            _,b,c,d = abcd(mesh.p,nds,nd)

            H_vectorField[1,k] -= u[nd]*b;
            H_vectorField[2,k] -= u[nd]*c;
            H_vectorField[3,k] -= u[nd]*d;
        end
    end

    # Magnetic field intensity
    H = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = norm(H_vectorField[:,k])
    end

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Magnetic field H")
    
    scatterPlot = scatter!(ax, 
        centroids[1,:],
        centroids[2,:],
        centroids[3,:], 
        color = H, 
        colormap=:rainbow, 
        markersize=20 .* mesh.VE./maximum(mesh.VE))

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    save("c++.png",fig)
    

end # end of main

meshSize = 4
showGmsh = false
saveMesh = false

main(meshSize,showGmsh,saveMesh)

