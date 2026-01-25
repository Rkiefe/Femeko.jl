#=
    Solves the scalar potential for a linear magnetic material 
    inside a uniform magnetic field  using the Conjugate Gradient method

    First with C++ then with Julia, and plots both outputs side by side
        C++ is about 5x faster
=#

include("../../src/Femeko.jl")

using IterativeSolvers
using GLMakie

function main(meshSize::Float64=0.0
            , localSize::Float64=0.0
            , showGmsh=false)

    gmsh.initialize()

    mu0 = pi*4e-7                       # vacuum magnetic permeability
    Hext::Vector{Float64} = [1.0,0,0]   # T
    
    cells = [] # List of volume cells (dim, tag)

    addCuboid([0,0,0], [1.0, 1.0, 1.0], cells) # position, dimensions, cell list, update cell list
    box = addSphere([0,0,0], 5.0) # Create a bounding shell

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id, box = unifyModel(cells, box)

    # Generate Mesh
    extendLocalRefinement(0) # Keep local refinement on the boundary
    mesh = Mesh(cells, meshSize, localSize)

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

    if showGmsh
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()

    # Adjust to 0 indexing for C++
    t::Matrix{Int32} = mesh.t .- 1
    surfaceT::Matrix{Int32} = mesh.surfaceT .- 1
    shell::Int32 = shell_id[1] - 1

    # Permeability
    mu::Vector{Float64} = zeros(mesh.nt) .+ mu0
    mu[mesh.InsideElements] .= mu0*(1.0 + 2.0)

    # Prepare the output
    u::Vector{Float64} = zeros(mesh.nv)

    println("\nSolving the scalar potential with C++")
    @time begin # Solve with C++

        # Run C++ version of FEM magnetostatic solver
        # Updates the input scalar potential 'u'
        @ccall "./magnetostatics.so".scalarPotential(
            u::Ptr{Float64},
            mesh.p::Ptr{Float64},       # Node coordinates, 3 by nv
            t::Ptr{Int32},              # Node connectivity, 4 by nt
            surfaceT::Ptr{Int32},       # Surface connectivity, 3 by ns
            mesh.normal::Ptr{Float64},  # Surface normals, 3 by ns
            mesh.nv::Int32,         # Number of nodes
            mesh.nt::Int32,         # Number of elements
            mesh.ns::Int32,         # Number of surface elements
            mesh.VE::Ptr{Float64},  # Volume of each element
            mu::Ptr{Float64},       # Permeability of each element
            Hext::Ptr{Float64},     # External field (3D vector)
            shell::Int32            # ID of the bounding shell
        )::Cvoid
    
    end

    # Magnetic field
    println("Calculating the magnetic field")
    Hfield::Matrix{Float64} = zeros(3, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k];

        # Sum the contributions
        for nd in nds
            # obtain the element parameters
            _,b,c,d = abcd(mesh.p,nds,nd)

            Hfield[1, k] -= u[nd]*b;
            Hfield[2, k] -= u[nd]*c;
            Hfield[2, k] -= u[nd]*d;
        end
        
    end # 3D vector field

    # Magnetic field intensity
    H::Vector{Float64} = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = norm(Hfield[:, k])
    end

    # -- Plot the result --
    centroids::Matrix{Float64} = zeros(3, mesh.nt) # Element centroids
    for k in 1:mesh.nt
        nds = mesh.t[:, k]
        centroids[:, k] = mean(mesh.p[:, nds], 2)
    end

    println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="C++ with Eigen")

    graph = arrows3d!(  ax
                        , centroids[1, :]
                        , centroids[2, :]
                        , centroids[3, :]
                        , mu0*Hfield[1, :]
                        , mu0*Hfield[2, :]
                        , mu0*Hfield[3, :]
                        , color = H
                        , lengthscale = 0.1
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    # Add a colorbar
    Colorbar(fig[2, 1], graph, vertical=false)
    

    # Solve with Julia
    println("\nSolving the scalar potential with Julia")
    @time begin
        # Stiffness matrix
        A = stiffnessMatrix(mesh, mu)
        RHS = BoundaryIntegral(mesh, Hext, shell_id)

        u = cg(A, -RHS)
    end

    println("Calculating the magnetic field")
    Hfield .= 0.0
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k];

        # Sum the contributions
        for nd in nds
            # obtain the element parameters
            _,b,c,d = abcd(mesh.p,nds,nd)

            Hfield[1, k] -= u[nd]*b;
            Hfield[2, k] -= u[nd]*c;
            Hfield[2, k] -= u[nd]*d;
        end
        
    end # 3D vector field

    # Magnetic field intensity
    H .= 0.0
    for k in 1:mesh.nt
        H[k] = norm(Hfield[:, k])
    end

    println("Including Julia solution to the plot")
    ax2 = Axis3(fig[1, 2], aspect = :data, title="Julia with 'IterativeSolvers'")
    graph = arrows3d!(  ax2
                        , centroids[1, :]
                        , centroids[2, :]
                        , centroids[3, :]
                        , mu0*Hfield[1, :]
                        , mu0*Hfield[2, :]
                        , mu0*Hfield[3, :]
                        , color = H
                        , lengthscale = 0.1
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    # Add a colorbar
    Colorbar(fig[2, 2], graph, vertical=false)
    wait(display(fig))


end # end of main

main(2.0, 0.1, false)
