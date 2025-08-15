#=
    A magnetic field simulation example in 2D
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)
	gmsh.initialize()

	cells = []

	# Add a 2D rectangle
	id = addRectangle([0,0,0], [1, 1], cells)
    # id = addDisk([0,0,0], 1, cells)

	# Add a container
	# box = addRectangle([0,0,0], [2, 4])
    box = addDisk([0,0,0], 4)

	# Combine the geometries
	gmsh.model.occ.fragment(vcat(cells,[(2,box)]), [])
	gmsh.model.occ.synchronize()

    # Generate mesh
	mesh::MESH = Mesh2D(cells, meshSize, localSize)

    println("\nNumber of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

	# Run Gmsh GUI
    if showGmsh
	   gmsh.fltk.run()
    end
	gmsh.fltk.finalize()


    # Element centroids
    centroids::Matrix{Float64} = zeros(2,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[1,k] = sum(mesh.p[1,nds])/3
        centroids[2,k] = sum(mesh.p[2,nds])/3
    end

    # Magnetostatic simulation
    
    # Lagrange multiplier technique
    Lag::Vector{Float64} = zeros(mesh.nv)
    for k = 1:mesh.nt
        nds = @view mesh.t[:,k]
        Lag[nds] .+= mesh.VE[k]/3
    end

    # Boundary conditions
    Hext::Vector{Float64} = [1,0]
    bc::Vector{Float64} = zeros(mesh.ne)
    for e in 1:mesh.ne
        if mesh.surfaceT[3,e] == 1 # Only the container boundary
            bc[e] = dot(Hext, mesh.normal[:,e])
        end
    end

    # Magnetic permeability
    mu::Vector{Float64} = ones(mesh.nt)
    mu[mesh.InsideElements] .= 2

    # Global stiffness matrix
    A = stiffnessMatrix2D(mesh, mu)

    # # Stiffness matrix
    # A = spzeros(mesh.nv,mesh.nv)

    # # Local stiffness matrix
    # Ak::Matrix{Float64} = zeros(9, mesh.nt)
    # b::Vector{Float64} = [0,0,0]
    # c::Vector{Float64} = [0,0,0]
    # for k in 1:mesh.nt
    #     nds = @view mesh.t[:,k]
    #     for i in 1:3
    #         _, b[i], c[i] = abc(mesh.p, nds, nds[i])
    #     end

    #     aux = mesh.VE[k]*mu[k]*(b*b' + c*c')
    #     Ak[:,k] = aux[:]
    # end

    # # Update sparse global matrix
    # n = 0
    # for i in 1:3
    #     for j in 1:3
    #         n += 1
    #         A += sparse(mesh.t[i,:],mesh.t[j,:],Ak[n,:],mesh.nv,mesh.nv)
    #     end
    # end

    # Vector of the border conditions
    RHS::Vector{Float64} = zeros(mesh.nv)
    for e in 1:mesh.ne
        nds = @view mesh.surfaceT[:, e]

        # Length of the edge
        l::Float64 = norm(mesh.p[:,nds[2]] - mesh.p[:,nds[1]])

        RHS[nds[1]] = RHS[nds[1]] + 0.5*l*bc[e]
        RHS[nds[2]] = RHS[nds[2]] + 0.5*l*bc[e]
    end

    # Magnetic scalar potential
    u::Vector{Float64} = [A Lag;Lag' 0]\[-RHS;0]
    u = u[1:mesh.nv]

    # Magnetic field
    H_vectorField::Matrix{Float64} = zeros(mesh.nt,2)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k];

        # Sum the contributions
        for nd in nds
            # obtain the element parameters
            _, bi, ci = abc(mesh.p, nds, nd)

            H_vectorField[k,1] -= u[nd]*bi;
            H_vectorField[k,2] -= u[nd]*ci;
        end
    end

    # Magnetic field intensity
    H::Vector{Float64} = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = norm(H_vectorField[k,:])
    end

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="Magnetic field H")

    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        color = H[mesh.InsideElements], 
        colormap=:rainbow, 
        markersize=20) # 20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements])

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl

end

main(1.0, 0.05, true)