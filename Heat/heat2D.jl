#=
    2D heat equation solver
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Simulation settings
    timeStep::Float64  = 0.001
    totalTime::Float64 = 0.1

    # Create model
	gmsh.initialize()
	
    # Add a 2D rectangle
	cells = []
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

    # Heat conduction
    thermalCond::Vector{Float64} = zeros(mesh.nt) .+ 1
    thermalCond[mesh.InsideElements] .= 1

    # Initial temperature
    T::Vector{Float64} = zeros(mesh.nv)
    T .= 0
    T[mesh.InsideNodes] .= 10

    # Stiffness matrix
    A = spzeros(mesh.nv,mesh.nv)

    # Local stiffness matrix
    Ak::Matrix{Float64} = zeros(9, mesh.nt)
    b::Vector{Float64} = [0,0,0]
    c::Vector{Float64} = [0,0,0]
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        for i in 1:3
            _, b[i], c[i] = abc(mesh.p, nds, nds[i])
        end

        aux = mesh.VE[k]*thermalCond[k]*(b*b' + c*c')
        Ak[:,k] = aux[:]
    end

    # Update sparse global matrix
    n = 0
    for i in 1:3
        for j in 1:3
            n += 1
            A += sparse(mesh.t[i,:],mesh.t[j,:],Ak[n,:],mesh.nv,mesh.nv)
        end
    end

    # Mass matrix
    Mlocal::Matrix{Float64} = 1/12 *[2 1 1;
                                     1 2 1;
                                     1 1 2]

    M = spzeros(mesh.nv,mesh.nv)
    Mk::Matrix{Float64} = zeros(9, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        Mk[:,k] = mesh.VE[k]*Mlocal[:];
    end

    # Update sparse global matrix
    n = 0
    for i in 1:3
        for j in 1:3
            n += 1
            M += sparse(mesh.t[i,:],mesh.t[j,:],Mk[n,:],mesh.nv,mesh.nv)
        end
    end

    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="Heat simulation")
    scatterPlot = scatter!(ax, 
        mesh.p[1,:],
        mesh.p[2,:],
        color = T, 
        colormap=:thermal, 
        colorrange = (minimum(T), maximum(T)),
        markersize=5) 

    Colorbar(fig[1, 2], scatterPlot, label="Temperature")
    
    # Display the figure (this will open an interactive window)
    display(fig) # This is required only if runing outside the repl

    # Time step
    time::Float64 = 0.0
    while time < totalTime
        time += timeStep
        T = (M + timeStep*A)\(M*T)

        ax.title = string(time)*" s"
        scatterPlot.color = T
        # display(fig)
        sleep(0.1)
    end
    wait() # Wait before closing the figure

end

main(1.0, 0.05, false)