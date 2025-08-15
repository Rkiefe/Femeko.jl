#=
    2D heat equation solver
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Simulation settings
    timeStep::Float64  = 0.001
    totalTime::Float64 = 0.05
    nSteps::Int32 = totalTime/timeStep + 1

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
    thermalCond[mesh.InsideElements] .= 2

    # Initial temperature
    T::Vector{Float64} = zeros(mesh.nv)
    T .= 0
    T[mesh.InsideNodes] .= 10

    # Interpolate the result in this coordinate
    xq = 0.75
    yq = 0.0

    # Stiffness matrix
    A = stiffnessMatrix2D(mesh, thermalCond)

    # Mass matrix
    M = massMatrix2D(mesh::MESH)

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
    Tq::Vector{Float64} = zeros(nSteps)
    time::Float64 = 0.0
    for i in 1:nSteps
        
        time += timeStep

        # New temperature
        T = (M + timeStep*A)\(M*T)

        # Interpolate the result in a target coordinate
        Tq[i] = interp2Dmesh(mesh, xq, yq, T)

        # Update plot
        ax.title = string(time)*" s"
        scatterPlot.color = T

        # Pause to let user see the change gradually
        sleep(0.1)
    end

    println("Simulation finished")
    wait(display(fig)) # Wait before closing the figure
    
    # Plot the interpolated result
    fig, ax = scatter(0.0:timeStep:totalTime, Tq)
    ax.ylabel = "Temperature"
    ax.xlabel = "Time"

    wait(display(fig)) # Wait before closing the figure


end

showGmsh = false
main(4.0, 0.1, showGmsh) # 1.0, 0.05