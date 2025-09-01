# 3D heat equation solver

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

using GLMakie

function main(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)

    # Simulation settings
    timeStep::Float64  = 0.001
    totalTime::Float64 = 0.1

    # Create a model
    gmsh.initialize()

    # List of cells inside the container
    cells = []

    # Add a cuboid
    addCuboid([0,0,0], [1,1,1], cells, true)

    # Create a bounding shell
    box = addSphere([0,0,0], 5)

    # Unify the volumes for a single geometry and get the bounding shell
    unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize, saveMesh)

    if showGmsh # Show GUI
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))


    # Element centroids
    centroids::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[:,k] = mean(mesh.p[:,nds],2)
    end

    # Heat conduction
    thermalCond::Vector{Float64} = zeros(mesh.nt) .+ 1
    thermalCond[mesh.InsideElements] .= 2

    # Initial Temperature
    T::Vector{Float64} = zeros(mesh.nv)
    T .= 0
    T[mesh.InsideNodes] .= 10

    # Stiffness matrix
    A = stiffnessMatrix(mesh, thermalCond)

    # Mass matrix
    Mlocal::Matrix{Float64} = 1/20 * [2 1 1 1;
                                      1 2 1 1;
                                      1 1 2 1;
                                      1 1 1 2]; # 3D local mass matrix

    M = spzeros(mesh.nv,mesh.nv)
    Mk::Matrix{Float64} = zeros(16, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        Mk[:,k] = mesh.VE[k]*Mlocal[:];
    end

    # Update sparse global matrix
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            M += sparse(mesh.t[i,:],mesh.t[j,:],Mk[n,:],mesh.nv,mesh.nv)
        end
    end

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Heat simulation")

    scatterPlot = scatter!(ax, 
        mesh.p[1,:],
        mesh.p[2,:],
        mesh.p[3,:],
        color = T, 
        colormap=:thermal, 
        colorrange = (minimum(T), maximum(T)),
        markersize=5) 

    Colorbar(fig[1, 2], scatterPlot, label="Temperature")

    display(fig)

    # Time step
    time::Float64 = 0.0
    while time < totalTime
        time += timeStep
        T = (M + timeStep*A)\(M*T)

        ax.title = string(time)*" s"
        scatterPlot.color = T
        sleep(0.1)
    end
    wait(display(fig)) # Wait before closing the figure

end

main(5.0, 0.1, true, false)