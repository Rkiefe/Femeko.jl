#  Heat equation with quadratic lagrange elements

include("../Fluids/fluids2D.jl")
using GLMakie

# Mesh settings
meshSize = 5.0
localSize = 0.2
showGmsh = false

mutable struct DATA

    # Mass density (g/cm3)
    density::Float64

    # Thermal conductivity W/(m K)
    k::Float64

    # Specific heat J/(kg K)
    Cp::Float64

    # Fluid viscosity mPa s or kg/(m s)
    mu::Float64

    # Constructor
    DATA() = new(1.0, 1.0, 1.0, 1.0)
end

function main(meshSize=0.0, localSize=0.0, showGmsh=false)
    gmsh.initialize()

    # Simulation settings
    timeStep::Float64 = 1e-3
    totalTime::Float64 = 1e-1
    maxSteps::Int32 = floor(totalTime/timeStep) + 1

    # List of materials
    materialProperties = Dict("blank" => DATA(),
                              "water" => DATA())

    materialProperties["blank"].mu = 1e3

    # Create model
    cells = []      # Cell IDs (dim, tag)
    cellLabels = [] # Tag of cell property

    # Add an obstacle
    # id = addDisk([0,0,0], 1.0, cells)
    id = addRectangle([0,0,0], [1.0, 1.0], cells)
    push!(cellLabels, "blank")
    
    box = addRectangle([0,0,0], [8, 8], cells) # Add tube
    push!(cellLabels, "water")

    # Combine model to create a conforming mesh
    unifyModel([(2, id)], box)

    # Generate 2nd order mesh
    # extendLocalRefinement()
    mesh::MESH = Mesh2D([(2, id)], meshSize, localSize, 2)

    println("\nNumber of elements ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of nodes ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ns)
    println("Mesh Order: ", mesh.order, "\n")

    # In Flow | Curve id: 2
    inFlow::Int32 = 2

    # Walls
    walls::Vector{Int32} = [1, 4, 5]

    # Initial temperature
    T::Vector{Float64} = zeros(mesh.nv) .+ 10.0
    T[mesh.InsideNodes] .= 1.0

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()
 
    # Define the diffusivity on the domain
    epsi::Vector{Float64} = zeros(mesh.nt)

    for i in 1:length(cells)

        id = cells[i][2] # Cell ID
        key = cellLabels[i] # Get the data set of current cell ID

        println("id: ", id, " ; key: ", key)

        # Find all elements of current cell ID
        elements = findall(x -> x==id, mesh.elementID)

        # Update viscosity value on this cell elements
        epsi[elements] .= materialProperties[key].k/(materialProperties[key].Cp * materialProperties[key].density)

    end

    # Prepare the heat transfer simulation
    println("Calculating the 2nd order mass matrix")
    M = quadraticMassMatrix2D(mesh)

    println("Calculating the 2nd order stiffness matrix")
    A = quadraticStiffnessMatrix2D(mesh, epsi)
    
    # Plot the heat transfer in real-time
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="0.0 s")
    graph = scatter!( ax
                    , mesh.p[1, :]
                    , mesh.p[2, :]
                    , color = T 
                    , colormap=:thermal 
                    , colorrange = (minimum(T), maximum(T))
                    # , markersize=5
                    )
    Colorbar(fig[1, 2], graph, label="T")
    display(fig)

    # Time iterations
    LM = M + timeStep*A # Backward Euler
    for frame in 1:maxSteps

        # Get the new temperature
        T = LM\(M*T)

        # Update plot
        # round(frame*timeStep*100.0)/100.0
        ax.title = string(frame*timeStep)*" s"
        graph.color = T

        totalHeat = sum(M * T)  # or ones(size(T))' * (M * T)
        println("Total heat at step $frame: $totalHeat")

        sleep(0.1)
    end

    wait(fig.scene)

end
main(meshSize, localSize, showGmsh)