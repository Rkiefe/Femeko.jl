#  Full 2D heat equation with convection to a passing fluid

include("../Fluids/fluids2D.jl")
using GLMakie

# Mesh settings
meshSize = 5.0
localSize = 0.2
showGmsh = true

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
    velocity::Vector{Float64} = [1.0, 0.0] # In-flow velocity
    viscosity::Float64 = 1.0

    # List of materials
    materialProperties = Dict("blank" => DATA(),
                              "water" => DATA())

    materialProperties["blank"].mu = 1e3

    # Create model
    cells = []      # Cell IDs (dim, tag)
    cellLabels = [] # Tag of cell property

    # Add an obstacle
    id = addDisk([-5,0,0], 1.0, cells)
    push!(cellLabels, "blank")
    
    box = addRectangle([0,0,0], [20, 5], cells) # Add tube
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

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()
 
    # Define the viscosity and the diffusivity on the domain
    mu::Vector{Float64} = zeros(mesh.nt)
    epsi::Vector{Float64} = zeros(mesh.nt)

    for i in 1:length(cells)

        id = cells[i][2] # Cell ID
        key = cellLabels[i] # Get the data set of current cell ID

        println("id: ", id, " ; key: ", key)

        # Find all elements of current cell ID
        elements = findall(x -> x==id, mesh.elementID)

        # Update viscosity value on this cell elements
        mu[elements] .= materialProperties[key].mu
        epsi[elements] .= materialProperties[key].k/(materialProperties[key].Cp * materialProperties[key].density)

    end

    # Run fluid simulation
    u::Matrix{Float64}, 
    p::Vector{Float64}, 
    vertexID::Vector{Int32},
    nVertices::Int32 = fluid2D(mesh, 
                               velocity,         # Intake fluid velocity
                               mu,               # Viscosity
                               inFlow, walls)    # Boundary IDs

    velocityNorm::Vector{Float64} = zeros(mesh.nv)
    for i in 1:mesh.nv
        velocityNorm[i] = sqrt(sum(u[i,:].^2))   
    end


    # Heat matrices
    println("Calculating the 2nd order mass matrix")
    M = quadraticMassMatrix2D(mesh)

    println("Calculating the 2nd order stiffness matrix")
    A = quadraticStiffnessMatrix2D(mesh, epsi)

    # return

    # ----- Plot results ------

    # Sort the coordinates
    x::Vector{Float64} = zeros(mesh.nv)
    x[vertexID[1:mesh.nv]] .= mesh.p[1,:]

    y::Vector{Float64} = zeros(mesh.nv)
    y[vertexID[1:mesh.nv]] .= mesh.p[2,:]

    println("Generating plots")

    # Plot results
    fig = Figure()
    Axis(fig[1, 1], aspect = DataAspect(), title="Velocity field")
    velocity_plot = arrows2d!(x, y, u[:,1], u[:,2], 
                              lengthscale = 0.5,
                              color = velocityNorm,
                              colormap = :thermal)

    Colorbar(  fig[2, 1], velocity_plot 
             , label = "Velocity"
             , vertical = false
             )

    # Add Pressure plot
    xPressure::Vector{Float64} = zeros(nVertices)
    yPressure::Vector{Float64} = zeros(nVertices)

    # Vertices
    vertices::Vector{Int32} = unique(vec(mesh.t[1:3,:]))
    xPressure[vertexID[vertices]] .= mesh.p[1,vertices]
    yPressure[vertexID[vertices]] .= mesh.p[2,vertices]

    ax2 = Axis(fig[1, 2], aspect = DataAspect(), title="Pressure")
    sc2 = scatter!( ax2, 
                    xPressure, yPressure 
                  , color = p
                  , colormap = :batlow
                  , markersize = 10
                  )
    
    Colorbar(fig[2, 2], sc2
             , label = "Pressure"
             , vertical = false
             )
    wait(display(fig))


end
main(meshSize, localSize, showGmsh)