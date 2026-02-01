
#  Full 2D heat equation with convection to a passing fluid

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
    velocity::Vector{Float64} = [10.0, 0.0] # Intake fluid velocity
    viscosity::Float64 = 1.0
    timeStep::Float64 = 1e-2
    totalTime::Float64 = 2.0
    maxSteps::Int32 = floor(totalTime/timeStep) + 1

    # List of materials
    materialProperties = Dict("blank" => DATA(),
                              "water" => DATA())

    materialProperties["blank"].mu = 1e3
    materialProperties["blank"].Cp = 3.0


    # Create model
    cells = []      # Cell IDs (dim, tag)
    cellLabels = [] # Tag of cell property

    # Add an obstacle
    id = addDisk([-2.5,0,0], 1.0, cells)
    push!(cellLabels, "blank")
    
    box = addRectangle([0,0,0], [10, 5], cells) # Add tube
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

    # Boundary ID
    inFlow::Int32 = 2  # (intake)
    outFlow::Int32 = 3 # (exhaust)

    # Walls
    walls::Vector{Int32} = [1, 4, 5]

    # Initial temperature
    T::Vector{Float64} = zeros(mesh.nv) .+ 0.0
    T[mesh.InsideNodes] .= 20.0

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
    u::Matrix{Float64}, # 2 by nv
    p::Vector{Float64}, # size = nVertices
    vertexID::Vector{Int32},
    nVertices::Int32,
    vertices::Vector{Int32} = fluid2D(mesh, 
                                      velocity, # Intake fluid velocity
                                      mu,       # Viscosity
                                      inFlow,
                                      walls)    # Boundary IDs

    velocityNorm::Vector{Float64} = zeros(mesh.nv)
    for i in 1:mesh.nv
        velocityNorm[i] = norm(u[:, i])   
    end

    # Plot velocity field and pressure
    println("Generating plots...")
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="Velocity field")
    velocity_plot = arrows2d!(  ax
                              , mesh.p[1, :]
                              , mesh.p[2, :]
                              , u[1, :]
                              , u[2, :]
                              , lengthscale = 1/norm(velocity)
                              , color = velocityNorm
                              , colormap = :thermal)

    Colorbar(  fig[1, 2], velocity_plot 
             , label = "Fluid velocity field"
             # , vertical = false
             )

    ax2 = Axis(fig[2, 1], aspect = DataAspect(), title="Pressure")
    sc2 = scatter!( ax2 
                  , mesh.p[1, vertices] # xPressure
                  , mesh.p[2, vertices] # yPressure 
                  , color = p[vertices]
                  , colormap = :batlow
                  , markersize = 10
                  )
    
    Colorbar(fig[2, 2], sc2
             , label = "Pressure"
             # , vertical = false
             )
    # wait(display(fig))
    display(GLMakie.Screen(), fig)
    # return

    # Define the boundary conditions at the intake and exhaust
    intakeBC::Vector{Float64} = zeros(mesh.ns)
    exhaustBC::Vector{Float64} = zeros(mesh.ns)
    boundaryConditions::Vector{Float64} = zeros(mesh.ns)
    for s in 1:mesh.ns
        
        ID = mesh.surfaceT[4, s] # ID of the boundary of current edge
        nds = @view mesh.surfaceT[:, s] # Nodes of current edge (and the boundary ID)
        
        if ID == inFlow

            # Define the diffusivity at the in-flow boundary as infinite
            # to approximate Dirichlet boundary conditions with Robin b.c
            intakeBC[s] = 1e6

            boundaryConditions[s] = 1e6

        elseif ID == outFlow 

            # Use the velocity on the edge midpoint as the normal velocity
            exhaustBC[s] = dot(u[:, nds[3]], mesh.normal[:, s])

            boundaryConditions[s] = dot(u[:, nds[3]], mesh.normal[:, s])

        end


    end # Boundary conditions


    # Prepare the heat transfer simulation
    println("Calculating the 2nd order mass matrix")
    M = @time quadraticMassMatrix2D(mesh)

    println("Calculating the 2nd order stiffness matrix")
    A = @time quadraticStiffnessMatrix2D(mesh, epsi)

    println("Calculating the 2nd order convection matrix")
    C = @time quadraticConvectionMatrix2D(mesh, u)

    println("Calculating the surface integral matrix (1D quadratic mass matrix)")
    R = @time quadraticMassMatrix1D(mesh, intakeBC) 

    # Plot the heat transfer in real-time
    println("Creating a new Makie.jl window")
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

    println("Running heat transport simulation")
    # Time iterations
    LM = M + timeStep*(A+R+C) # Backward Euler
    for frame in 1:maxSteps

        # Get the new temperature
        T = LM\(M*T)

        # Update plot
        # round(frame*timeStep*100.0)/100.0
        ax.title = string(frame*timeStep)*" s"
        graph.color = T

        sleep(timeStep)
    end

    wait(fig.scene)

end
main(meshSize, localSize, showGmsh)