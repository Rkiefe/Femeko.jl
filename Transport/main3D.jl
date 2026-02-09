# 3D Full heat equation with convection
include("../Fluids/fluids3D.jl")

using GLMakie

meshSize = 0.0
localSize = 0.0
showGmsh = true

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    @warn "This implementation is not fully confirmed to be working"

    # Setup
    viscosity = 1.0                   # Fluid viscosity
    velocity::Vector{Float64} = [0.0, 1.0, 0.0] # Intake fluid velocity
    timeStep::Float64 = 1e-4
    totalTime::Float64 = 0.2
    maxSteps::Int32 = floor(totalTime/timeStep) + 1

    # List of materials
    materialProperties = Dict("blank" => DATA(),
                           "water" => DATA())

    materialProperties["blank"].mu = 1e3
    materialProperties["blank"].Cp = 3.0


    # Create 3D model
    gmsh.initialize()
    cells = [] # Store the obstacles cell ID
    cellLabels = []

    id = addSphere([0,-2.5,0], 1.0, cells)              # Add obstacle
    push!(cellLabels, "blank")

    box = addCylinder([0,-5,0], [0, 10, 0], 2.5, cells)  # Add tube
    push!(cellLabels, "water")

    # Unify the volumes for a single geometry
    _, box = unifyModel([(3, id)], box)

    # Generate Mesh
    # extendLocalRefinement()
    mesh = Mesh([(3, id)], meshSize, localSize, false, 2)

    println("\nNumber of elements: ", mesh.nt)
    println("Number of vertices: ", mesh.nv)
    println("Number of surface elements ", mesh.ns)
    println("Number of edges: ", mesh.ne)
    println("")

    if showGmsh
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end

    # Boundary surface IDs (Check in Gmsh GUI)
    inFlow = 4
    outFlow = 3
    walls = [1, 2]

    # Initial temperature
    T::Vector{Float64} = zeros(mesh.nv) .+ 0.0
    T[mesh.InsideNodes] .= 10.0

    # Define the viscosity and the diffusivity on the domain
    mu::Vector{Float64} = zeros(mesh.nt)
    rhoCp::Vector{Float64} = zeros(mesh.nt)
    conductivity::Vector{Float64} = zeros(mesh.nt)
    for i in 1:length(cells)

        id = cells[i][2] # Cell ID
        key = cellLabels[i] # Get the data set of current cell ID

        println("id: ", id, " ; key: ", key)

        # Find all elements of current cell ID
        elements = findall(x -> x==id, mesh.elementID)

        # Update viscosity value on this cell elements
        mu[elements] .= materialProperties[key].mu
        conductivity[elements] .= materialProperties[key].k
        rhoCp[elements] .= materialProperties[key].Cp * materialProperties[key].density
    end

    # Run the fluid simulation
    u, P, vertices, S = fluid3D( mesh, velocity, mu, inFlow, walls)
    gmsh.finalize()

    # Norm of velocity (defined on global node IDs)
    uNorm::Vector{Float64} = zeros(mesh.nv) 
    for i in 1:mesh.nv
        uNorm[i] = norm(u[:, i])
    end

    # Plot result
    println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data
                , title="Velocity field"
              )

    # Velocity vector field
    graph = arrows3d!(  ax
                        , mesh.p[1, :]
                        , mesh.p[2, :]
                        , mesh.p[3, :]
                        , u[1, :]
                        , u[2, :]
                        , u[3, :]
                        , color = uNorm
                        , lengthscale = 1.0/maximum(uNorm)
                        , colormap = :rainbow  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    Colorbar(fig[2, 1], graph, 
             label = "Velocity", vertical = false)

    # Pressure plot
    ax2 = Axis3(fig[1,2], aspect=:data
                , title="Velocity field"
              )

    graph = scatter!(ax2
                     , mesh.p[1, vertices] # x
                     , mesh.p[2, vertices] # y
                     , mesh.p[3, vertices] # z
                     , color = P[vertices] # Pressure
                     , markersize = 20
                     , colormap = :rainbow  # :CMRmap :viridis :redsblues :turbo :rainbow
                    )
    
    Colorbar(fig[2, 2], graph, 
             label = "Pressure", vertical = false)
    
    display(GLMakie.Screen(), fig)
    # display(fig)
    # wait(display(fig))

    # Define the boundary conditions at the intake
    intakeBC::Vector{Float64} = zeros(mesh.ns)
    for s in 1:mesh.ns        
        if inFlow == mesh.surfaceT[7, s] # ID of the boundary of current surface triangle (6 nodes + ID)
            # Approximate Dirichlet boundary conditions by penalty
            intakeBC[s] = 1e6
            T[mesh.surfaceT[1:6, s]] .= 0.0 # Temperature of the intake fluid            
        end
    end # Boundary conditions

    # Prepare the heat transfer simulation
    println("Calculating the 2nd order mass matrix")
    M = @time quadraticMassMatrix(mesh, S, rhoCp)

    println("Calculating the 2nd order stiffness matrix")
    A = @time quadraticStiffnessMatrix(mesh, S, conductivity)

    # println("Calculating the 2nd order convection matrix")
    C = @time quadraticConvectionMatrix(mesh, S, u, rhoCp)

    # println("Calculating the surface integral matrix (1D quadratic mass matrix)")
    R = @time quadraticBoundaryMassMatrix(mesh, intakeBC) 

    # Plot the heat transfer in real-time
    println("Creating a new Makie.jl window")
    fig = Figure()
    ax3D = Axis3(fig[1, 1], aspect = :data, title="0.0 s")
    graph3D = scatter!( ax3D
                    , mesh.p[1, :]
                    , mesh.p[2, :]
                    , mesh.p[3, :]
                    , color = T 
                    , colormap = :rainbow 
                    , colorrange = (minimum(T), maximum(T))
                    # , markersize=5
                    )

    Colorbar(fig[2, 1], graph3D, label="T", vertical=false)

    # XoZ plane
    X, Y, Z = plane([0,1,0], [0,0,1], [0,-1.5,0], 1.0, 15) # direction 1, direction 2, origin, radius, grid points
    Tq = zeros(size(X))
    for i in 1:size(X, 1)
        for j in 1:size(X, 2)
            Tq[i, j] = interp3Dmesh(mesh, X[i, j], Y[i, j], Z[i, j], T)
        end
    end

    # Plot slice view
    println("Adding slice view to the plot")
    ax = Axis(fig[1, 2], aspect = DataAspect(), title="Slice view")
    graph = scatter!(  ax
                     , Y[:]
                     , Z[:]
                     , color = Tq[:]
                     , colormap = :rainbow  # :CMRmap :viridis :redsblues :turbo :rainbow :thermal
                     , colorrange = (minimum(T), maximum(T))
                        # , markersize = 5
                      )

    # Add a colorbar
    Colorbar(fig[2, 2], graph, vertical=false)
    display(fig)

    println("Running heat transport simulation")
    # Time iterations
    LM = M + timeStep*(A+R+C) # Backward Euler
    for frame in 1:maxSteps-1

        # Get the new temperature
        T = LM\(M*T)

        # Update plot
        ax3D.title = string(frame*timeStep)*" s"
        graph3D.color = T

        # Create a slice view
        Tq = zeros(size(X))
        for i in 1:size(X, 1)
            for j in 1:size(X, 2)
                Tq[i, j] = interp3Dmesh(mesh, X[i, j], Y[i, j], Z[i, j], T)
            end
        end
        graph.color = Tq[:]

        sleep(1/24)
    end

    println("Simulation finished")
    wait(fig.scene)

end # main()

main(meshSize, localSize, showGmsh)