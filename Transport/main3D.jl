# 3D Full heat equation with convection
include("../Fluids/fluids3D.jl")

using GLMakie

meshSize = 0.0
localSize = 0.0
showGmsh = true

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Setup
    viscosity = 1.0                   # Fluid viscosity
    velocity::Vector{Float64} = [0.0, 1.0, 0.0] # Intake fluid velocity
    timeStep::Float64 = 1e-2
    totalTime::Float64 = 2.0
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

    id = addSphere([0,-2.5,0], 0.5, cells)              # Add obstacle
    push!(cellLabels, "blank")

    box = addCylinder([0,-5,0], [0, 10, 0], 2.0, cells)  # Add tube
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
    T[mesh.InsideNodes] .= 20.0

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

    # Run the fluid simulation
    u, P, vertices = fluid3D( mesh, velocity, mu, inFlow, walls)
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
                        , colormap = :redsblues,  # :CMRmap :viridis :redsblues :turbo :rainbow
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
                     , colormap = :CMRmap  # :CMRmap :viridis :redsblues :turbo :rainbow
                    )
    
    Colorbar(fig[2, 2], graph, 
             label = "Pressure", vertical = false)
    
    wait(display(fig))

end # main()

main(meshSize, localSize, showGmsh)