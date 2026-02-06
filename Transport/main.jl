#=  3D Full heat equation with convection
    
    Solves the fluid velocity and pressure using a quadratic mesh
    then solves the heat equation to the passing fluid using the linear mesh.

    This makes each time iteration much faster
=#
include("../Fluids/fluids3D.jl")

using GLMakie

meshSize = 0.0
localSize = 0.0
showGmsh = true

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    @warn "This implementation is not fully confirmed to be working"

    # Setup
    viscosity = 1.0                   # Fluid viscosity
    velocity::Vector{Float64} = [0.0, 10.0, 0.0] # Intake fluid velocity
    timeStep::Float64 = 1e-3
    totalTime::Float64 = 1.0
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

    # Map the global Gmsh node IDs to a local ordered ID
    localNodeID::Vector{Int32} = zeros(mesh.nv)
    nVertices = length(vertices) 

    # Map global 1st order mesh nodes to a local ID
    for (i, ID) in enumerate(vertices) # (local node ID, Global node ID)
        localNodeID[ID] = i
    end

    # Map the global 2nd order mesh nodes (edges) to a local ID (starting after nVertices)
    for (i, ID) in enumerate(mesh.edges) # (local edge ID, Global edge ID)
        localNodeID[ID] = nVertices + i
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

    # # Define the boundary conditions at the intake
    # intakeBC::Vector{Float64} = zeros(mesh.ns)
    # for s in 1:mesh.ns        
    #     if inFlow == mesh.surfaceT[7, s] # ID of the boundary of current surface triangle (6 nodes + ID)
    #         # Define the diffusivity at the in-flow boundary as infinite
    #         # to approximate Dirichlet boundary conditions with Robin b.c
    #         intakeBC[s] = 1e6
    #     end
    # end # Boundary conditions

    # Prepare the heat transfer simulation
    println("Computing the linear basis function on each node of each element")
    a::Matrix{Float64} = zeros(4, mesh.nt)
    b::Matrix{Float64} = zeros(4, mesh.nt)
    c::Matrix{Float64} = zeros(4, mesh.nt)
    d::Matrix{Float64} = zeros(4, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k]
        for i in 1:4
            a[i, k], b[i, k], c[i, k], d[i, k] = abcd(mesh.p, nds, nds[i])
        end
    end

    println("Calculating the 1nd order mass matrix")
    M = spzeros(nVertices, nVertices)
    Mk::Matrix{Float64} = 1/20 *[2 1 1 1;
                                 1 2 1 1;
                                 1 1 2 1;
                                 1 1 1 2]

    @time for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k]
        dV::Float64 = rhoCp[k]*mesh.VE[k]
        for i in 1:4
            nd1 = localNodeID[nds[i]]
            for j in 1:4
                nd2 = localNodeID[nds[j]]
                M[nd1, nd2] += Mk[i, j]*dV
            end
        end
    end

    println("Calculating the 1nd order stiffness matrix")
    A = spzeros(nVertices, nVertices)
    @time for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k]
        dV::Float64 = conductivity[k]*mesh.VE[k]

        for i in 1:4
            nd1 = localNodeID[nds[i]]
            for j in 1:4
                nd2 = localNodeID[nds[j]]
                A[nd1, nd2] += (b[i,k]*b[j,k] +
                                c[i,k]*c[j,k] +
                                d[i,k]*d[j,k]) * dV
            end
        end
    end

    println("Calculating the 1st order convection matrix")
    C = spzeros(nVertices, nVertices)
    @time for k in 1:mesh.nt
        
        dV::Float64 = rhoCp[k]*mesh.VE[k]
        nds = @view mesh.t[1:4, k]

        u_avg = mean(u[:, nds], 2) # Average velocity field on the element
                                   # note: grad phi is constant on each element
        for j in 1:4
            nd2 = localNodeID[nds[j]]
            gradPhi::Vector{Float64} = [b[j,k], c[j,k], d[j,k]]
            for i in 1:4
                nd1 = localNodeID[nds[i]]

                C[nd1, nd2] += dot(u_avg, gradPhi) *0.25*dV
            end
        end

    end # Loop over the mesh elements
 
        
    println("Calculating the boundary matrix")
    R = spzeros(nVertices, nVertices)
    Rk::Matrix{Float64} = 1/12 *[2 1 1;
                                 1 2 1;
                                 1 1 2]
    for s in 1:mesh.ns

        ID = mesh.surfaceT[7, s]
        if ID != inFlow
            continue
        end

        nds = @view mesh.surfaceT[1:3, s]
        T[nds] .= 0.0
        
        ds::Float64 = mesh.AE[s]*1e6 # Penalty method to set Dirichlet boundary conditions
        for i in 1:3
            nd1 = localNodeID[nds[i]]
            for j in 1:3
                nd2 = localNodeID[nds[j]]

                R[nd1, nd2] += Rk[i, j]*ds
            end
        end

    end


    # Plot the heat transfer in real-time
    T = T[vertices]

    println("Creating a new Makie.jl window")
    fig = Figure()
    ax3D = Axis3(fig[1, 1], aspect = :data, title="0.0 s")
    graph3D = scatter!( ax3D
                    , mesh.p[1, vertices]
                    , mesh.p[2, vertices]
                    , mesh.p[3, vertices]
                    , color = T 
                    , colormap = :rainbow 
                    , colorrange = (minimum(T), maximum(T))
                    # , markersize=5
                    )

    Colorbar(fig[2, 1], graph3D, label="T", vertical=false)
    display(fig)

    println("Running heat transport simulation")
    # Time iterations
    LM = M + timeStep*(A + C + R) # C+R # Backward Euler
    for frame in 1:maxSteps
        # println(frame, " | ", timeStep*frame, " s")

        # Get the new temperature
        T = LM\(M*T)

        # Update plot
        ax3D.title = string(frame*timeStep)*" s"
        graph3D.color = T

        sleep(1/24)
    end

    println("\nSimulation finished")
    wait(fig.scene)

end # main()

main(meshSize, localSize, showGmsh)