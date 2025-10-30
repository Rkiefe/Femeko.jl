#=
    2D viscous fluid simulation

    Gets the static velocity and pressure of a viscous fluid
    flowing around an obstacle, following the stokes equation
=#

include("../src/Femeko.jl")

using GLMakie

# Run the matrix assemblies and solve the system for pressure and velocity
function fluid2D(mesh::MESH, velocity::Vector{Float64}, mu::Vector{Float64}, inFlow, walls)

    # Gmsh orders the nodes arbitrarily
    # So I have to re-label the vertices and the midpoints
    vertexID::Vector{Int32}, nVertices::Int32 = sortMeshNodes2D(mesh)

    println("Building stiffness matrix")

    # Global Stiffness matrix
    A = spzeros(mesh.nv, mesh.nv)

    # Local Stiffness matrix
    Ak::Matrix{Float64} = quadraticLocalStiffnessMatrix2D(mesh, mu)

    # Update sparse global matrix
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            A += sparse(vertexID[mesh.t[i,:]],      # Convert original node ID to sorted node ID
                        vertexID[mesh.t[j,:]],      # Convert original node ID to sorted node ID
                        Ak[n,:], mesh.nv, mesh.nv)
        end
    end

    println("Building divergence matrix")

    # Pressure matrix
    B1::Matrix{Float64}, B2::Matrix{Float64} = divergenceMatrix2D(mesh, 
                                                                  vertexID, 
                                                                  nVertices)
    # Full matrix
    LHS =  [A spzeros(mesh.nv, mesh.nv) B1'; 
            spzeros(mesh.nv, mesh.nv) A B2';
            B1 B2 spzeros(nVertices, nVertices)]

    # Boundary conditions

    # Fluid intake
    inFlowNodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(1, inFlow)

    # No-slip boundary conditions
    wallNodes::Vector{Int32} = Int32[]
    for id in walls
        nodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(1, id)
        append!(wallNodes, nodes)
    end

    # Nodes with boundary conditions
    fixed = [
             vertexID[inFlowNodes]; 
             vertexID[wallNodes];
             mesh.nv .+ vertexID[inFlowNodes];
             mesh.nv .+ vertexID[wallNodes];
            ]

    free = setdiff(1:(2*mesh.nv + nVertices), fixed)

    gD::Vector{Float64} = zeros(2*mesh.nv + nVertices)
    gD[vertexID[inFlowNodes]] .= velocity[1]
    gD[mesh.nv .+ vertexID[inFlowNodes]] .= velocity[2]

    RHS = -LHS[free,fixed]*gD[fixed]

    println("Solving matrix equation")

    # Velocity and pressure (u and p)
    UP::Vector{Float64} = zeros(2*mesh.nv + nVertices)
    UP[fixed] .= gD[fixed]
    UP[free] = LHS[free,free]\RHS    
    
    # Velocity
    u::Matrix{Float64} = zeros(mesh.nv, 2)
    u[:, 1] .= UP[1:mesh.nv]
    u[:, 2] .= UP[mesh.nv+1:2*mesh.nv]

    # Pressure
    p::Vector{Float64} = UP[2*mesh.nv.+(1:nVertices)]

    return u, p, vertexID, nVertices
end


function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Simulation settings
    velocity::Vector{Float64} = [1.0, 0.0]
    viscosity::Float64 = 1.0

    # Create model
    gmsh.initialize()
    
    # Add an obstacle
    cells = []
    # id = addRectangle([0,0,0], [5, 3], cells)
    id = addDisk([-5,0,0], 1.0, cells)

    # Add a container
    box = addRectangle([0,0,0], [20, 5])
    # box = addDisk([0,0,0], 4)

    # Combine the geometries
    gmsh.model.occ.fragment(vcat(cells,[(2,box)]), [])
    gmsh.model.occ.synchronize()

    # Generate mesh
    mesh::MESH = Mesh2D(cells, meshSize, localSize, 2)

    println("\nNumber of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("Mesh Order: ", mesh.order)

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()
 
    # Define the viscosity on the domain
    mu::Vector{Float64} = zeros(mesh.nt) .+ viscosity

    # Define the obstacle as an extremely viscous object to set the internal velocity to zero
    mu[mesh.InsideElements] .= 1e3 * viscosity

    # Apply boundary conditions
    # NOTE! The boundary ID must be checked manually from Gmsh GUI directly

    # In Flow | Curve id: 2
    inFlow::Int32 = 2

    # Walls
    walls::Vector{Int32} = [1, 4, 5]
    # note: the missing wall is an outflow




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





    # ----- Plot results ------

    # Sort the coordinates
    x::Vector{Float64} = zeros(mesh.nv)
    x[vertexID[1:mesh.nv]] .= mesh.p[1,:]

    y::Vector{Float64} = zeros(mesh.nv)
    y[vertexID[1:mesh.nv]] .= mesh.p[2,:]

    println("Printing plots")

    # Plot results
    fig = Figure()

    # Add vector field
    Axis(fig[1, 1], aspect = DataAspect(), title="Velocity field")
    velocity_plot = arrows2d!(x, y, u[:,1], u[:,2], 
                              lengthscale = 0.5,
                              color = velocityNorm,
                              colormap = :thermal)

    Colorbar(fig[2, 1], velocity_plot, 
             label = "Velocity", vertical = false)

    # Add Pressure plot
    xPressure::Vector{Float64} = zeros(nVertices)
    yPressure::Vector{Float64} = zeros(nVertices)

    # Vertices
    vertices::Vector{Int32} = unique(vec(mesh.t[1:3,:]))
    xPressure[vertexID[vertices]] .= mesh.p[1,vertices]
    yPressure[vertexID[vertices]] .= mesh.p[2,vertices]


    ax2 = Axis(fig[1, 2], aspect = DataAspect(), title="Pressure")
    sc2 = scatter!( ax2, 
                    xPressure, yPressure, 
                    color=p, 
                    colormap=:batlow,
                    markersize=10)
    
    Colorbar(fig[2, 2], sc2,
             label = "Pressure", vertical = false)

    wait(display(fig))


end

showGmsh = true
main(2.5, 0.25, showGmsh)