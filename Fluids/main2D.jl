include("fluids2D.jl")

using GLMakie
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
    unifyModel(cells, box)

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
    u::Matrix{Float64}, # 2 by nv
    p::Vector{Float64}, # size = nVertices
    vertexID::Vector{Int32},
    nVertices::Int32,
    vertices::Vector{Int32} = fluid2D(mesh, 
                                      velocity, # Intake fluid velocity
                                      mu,       # Viscosity
                                      inFlow,
                                      walls)    # Boundary IDs
    gmsh.fltk.finalize()

    
    velocityNorm::Vector{Float64} = zeros(mesh.nv)
    for i in 1:mesh.nv
        velocityNorm[i] = norm(u[:, i])   
    end

    # Plot results
    println("Generating plots...")
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="Velocity field")
    velocity_plot = arrows2d!(  ax
                              , mesh.p[1, :]
                              , mesh.p[2, :]
                              , u[1, :]
                              , u[2, :]
                              , lengthscale = 0.5
                              , color = velocityNorm
                              , colormap = :thermal)

    Colorbar(  fig[2, 1], velocity_plot 
             , label = "Fluid velocity field"
             , vertical = false
             )

    ax2 = Axis(fig[1, 2], aspect = DataAspect(), title="Pressure")
    sc2 = scatter!( ax2 
                  , mesh.p[1, vertices] # xPressure
                  , mesh.p[2, vertices] # yPressure 
                  , color = p[vertices]
                  , colormap = :batlow
                  , markersize = 10
                  )
    
    Colorbar(fig[2, 2], sc2
             , label = "Pressure"
             , vertical = false
             )
    wait(display(fig))


end

showGmsh = true
main(2.5, 0.25, showGmsh)