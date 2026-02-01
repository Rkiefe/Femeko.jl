# 3D Viscous fluid simulation (pressure and velocity)
include("fluids3D.jl")
using GLMakie

meshSize = 0.0
localSize = 0.0
showGmsh = false

function main(meshSize=0.0, localSize=0.0, showGmsh=false)
    gmsh.initialize()

    # Setup
    viscosity = 1.0                   # Fluid viscosity
    velocity::Vector{Float64} = [0.0, 
                                 1.0,
                                 0.0] # Inflow

    # Boundary surface IDs 
    inFlow = 4
    walls = [1, 2]
    outFlow = 3

    # Create 3D model
    cells = [] # Store the obstacles cell ID
    
    addSphere([0,-2.5,0], 0.5, cells)              # Add obstacle
    box = addCylinder([0,-5,0], [0, 10, 0], 2.0)  # Add tube

    # Unify the volumes for a single geometry
    _, box = unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize, false, 2)

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

    # Define the viscosity on each element of the mesh
    mu::Vector{Float64} = zeros(mesh.nt) .+ viscosity
    # mu[mesh.InsideElements] .= 1e3*viscosity

    # Run the fluid simulation
    u, P, vertices, _ = fluid3D( mesh, velocity, mu, inFlow, walls)
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