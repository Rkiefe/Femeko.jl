#=
        Showcase of Femeko data-structure of second order tetrahedral mesh elements (10 nodes, 4 vertices and 6 edge midpoints)
    
    The mesh is generated manually, instead of with 'Mesh()', but here I demonstrate how to get
    a 'Femeko-like' data strcuture
=#

include("src/gmsh_wrapper.jl")
using GLMakie

function SecondOrder(meshSize = 0.0, showGmsh=true)
    #=
        Make your own model with cuboids and spheres
    =#
    gmsh.initialize()

    # Add a cube
    cell = addCuboid([0,0,0], [1.0, 1.0, 1.0])

    # Generate mesh
    mesh = Mesh([], meshSize, 0.0, false, 2)

    # Edges (mid-points)
    edges = unique(mesh.t[5:10, :])
    ne = length(edges)

    println("\nNumber of elements: ", mesh.nt)
    println("Number of vertices: ", mesh.nv)
    println("Number of surface elements ", size(mesh.surfaceT, 2))
    println("Number of edges: ", ne)
    println("")

    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    # Plot a single volume element
    nds = mesh.t[:,1]

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)

    # All nodes of the element
    scatter!(ax, mesh.p[1,nds], mesh.p[2,nds], mesh.p[3,nds])
    
    # Add the first 4 nodes | the nodes of the linear mesh
    scatter!(ax, mesh.p[1,nds[1:4]], mesh.p[2,nds[1:4]], mesh.p[3,nds[1:4]])
    wait(display(fig))

end

SecondOrder(0.0, true)

