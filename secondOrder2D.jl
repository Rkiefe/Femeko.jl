#=
    An example on how to generate a 2nd order 2D triangular mesh
=#

include("src/gmsh_wrapper.jl")

using GLMakie

function test(showGmsh=true)
    
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.ElementOrder", 2)       # Set to quadratic
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)  # Dont conform at the boundary

    # Add object
    cell = addDisk([0.0, 0.0], 0.5)

    # Mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5) # Max element size
    gmsh.model.mesh.generate(2) # generate 2D mesh

    # Create a FEMEKO mesh data structure
    mesh = MESH()

    # Node connectivity
    t_tags::Vector{Int32}, t::Vector{Int32} = gmsh.model.mesh.getElementsByType(9) # 2 for 1st order triangle ; 9 for 2nd order 6 node triangles
    mesh.t = reshape(t, 6, Int(length(t)/6))
    mesh.nt = size(mesh.t, 2)

    # Node coordinates
    _, p, _ = gmsh.model.mesh.getNodes()
    mesh.p = reshape(p, 3, Int(length(p)/3))
    mesh.nv = size(mesh.p, 2)

    println("")
    println("Number of elements ", mesh.nt)
    println("Number of nodes ", mesh.nv)
    # println("Number of surface elements ",size(mesh.surfaceT,2))
    println("")

    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    # Plot a single volume element
    nds = mesh.t[:, 1]
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    scatter!(ax, 
             mesh.p[1, nds], 
             mesh.p[2, nds], 
             mesh.p[3, nds])
    
    # The first 3 nodes are the nodes of the linear triangle
    scatter!(ax, 
             mesh.p[1, nds[1:3]], 
             mesh.p[2, nds[1:3]], 
             mesh.p[3, nds[1:3]])
    
    wait(display(fig))

end