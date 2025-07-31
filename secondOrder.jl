#=
    Testing second order tetrahedral mesh elements (10 nodes, 4 vertices and 6 edge midpoints)
=#

include("src/gmsh_wrapper.jl")
using GLMakie

function SecondOrder(showGmsh=true)
    #=
        Make your own model with cuboids and spheres
    =#
    gmsh.initialize()

    # Cuboid dimensions
    L::Vector{Float64} = [1.0, 1.0, 1.0]

    cell = addCuboid([0,0,0], L)

    # Mesh
    # gmsh.option.setNumber("Mesh.MeshSizeMax", 20.0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 15.0)

    gmsh.option.setNumber("Mesh.ElementOrder", 2) # Set to quadratic
    gmsh.model.mesh.generate(3)

    mesh = MESH()

    # Node connectivity
    t_tags, t = gmsh.model.mesh.getElementsByType(11) # 11 for 2nd order tetrahedral ; 4 for linear
    mesh.t = reshape(t, 10, Int(size(t,1)/10))
    mesh.nt = size(mesh.t, 2)

    # Node coordinates
    _,p,_ = gmsh.model.mesh.getNodes()
    mesh.p = reshape(p, 3, Int(size(p,1)/3))
    mesh.nv = size(mesh.p, 2)

    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    # Plot a single volume element
    nds = mesh.t[:,1]
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)
    scatter!(ax, mesh.p[1,nds], mesh.p[2,nds], mesh.p[3,nds])
    
    # The first 4 nodes are the nodes of the linear mesh
    scatter!(ax, mesh.p[1,nds[1:4]], mesh.p[2,nds[1:4]], mesh.p[3,nds[1:4]])
    wait(display(fig))
end

SecondOrder(false)

