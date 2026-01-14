#=
    3D viscous fluid simulation

    Gets the static velocity and pressure of a viscous fluid
    flowing around an obstacle, following the stokes equation

    The implementation considers mixed-elements. The pressure is defined
    over linear lagrange elements, and the velocity is defined over
    quadratic lagrange elements.
=#


include("../src/Femeko.jl")
# using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    gmsh.initialize()

    cells = []
    
    # Add obstacle
    addSphere([0,0,0], 1.0, cells)

    # Add tube
    box = addCuboid([0,0,0], [5.0, 20.0, 5.0])

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id, box = unifyModel(cells, box)

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
    gmsh.finalize()


    # mesh.edges # Global edge ID (mid-points)

    # Sort the quadratic mesh vertices and edge midpoints
    # Vertices must start from 1 to 'nVertices'. Edge midpoints  must 
    # start from 'nVertices'+1 to mesh.nv

    # 1st order element nodes
    vertices::Vector{Int32} = unique(mesh.t[1:4, :])
    nVertices::Int32 = length(vertices)
    
    # 2nd order element nodes -> mesh.edges::Vector{Int32}
    
    # Map the global Gmsh node IDs to a local ordered ID
    localNodeID::Vector{Int32} = zeros(mesh.nv)

    # Map global 1st order mesh nodes to a local ID
    for (i, ID) in enumerate(vertices) # (local node ID, Global node ID)
        localNodeID[ID] = i
    end

    # Map the global 2nd order mesh nodes (edges) to a local ID (starting after nVertices)
    for (i, ID) in enumerate(mesh.edges) # (local edge ID, Global edge ID)
        localNodeID[ID] = nVertices + i
    end


    # Testing quadratic basis function
    nds = @view mesh.t[:, 1] # All nodes of element
        
    # Target first node
    S = quadraticBasis(mesh, nds, nds[1])

    # Evaluate at target point:
    for i in 1:10
        xt = mesh.p[1, nds[i]]
        yt = mesh.p[2, nds[i]]
        zt = mesh.p[3, nds[i]]

        u = S[1] + S[2]*xt + S[3]*yt + S[4]*zt + 
            S[5]*xt^2 + S[6]*xt*yt + S[7]*xt*zt + 
            S[8]*yt^2 + S[9]*yt*zt + S[10]*zt^2

        println(u)
    end

    # Passed the test


end # main()



meshSize = 0.0
localSize = 0.0
showGmsh = false

main(meshSize, localSize, showGmsh)