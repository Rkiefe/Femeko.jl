#=
    An example on how to calculate the lowest order nedelec shape function
    on any edge of the mesh and how to attribute a local edge label 
    to the Gmsh edges using Femeko datastructures
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")
include("../src/magneticProperties.jl")

using GLMakie

# Get the global nodes of a local edge label 'ie'
# of a global element label 'k'
function NodesFromLocalEdge( mesh::MESH, 
                             k::Int, # Element index
                             ie::Int # Local edge index (1 to 3)
                            )

    # Triangle nodes
    nds = @view mesh.t[1:3, k]
    
    # Global edge label
    edge = mesh.t[3+ie, k]

    # The edge nodes must be sorted
    i = ie
    j = ie + 1
    if j > 3
        j = 1
    end

    if nds[i] > nds[j]
        aux = i
        i = j
        j = aux
    end

    edge_nds = [nds[i], nds[j]]

    return edge_nds
end


function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Create model
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.ElementOrder", 2) # Set to quadratic
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1) # Dont conform at the boundary

    cells = []

    # Add magnetic geometry
    # id = addRectangle([0,0,0], [1.0, 1.0], cells)
    id = addDisk([0,0,0], 0.5, cells)

    # Add a container
    box = addDisk([0,0,0], 5.0)

    # Combine the geometries
    unifyModel(cells, box)

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", meshSize) # Max element size
    gmsh.model.mesh.generate(2) # 2D mesh
    mesh = MESH()

    # Node connectivity
    t_tags::Vector{Int32}, t::Vector{Int32} = gmsh.model.mesh.getElementsByType(9) # 2 for 1st order triangle ; 9 for 2nd order 6 node triangles
    mesh.t = reshape(t, 6, Int(length(t)/6))
    mesh.nt = size(mesh.t, 2)

    # Node coordinates
    _, p, _ = gmsh.model.mesh.getNodes()
    mesh.p = reshape(p, 3, Int(length(p)/3))
    mesh.nv = size(mesh.p, 2)

    println("\nNumber of elements: ", mesh.nt)
    println("Number of vertices: ", mesh.nv)
    println("")

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()

    # Edges (mid-points)
    edges = unique(mesh.t[4:6, :])
    ne = length(edges)

    # Map the global edge label to an ordered, local edge label
    # (the global edge label is not ordered with the node label)
    global2local_edge::Vector{Int32} = zeros(maximum(edges))
    for e in 1:ne
        edge = edges[e] # Global edge ID (its unique)
        global2local_edge[edge] = e
    end

    # Testing Nedelec basis function
    k = 1   # Global element label (any from 1 to mesh.nt)
    ie = 3  # Local edge label (from 1 to 3)
    
    # Global node labels of the edge
    edge_nds = NodesFromLocalEdge(mesh, k, ie)

    # Edge length
    edgeLength = norm(mesh.p[1:2, edge_nds[2]]-mesh.p[1:2, edge_nds[1]])

    # Nodes of the linear triangle
    nds = @view mesh.t[1:3, k]

    # 1st order Lagrange basis function (hat function)
    ai, bi, ci = abc(mesh.p, nds, edge_nds[1])
    aj, bj, cj = abc(mesh.p, nds, edge_nds[2])

    # Nedelec shape function on each node of the 2nd order triangle
    N = zeros(2, 6)::Matrix{Float64}
    normN = zeros(6)::Vector{Float64}

    x = zeros(6)
    y = zeros(6)
    for i in 1:6
        x[i] = mesh.p[1, mesh.t[i, k]]
        y[i] = mesh.p[2, mesh.t[i, k]]
        
        N[1, i] = (ai + bi*x[i] + ci*y[i])*bj - (aj + bj*x[i] + cj*y[i])*bi
        N[2, i] = (ai + bi*x[i] + ci*y[i])*cj - (aj + bj*x[i] + cj*y[i])*ci

        # Multiply by the length of the edge
        N[:, i] .*= edgeLength

        # Norm of the Nedelec shape function on (x, y)
        normN[i] = norm(N[:, i])
    end

    # Plot
    println("Generating plots...")
    fig = Figure()
    ax = Axis(fig[1,1], aspect = DataAspect())

    scatter!(ax, mesh.p[1, nds], mesh.p[2, nds])

    graph = arrows2d!(  ax,
                        x,
                        y,
                        N[1,:],
                        N[2,:]
                        , color = normN
                        , lengthscale = 0.1
                        , colormap = :turbo
                    )

    wait(display(fig))

end


main(0.0, 0.0, false)
