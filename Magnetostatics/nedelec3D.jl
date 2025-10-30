#=
        Calculate the magnetostatic interaction between non-linear magnetic materials and a source field
    
    Using: 
        - Nedelec shape functions (instead of linear Lagrange)
        - The magnetostatic vector potential (instead of the scalar potential)
        - Magnetic Reluctance (instead of permeability)

    Why:
        - Should be much more stable for non-linear media
        - Can include Surface and Volume electric currents
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")
include("../src/magneticProperties.jl")

using GLMakie

# Get the global nodes of a local edge label 'ie'
# of a global element label 'k'
function NodesFromLocalEdge( mesh::MESH, 
                             k, # Element index
                             ie # Local edge index (1 to 3)
                            )

    # Triangle nodes
    nds = @view mesh.t[1:4, k]
    

    if ie == 1 || ie == 2
        i = ie
        j = ie+1
    
    elseif ie == 3 || ie == 4
        i = ie
        j = 1
    
    elseif ie == 5
        i = 3
        j = 4
    
    elseif ie == 6
        i = 2
        j = 4
    end
    
    # The edge nodes must be sorted
    if nds[i] > nds[j]
        aux = i
        i = j
        j = aux
    end

    edge_nds = [nds[i], nds[j]]

    return edge_nds
end

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # vacuum magnetic permeability
    mu0::Float64 = pi*4e-7

    # Temperature
    T::Float64 = 293.0

    # Applied field | Tesla
    Bext::Vector{Float64} = [1.0, 
                             0.0,
                             0.0]

    # Convergence criteria
    picardDeviation::Float64 = 1e-4
    maxAtt::Int32 = 100

    # Data of magnetic materials
    data = DATA()
    loadMaterial( data,
                 "Materials", # Folder with materials
                 "Gd_MFT",    # Data folder of target material
                 "Gd",        # Material name
                 7.9,
                 T)

    spl = Spline1D(data.B, mu0.*(data.HofM./data.B))

    # Create model
    gmsh.initialize()
    cells = []

    # Add magnetic geometry
    # id = addCuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], cells)
    id = addSphere([0.0, 0.0, 0.0], 0.5, cells)

    # Add a container
    box = addSphere([0.0, 0.0, 0.0], 2.5)

    # Combine the geometries
    shell_id = unifyModel(cells, box)

    # Generate mesh
    mesh = Mesh(cells, meshSize, localSize, false, 2)

    # Edges (mid-points)
    edges = unique(mesh.t[5:10, :])
    ne = length(edges)

    # Map the global edge label to an ordered, local edge label
    # (the global edge label is not ordered with the node label)
    global2local_edge::Vector{Int32} = zeros(maximum(edges))
    for e in 1:ne
        edge = edges[e] # Global edge ID (its unique)
        global2local_edge[edge] = e
    end

    # Get the outer boundary edges
    id = findall(x -> x in shell_id, mesh.surfaceT[end, :])
    surfaceEdges = unique(mesh.surfaceT[4:6, id])
    localSurfaceEdges = global2local_edge[surfaceEdges]


    println("\nNumber of elements: ", mesh.nt)
    println("Number of Inside elements ", length(mesh.InsideElements))
    println("Number of edges: ", ne)
    println("Number of vertices: ", mesh.nv)
    println("Number of Inside nodes ", length(mesh.InsideNodes))
    println("Number of surface elements ", size(mesh.surfaceT, 2))
    println("")

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()

    # Relative magnetic reluctance
    nu::Vector{Float64} = ones(mesh.nt)

    # Prepare output
    Bfield::Matrix{Float64} = zeros(3, mesh.nt)
    B::Vector{Float64} = zeros(mesh.nt)
    Bold::Vector{Float64} = zeros(mesh.nt)

    att::Int32 = 0
    div::Float64 = Inf
    while att < maxAtt && div > picardDeviation
        att += 1
        Bold .= B

        # Global stiffness matrix
        A = spzeros(ne, ne)

        # Local stiffness matrix
        Ak::Matrix{Float64} = zeros(36, mesh.nt) # 6 edges per volume element
        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k] # Nodes of the linear volume element

            n::Int32 = 0
            for ie in 1:6 # For each edge of the tetrahedron

                # Global edge label
                # edge = mesh.t[4+ie, k]
                # i = global2local_edge[edge] # Local edge label
                
                # Global node labels of the edge
                edge_nds = NodesFromLocalEdge(mesh, k, ie)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, edge_nds[2]] - mesh.p[1:3, edge_nds[1]])

                # 1st order Lagrange basis function (hat function)
                _, bi, ci, di = abcd(mesh.p, nds, edge_nds[1])
                _, bj, cj, dj = abcd(mesh.p, nds, edge_nds[2])

                curlN_i = 2.0*edgeLength*cross([bi, ci, di],[bj, cj, dj])

                for je in 1:6
                    n += 1

                    # Global edge label
                    # edge = mesh.t[4+je, k]
                    # j = global2local_edge[edge] # Local edge label

                    # Global node labels of the edge
                    edge_nds = NodesFromLocalEdge(mesh, k, je)

                    # Length of edge
                    edgeLength = norm(mesh.p[1:3, edge_nds[2]] - mesh.p[1:3, edge_nds[1]])

                    # 1st order Lagrange basis function (hat function)
                    _, bi, ci, di = abcd(mesh.p, nds, edge_nds[1])
                    _, bj, cj, dj = abcd(mesh.p, nds, edge_nds[2])

                    curlN_j = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

                    # Update the stiffness matrix
                    Ak[n, k] = mesh.VE[k]*dot(curlN_i, curlN_j)*nu[k]
                
                end # Loop over the edges of the element
            end     # Loop over the edges of the element
            
        end # Loop over the volume element labels

        # Update global stiffness matrix
        n = 0
        for i in 1:6
            edge1 = global2local_edge[mesh.t[4+i, :]]
            
            for j in 1:6
                n += 1

                edge2 = global2local_edge[mesh.t[4+j, :]]
                A += sparse(edge1, edge2, Ak[n,:], ne, ne)
            end
        end

        # Load vector
        RHS::Vector{Float64} = zeros(ne)
        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k] # Nodes of the linear volume element
            
            for ie in 1:6 # For each edge of the tetrahedron
                
                # Global edge label
                edge = mesh.t[4+ie, k]
                i = global2local_edge[edge] # Local edge label
                
                # Global node labels of the edge
                edge_nds = NodesFromLocalEdge(mesh, k, ie)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, edge_nds[2]] - mesh.p[1:3, edge_nds[1]])

                # 1st order Lagrange basis function (hat function)
                _, bi, ci, di = abcd(mesh.p, nds, edge_nds[1])
                _, bj, cj, dj = abcd(mesh.p, nds, edge_nds[2])

                curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

                RHS[i] += mesh.VE[k]*dot(Bext, curlN)
            end
        
        end # Loop over the volume element labels

        # Solve the magnetostatic vector potential
        u = A\RHS

        # Get the magnetic flux B
        Bfield .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k] # Nodes of the linear volume element
            
            for ie in 1:6 # For each edge of the tetrahedron
                
                # Global edge label
                edge = mesh.t[4+ie, k]
                i = global2local_edge[edge] # Local edge label
                
                # Global node labels of the edge
                edge_nds = NodesFromLocalEdge(mesh, k, ie)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, edge_nds[2]] - mesh.p[1:3, edge_nds[1]])

                # 1st order Lagrange basis function (hat function)
                _, bi, ci, di = abcd(mesh.p, nds, edge_nds[1])
                _, bj, cj, dj = abcd(mesh.p, nds, edge_nds[2])

                curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

                Bfield[:, k] += u[i].*curlN
            end
        end

        # Norm of flux field B
        B .= 0.0
        for k in 1:mesh.nt
            B[k] = norm(Bfield[:, k])
        end

        # Update reluctance
        nu[mesh.InsideElements] = spl(B[mesh.InsideElements])

        div = maximum(abs.(B .- Bold))
        println(att, " , |B(n)-B(n-1)| = ", div)

    end

    Hfield = zeros(3, mesh.nt)
    H = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = nu[k]/mu0 * B[k]
        Hfield[:, k] = nu[k]/mu0 .* Bfield[:, k]
    end

    println("\nGenerating plots...")
    elements =  
                1:mesh.nt
                # mesh.InsideElements

    x = zeros(length(elements))
    y = zeros(length(elements))
    z = zeros(length(elements))

    for (i, k) in enumerate(elements)
        nds = @view mesh.t[1:4, k]
        
        x[i] = mean(mesh.p[1, nds])
        y[i] = mean(mesh.p[2, nds])
        z[i] = mean(mesh.p[3, nds])
    end
    
    fig = Figure()
    ax = Axis3(fig[1,1], aspect = :data)
    
    graph = arrows3d!(  ax,
                        x, 
                        y, 
                        z,
                        Bfield[1, elements], 
                        Bfield[2, elements], 
                        Bfield[3, elements]
                        , color = B[elements]
                        , lengthscale = 0.2
                        , colormap = :turbo
                    )

    # graph = scatter!(ax,
    #                  x,
    #                  y,
    #                  z
    #                  , color = B
    #                  , colormap = :turbo)

    Colorbar( fig[1, 2], graph
             , label = "B (T)"
             # , vertical = false
             )

    wait(display(fig))

end




main(5.0, 0.5, false)
