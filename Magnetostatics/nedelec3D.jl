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

include("../src/Femeko.jl")
include("../src/magneticProperties.jl")

using IterativeSolvers
using GLMakie

# Local stiffness matrix with Nedelec shape elements
function nedelecStiffness(mesh::MESH)

    # Local stiffness matrix
    Ak::Matrix{Float64} = zeros(36, mesh.nt) # 6x6 edges per volume element
    for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k] # Nodes of the linear volume element

        # Hat shape element for each of the 4 nodes
        hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
        for i in 1:4
            hat[1, i], 
            hat[2, i], 
            hat[3, i], 
            hat[4, i] = abcd(mesh.p, nds, nds[i])
        end

        curlN::Matrix{Float64} = zeros(3, 6)
        for ie in 1:6
            # Global node labels of the edge
            ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

            # Length of edge
            edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])
            
            # Curl of Nedelec shape element
            curlN[:, ie] = 2.0*edgeLength*cross(hat[2:4, ndi], hat[2:4, ndj])
        end

        n = 0
        # Update the stiffness matrix
        for ie in 1:6 # For each edge of the tetrahedron
            for je in 1:6
                n += 1
                Ak[n, k] = mesh.VE[k]*dot(curlN[:, ie], curlN[:, je])
            end
        end     # Loop over the edges of the element
        
    end # Loop over the volume element labels

    return Ak
end # Local stiffness matrix with Nedelec shape elements


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

    # Load vector
    RHS::Vector{Float64} = zeros(ne)
    for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k] # Nodes of the linear volume element
        
        # Hat shape element for each of the 4 nodes
        hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
        for i in 1:4
            hat[1, i], 
            hat[2, i], 
            hat[3, i], 
            hat[4, i] = abcd(mesh.p, nds, nds[i])
        end

        for ie in 1:6 # For each edge of the tetrahedron
            
            # Global edge label
            edge = mesh.t[4+ie, k]
            i = global2local_edge[edge] # Local edge label
            
            # Global node labels of the edge
            ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

            # Length of edge
            edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])

            # 1st order Lagrange basis function (hat function)
            _, bi, ci, di = hat[:, ndi]
            _, bj, cj, dj = hat[:, ndj]

            curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

            RHS[i] += mesh.VE[k]*dot(Bext, curlN)
        end
    
    end # Loop over the volume element labels

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
        Ak::Matrix{Float64} = zeros(36, mesh.nt) # 6x6 edges per volume element
        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k] # Nodes of the linear volume element

            # Hat shape element for each of the 4 nodes
            hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
            for i in 1:4
                hat[1, i], 
                hat[2, i], 
                hat[3, i], 
                hat[4, i] = abcd(mesh.p, nds, nds[i])
            end

            curlN::Matrix{Float64} = zeros(3, 6)
            for ie in 1:6
                # Global node labels of the edge
                ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])
                
                # Curl of Nedelec shape element
                curlN[:, ie] = 2.0*edgeLength*cross(hat[2:4, ndi], hat[2:4, ndj])
            end

            n::Int32 = 0
            # Update the stiffness matrix
            for ie in 1:6 # For each edge of the tetrahedron
                for je in 1:6
                    n += 1
                    Ak[n, k] = mesh.VE[k]*dot(curlN[:, ie], curlN[:, je])*nu[k]
                end
            end     # Loop over the edges of the element
            
        end # Loop over the volume element labels

        # Update sparse global stiffness matrix
        n = 0
        for i in 1:6
            edge1 = global2local_edge[mesh.t[4+i, :]]
            
            for j in 1:6
                n += 1

                edge2 = global2local_edge[mesh.t[4+j, :]]
                A += sparse(edge1, edge2, Ak[n,:], ne, ne)
            end
        end

        # Solve the magnetostatic vector potential
        u = cg(A, RHS)

        # Get the magnetic flux B
        Bfield .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k] # Nodes of the linear volume element
            
            # Hat shape element for each of the 4 nodes
            hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
            for i in 1:4
                hat[1, i], 
                hat[2, i], 
                hat[3, i], 
                hat[4, i] = abcd(mesh.p, nds, nds[i])
            end

            for ie in 1:6 # For each edge of the tetrahedron
                
                # Global edge label
                edge = mesh.t[4+ie, k]
                i = global2local_edge[edge] # Local edge label
                
                # Global node labels of the edge
                ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])

                # 1st order Lagrange basis function (hat function)
                _, bi, ci, di = hat[:, ndi]
                _, bj, cj, dj = hat[:, ndj]

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




main(5.0, 0.1, false)
