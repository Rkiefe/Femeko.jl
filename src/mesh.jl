#=
        Femeko mesh data structure
    
    This file holds the logic to automatically create a mesh with Gmsh with local 
    refinement on target cells and with 1st or 2nd order elements for 2D and 3D
=#

# Holds the mesh information needed for FEM simulations
mutable struct MESH

    # Order | 1 - Linear, 2 - Quadratic
    order::Int

    # Node coordinates
    p::Matrix{Float64}

    # Connectivity list
    t::Matrix{Int32}
    
    # Surface triangles
    surfaceT::Matrix{Int32}

    # Edges (updated when doing second order meshes)
    edges::Vector{Int32} # Gmsh label of each edge midpoint (2nd order mesh node)
    
    # Map the global edge label to a local edge from 1 to 'ne'
    edge2localMap::Vector{Int32} 

    # Map mesh volume elements to the parent Cell/Volume ID
    elementID::Vector{Int32}

    # Elements inside material
    InsideElements::Vector{Int32}
    
    # Nodes inside material
    InsideNodes::Vector{Int32}
    
    # Volume of each mesh element
    VE::Vector{Float64}
    
    # Area of each surface mesh element
    AE::Vector{Float64} 

    # Normal of each surface triangle
    normal::Matrix{Float64}
    
    nv::Int32           # Number of nodes
    nt::Int32           # Number of elements
    ns::Int32           # Number of surface elements
    ne::Int32           # Number of edges

    nInside::Int32      # Number of elements for the inside cells
    nInsideNodes::Int32 # Numberf of nodes for the inside cells

    shell_id # Id of the surface boundaries of the bounding shell

    # Constructor
    MESH() = new()
end

# 3D
function Mesh(cells, meshSize=0.0, localSize=0.0, saveMesh::Bool=false, order=1)
    #=
        Inputs
            cells       -> target geometries for local refinement (localSize) 
            meshSize    -> Maximum mesh element size
            localSize   -> Local mesh refinement
    
        note: you can use extendLocalRefinement(0) before Mesh() 
        to keep the local refinement near the boundaries
    =#

    # Check if the model is 3D
    if !isempty(cells) && cells[1][1] < 3 # Dimension of the cell

        return Mesh2D(cells, meshSize, localSize, order, saveMesh)

    end # If no cell is provided, assume its 3D

    if order > 1
        gmsh.option.setNumber("Mesh.ElementOrder", 2)       # Set to quadratic
        gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)  # Dont conform at the boundary (keeps elements regular)
    end

    # Make a mesh object
    mesh = MESH()

    # Get the volume IDs of the cells inside the container
    volumeID = []
    if !isempty(cells)
        for i in cells
            append!(volumeID, i[2])
        end
    end
    
    # Set local mesh size
    if localSize>0 && !isempty(cells)
        refineCell(cells, localSize, meshSize) # Set local refinement on the sphere Cell
    end

    # Set maximum element size
    if meshSize > 0
        gmsh.option.setNumber("Mesh.MeshSizeMax", meshSize)
    end

    # Generate mesh
    gmsh.model.mesh.generate(3)
    
    # -- Optimize the mesh --
    # "" (default is empty), "NetGen", "HighOrder", 
    # "Relocate3D", "HighOrderElastic", "UntangleMeshGeometry"
    if order < 2
        gmsh.model.mesh.optimize()
    else
        gmsh.model.mesh.optimize("HighOrder")   # Other wise the mesh MIGHT not be proper.
                                                # Resulting in some node labels = 0 (which should never happen since indexing starts at 1)
    end
    
    # Get node coordinates
    _,p,_ = gmsh.model.mesh.getNodes()
    mesh.p = reshape(p, 3, Int(size(p,1)/3))
    mesh.nv = size(mesh.p, 2)

    # Get all tetrahedral element GMSH tags
    # 4 - linear tetra. ; 11 - 2nd order tetra. (10 nodes)
    if order < 2
        t_tags, t = gmsh.model.mesh.getElementsByType(4)
        mesh.t = reshape(t, 4, Int(length(t)/4))
    
    elseif order > 1 # Assume quadratic order
        t_tags, t = gmsh.model.mesh.getElementsByType(11)
        mesh.t = reshape(t, 10, Int(length(t)/10))

    end
    mesh.nt = size(mesh.t, 2)

    # Set a cell ID to each mesh element
    mesh.elementID = zeros(mesh.nt)
    for k in 1:mesh.nt
        # element type , nodes of the element , dimension , id
        _, _, _, id = gmsh.model.mesh.getElement(t_tags[k]) # Cell ID of this element
        mesh.elementID[k] = id
    end

    # Get all surface triangles
    if order < 2
        surfaceT_tags, surfaceT = gmsh.model.mesh.getElementsByType(2)
        mesh.surfaceT = reshape(surfaceT, 3, Int(length(surfaceT)/3))
    
    elseif order > 1 # 2nd order triangles (6 nodes)
        surfaceT_tags, surfaceT = gmsh.model.mesh.getElementsByType(9) 
        mesh.surfaceT = reshape(surfaceT, 6, Int(length(surfaceT)/6))
    
    end

    # Expand surface triangles to include boundary id
    mesh.surfaceT = vcat( mesh.surfaceT, 
                          zeros(1, Int(size(mesh.surfaceT, 2)) ) 
                        )

    for i in 1:size(mesh.surfaceT, 2)
        # Get ID of the element of the current surface triangle
        _,_,_, id = gmsh.model.mesh.getElement(surfaceT_tags[i])
        mesh.surfaceT[end, i] = id
    end
    mesh.ns = size(mesh.surfaceT, 2) # Number of surface elements

    # Mesh edges (mid-points in a 2nd order mesh)
    if order > 1 # assume 2nd order
        mesh.edges = unique(mesh.t[5:10, :]) # Edges (mid-points)
        mesh.ne = length(mesh.edges)

        # Map the global edge label to an ordered, local edge label
        # (Gmsh does not order the 2nd order mesh nodes after the 1st order nodes)
        mesh.edge2localMap = zeros(maximum(mesh.edges))
        for e in 1:mesh.ne # For each edge
            
            # Get the global edge label (its unique)
            edge = mesh.edges[e]

            # Map the global edge label 'edge' to the local edge label 'e'
            mesh.edge2localMap[edge] = e
        end # Map global to local edge labels
    
    end # Manage mesh edges

    # Mesh elements inside the container
    if !isempty(cells)

        mesh.InsideElements = zeros(mesh.nt)
        mesh.nInside = 0
        for k::Int32 in 1:mesh.nt
            # element type , nodes of the element , dimension , id
            _, _, _, id = gmsh.model.mesh.getElement(t_tags[k])

            if id in volumeID # Current cell is part of 'cells'
                mesh.nInside += 1
                mesh.InsideElements[mesh.nInside] = k
            end
        end

        # Remove non-zeros
        mesh.InsideElements = mesh.InsideElements[1:mesh.nInside]

    else # No cells inside a container
        mesh.InsideElements = []
        mesh.nInside = 0
    end

    # Inside nodes
    aux = mesh.t[:, mesh.InsideElements]
    mesh.InsideNodes = unique(vec(aux))
    mesh.nInsideNodes = length(mesh.InsideNodes)

    # Element volumes
    mesh.VE = zeros(mesh.nt)
    for k in 1:mesh.nt
        mesh.VE[k] = elementVolume(mesh.p ,mesh.t[1:4, k])
    end

    # List of all surface triangle normals are areas
    mesh.normal = zeros(3, mesh.ns)
    mesh.AE = zeros(mesh.ns)
    for i in 1:mesh.ns
        nds = @view mesh.surfaceT[1:3, i]
        mesh.normal[:, i] = normalSurface(mesh.p, nds);
        mesh.AE[i] = areaTriangle(mesh.p[1, nds], 
                                  mesh.p[2, nds],
                                  mesh.p[3, nds])
    end

    # Save mesh 
    if saveMesh
        save2file("t.txt",mesh.t) # Save connectivity list to a .txt file
        save2file("p.txt",mesh.p)
        save2file("InsideNodes.txt",mesh.InsideNodes)
        save2file("surfaceT.txt",mesh.surfaceT)
        save2file("InsideElements.txt",mesh.InsideElements)
        save2file("VE.txt",mesh.VE)
        save2file("AE.txt",mesh.AE)
        save2file("normal.txt",mesh.normal)
    end

    return mesh
end


# 2D
function Mesh2D(cells, meshSize=0.0, localSize=0.0, order::Int=1, saveMesh=false)
    #=
        Generates a 2d triangle mesh considering that the model is made of 
        1 container and every other volume beyond the container is listed in the 'cells'

        Inputs
            cells       -> geometries that are inside the container
            meshSize    -> overall target mesh size
            localSize   -> Local mesh refinement
    =#

    mesh = MESH()

    # Store the order
    mesh.order = order

    # Get the volume IDs of the cells inside the container
    volumeID = []
    if !isempty(cells)
        for i in cells
            append!(volumeID,i[2])
        end
    end

    # Set local mesh size
    if localSize>0.0 && !isempty(cells)
        refineCell(cells, localSize, meshSize) # Set local refinement on the sphere Cell
    end

    # Set maximum element size
    if meshSize > 0
        gmsh.option.setNumber("Mesh.MeshSizeMax", meshSize)
    end

    # Set mesh order
    gmsh.option.setNumber("Mesh.ElementOrder", order) # Set to quadratic

    # Generate 2D mesh
    gmsh.model.mesh.generate(2)

    # Get element connectivity

    if order == 1
        t_tags, t = gmsh.model.mesh.getElementsByType(2) # 3 node first order triangle
        mesh.t = reshape(t,3,Int(size(t,1)/3))
    
    elseif order == 2
        t_tags, t = gmsh.model.mesh.getElementsByType(9) # 6 node second order triangle
        mesh.t = reshape(t, 6, Int(size(t, 1)/6))
    
    end

    mesh.nt = size(mesh.t,2)

    # Set a cell ID to each mesh element
    mesh.elementID = zeros(mesh.nt)
    for k in 1:mesh.nt
        # element type , nodes of the element , dimension , id
        _, _, _, id = gmsh.model.mesh.getElement(t_tags[k]) # Cell ID of this element
        mesh.elementID[k] = id
    end

    # Get node coordinates
    _,p,_ = gmsh.model.mesh.getNodes()
    mesh.p = reshape(p, 3, Int(size(p,1)/3))
    mesh.nv = size(mesh.p,2)

    # Get all boundary edges
    if order == 1
        edge_tags, edges = gmsh.model.mesh.getElementsByType(1) # 2 node line
        edges = reshape(edges, 2, Int(size(edges,1)/2))

    elseif order == 2
        edge_tags, edges = gmsh.model.mesh.getElementsByType(8) # 3 node second order line (2 nodes - vertices, 1 node - midpoint)
        edges = reshape(edges, 3, Int(size(edges,1)/3))
    end

    # Expand boundary edges to include boundary id
    edges = [edges; UInt.(zeros(1,size(edges,2)))]
    for i in 1:size(edges,2)
        # Get ID of the element of the current surface triangle
        _,_,_, id = gmsh.model.mesh.getElement(edge_tags[i]) # Element type, nodeTags, dim, tag
        edges[end, i] = id; # Set the boundary edge id
    end

    mesh.surfaceT = edges
    mesh.ns = size(edges, 2) # Its a surface element and an edge

    # Mesh elements inside the container
    if !isempty(cells)

        mesh.InsideElements = zeros(mesh.nt)
        mesh.nInside = 0
        for k in 1:mesh.nt
            # element type , nodes of the element , dimension , id
            _, _, _, id = gmsh.model.mesh.getElement(t_tags[k])

            if id in volumeID
                mesh.nInside += 1
                mesh.InsideElements[mesh.nInside] = k
            end
        end

        # Remove non-zeros
        mesh.InsideElements = mesh.InsideElements[1:mesh.nInside]

    else
        mesh.InsideElements = []
        mesh.nInside = 0
    end

    # Inside nodes
    aux::Matrix{Int32} = mesh.t[:,mesh.InsideElements]
    mesh.InsideNodes = unique(vec(aux))
    mesh.nInsideNodes = length(mesh.InsideNodes)

    # Area of each element
    mesh.VE = zeros(mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[1:3,k]
        mesh.VE[k] = areaTriangle(  mesh.p[1,nds],
                                    mesh.p[2,nds],
                                    mesh.p[3,nds])
    end

    # Edge normal and length
    mesh.normal = zeros(2, mesh.ns)
    mesh.AE = zeros(mesh.ns)
    for e in 1:mesh.ns
        nds = @view mesh.surfaceT[1:2, e]
        mesh.normal[:, e] = normalEdge(mesh.p, nds)

        mesh.AE[e] = norm(mesh.p[1:2, nds[2]] - mesh.p[1:2, nds[1]])
    end

    # Save mesh 
    if saveMesh
        save2file("t.txt",mesh.t) # Save connectivity list to a .txt file
        save2file("p.txt",mesh.p)
        save2file("InsideNodes.txt",mesh.InsideNodes)
        save2file("surfaceT.txt",mesh.surfaceT)
        save2file("InsideElements.txt",mesh.InsideElements)
        save2file("VE.txt",mesh.VE)
        save2file("AE.txt",mesh.AE)
        save2file("normal.txt",mesh.normal)
    end

    return mesh
end # 2D mesh generation

# Sort the quadratic mesh vertices and edge midpoints
function sortMeshNodes2D(mesh::MESH)
    # Vertices must start from 1 to 'nVertices'. Edge midpoints  must 
    # start from 'nVertices'+1 to mesh.nv

    vertices::Vector{Int32} = unique(vec(mesh.t[1:3,:]))
    nVertices::Int32 = length(vertices)

    edgeMidPoints::Vector{Int32} = unique(vec(mesh.t[4:6,:]))
    nEdges::Int32 = length(edgeMidPoints)
    
    # Map the mesh nodes to the ordered array of nodes + edge midpoints
    vertexID::Vector{Int32} = zeros(mesh.nv)
    for i in 1:nVertices
        vertexID[vertices[i]] = i

    end

    # Map the edge midpoints to the ordered array of nodes + edge midpoints
    for i in 1:nEdges
        vertexID[edgeMidPoints[i]] = nVertices + i
    end

    return vertexID, nVertices
end # Sort 2D quadratic mesh

# Get the normal vector to the boundary edge
function normalEdge(p, nds)
    p1 = p[:, nds[1]]
    p2 = p[:, nds[2]]
    v = p2-p1
    n = [-v[2], v[1]] # Normal to edge

    return n./norm(n)
end

# Normal to surface triangle
function normalSurface(p, nds)
    # Reshape coords into 3 points (x,y,z)
    p1 = p[:,nds[1]]
    p2 = p[:,nds[2]]
    p3 = p[:,nds[3]]
    # Edge vectors
    v1 = p2 - p1
    v2 = p3 - p1
    # Cross product (normal vector)
    n = [v1[2]*v2[3] - v1[3]*v2[2],
         v1[3]*v2[1] - v1[1]*v2[3],
         v1[1]*v2[2] - v1[2]*v2[1]]

    return n ./ sqrt(n[1]^2 + n[2]^2 + n[3]^2)
end # Normal to surface triangle

# Sets a volume element to each surface element
function surface2volume(mesh::MESH)
    surface2element::Vector{Int32} = zeros(mesh.ns)

    for s in 1:mesh.ns
        surface_nds = @view mesh.surfaceT[1:3, s]

        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k]
            overlap = intersect(nds, surface_nds)

            if length(overlap) > 2 # Found the surface element
                surface2element[s] = k
                break
            end
        end
    end

    return surface2element
end

# Area of the 3D triangle
function areaTriangle(xt::AbstractVector{Float64},
                      yt::AbstractVector{Float64},
                      zt::AbstractVector{Float64})

    Atr::Float64 = 0.5*sqrt(det([xt';yt';[1 1 1]])^2 + det([yt';zt';[1 1 1]])^2 + det([zt';xt';[1 1 1]])^2)

    return Atr
end # Area of the 3D triangle

# Mesh element volume
function elementVolume(p,nds)
    # Extract the four nodes (columns of p)
    A = p[:, nds[1]]
    B = p[:, nds[2]]
    C = p[:, nds[3]]
    D = p[:, nds[4]]

    # Compute vectors AB, AC, AD
    AB = B - A
    AC = C - A
    AD = D - A

    # Compute the scalar triple product (AB ⋅ (AC × AD))
    cross_AC_AD = [AC[2]*AD[3] - AC[3]*AD[2],
                   AC[3]*AD[1] - AC[1]*AD[3],
                   AC[1]*AD[2] - AC[2]*AD[1]]
    triple_product = AB[1] * cross_AC_AD[1] + AB[2] * cross_AC_AD[2] + AB[3] * cross_AC_AD[3]

    # Volume = (1/6) * |triple_product|
    volume = abs(triple_product) / 6.0
    return volume
end # Mesh element volume

# Set mesh extension to on/off (default off)
function extendLocalRefinement(value = 0)
    # 0 -> Dont extend mesh size, to prevent over-refinement inside entity
    # 1 -> Extend mesh refinement
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", value)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", value)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", value)
end

# Local mesh refinement on target cell
function refineCell(cell,localSize,meshSize)
    # Sets every volume in 'cell' to be locally refined with target 'localSize' 

    # Get the boundary of the cell 
    cell_boundary = gmsh.model.getBoundary(cell, false, false, false)

    # Create a distance field for local refinement
    distance_field = gmsh.model.mesh.field.add("Distance")

    if cell_boundary[1][1] < 2 # 1 -> curves
        gmsh.model.mesh.field.setNumbers(distance_field, "CurvesList", [s[2] for s in cell_boundary])
    else # 2 -> faces
        gmsh.model.mesh.field.setNumbers(distance_field, "FacesList", [s[2] for s in cell_boundary])
    end

    # Set number of sampling points over the surfaces
    gmsh.model.mesh.field.setNumber(distance_field, "Sampling", 500) # default is 100

    # Enforce the distance field by a threshold field
    setDistanceField(distance_field, meshSize, localSize, 2*meshSize)

end # Local mesh refinement on target cell

function setDistanceField(distance_field, meshSize, localSize, refineRange)

    # Create a threshold field that defines the refinement region
    threshold_field = gmsh.model.mesh.field.add("Threshold")

    gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", localSize)
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", meshSize)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", localSize)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", refineRange)

    # gmsh.model.mesh.field.setNumber(threshold_field, "Sigmoid", true)
    gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)
end

function findNodes(mesh::MESH, region::String, id)
    # region    - face | volume
    # id        - Int or vector of Int 

    nodesFound::Vector{Int32} = zeros(mesh.nv)
    n::Int32 = 0 # Number of nodes found
    
    if region == "face" || region == "Face" # added "Face" to handle variations

        # Go to each surface triangle
        for s in 1:mesh.ns
            # Get the surface id of current triangle
            current_Id::Int32 = mesh.surfaceT[end,s]
            if current_Id in id
                # Nodes of current triangle
                nds = @view mesh.surfaceT[1:3,s]

                # Update number of nodes found
                for nd in nds
                    if nodesFound[nd] < 1   # Only count those who were not found yet
                        n += 1                  # count it
                    end
                end

                # Update the nodes found with desired face ID
                nodesFound[nds] .= 1
            end
        end

        # Prepare the output
        nodes::Vector{Int32} = zeros(n)
        j::Int32 = 0
        for nd in 1:mesh.nv
            if nodesFound[nd] > 0
                j += 1
                nodes[j] = nd
            end
        end

    else # volume
        # not implemented yet
    end

    return nodes
end

# Get the global nodes of a local edge index 'ie'
# of a global element index 'k'
function NodesFromLocalEdge( mesh::MESH, 
                             k, # Global element index
                             ie # Local edge index (1 to 6)
                            )
# Follows the schematics from the GMSH reference manual

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

    return i, j
end

# Find the element (triangle) that contains the node (xq, yq, 0.0)
function findElement2D(mesh::MESH, xq::Float64, yq::Float64)
    # Check what triangle has the point P inside it. Algorithm from
    # https://blackpawn.com/texts/pointinpoly/#:~:text=//%20Compute%20vectors%20v0%20=%20C%20%2D,0)%20&&%20(u%20+%20v%20%3C%201)
    
    tag::Int32 = 0 # Tag of the element that contains P

    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        AC = mesh.p[:, nds[3]] - mesh.p[:, nds[1]]
        AB = mesh.p[:, nds[2]] - mesh.p[:, nds[1]]
        AP = [xq, yq, 0.0] - mesh.p[:, nds[1]]

        ACAC = dot(AC, AC) # 00
        ACAB = dot(AC, AB) # 01
        ACAP = dot(AC, AP) # 02
        ABAB = dot(AB, AB) # 11
        ABAP = dot(AB, AP) # 12

        # Barycentric coordinates
        invDenom = 1.0 / (ACAC*ABAB - ACAB*ACAB)
        u = (ABAB*ACAP - ACAB*ABAP) * invDenom
        v = (ACAC*ABAP - ACAB*ACAP) * invDenom

        # Check if P is inside
        if u >= 0 && v >= 0 && (u + v < 1) # Its inside
            tag = k
            break
        end # Check if P is inside

    end # Search each mesh element (triangle)

    return tag
end # Find the element that contains the coordinate P inside it

# Interpolation over a scattered mesh of nodes
function interp2Dmesh(mesh::MESH, xq::Float64, yq::Float64, T::Vector{Float64})

    # Find the mesh element closest/containing the target node
    tag = findElement2D(mesh, xq, yq)

    nds = @view mesh.t[1:3, tag]
    P1 = mesh.p[:, nds[1]]
    P2 = mesh.p[:, nds[2]]
    P3 = mesh.p[:, nds[3]]

    # Set the Z value to the input values
    if length(T) == mesh.nt # T is expected to be defined over the elements of the mesh
        P1[3] = T[tag]
        P2[3] = T[tag]
        P3[3] = T[tag]
    
    else  # T is defined over the nodes of the mesh
        P1[3] = T[nds[1]]
        P2[3] = T[nds[2]]
        P3[3] = T[nds[3]]
    end

    # Create a plane for interpolation
    AB = P2-P1
    AC = P3-P1
    n = cross(AB, AC)
    a, b, c = n
    d = -dot(n, P1)

    # Get the interpolation of the solution on the mesh element
    zq = -(a*xq + b*yq + d)/c

    return zq
end # Interpolate value over the mesh

# Find the element (triangle) that contains the node (xq, yq, 0.0)
function findElement3D(mesh::MESH, xq::Float64, yq::Float64, zq::Float64)
    # Check what tetrahedron has the point P inside it. Algorithm extended from 2D to 3D, based on
    # https://blackpawn.com/texts/pointinpoly/#:~:text=//%20Compute%20vectors%20v0%20=%20C%20%2D,0)%20&&%20(u%20+%20v%20%3C%201)
    
    tag::Int32 = 0 # Tag of the element that contains P

    for k in 1:mesh.nt
        nds = @view mesh.t[:, k] # Nodes of the tetrahedron

        A = mesh.p[:, nds[1]]
        B = mesh.p[:, nds[2]]
        C = mesh.p[:, nds[3]]
        D = mesh.p[:, nds[4]]

        AB = B - A
        AC = C - A
        AD = D - A
        AP = [xq, yq, zq] - A

        M = [AB AC AD] # 3x3 matrix
        detM = dot(AB, cross(AC, AD))

        if abs(detM) < 1e-12 # Degenerate tetrahedron
            @warn "findElement3D(), skipping degenerate tetrahedron"
            continue
        end

        invDet = 1.0 / detM

        # Barycentric coordinates
        v = dot(cross(AC, AD), AP) * invDet
        u = dot(cross(AP, AD), AB) * invDet  # Note: cross(AB, AP) would give negative
        w = dot(cross(AB, AC), AP) * invDet

        # Check if point is inside tetrahedron
        # All barycentric coordinates should be >= 0 and sum <= 1
        if v >= 0 && u >= 0 && w >= 0 && (v + u + w) <= 1 # Its inside
           tag = k
           break
        end # Check if P is inside

    end # Search each mesh element (triangle)

    return tag
end # Find the element that contains the coordinate P inside it

# Interpolation over a scattered mesh of nodes
function interp3Dmesh(mesh::MESH
                    , xq::Float64, yq::Float64, zq::Float64 # Input coordinate query
                    , T::Vector{Float64} # Value of the result over the mesh  
                     )
# Note: T can be on either the nodes or the elements
# Note on performance: This interpolation method 
# is about 1000x faster than calculating the linear basis function of each node

    # Find the mesh element closest/containing the target node
    tag = findElement3D(mesh, xq, yq, zq)

    nds = @view mesh.t[1:4, tag]
    P1 = [mesh.p[:, nds[1]]; 0.0] # x,y,z, W
    P2 = [mesh.p[:, nds[2]]; 0.0] # ...
    P3 = [mesh.p[:, nds[3]]; 0.0]
    P4 = [mesh.p[:, nds[4]]; 0.0]

    # Set the W value to the result
    if length(T) == mesh.nt # T is expected to be defined over the elements of the mesh
        P1[4] = T[tag]
        P2[4] = T[tag]
        P3[4] = T[tag]
        P4[4] = T[tag]

    else  # T is defined over the nodes of the mesh
        P1[4] = T[nds[1]]
        P2[4] = T[nds[2]]
        P3[4] = T[nds[3]]
        P4[4] = T[nds[4]]

    end

    # Create a plane for interpolation
    AB = P2-P1
    AC = P3-P1
    AD = P4-P1

    M = [[P1 P2 P3 P4]' ones(4)] # 4 by 5 matrix, underdetermined

    # Solution is 1D nullspace of M
    NS = nullspace(M)
    if size(NS, 2) < 1
        @error "interp3Dmesh() , No nontrivial solution found for the hyperplane interpolation"
    end

    # Take the first basis vector as the normal to the hyperplane
    n = NS[:, 1]  # ax + by + cz + dw + e = 0 <-> dot(n, [x,y,z,w]) + e = 0
    n ./= norm(n)

    # Get the interpolation of the solution on the mesh element
    w = -(n[1]*xq + n[2]*yq + n[3]*zq + n[5])/n[4]
    return w

end # Interpolate value over the mesh

# View the mesh using Makie
function viewMesh(mesh::MESH)
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect=:data, title="")
    
    # Convert surface triangles to Makie format
    faces = [GLMakie.GLTriangleFace(mesh.surfaceT[1,i], 
                                    mesh.surfaceT[2,i], 
                                    mesh.surfaceT[3,i]) for i in 1:size(mesh.surfaceT,2)]
    
    # Create mesh plot using surface triangles
    mesh!(ax, mesh.p', faces,
            color=:lightblue,
            transparency=true,
            alpha=0.3)

    wait(display(fig))
end # View the mesh using Makie