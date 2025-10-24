#=
    A wrapper to simplify gmsh interactions. 
    It helps users with:
        . Create simple geometries
        . Combine them into a unified model
        . Create a local refinement field using a single scalar value (localSize)
        . Convert the gmsh mesh to a more user friendly MESH() struct
        . Make a bounding shell for open boundary problems
        . Manage surface triangles and surface ID's
=#

using Gmsh, LinearAlgebra

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
    ne::Int32           # Number of surface elements
    nInside::Int32      # Number of elements for the inside cells
    nInsideNodes::Int32 # Numberf of nodes for the inside cells

    shell_id # Id of the surface boundaries of the bounding shell

    # Constructor
    MESH() = new()
end

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

    setDistanceField(distance_field, meshSize, localSize)

end # Local mesh refinement on target cell

function setDistanceField(distance_field, meshSize, localSize)

    # Create a threshold field that defines the refinement region
    threshold_field = gmsh.model.mesh.field.add("Threshold")

    gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", localSize)
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", meshSize)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", max(0.1*meshSize,10*localSize))

    gmsh.model.mesh.field.setNumber(threshold_field, "Sigmoid", true)
    gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

end


# Add 2d rectangle from its center
function addRectangle(position, dimensions, cells=[])

    x = position[1] - dimensions[1]/2
    y = position[2] - dimensions[2]/2

    id = gmsh.model.occ.addRectangle(x, y, position[3], 
                                    dimensions[1], dimensions[2])

    push!(cells,(2, id))

    # Sync kernel before exiting
    gmsh.model.occ.synchronize()
    
    return id
end

# Add 2d disk (circular)
function addDisk(position, radius, cells=[])
    id = gmsh.model.occ.addDisk(position[1],
                                position[2],
                                position[3],
                                radius,
                                radius)
    push!(cells,(2, id))
    
    # Sync kernel before exiting
    gmsh.model.occ.synchronize()
    
    return id
end

# Add a cuboid based on its center
function addCuboid(position,dimensions,cells=[],updateCells::Bool=false)
    #=
        Makes a cuboid based on its centroid position
        Updates the cells list in case this cuboid is not meant to be
        the container of the simulation
    =#

    r = position - dimensions/2
    box = gmsh.model.occ.addBox(r[1], r[2], r[3], dimensions[1], dimensions[2], dimensions[3])
    
    if updateCells
        cells = append!(cells,[(3,box)])
    end

    # Sync kernel before exiting
    gmsh.model.occ.synchronize()

    return box
end # Make a cuboid based on its center

# Add a sphere
function addSphere(position,radius,cells=[],updateCells::Bool=false)
    #=
        Inputs:
            Position vector
            radius value
            cells <- a list of volumes (cells) that are inside a container
    =# 
    
    # Add a sphere to the current model
    sphere = gmsh.model.occ.addSphere(position[1],position[2],position[3],radius)

    # If sphere is not the container
    if updateCells    
        cells = append!(cells,[(3,sphere)])
    end

    # Sync kernel before exiting
    gmsh.model.occ.synchronize()
    
    return sphere
end # Make a sphere

# Create container based on current model surface (3D only)
function makeContainer(scale::Float64=5.0)

    # Get all surface entities
    surfaces = gmsh.model.getEntities(2)

    # Initialize min/max coordinates
    x_min = Inf
    y_min = Inf
    z_min = Inf
    x_max = -Inf
    y_max = -Inf
    z_max = -Inf

    # Find global bounding box of the STL
    for s in surfaces
        bb = gmsh.model.getBoundingBox(s[1], s[2])
        x_min = min(x_min, bb[1])
        y_min = min(y_min, bb[2])
        z_min = min(z_min, bb[3])
        x_max = max(x_max, bb[4])
        y_max = max(y_max, bb[5])
        z_max = max(z_max, bb[6])
    end

    Lx = x_max-x_min
    Ly = y_max-y_min
    Lz = z_max-z_min
    L = maximum([Lx,Ly,Lz])

    # Container position and dimensions
    center = [(x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2]
    dimensions = scale*L

    box = addSphere(center,dimensions)

    return box
end # Create container based on current model surface

# Import cad geometry file
function importCAD(file::String, cells, setContainer::Bool=true, scale::Float64=5.0)
    #= 
        Import cad geometry file and create a container
        if there is none
    =#

    # Import CAD file (BREP, STEP or IGES)
    gmsh.model.occ.importShapes(file)
    volume = gmsh.model.occ.healShapes()
    _,volume = volume[1]

    gmsh.model.occ.synchronize()
    cells = append!(cells,[(3,volume)])

    # Make a container for the stl file
    if setContainer
        box = makeContainer(scale)
        return box
    end

end # Import cad geometry file

# Boolean subtract the geometries
function occCut(cell1, cell2)

    outDimTags = gmsh.model.occ.cut(cell1, cell2)[1]
    gmsh.model.occ.synchronize()

    return outDimTags
end

function findNodes(mesh::MESH, region::String, id)
    # region    - face | volume
    # id        - Int or vector of Int 

    nodesFound::Vector{Int32} = zeros(mesh.nv)
    n::Int32 = 0 # Number of nodes found
    
    if region == "face" || region == "Face" # added "Face" to handle variations

        # Go to each surface triangle
        for s in 1:mesh.ne
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

# Unify the volumes and get the surface IDs of the bounding shell
function unifyModel(cells, box=-1)
    # Works for both 2D and 3D geometries
    
    dim = cells[1][1] # Get dimension of the model

    inputBox::Bool = box > 0 # Check if a bounding cell was provided

    # Fragment to make a unified geometry
    if box > 0
        fragments, _ = gmsh.model.occ.fragment(vcat(cells,[(dim, box)]), [])
    else
        fragments, _ = gmsh.model.occ.fragment(cells, [])
    end

    gmsh.model.occ.synchronize()

    # Update cell ids
    cells .= box > 0 ? fragments[1:end-1] : fragments[1:end]
    
    # Set the box to the last volume
    box = fragments[end][2]

    # Get bounding shell surface id
    shell_id = gmsh.model.getBoundary([(dim, box)], false, false, false) # (dim, tag)
    shell_id = [s[2] for s in shell_id] # tag

    # Volume surface ids
    if inputBox
        internal_surfaces = gmsh.model.getBoundary(cells, false, false, false) # (dim, tag)
    
    else # The last cell is actually a box
        internal_surfaces = gmsh.model.getBoundary(cells[1:end-1], false, false, false) # (dim, tag)
    
    end
    internal_surfaces = [s[2] for s in internal_surfaces] # tag


    shell_id = setdiff(shell_id, internal_surfaces) # Only the outer surfaces

    # Must separate return logic to not mess the logic of older code
    # otherwise shell_id will be a tuple in some cases
    if inputBox
        return shell_id
    else
        return shell_id, box
    end

end # Unify the volumes

# 3D
function Mesh(cells, meshSize=0, localSize=0, saveMesh::Bool=false)
    #=
        Generates a 3d tetrahedral mesh considering that the model is made of 
        1 container and every other volume beyond the container is listed in the 'cells'

        Inputs
            cells       -> geometries that are inside the container
            meshSize    -> overall target mesh size
            localSize   -> Local mesh refinement
    =#

    # Make a mesh object
    mesh = MESH()

    # Get the volume IDs of the cells inside the container
    volumeID = []
    if !isempty(cells)
        for i in cells
            append!(volumeID,i[2])
        end
    end

    # >> Mesh settings
    
    # Set local mesh size
    if localSize>0 && !isempty(cells)
        refineCell(cells,localSize,meshSize) # Set local refinement on the sphere Cell
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
    gmsh.model.mesh.optimize()
    
    # Get all tetrahedral elements (4 - tetrahedrons)
    t_tags, t = gmsh.model.mesh.getElementsByType(4)
    mesh.t = reshape(t,4,Int(size(t,1)/4))
    mesh.nt = size(mesh.t,2)

    # Get node coordinates
    _,p,_ = gmsh.model.mesh.getNodes()
    mesh.p = reshape(p, 3, Int(size(p,1)/3))
    mesh.nv = size(mesh.p,2)

    # Get all surface triangles
    surfaceT_tags, surfaceT = gmsh.model.mesh.getElementsByType(2)
    surfaceT = reshape(surfaceT, 3, Int(size(surfaceT,1)/3))

    # Expand surface triangles to include boundary id
    surfaceT = [surfaceT;UInt.(zeros(1,size(surfaceT,2)))]
    for i in 1:size(surfaceT,2)
        # Get ID of the element of the current surface triangle
        _,_,_, id = gmsh.model.mesh.getElement(surfaceT_tags[i])
        surfaceT[end,i] = id; # Set the surface triangle boundary id
    end

    mesh.surfaceT = surfaceT
    mesh.ne = size(mesh.surfaceT,2)

    # Mesh elements inside the container
    if !isempty(cells)

        mesh.InsideElements = zeros(mesh.nt)
        mesh.nInside = 0
        for k::Int32 in 1:mesh.nt
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

    # Element volumes
    mesh.VE = zeros(mesh.nt)
    for k in 1:mesh.nt
        mesh.VE[k] = elementVolume(mesh.p,mesh.t[:,k])
    end

    # List of all surface triangle normals are areas
    mesh.normal = zeros(3,mesh.ne)
    mesh.AE = zeros(mesh.ne)
    for i in 1:mesh.ne
        nds = @view mesh.surfaceT[1:3,i]
        mesh.normal[:,i] = normalSurface(mesh.p,nds);
        mesh.AE[i] = areaTriangle(mesh.p[1,nds],mesh.p[2,nds],mesh.p[3,nds])
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
        _,_,_, id = gmsh.model.mesh.getElement(edge_tags[i])
        edges[end, i] = id; # Set the boundary edge id
    end

    mesh.surfaceT = edges
    mesh.ne = size(edges, 2)

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

    # Normal to each edge
    mesh.normal = zeros(2, mesh.ne)
    for e in 1:mesh.ne
        nds = @view mesh.surfaceT[1:2, e]
        mesh.normal[:, e] = normalEdge(mesh.p, nds)
    end

    # Save mesh 
    if saveMesh
        save2file("t.txt",mesh.t) # Save connectivity list to a .txt file
        save2file("p.txt",mesh.p)
        save2file("InsideNodes.txt",mesh.InsideNodes)
        save2file("surfaceT.txt",mesh.surfaceT)
        save2file("InsideElements.txt",mesh.InsideElements)
        save2file("VE.txt",mesh.VE)
        # save2file("AE.txt",mesh.AE)
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

# Save variable to .txt
function save2file(fileName,input)
    # Saves matrix to a .txt file
    open(fileName, "w") do io
        for row in eachrow(input)
            println(io, join(row, " , "))  # Space-separated
        end
    end
end # Save matrix to .txt file

# Interpolation over a scattered mesh of nodes
function interp2Dmesh(mesh::MESH, xq::Float64, yq::Float64, T::Vector{Float64})

    # Find closest element to target coordinate
    tag = 0
    r = Inf
    for i in 1:mesh.nt
        nds = @view mesh.t[:,i]
        center = mean(mesh.p[:,nds],2)

        new_r = sqrt((center[1] - xq)^2 + (center[2] - yq)^2)
        if new_r < r
            r = new_r
            tag = i
        end
    end

    nds = @view mesh.t[:, tag]
    P1 = mesh.p[:, nds[1]]
    P2 = mesh.p[:, nds[2]]
    P3 = mesh.p[:, nds[3]]

    # Set the Z value to the temperature
    if length(T) == mesh.nv
        P1[3] = T[nds[1]]
        P2[3] = T[nds[2]]
        P3[3] = T[nds[3]]
    
    else
        P1[3] = T[tag]
        P2[3] = T[tag]
        P3[3] = T[tag]
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
end # Interpolate value over the mesh nodes

# Sets a volume element to each surface element
function surface2volume(mesh::MESH)
    surface2element::Vector{Int32} = zeros(mesh.ne)

    for s in 1:mesh.ne
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

# Calculates the size of the variable and prints if verbose
function memory(var, verbose::Bool=true)
    mem::Float64 = Base.summarysize(var)/(1024*1024)
    
    if verbose
        println("Var size (MB): ", mem)
    end

    return mem
end