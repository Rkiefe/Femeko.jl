#=
    2D viscous fluid simulation

    Gets the static velocity and pressure of a viscous fluid
    flowing around an obstacle, following the stokes equation

    The implementation considers mixed-elements. The pressure is defined
    over linear lagrange elements, and the velocity is defined over
    quadratic lagrange elements.
=#

include("../src/Femeko.jl")

# Run the matrix assemblies and solve the system for pressure and velocity
function fluid2D(mesh::MESH, velocity::Vector{Float64}, mu::Vector{Float64}, inFlow, walls)

    # Gmsh orders the nodes arbitrarily
    # So I have to re-label the vertices and the midpoints
    # because the velocity and pressure are solved in different basis
    vertexID::Vector{Int32}, nVertices::Int32 = sortMeshNodes2D(mesh)

    println("Building stiffness matrix")

    # Local Stiffness matrix
    Ak::Matrix{Float64} = quadraticLocalStiffnessMatrix2D(mesh)

    # Global Stiffness matrix
    A = spzeros(mesh.nv, mesh.nv)
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            A += sparse(  vertexID[mesh.t[i,:]]      # Convert original node ID to sorted node ID
                        , vertexID[mesh.t[j,:]]      # Convert original node ID to sorted node ID
                        , Ak[n,:].*mu
                        , mesh.nv, mesh.nv)
        end
    end

    println("Building divergence matrix")

    # Pressure matrix
    B1::Matrix{Float64}, B2::Matrix{Float64} = divergenceMatrix2D(mesh, 
                                                                  vertexID, 
                                                                  nVertices)
    # Full matrix
    LHS =  [A spzeros(mesh.nv, mesh.nv) B1'; 
            spzeros(mesh.nv, mesh.nv) A B2';
            B1 B2 spzeros(nVertices, nVertices)]

    # Boundary conditions

    # Fluid intake
    inFlowNodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(1, inFlow)

    # No-slip boundary conditions
    wallNodes::Vector{Int32} = Int32[]
    for id in walls
        nodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(1, id, true)
        append!(wallNodes, nodes)
    end

    # Nodes with boundary conditions
    fixed = [
             vertexID[inFlowNodes]; 
             vertexID[wallNodes];
             mesh.nv .+ vertexID[inFlowNodes];
             mesh.nv .+ vertexID[wallNodes];
            ]

    free = setdiff(1:(2*mesh.nv + nVertices), fixed)

    gD::Vector{Float64} = zeros(2*mesh.nv + nVertices)
    gD[vertexID[inFlowNodes]] .= velocity[1]
    gD[mesh.nv .+ vertexID[inFlowNodes]] .= velocity[2]

    RHS = -LHS[free,fixed]*gD[fixed]

    println("Solving matrix equation")

    # Velocity and pressure (u and p)
    UP::Vector{Float64} = zeros(2*mesh.nv + nVertices)
    UP[fixed] .= gD[fixed]
    UP[free] = LHS[free,free]\RHS    
    
    # Velocity
    u::Matrix{Float64} = zeros(2, mesh.nv)
    u[1, :] .= UP[1:mesh.nv]
    u[2, :] .= UP[mesh.nv+1:2*mesh.nv]

    # Pressure
    p::Vector{Float64} = UP[2*mesh.nv.+(1:nVertices)]

    # Remap the output to the original Gmsh node labels
    println("Remapping the velocity and pressure field to the original node labels")
    uOrdered::Matrix{Float64} = zeros(2, mesh.nv) 
    for i in 1:mesh.nv
        uOrdered[:, i] .= u[:, vertexID[i]]
    end

    # Get the 1st order mesh nodes (vertices)
    vertices::Vector{Int32} = unique(vec(mesh.t[1:3,:]))

    pOrdered::Vector{Float64} = zeros(mesh.nv)
    for nd in vertices # Global node labels
        i = vertexID[nd] # Local node label of current vertex
        pOrdered[nd] = p[i]
    end
    println("Fluid simulation concluded")

    return uOrdered, pOrdered, vertexID, nVertices, vertices
end

