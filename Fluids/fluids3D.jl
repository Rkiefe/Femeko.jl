#=
    3D viscous fluid simulation

    Gets the static velocity and pressure of a viscous fluid
    flowing around an obstacle, following the stokes equation

    The implementation considers mixed-elements. The pressure is defined
    over linear lagrange elements, and the velocity is defined over
    quadratic lagrange elements.
=#

include("../src/Femeko.jl")
using IterativeSolvers

function fluid3D( mesh::MESH
                , velocity::Vector{Float64}
                , mu::Vector{Float64}
                , inFlow
                , walls)
    
    # Apply boundary conditions at ...
    
    # Intake
    inFlowNodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(2, inFlow) # 2 is for 3-node triangles

    # Walls (no-slip)
    wallNodes::Vector{Int32} = Int32[]
    for id in walls
        nodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(2, id, true) # 2 is for 3-node triangles, true -> include the boundary
        append!(wallNodes, nodes)
    end

    # Remove the overlapping nodes
    inFlowNodes = setdiff(inFlowNodes, wallNodes)


    # Sort the quadratic mesh vertices and edge midpoints
        # Vertices must start from 1 to 'nVertices'. Edge midpoints  must 
        # start from 'nVertices'+1 to mesh.nv

    # 1st order nodes
    vertices::Vector{Int32} = unique(mesh.t[1:4, :])
    nVertices::Int32 = length(vertices)

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

    # This way the P1 nodes are from 1 to nVertices, 
    # the P2 nodes are from nVertices+1 to mesh.nv

    # Pre compute the quadratic basis function for the entire mesh
    println("Computing the quadratic basis function on each node of each element")
    S = zeros(10, 10, mesh.nt)
    @time for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        for i in 1:10
            # a, b, c, ... of current node 'i' and element 'k'
            S[:, i, k] = quadraticBasis(mesh, nds, nds[i])
        end

    end # Quadratic basis function for every node and element

    # Quadratic Local Stiffness matrix
    println("Computing the quadratic stiffness matrix...")
    Ak = @time quadraticLocalStiffnessMatrix(mesh, S, mu)
    
    # Global stiffness matrix
    A = spzeros(mesh.nv, mesh.nv)
    n = 0
    for i in 1:10
        iLocal = localNodeID[mesh.t[i, :]]
        for j in 1:10
            jLocal = localNodeID[mesh.t[j, :]]
            n += 1
            A += sparse(iLocal, jLocal, Ak[n, :], mesh.nv, mesh.nv)
        end
    end    

    # Local 3D divergence matrix
    println("Computing the local Divergence xyz matrices...")
    B1k, B2k, B3k = @time localDivergenceMatrix(mesh, S)

    # Global Divergence matrix
    B1 = spzeros(nVertices, mesh.nv) # Vertices by Nodes
    B2 = spzeros(nVertices, mesh.nv) # ...
    B3 = spzeros(nVertices, mesh.nv) # ...

    n = 0
    for i in 1:4
        iLocal = localNodeID[mesh.t[i, :]]
        for j in 1:10
            n += 1

            jLocal = localNodeID[mesh.t[j, :]]
            B1 += sparse(iLocal, jLocal, B1k[n, :], nVertices, mesh.nv)
            B2 += sparse(iLocal, jLocal, B2k[n, :], nVertices, mesh.nv)
            B3 += sparse(iLocal, jLocal, B3k[n, :], nVertices, mesh.nv)
        end
    end # Update global divergence matrix

    # Full matrix
    LHS =  [A spzeros(mesh.nv, mesh.nv) spzeros(mesh.nv, mesh.nv) B1'; 
            spzeros(mesh.nv, mesh.nv) A spzeros(mesh.nv, mesh.nv) B2';
            spzeros(mesh.nv, mesh.nv) spzeros(mesh.nv, mesh.nv) A B3';
            B1 B2 B3 spzeros(nVertices, nVertices)]

    # Schematic of the 2D equation to help build the 3D case
    # |A  0  B1| |ux| = f
    # |0  A  B2| |uy|
    # |B1 B2 0 | |p |

    # Nodes with boundary conditions
    fixed = [

            # Velocity field
             localNodeID[inFlowNodes];              # .x              
             localNodeID[wallNodes];                # .x

             mesh.nv .+ localNodeID[inFlowNodes];   # .y
             mesh.nv .+ localNodeID[wallNodes];     # .y

             2*mesh.nv .+ localNodeID[inFlowNodes]; # .z
             2*mesh.nv .+ localNodeID[wallNodes];   # .z
            ]

    DOF = 3*mesh.nv + nVertices # 3D vector field + Scalar field
    free = setdiff(1:DOF, fixed)

    # Define boundary conditions 
    gD::Vector{Float64} = zeros(DOF)

    gD[localNodeID[inFlowNodes]]              .= velocity[1]
    gD[mesh.nv   .+ localNodeID[inFlowNodes]] .= velocity[2]
    gD[2*mesh.nv .+ localNodeID[inFlowNodes]] .= velocity[3]

    # Apply the boundary conditions
    RHS = -LHS[free, fixed]*gD[fixed]

    # Solution (Velocity xyz then pressure)
    UP::Vector{Float64} = zeros(DOF)
    UP[fixed] .= gD[fixed] # Dirichlet conditions
    
    # Solve for the velocity and pressure (u and p)
    print("Solving matrix equation for the fluid velocity and pressure with the Conjugate Gradient method...")
    # UP[free] = LHS[free, free]\RHS
    UP[free] = cg(LHS[free, free], RHS)
    println(" Done")

    # Velocity (defined on the local node IDs)
    u::Matrix{Float64} = zeros(3, mesh.nv)
    u[1, :] .= UP[1:mesh.nv]
    u[2, :] .= UP[mesh.nv+1:2*mesh.nv]
    u[3, :] .= UP[2*mesh.nv+1:3*mesh.nv]

    # Pressure (defined on the local node IDs)
    p::Vector{Float64} = UP[3*mesh.nv .+ (1:nVertices)]

    # Remap the output to the original Gmsh node labels
    println("Remapping the velocity and pressure field to the original node labels")
    uOrdered::Matrix{Float64} = zeros(3, mesh.nv) 
    for i in 1:mesh.nv
        uOrdered[:, i] .= u[:, localNodeID[i]]
    end

    pOrdered::Vector{Float64} = zeros(mesh.nv)
    for nd in vertices # Global node labels of 1st order element
        i = localNodeID[nd] # Local node label of current vertex
        pOrdered[nd] = p[i]
    end
    println("Fluid simulation concluded")

    return uOrdered, pOrdered, vertices
end

