#=
    2D viscous fluid simulation

    Gets the static velocity and pressure of a viscous fluid
    flowing around an obstacle, following the stokes equation
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

using GLMakie

# Sort the quadratic mesh vertices and edge midpoints
function sortMeshNodes(mesh::MESH)
    # Vertices must start from 1 to 'nVertices'. Edge midpoints  must 
    # start from 'nVertices'+1 to mesh.nv

    vertices::Vector{Int32} = unique(vec(mesh.t[1:3,:]))
    nVertices = length(vertices)

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
end

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Simulation settings
    velocity::Vector{Float64} = [1.0, 0.0]
    viscosity::Float64 = 1.0

    # Create model
    gmsh.initialize()
    
    # Add an obstacle
    cells = []
    # id = addRectangle([0,0,0], [5, 3], cells)
    id = addDisk([-5,0,0], 1.0, cells)

    # Add a container
    box = addRectangle([0,0,0], [20, 5])
    # box = addDisk([0,0,0], 4)

    # Combine the geometries
    gmsh.model.occ.fragment(vcat(cells,[(2,box)]), [])
    gmsh.model.occ.synchronize()

    # Generate mesh
    mesh::MESH = Mesh2D(cells, meshSize, localSize, 2)

    println("\nNumber of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("Mesh Order: ", mesh.order)

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()

    # Gmsh orders the nodes arbitrarily
    # So I have to re-label the vertices and the midpoints
    vertexID::Vector{Int32}, nVertices::Int32 = sortMeshNodes(mesh)

    # Define the viscosity on the domain
    mu::Vector{Float64} = zeros(mesh.nt) .+ viscosity

    # Define the obstacle as an extremely viscous object to set the internal velocity to zero
    mu[mesh.InsideElements] .= 1e3 * viscosity

    println("Building stiffness matrix")

    # Global Stiffness matrix
    A = spzeros(mesh.nv, mesh.nv)

    # Local stiffness matrix
    Ak::Matrix{Float64} = zeros(36, mesh.nt) # 6 x 6
    temp::Matrix{Float64} = zeros(6,6)

    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]

        temp .= 0
        for i in 1:length(nds)
            Si = quadraticBasis2D(mesh.p, nds, nds[i])
            
            for j in i:length(nds)
                Sj = quadraticBasis2D(mesh.p, nds, nds[j])

                # 6 node quadrature
                aux::Float64 = 0.0
                for n in 1:6
                    dxi::Float64 = Si[2] + 2*Si[4]*mesh.p[1,nds[n]] + Si[5]*mesh.p[2,nds[n]]
                    dyi::Float64 = Si[3] + Si[5]*mesh.p[1,nds[n]] + 2*Si[6]*mesh.p[2,nds[n]] 
                    dxj::Float64 = Sj[2] + 2*Sj[4]*mesh.p[1,nds[n]] + Sj[5]*mesh.p[2,nds[n]]
                    dyj::Float64 = Sj[3] + Sj[5]*mesh.p[1,nds[n]] + 2*Sj[6]*mesh.p[2,nds[n]]
                    aux += dxi*dxj + dyi*dyj
                end 
                aux /= 6

                temp[i,j] = mu[k]*aux*mesh.VE[k]
                temp[j,i] = temp[i,j] # It is symmetric
            end
        end

        Ak[:,k] = vec(temp)

    end # Local stiffness matrix

    # Update sparse global matrix
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            A += sparse(vertexID[mesh.t[i,:]],      # Convert original node ID to sorted node ID
                        vertexID[mesh.t[j,:]],      # Convert original node ID to sorted node ID
                        Ak[n,:], mesh.nv, mesh.nv)
        end
    end

    println("Building divergence matrix")

    # Pressure matrix
    B1::Matrix{Float64} = zeros(nVertices, mesh.nv) # Vertices x Nodes
    B2::Matrix{Float64} = zeros(nVertices, mesh.nv) # Vertices x Nodes
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        for i in 1:3
            a::Float64, b::Float64, c::Float64 = abc(mesh.p, nds[1:3], nds[i])
            for j in 1:6
                S = quadraticBasis2D(mesh.p, nds, nds[j])
                
                # 6 Node quadrature
                b1::Float64 = 0.0
                b2::Float64 = 0.0
                for n in 1:6
                    b1 -= (a + b*mesh.p[1,nds[n]] + c*mesh.p[2,nds[n]])*            # Linear
                          (S[2] + 2*S[4]*mesh.p[1,nds[n]] + S[5]*mesh.p[2,nds[n]])  # Quadratic
                    
                    b2 -= (a + b*mesh.p[1,nds[n]] + c*mesh.p[2,nds[n]])*            # Linear
                          (S[3] + S[5]*mesh.p[1,nds[n]] + 2*S[6]*mesh.p[2,nds[n]])  # Quadratic
                end # 6 node quadrature (quadratic nodes)

                # Convert original node ID to sorted node ID and update matrix
                B1[vertexID[nds[i]], vertexID[nds[j]]] += mesh.VE[k]*b1/6
                B2[vertexID[nds[i]], vertexID[nds[j]]] += mesh.VE[k]*b2/6
            
            end # Quadratic nodes loop
        end # Linear nodes loop
    end # Element loop

    # Full matrix
    LHS =  [A zeros(mesh.nv, mesh.nv) B1'; 
            zeros(mesh.nv, mesh.nv) A B2';
            B1 B2 zeros(nVertices, nVertices)]

    # Apply boundary conditions
    # NOTE! The boundary ID must be checked manually from Gmsh GUI directly

    # In Flow | Curve id: 2
    inFlow::Int32 = 2
    inFlowNodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(1, inFlow)

    # No-slip boundary conditions
    walls::Vector{Int32} = [1, 4, 5]
    wallNodes::Vector{Int32} = Int32[]
    for id in walls
        nodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(1, id)
        append!(wallNodes, nodes)
    end

    # Force velocity on the obstacle to zero
    # println(mesh.InsideNodes)

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
    u::Matrix{Float64} = zeros(mesh.nv, 2)
    u[:, 1] .= UP[1:mesh.nv]
    u[:, 2] .= UP[mesh.nv+1:2*mesh.nv]

    velocityNorm::Vector{Float64} = zeros(mesh.nv)
    for i in 1:mesh.nv
        velocityNorm[i] = sqrt(sum(u[i,:].^2))   
    end

    # Pressure
    p = UP[2*mesh.nv.+(1:nVertices)]

    # Sort the coordinates
    x::Vector{Float64} = zeros(mesh.nv)
    x[vertexID[1:mesh.nv]] .= mesh.p[1,:]

    y::Vector{Float64} = zeros(mesh.nv)
    y[vertexID[1:mesh.nv]] .= mesh.p[2,:]

    println("Printing plots")

    # Plot results
    fig = Figure()

    # Add vector field
    Axis(fig[1, 1], aspect = DataAspect(), title="Velocity field")
    velocity_plot = arrows2d!(x, y, u[:,1], u[:,2], 
                              lengthscale = 0.5,
                              color = velocityNorm,
                              colormap = :thermal)

    Colorbar(fig[2, 1], velocity_plot, 
             label = "Velocity", vertical = false)

    # Add Pressure plot
    xPressure::Vector{Float64} = zeros(nVertices)
    yPressure::Vector{Float64} = zeros(nVertices)

    # Vertices
    vertices::Vector{Int32} = unique(vec(mesh.t[1:3,:]))
    xPressure[vertexID[vertices]] .= mesh.p[1,vertices]
    yPressure[vertexID[vertices]] .= mesh.p[2,vertices]


    ax2 = Axis(fig[1, 2], aspect = DataAspect(), title="Pressure")
    sc2 = scatter!( ax2, 
                    xPressure, yPressure, 
                    color=p, 
                    colormap=:batlow,
                    markersize=20)
    
    Colorbar(fig[2, 2], sc2,
             label = "Pressure", vertical = false)

    wait(display(fig))


end

showGmsh = true
main(2.5, 0.25, showGmsh)