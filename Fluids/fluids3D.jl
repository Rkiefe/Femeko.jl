#=
    3D viscous fluid simulation

    Gets the static velocity and pressure of a viscous fluid
    flowing around an obstacle, following the stokes equation

    The implementation considers mixed-elements. The pressure is defined
    over linear lagrange elements, and the velocity is defined over
    quadratic lagrange elements.
=#


include("../src/Femeko.jl")
using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Setup
    viscosity = 1.0                   # Fluid viscosity
    velocity::Vector{Float64} = [0.0, 
                                 1.0,
                                 0.0] # Inflow

    # Boundary surface IDs 
    inFlow = 4
    walls = [1, 2]
    outFlow = 3

    # Create 3D model
    gmsh.initialize()
    cells = [] # Store the obstacles cell ID
    
    addSphere([0,-2.5,0], 0.5, cells)              # Add obstacle
    box = addCylinder([0,-5,0], [0, 10, 0], 2.0)  # Add tube

    # Unify the volumes for a single geometry
    _, box = unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize, false, 2)

    println("\nNumber of elements: ", mesh.nt)
    println("Number of vertices: ", mesh.nv)
    println("Number of surface elements ", mesh.ns)
    println("Number of edges: ", mesh.ne)
    println("")

    if showGmsh
        # gmsh.option.setNumber("Mesh.Clip", 1)
        # gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        # gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end

    # Apply boundary conditions at ...
    # Fluid intake
    inFlowNodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(2, inFlow) # 2 is for 3-node triangles

    # Walls (no-slip)
    wallNodes::Vector{Int32} = Int32[]
    for id in walls
        nodes::Vector{Int32}, _, _ = gmsh.model.mesh.getNodes(2, id, true) # 2 is for 3-node triangles, true -> include the boundary
        append!(wallNodes, nodes)
    end

    # Remove the overlapping nodes
    inFlowNodes = setdiff(inFlowNodes, wallNodes)

    gmsh.finalize() # No more need for Gmsh

    # # Test the wallNodes
    # fig, ax, _ = scatter(mesh.p[1,:], mesh.p[2,:], mesh.p[3,:])
    # scatter!(ax, mesh.p[1, wallNodes], mesh.p[2, wallNodes], mesh.p[3, wallNodes])
    # wait(display(fig))
    # # return

    # # Test the inFlowNodes
    # fig, ax, _ = scatter(mesh.p[1,:], mesh.p[2,:], mesh.p[3,:])
    # scatter!(ax, mesh.p[1, inFlowNodes], mesh.p[2, inFlowNodes], mesh.p[3, inFlowNodes])
    # wait(display(fig))
    # return



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

    # This way the P1 nodes are from 1 to nVertices, 
    # the P2 nodes are from nVertices+1 to mesh.nv

    # Pre compute the quadratic basis function for the entire mesh
    S = zeros(10, 10, mesh.nt)
    @time for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        for i in 1:10
            # a, b, c, ... of current node and element
            S[:, i, k] = quadraticBasis(mesh, nds, nds[i])
        end

    end # Quadratic basis function for every node and element

    # Define the viscosity on each element of the mesh
    mu::Vector{Float64} = zeros(mesh.nt) .+ viscosity
    # mu[mesh.InsideElements] .= 1e3*viscosity

    # Quadratic Stiffness matrix (all nodes)
    A = spzeros(mesh.nv, mesh.nv)
    temp = zeros(10, 10)     # Local element wise stiffness matrix

    @time for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        for i in 1:10
            Si = @view S[:, i, k]
            
            for j in 1:10
                Sj = @view S[:, j, k]

                # 10 node quadrature
                aux::Float64 = 0.0
                for n in 1:10

                    x = mesh.p[1, nds[n]] 
                    y = mesh.p[2, nds[n]] 
                    z = mesh.p[3, nds[n]]
                    
                    dxi::Float64 = Si[2] + 2*Si[5] *x + Si[6]*y + Si[7]*z
                    dyi::Float64 = Si[3] + 2*Si[8] *y + Si[6]*x + Si[9]*z
                    dzi::Float64 = Si[4] + 2*Si[10]*z + Si[7]*x + Si[9]*y

                    dxj::Float64 = Sj[2] + 2*Sj[5] *x + Sj[6]*y + Sj[7]*z
                    dyj::Float64 = Sj[3] + 2*Sj[8] *y + Sj[6]*x + Sj[9]*z
                    dzj::Float64 = Sj[4] + 2*Sj[10]*z + Sj[7]*x + Sj[9]*y

                    aux += dxi*dxj + dyi*dyj + dzi*dzj
                end 
                aux /= 10

                temp[i,j] = mesh.VE[k]*mu[k]*aux
                # temp[j,i] = temp[i,j] # Stiffness matrix is symmetric
            end
        end # Local element wise stiffness matrix

        # Update global stiffness matrix
        A[localNodeID[nds], localNodeID[nds]] += temp

    end # Local stiffness matrix

    # 3D Divergence matrix
        
    # Global Divergence matrix
    B1 = zeros(nVertices, mesh.nv) # Vertices by Nodes
    B2 = zeros(nVertices, mesh.nv) # ...
    B3 = zeros(nVertices, mesh.nv) # ...

    # Local Divergence matrix
    B1temp::Matrix{Float64} = zeros(4, 10) # Element wise matrix
    B2temp::Matrix{Float64} = zeros(4, 10) # ...
    B3temp::Matrix{Float64} = zeros(4, 10) # ...

    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        for i in 1:4
            
            a, b, c, d = abcd(mesh.p, nds[1:4], nds[i]) # Linear basis function
            
            for j in 1:10
                Sj = @view S[:, j, k]
                
                # 10 Node quadrature
                b1::Float64 = 0.0
                b2::Float64 = 0.0
                b3::Float64 = 0.0
                for n in 1:10

                    x = mesh.p[1, nds[n]]
                    y = mesh.p[2, nds[n]]
                    z = mesh.p[3, nds[n]]

                    b1 -= (a + b*x + c*y + d*z)*                    # Linear
                          (Sj[2] + 2*Sj[5] *x + Sj[6]*y + Sj[7]*z)  # Quadratic
                    
                    b2 -= (a + b*x + c*y + d*z)*                    # Linear
                          (Sj[3] + 2*Sj[8] *y + Sj[6]*x + Sj[9]*z)  # Quadratic

                    b3 -= (a + b*x + c*y + d*z)*                    # Linear
                          (Sj[4] + 2*Sj[10]*z + Sj[7]*x + Sj[9]*y)  # Quadratic


                end # 10 node quadrature (quadratic nodes)

                # Element wise divergence matrix
                B1temp[i, j] = mesh.VE[k]*b1/10
                B2temp[i, j] = mesh.VE[k]*b2/10
                B3temp[i, j] = mesh.VE[k]*b3/10

            end # Quadratic nodes loop
        end # Linear nodes loop

        # Update global divergence matrix
        B1[localNodeID[nds[1:4]], localNodeID[nds[1:10]]] += B1temp
        B2[localNodeID[nds[1:4]], localNodeID[nds[1:10]]] += B2temp
        B3[localNodeID[nds[1:4]], localNodeID[nds[1:10]]] += B3temp

    end # Local divergence matrix


    # # Full matrix
    LHS =  [A zeros(mesh.nv, mesh.nv) zeros(mesh.nv, mesh.nv) B1'; 
            zeros(mesh.nv, mesh.nv) A zeros(mesh.nv, mesh.nv) B2';
            zeros(mesh.nv, mesh.nv) zeros(mesh.nv, mesh.nv) A B3';
            B1 B2 B3 zeros(nVertices, nVertices)]

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

    # Solve for the velocity and pressure (u and p)
    UP::Vector{Float64} = zeros(DOF)
    UP[fixed] .= gD[fixed]
    UP[free] = LHS[free, free]\RHS    

    # Velocity (defined on the local node IDs)
    u::Matrix{Float64} = zeros(3, mesh.nv)
    u[1, :] .= UP[1:mesh.nv]
    u[2, :] .= UP[mesh.nv+1:2*mesh.nv]
    u[3, :] .= UP[2*mesh.nv+1:3*mesh.nv]

    # Pressure (defined on the local node IDs)
    p::Vector{Float64} = UP[3*mesh.nv .+ (1:nVertices)]

    # Norm of velocity (defined on the local node IDs)
    uNorm::Vector{Float64} = zeros(mesh.nv)
    for i in 1:mesh.nv
        uNorm[i] = norm(u[:, i])
    end

    # Convert x,y,z coordinates to local node IDs

    # Plot result
    println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1,1], aspect=:data
                , title="Velocity field"
              )

    # Velocity vector field
    x = zeros(mesh.nv)
    y = zeros(mesh.nv)
    z = zeros(mesh.nv)
    
    x[localNodeID[1:mesh.nv]] .= mesh.p[1, :]
    y[localNodeID[1:mesh.nv]] .= mesh.p[2, :]
    z[localNodeID[1:mesh.nv]] .= mesh.p[3, :]

    graph = arrows3d!(  ax
                        , x, y, z
                        , u[1, :]
                        , u[2, :]
                        , u[3, :]
                        , color = uNorm
                        , lengthscale = 0.5
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )
    Colorbar(fig[2, 1], graph, 
             label = "Velocity", vertical = false)

    # Pressure plot
    ax2 = Axis3(fig[1,2], aspect=:data
                , title="Velocity field"
              )
    x = zeros(nVertices)
    y = zeros(nVertices)
    z = zeros(nVertices)

    x[localNodeID[vertices]] .= mesh.p[1, vertices]
    y[localNodeID[vertices]] .= mesh.p[2, vertices]
    z[localNodeID[vertices]] .= mesh.p[3, vertices]

    graph = scatter!(ax2, x, y, z
                     , color = p
                     , colormap = :CMRmap  # :CMRmap :viridis :redsblues :turbo :rainbow
                    )
    
    Colorbar(fig[2, 2], graph, 
             label = "Pressure", vertical = false)
    
    wait(display(fig))

end # main()

meshSize = 0.0
localSize = 0.0
showGmsh = false

main(meshSize, localSize, showGmsh)