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

    # Setup
    viscosity = 1.0                   # Fluid viscosity
    velocity::Vector{Float64} = [1.0, 
                                 0.0,
                                 0.0] # Inflow

    # Create 3D model
    gmsh.initialize()
    cells = [] # Store the obstacles cell ID
    
    addSphere([0,0,0], 1.0, cells)              # Add obstacle
    box = addCuboid([0,0,0], [5.0, 20.0, 5.0])  # Add tube

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
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()

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

    # # Target first node
    # nds = @view mesh.t[:, 1] # All nodes of element
    # S = quadraticBasis(mesh, nds, nds[1])
    # u = S[1] + S[2]*xt + S[3]*yt + S[4]*zt + 
    # S[5]*xt^2 + S[6]*xt*yt + S[7]*xt*zt + 
    # S[8]*yt^2 + S[9]*yt*zt + S[10]*zt^2

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
    mu[mesh.InsideElements] .= 1e3*viscosity

    # Quadratic Stiffness matrix (all nodes)
    A = zeros(mesh.nv, mesh.nv)
    temp = zeros(10, 10) # Local element wise stiffness matrix

    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        for i in 1:10
            Si = @view S[:, i, k]
            
            for j in 1:10
                Sj = @view S[:, j, k]

                # 10 node quadrature
                aux::Float64 = 0.0
                for n in 1:10

                    x::Float64 = mesh.p[1, nds[n]] 
                    y::Float64 = mesh.p[2, nds[n]] 
                    z::Float64 = mesh.p[3, nds[n]]
                    
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
                temp[j,i] = temp[i,j] # Stiffness matrix is symmetric
            end
        end # Local element wise stiffness matrix

        # Update global stiffness matrix
        A[nds, nds] += temp

    end # Local stiffness matrix



end # main()



meshSize = 0.0
localSize = 0.0
showGmsh = false

main(meshSize, localSize, showGmsh)