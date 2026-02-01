
# Missing the outflow boundary condition

#  Full 2D heat equation with convection to a passing fluid

include("../Fluids/fluids2D.jl")
using GLMakie

# Mesh settings
meshSize = 5.0
localSize = 0.2
showGmsh = false


function quadraticConvectionMatrix2D(mesh::MESH, u::Matrix{Float64})

    # Quadratic shape function times its gradient -> Order 3 integrand
    weights::Vector{Float64}, 
    points::Matrix{Float64} = GaussQuadrature2D(3)

    # Local convection matrix
    Clocal = zeros(36, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]
        vertices = mesh.p[:, nds[1:3]]  # Triangle vertices
        
        # Jacobian for affine transformation from reference to physical triangle
        J11 = vertices[1, 2] - vertices[1, 1]
        J12 = vertices[1, 3] - vertices[1, 1]
        J21 = vertices[2, 2] - vertices[2, 1]
        J22 = vertices[2, 3] - vertices[2, 1]
        
        detJ = abs(J11 * J22 - J12 * J21)

        # Precompute basis coefficients for all 6 nodes
        S::Matrix{Float64} = zeros(6 ,6)  
        for i in 1:6
            S[:, i] = quadraticBasis2D(mesh.p, nds, nds[i])
        end
        
        # Initialize element-wise matrix
        Ck = zeros(6, 6)

        # Loop over quadrature points
        for q in 1:length(weights)
            xi, eta = points[q, 1:2]
            
            # Transform from reference to physical coordinates
            x = vertices[1, 1] + J11 * xi + J12 * eta
            y = vertices[2, 1] + J21 * xi + J22 * eta
            
            # Evaluate on the current quadrature point the
                #  2nd order Lagrange shape function
                #  The velocity field
                #  The gradient of the 2nd order Lagrange shape function
            
            phi::Vector{Float64} = zeros(6)         # Shape function 
            ux::Float64 = 0.0; uy::Float64 = 0.0    # Velocity field
            gradPhi::Matrix{Float64} = zeros(2, 6)  # Gradient of shape function
            
            for i in 1:6 # All nodes of the element
            
                # Basis function on the quadrature point
                phi[i] = S[1, i] + S[2, i]*x + S[3, i]*y + S[4, i]*x^2 + S[5, i]*x*y + S[6, i]*y^2
                
                # Velocity on the quadrature point
                ux += u[1, nds[i]]*phi[i]
                uy += u[2, nds[i]]*phi[i]

                # Gradient of basis function on quadrature point
                gradPhi[1, i] = S[2, i] + 2*S[4, i]*x + S[5, i]*y
                gradPhi[2, i] = S[3, i] + 2*S[6, i]*y + S[5, i]*x
            end

            # Accumulate to local matrix
            w = weights[q] * detJ
            for i in 1:6
                grad_i = gradPhi[1, i] * ux + gradPhi[2, i] * uy
                for j in 1:6
                    Ck[i, j] += w * grad_i * phi[j]
                end
            end
        end
        
        Clocal[:, k] = vec(Ck)
    end
    
    # Sparse global convection matrix
    C = spzeros(mesh.nv, mesh.nv)
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            C += sparse(  mesh.t[i, :]
                        , mesh.t[j, :]
                        , Clocal[n, :]
                        , mesh.nv, mesh.nv)
        end
    end

    return C
end

mutable struct DATA

    # Mass density (g/cm3)
    density::Float64

    # Thermal conductivity W/(m K)
    k::Float64

    # Specific heat J/(kg K)
    Cp::Float64

    # Fluid viscosity mPa s or kg/(m s)
    mu::Float64

    # Constructor
    DATA() = new(1.0, 1.0, 1.0, 1.0)
end

function main(meshSize=0.0, localSize=0.0, showGmsh=false)
    gmsh.initialize()

    # Simulation settings
    velocity::Vector{Float64} = [1.0, 0.0] # In-flow velocity
    viscosity::Float64 = 1.0
    timeStep::Float64 = 1e-3
    totalTime::Float64 = 1.0
    maxSteps::Int32 = floor(totalTime/timeStep) + 1

    # List of materials
    materialProperties = Dict("blank" => DATA(),
                              "water" => DATA())

    materialProperties["blank"].mu = 1e3

    # Create model
    cells = []      # Cell IDs (dim, tag)
    cellLabels = [] # Tag of cell property

    # Add an obstacle
    # id = addDisk([0,0,0], 1.0, cells)
    id = addRectangle([0,0,0], [1.0, 1.0], cells)
    push!(cellLabels, "blank")
    
    box = addRectangle([0,0,0], [8, 8], cells) # Add tube
    push!(cellLabels, "water")

    # Combine model to create a conforming mesh
    unifyModel([(2, id)], box)

    # Generate 2nd order mesh
    # extendLocalRefinement()
    mesh::MESH = Mesh2D([(2, id)], meshSize, localSize, 2)

    println("\nNumber of elements ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of nodes ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ns)
    println("Mesh Order: ", mesh.order, "\n")

    # In Flow | Curve id: 2
    inFlow::Int32 = 2

    # Walls
    walls::Vector{Int32} = [1, 4, 5]

    # Initial temperature
    T::Vector{Float64} = zeros(mesh.nv) .+ 0.0
    T[mesh.InsideNodes] .= 10.0

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()
 
    # Define the viscosity and the diffusivity on the domain
    mu::Vector{Float64} = zeros(mesh.nt)
    epsi::Vector{Float64} = zeros(mesh.nt)

    for i in 1:length(cells)

        id = cells[i][2] # Cell ID
        key = cellLabels[i] # Get the data set of current cell ID

        println("id: ", id, " ; key: ", key)

        # Find all elements of current cell ID
        elements = findall(x -> x==id, mesh.elementID)

        # Update viscosity value on this cell elements
        mu[elements] .= materialProperties[key].mu
        epsi[elements] .= materialProperties[key].k/(materialProperties[key].Cp * materialProperties[key].density)

    end

    # Run fluid simulation
    u::Matrix{Float64}, # 2 by nv
    p::Vector{Float64}, # size = nVertices
    vertexID::Vector{Int32},
    nVertices::Int32,
    vertices::Vector{Int32} = fluid2D(mesh, 
                                      velocity, # Intake fluid velocity
                                      mu,       # Viscosity
                                      inFlow,
                                      walls)    # Boundary IDs

    velocityNorm::Vector{Float64} = zeros(mesh.nv)
    for i in 1:mesh.nv
        velocityNorm[i] = norm(u[:, i])   
    end

    # # Plot velocity field and pressure
    # println("Generating plots...")
    # fig = Figure()
    # ax = Axis(fig[1, 1], aspect = DataAspect(), title="Velocity field")
    # velocity_plot = arrows2d!(  ax
    #                           , mesh.p[1, :]
    #                           , mesh.p[2, :]
    #                           , u[1, :]
    #                           , u[2, :]
    #                           , lengthscale = 0.5
    #                           , color = velocityNorm
    #                           , colormap = :thermal)

    # Colorbar(  fig[1, 2], velocity_plot 
    #          , label = "Fluid velocity field"
    #          # , vertical = false
    #          )

    # ax2 = Axis(fig[2, 1], aspect = DataAspect(), title="Pressure")
    # sc2 = scatter!( ax2 
    #               , mesh.p[1, vertices] # xPressure
    #               , mesh.p[2, vertices] # yPressure 
    #               , color = p[vertices]
    #               , colormap = :batlow
    #               , markersize = 10
    #               )
    
    # Colorbar(fig[2, 2], sc2
    #          , label = "Pressure"
    #          # , vertical = false
    #          )
    # # wait(display(fig))
    # display(GLMakie.Screen(), fig)
    # return


    # Prepare the heat transfer simulation
    println("Calculating the 2nd order mass matrix")
    M = quadraticMassMatrix2D(mesh)

    println("Calculating the 2nd order stiffness matrix")
    A = quadraticStiffnessMatrix2D(mesh, epsi)

    println("Calculating the 2nd order convection matrix")
    C = quadraticConvectionMatrix2D(mesh, u)
    
    return

    # Plot the heat transfer in real-time
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="0.0 s")
    graph = scatter!( ax
                    , mesh.p[1, :]
                    , mesh.p[2, :]
                    , color = T 
                    , colormap=:thermal 
                    , colorrange = (minimum(T), maximum(T))
                    # , markersize=5
                    )
    Colorbar(fig[1, 2], graph, label="T")
    display(fig)

    # Time iterations
    LM = M + timeStep*(A+C) # Backward Euler
    for frame in 1:maxSteps

        # Get the new temperature
        T = LM\(M*T)

        # Update plot
        # round(frame*timeStep*100.0)/100.0
        ax.title = string(frame*timeStep)*" s"
        graph.color = T

        sleep(1.0/24.0)
    end

    wait(fig.scene)

end
main(meshSize, localSize, showGmsh)