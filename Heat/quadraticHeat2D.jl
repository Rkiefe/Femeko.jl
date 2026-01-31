#  Heat equation with quadratic lagrange elements

include("../Fluids/fluids2D.jl")
using GLMakie

# Mesh settings
meshSize = 5.0
localSize = 0.2
showGmsh = true

function quadraticMassMatrix2D(mesh::MESH)
    # Proper 7-point Gaussian quadrature for triangles (degree of precision 5)
    # Coordinates and weights for reference triangle (0,0), (1,0), (0,1)
    quad_points = [
        (1.0/3.0, 1.0/3.0),
        (0.059715871789770, 0.470142064105115),
        (0.470142064105115, 0.059715871789770),
        (0.470142064105115, 0.470142064105115),
        (0.797426985353087, 0.101286507323456),
        (0.101286507323456, 0.797426985353087),
        (0.101286507323456, 0.101286507323456)
    ]
    
    quad_weights = [
        0.225000000000000,
        0.132394152788506,
        0.132394152788506,
        0.132394152788506,
        0.125939180544827,
        0.125939180544827,
        0.125939180544827
    ]
    
    # Reference triangle area is 0.5, so weights already sum to 0.5
    
    Mlocal = zeros(36, mesh.nt)
    
    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]
        vertices = mesh.p[:, nds[1:3]]  # Triangle vertices
        
        # Jacobian for affine transformation from reference to physical triangle
        # x = v1 + J11*ξ + J12*η
        # y = v2 + J21*ξ + J22*η
        J11 = vertices[1, 2] - vertices[1, 1]
        J12 = vertices[1, 3] - vertices[1, 1]
        J21 = vertices[2, 2] - vertices[2, 1]
        J22 = vertices[2, 3] - vertices[2, 1]
        
        detJ = J11 * J22 - J12 * J21  # Determinant of Jacobian
        abs_detJ = abs(detJ)
        
        # Precompute basis coefficients for all 6 nodes
        basis_coeffs = Vector{Vector{Float64}}(undef, 6)
        for i in 1:6
            basis_coeffs[i] = quadraticBasis2D(mesh.p, nds, nds[i])
        end
        
        # Initialize local mass matrix
        Mk = zeros(6, 6)
        
        # Loop over quadrature points
        for q in 1:length(quad_weights)
            xi, eta = quad_points[q]
            weight = quad_weights[q]
            
            # Transform from reference to physical coordinates
            x = vertices[1, 1] + J11 * xi + J12 * eta
            y = vertices[2, 1] + J21 * xi + J22 * eta
            
            # Evaluate all basis functions at this quadrature point
            phi_vals = zeros(6)
            for i in 1:6
                S = basis_coeffs[i]
                phi_vals[i] = S[1] + S[2]*x + S[3]*y + S[4]*x^2 + S[5]*x*y + S[6]*y^2
            end
            
            # Accumulate to mass matrix
            w = weight * abs_detJ
            for i in 1:6
                for j in i:6  # Only compute upper triangle
                    Mk[i, j] += w * phi_vals[i] * phi_vals[j]
                end
            end
        end
        
        # Fill lower triangle (symmetric)
        for i in 1:6
            for j in 1:(i-1)
                Mk[i, j] = Mk[j, i]
            end
        end
        
        Mlocal[:, k] = vec(Mk)
    end
    
    # Assemble global sparse mass matrix
    M = spzeros(mesh.nv, mesh.nv)
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            M += sparse(mesh.t[i, :], mesh.t[j, :], Mlocal[n, :], mesh.nv, mesh.nv)
        end
    end
    
    return M
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
    timeStep::Float64 = 1e-5
    totalTime::Float64 = 1e-2
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
    T::Vector{Float64} = zeros(mesh.nv) .+ 10.0
    T[mesh.InsideNodes] .= 1.0

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()
 
    # Define the diffusivity on the domain
    epsi::Vector{Float64} = zeros(mesh.nt)

    for i in 1:length(cells)

        id = cells[i][2] # Cell ID
        key = cellLabels[i] # Get the data set of current cell ID

        println("id: ", id, " ; key: ", key)

        # Find all elements of current cell ID
        elements = findall(x -> x==id, mesh.elementID)

        # Update viscosity value on this cell elements
        epsi[elements] .= materialProperties[key].k/(materialProperties[key].Cp * materialProperties[key].density)

    end

    # Prepare the heat transfer simulation
    println("Calculating the 2nd order mass matrix")
    M = quadraticMassMatrix2D(mesh)

    println("Calculating the 2nd order stiffness matrix")
    A = quadraticStiffnessMatrix2D(mesh, epsi)
    
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
    LM = M + timeStep*A # Backward Euler
    for frame in 1:maxSteps

        # Get the new temperature
        T = LM\(M*T)

        # Update plot
        # round(frame*timeStep*100.0)/100.0
        ax.title = string(frame*timeStep)*" s"
        graph.color = T

        totalHeat = sum(M * T)  # or ones(size(T))' * (M * T)
        println("Total heat at step $frame: $totalHeat")

        sleep(0.1)
    end

    wait(fig.scene)

end
# main(meshSize, localSize, showGmsh)