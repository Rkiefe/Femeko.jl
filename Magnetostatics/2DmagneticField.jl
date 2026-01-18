#=
    An example on how to calculate the magnetostatic energy
    of non linear materials in 2D
=#

include("../src/Femeko.jl")
include("../src/magneticProperties.jl")

using IterativeSolvers

using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # vacuum magnetic permeability
    mu0::Float64 = pi*4e-7

    # Temperature
    T::Float64 = 293.0

    # Applied field
    Hext::Vector{Float64} = 1.35.*[1,0]./mu0    # A/m

    # Scale of the model
    scale = 1e-4 # cm2
    depth = 1e-3 # 1 mm , same as FEMM default

    # Convergence criteria
    maxDeviation::Float64 = 1e-6
    maxAtt::Int32 = 10
    relax::Float64 = 1.0 # Relaxation factor for N-R ]0, 1.0]

    # Data of magnetic materials
    data = DATA()
    loadMaterial( data,
                 "Materials", # Folder with materials
                 "Gd_MFT",    # Data folder of target material
                 "Gd",        # Material name
                 7.9,
                 T)

    # Create a spline to interpolate the material properties
    spl = Spline1D(data.HofM,
                   data.mu
                   ) # ;bc="nearest") # nearest , extrapolate

    # Create model
    gmsh.initialize()

    cells = [] # Store the cell IDs (dim, tag)

    # Add a 2D rectangle
    # id = addRectangle([0,0,0], [0.25, 1], cells)
    id = addDisk([0,0,0], 1.0, cells)

    # Add a container
    box = addDisk([0,0,0], 5.0)

    # Combine the geometries to make a conforming mesh
    shell_id, box = unifyModel(cells, box)

    # Generate mesh
    extendLocalRefinement() # Don't extend the refinement to the volume
    mesh::MESH = Mesh2D(cells, meshSize, localSize)

    println("\nNumber of elements ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of nodes ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ns)

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()


    # Element centroids
    centroids::Matrix{Float64} = zeros(2, mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[:, k] = mean(mesh.p[1:2, nds], 2)
    end

    # Magnetostatic simulation
    
    # Boundary conditions
    bc::Vector{Float64} = zeros(mesh.ns)
    for e in 1:mesh.ns
        if !isempty(intersect(mesh.surfaceT[3, e], shell_id)) # Only the container boundary
            bc[e] = mu0*dot(Hext, mesh.normal[:, e])
        end
    end

    # Boundary integral
    RHS::Vector{Float64} = zeros(mesh.nv)
    for e in 1:mesh.ns
        nds = @view mesh.surfaceT[:, e]

        # Length of the edge
        l::Float64 = mesh.AE[e]

        RHS[nds[1]] = RHS[nds[1]] + 0.5*l*bc[e]
        RHS[nds[2]] = RHS[nds[2]] + 0.5*l*bc[e]
    end

    # Lagrange multiplier technique
    Lag::Vector{Float64} = zeros(mesh.nv)
    for k = 1:mesh.nt
        nds = @view mesh.t[:,k]
        Lag[nds] .+= mesh.VE[k]/3
    end

    # Magnetic permeability
    mu::Vector{Float64} = zeros(mesh.nt) .+ mu0

    # Magnetic scalar potential
    u::Vector{Float64} = zeros(mesh.nv)
    
    # Magnetic field
    Hfield::Matrix{Float64} = zeros(2, mesh.nt)
    H::Vector{Float64} = zeros(mesh.nt)
    Hold::Vector{Float64} = zeros(mesh.nt)
    
    # Picard iteration to handle non-linear material
    att::Int32 = 0
    div::Float64 = Inf
    while div > maxDeviation && att < maxAtt 

        att += 1
        Hold .= H

        # Stiffness matrix
        A = stiffnessMatrix2D(mesh, mu)

        # Magnetic scalar potential
        # u = [A Lag;Lag' 0]\[RHS;0]
        u = cg(A, RHS)

        # Magnetic field
        Hfield .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[:,k];

            # Sum the contributions
            for nd in nds
                # obtain the element parameters
                _, bi, ci = abc(mesh.p, nds, nd)

                Hfield[1, k] -= u[nd]*bi;
                Hfield[2, k] -= u[nd]*ci;
            end
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(Hfield[:, k])
        end

        # Update magnetic permeability
        mu[mesh.InsideElements] .= spl(H[mesh.InsideElements])

        # Check interpolation
        idx = findall(findErr -> !isfinite(findErr), mu)
        if !isempty(idx)
            println(idx)
            error("Nans or Infs in mu")
        end

        # Check deviation from previous result
        div = mu0*maximum(abs.(H[mesh.InsideElements].-Hold[mesh.InsideElements]))
        println(att, " | mu0 |H(n)-H(n-1)| = ", div)

    end # Picard Iteration


    # Magnetic flux density
    Bfield::Matrix{Float64} = mu' .*Hfield
    B::Vector{Float64} = mu.*H

    # Calculate the magnetostatic energy
    energy::Float64 = getEnergy(mesh, data, H, B)

    # Adjust the volume integral by the scale
    energy *= scale*depth

    println("Energy: ",energy, " J")

    # Plot result | Uncomment "using GLMakie"
    println("Generating plots...")
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="Magnetic field B")

    graph = arrows2d!(  ax,
                        centroids[1, :], centroids[2, :],
                        Bfield[1, :], Bfield[2, :]
                        , color = B
                        , lengthscale = 0.1
                        , colormap = :rainbow,  # :CMRmap :viridis :redsblues :turbo :rainbow
                    )

    Colorbar(fig[1, 2], graph, label="B field")
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl

end

meshSize = 1.0
localSize = 0.05
showGmsh = false

main(meshSize, localSize, showGmsh)