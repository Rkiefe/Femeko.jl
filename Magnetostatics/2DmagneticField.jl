#=
    An example on how to calculate the magnetostatic energy
    of non linear materials in 2D
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")
include("../src/magneticProperties.jl")

using GLMakie

# Wrapper to Fortran dierckx | Interpolation functions
using Dierckx

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
    picardDeviation::Float64 = 1e-4
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

    cells = []

    # Add a 2D rectangle
    # id = addRectangle([0,0,0], [0.25, 1], cells)
    id = addDisk([0,0,0], 1.0, cells)

    # Add a container
    # box = addRectangle([0,0,0], [2, 4])
    box = addDisk([0,0,0], 5.0)

    # Combine the geometries
    gmsh.model.occ.fragment(vcat(cells,[(2,box)]), [])
    gmsh.model.occ.synchronize()

    # Generate mesh
    mesh::MESH = Mesh2D(cells, meshSize, localSize)

    println("\nNumber of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()


    # Element centroids
    centroids::Matrix{Float64} = zeros(2,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[1,k] = sum(mesh.p[1,nds])/3
        centroids[2,k] = sum(mesh.p[2,nds])/3
    end

    # Magnetostatic simulation
    
    # Boundary conditions
    bc::Vector{Float64} = zeros(mesh.ne)
    for e in 1:mesh.ne
        if mesh.surfaceT[3,e] == 1 # Only the container boundary
            bc[e] = mu0*dot(Hext, mesh.normal[:,e])
        end
    end

    # Boundary integral
    RHS::Vector{Float64} = zeros(mesh.nv)
    for e in 1:mesh.ne
        nds = @view mesh.surfaceT[:, e]

        # Length of the edge
        l::Float64 = norm(mesh.p[:,nds[2]] - mesh.p[:,nds[1]])

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
    u::Vector{Float64} = zeros(mesh.nv+1)
    
    # Magnetic field
    H_field::Matrix{Float64} = zeros(mesh.nt, 2)
    H::Vector{Float64} = zeros(mesh.nt)
    Hold::Vector{Float64} = zeros(mesh.nt)
    
    # Picard iteration to handle non-linear material
    att::Int32 = 0
    div::Float64 = maxDeviation + 1.0
    while div > picardDeviation && att < maxAtt 

        att += 1
        Hold .= H

        # Stiffness matrix
        A = stiffnessMatrix2D(mesh, mu)

        # Magnetic scalar potential
        u = [A Lag;Lag' 0]\[-RHS;0]

        # Magnetic field
        H_field .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[:,k];

            # Sum the contributions
            for nd in nds
                # obtain the element parameters
                _, bi, ci = abc(mesh.p, nds, nd)

                H_field[k,1] -= u[nd]*bi;
                H_field[k,2] -= u[nd]*ci;
            end
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(H_field[k,:])
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
    B_field::Matrix{Float64} = mu.*H_field
    B::Vector{Float64} = mu.*H

    # Calculate the magnetostatic energy
    energy::Float64 = getEnergy(mesh, data, H, B)

    # Adjust the volume integral by the scale
    energy *= scale*depth

    println("Energy: ",energy, " J")


    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect(), title="Magnetic field H")

    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        color = mu0.*H[mesh.InsideElements], 
        colormap=:rainbow, 
        markersize=20) # 20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements])

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl

end

main(1.0, 0.05, true)