#=
    Newton Raphson (or Newton Galerkin) iteration method
    for magnetostatics    
=#


include("../src/Femeko.jl")
include("../src/magneticProperties.jl")

using IterativeSolvers
using GLMakie

function main(meshSize=0, localSize=0, showGmsh=true, verbose=false)
    #=
        Makes a model with cubes and spheres and refines the mesh on the spheres
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
    =#

    # vacuum magnetic permeability
    mu0 = pi*4e-7

    # Scale of model
    scale::Float64 = 1e-6 # cm^3 -> m^3

    # Temperature
    T::Float64 = 270.0

    # Applied field T
    Bext::Vector{Float64} = [0.5, 
                             0.0, 
                             0.0]

    # Convergence criteria
    maxDeviation::Float64 = 1e-4 # Tolerance in |H| (Tesla)
    maxAtt::Int32 = 100 # Maximum number of iterations allowed

    # Load Gd
    data = DATA()
    loadMaterial( data,
                  "Materials",  # Folder with materials
                  "Gd_MFT",     # Data folder of target material
                  "Gd",         # Material name
                  7.9,          # Density g/cm3
                  T)

    # Load Fe
    # data = DATA()
    # data.density = 7.874 # g/cm3
    # data.HofM = vec(readdlm("Materials/Pure_Iron_FEMM/H_Fe_extrap.dat"))  # A/m
    # data.B = vec(readdlm("Materials/Pure_Iron_FEMM/B_Fe_extrap.dat"))     # T
    # materialPermeability(data) # Get the permeability and its derivative

    # Spline for interpolation of the permeability
    spl = Spline1D(data.HofM, data.mu
                   ; bc="error" # nearest , extrapolate , error
                   )

    spl_dmu = Spline1D( data.HofM, data.dmu
                      ; bc="error" # nearest , extrapolate , error
                      ) 

    # 3D Model
    gmsh.initialize()
    
    # Array of volume cell IDs
    cells = []

    # Import cad file
    # box = importCAD("../STEP_Models/cube.STEP", cells, true)

    # Add material
    addCuboid([0,0,0], [1.0, 1.0, 1.0], cells)
    # addSphere([0,0,0], 0.5, cells)

    # Create bounding shell
    box = addSphere([0,0,0], 5.0)

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id, box = unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize)
    
    println("Number of elements ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of nodes ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ns)

    if showGmsh
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()

    # Element centroids
    centroids::Matrix{Float64} = zeros(3, mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:, k]
        centroids[:, k] = mean(mesh.p[:, nds], 2)
    end

    # Boundary conditions
    RHS::Vector{Float64} = - BoundaryIntegral(mesh, Bext, shell_id)

    # Local stiffness matrix
    Ak::Matrix{Float64} = localStiffnessMatrix(mesh)

    # Convert to a compressed sparse column format
    rowIDs::Vector{Int} = zeros(16*mesh.nt)
    colIDs::Vector{Int} = zeros(16*mesh.nt)
    Acsc = zeros(length(rowIDs)) # Compressed sparse column format of the stiffness matrix
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            for k in 1:mesh.nt
                rowIDs[(n-1)*mesh.nt + k] = mesh.t[i, k]
                colIDs[(n-1)*mesh.nt + k] = mesh.t[j, k]
            end
        end
    end # Compressed sparse stiffness matrix indices

    # Lagrange multiplier technique
    Lag::Vector{Float64} = lagrange(mesh)

    # Magnetic permeability
    mu::Vector{Float64} = zeros(mesh.nt) .+ mu0

    # Prepare output
    u::Vector{Float64} = zeros(mesh.nv) # Scalar potential
    
    Hfield::Matrix{Float64} = zeros(3, mesh.nt) # Vector field
    H::Vector{Float64} = zeros(mesh.nt)
    Hold::Vector{Float64} = zeros(mesh.nt)
    
    att::Int32 = 0
    div::Float64 = Inf

    # Newton-Raphson
    dmu::Vector{Float64} = zeros(mesh.nt)
    du::Vector{Float64} = zeros(mesh.nv+1)
    while div > maxDeviation && att < maxAtt # maxAtt

        att += 1
        Hold .= H

        # Updated the compressed stiffness matrix
        for n in 1:16
            for k in 1:mesh.nt
                Acsc[(n-1)*mesh.nt + k] = Ak[n, k] * mu[k]
            end
        end # Acsc

        # Update the global stiffness matrix
        A = sparse(rowIDs, colIDs, Acsc, mesh.nv, mesh.nv)

        # Tangential stiffness matrix
        if att > 1
            At = tangentialStiffnessMatrix(mesh, Hfield, dmu)
        else
            At = spzeros(mesh.nv, mesh.nv)
        end

        # Correction to the magnetic scalar potential
        du = [A+At Lag;Lag' 0]\[RHS-A*u;0]
        # du[1:mesh.nv] = cg(A+At, RHS-A*u)

        # Update the potential
        residue = lineSearch(u, du[1:mesh.nv], A, RHS) # residue = norm(RHS - A*u)

        # Magnetic field
        Hfield .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[:,k];

            # Sum the contributions
            for nd in nds
                # obtain the element parameters
                _,b,c,d = abcd(mesh.p,nds,nd)

                Hfield[1, k] -= u[nd]*b;
                Hfield[2, k] -= u[nd]*c;
                Hfield[2, k] -= u[nd]*d;
            end
            
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(Hfield[:, k])
        end

        # Check deviation from previous result
        div = mu0*maximum(abs.(H[mesh.InsideElements].-Hold[mesh.InsideElements]))
        verbose ? println(att, " | mu0 |H(n)-H(n-1)| = ", div, " , |y-Ax| = ", residue) : nothing

        # Update magnetic permeability
        mu[mesh.InsideElements] .= spl(H[mesh.InsideElements])
        
        # d/dH mu
        dmu[mesh.InsideElements] .= spl_dmu(H[mesh.InsideElements])

        # Check for nans
        if any(x -> !isfinite(x), mu) || any(x -> !isfinite(x), dmu)
            error("Nans/Infs in mu")
        end

    end # Newton iteration

    # Magnetic flux
    B::Vector{Float64} = mu.*H

    # Calculate the magnetostatic energy
    # energy::Float64 = getEnergy(mesh, data, H, B)
    # energy *= scale # Adjust the volume integral by the scale
    # println("Energy (J): ", energy)

    # Magnetization
    chi::Vector{Float64} = mu./mu0 .- 1.0
    M::Vector{Float64} = chi.*H

    Mfield::Matrix{Float64} = zeros(3, mesh.nt)
    Mfield[1, :] = chi.*Hfield[1, :]
    Mfield[2, :] = chi.*Hfield[2, :]
    Mfield[3, :] = chi.*Hfield[3, :]
    
    # Average magnetization
    M_avg::Float64 = 0.0
    volume::Float64 = 0.0
    for k in mesh.InsideElements
        volume += mesh.VE[k]
        M_avg += M[k]*mesh.VE[k]
    end # <H>
    M_avg /= volume
    println("<M> (emu/g) = ", M_avg/(data.density*1e3))

    # Output plot
    println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)

    graph = arrows3d!(  ax
                        , centroids[1, mesh.InsideElements]
                        , centroids[2, mesh.InsideElements]
                        , centroids[3, mesh.InsideElements]
                        , mu0*Mfield[1, mesh.InsideElements]
                        , mu0*Mfield[2, mesh.InsideElements]
                        , mu0*Mfield[3, mesh.InsideElements]
                        , color = M[mesh.InsideElements]./(data.density*1e3)
                        , lengthscale = 0.1
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    # Add a colorbar
    Colorbar(fig[1, 2], graph, label = "M (emu/g)")
    wait(display(fig))


    # -- Slice view of the vector field over a 2D grid --
    println("Interpreting the result over the slice view plane...")
    
    # XoZ plane
    X, Y, Z = plane([1,0,0], [0,0,1], [0,0,0], 1.0, 20)

    Hx = zeros(size(X))
    Hy = zeros(size(X))
    Hz = zeros(size(X))
    color = zeros(size(X)) # Norm of the H field in the slice plane
    
    @time for i in 1:size(X, 1)
        for j in 1:size(X, 2)
            
            xq = X[i, j]
            yq = Y[i, j]
            zq = Z[i, j]
                
            # Interpolate the vector field
            Hx[i, j] = interp3Dmesh(mesh, xq, yq, zq, Hfield[1, :])
            Hy[i, j] = interp3Dmesh(mesh, xq, yq, zq, Hfield[2, :])
            Hz[i, j] = interp3Dmesh(mesh, xq, yq, zq, Hfield[3, :])

            color[i, j] = sqrt(Hx[i, j]^2 + Hy[i, j]^2 + Hz[i, j]^2)
        end
    end

    # Plot slice view
    println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)

    graph = arrows3d!(  ax
                        , X[:]
                        , Y[:]
                        , Z[:]
                        , mu0*Hx[:]
                        , mu0*Hy[:]
                        , mu0*Hz[:]
                        , color = color[:]
                        , lengthscale = 0.2
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    # Add a colorbar
    # Colorbar(fig[1, 2], graph)
    wait(display(fig))

end # end of main

meshSize  = 5.0
localSize = 0.1
showGmsh = false
verbose = true

main(meshSize, localSize, showGmsh, verbose)