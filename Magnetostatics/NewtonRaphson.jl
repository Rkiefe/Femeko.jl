#=
    Newton Raphson (or Newton Galerkin) iteration method
    for magnetostatics    
=#


include("../src/Femeko.jl")
include("../src/magneticProperties.jl")

using GLMakie

function main(meshSize=0, localSize=0, showGmsh=true, saveMesh=false)
    #=
        Makes a model with cubes and spheres and refines the mesh on the spheres
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
            saveMesh  - Save mesh to a FEMCE compatible format 

    =#

    # vacuum magnetic permeability
    mu0 = pi*4e-7

    # Scale of model
    scale::Float64 = 1e-6 # cm^3 -> m^3

    # Temperature
    T::Float64 = 300.0

    # Applied field A/m
    Hext::Vector{Float64} = [1.0, 
                             0.0, 
                             0.0]*0.5/mu0

    # Convergence criteria
    picardDeviation::Float64 = 1e-4
    maxDeviation::Float64 = 1e-8 # Inf -> Don't run the N-R method
    maxAtt::Int32 = 100
    relax::Float64 = 1.0 # Relaxation factor for N-R ]0, 1.0]

    # Data of magnetic materials
    data = DATA()

    # Load Gd
    loadMaterial( data,
                  "Materials",  # Folder with materials
                  "Gd_MFT",     # Data folder of target material
                  "Gd",         # Material name
                  7.9,          # Density g/cm3
                  T)

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
    mesh = Mesh(cells,meshSize,localSize,saveMesh)
    
    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

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
    RHS::Vector{Float64} = BoundaryIntegral(mesh, mu0.*Hext, shell_id)

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

    # FEM
    u::Vector{Float64} = zeros(mesh.nv+1)
    
    Hfield::Matrix{Float64} = zeros(3, mesh.nt)
    H::Vector{Float64} = zeros(mesh.nt)
    Hold::Vector{Float64} = zeros(mesh.nt)
    
    att::Int32 = 0
    div::Float64 = Inf
    while div > picardDeviation && att < maxAtt 

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

        # Magnetic scalar potential
        u = [A Lag;Lag' 0]\[-RHS;0]
        
        # Check solution
        if any(x -> !isfinite(x), u)
            error("Nans/Infs in the scalar potential")
        end

        # Magnetic field
        Hfield .= 0.0
        for k in 1:mesh.nt
            nds = mesh.t[:,k];

            # Sum the contributions
            for nd in nds
                # obtain the element parameters
                _,b,c,d = abcd(mesh.p,nds,nd)

                Hfield[1,k] -= u[nd]*b;
                Hfield[2,k] -= u[nd]*c;
                Hfield[3,k] -= u[nd]*d;
            end
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(Hfield[:, k])
        end

        # Update magnetic permeability            
        mu[mesh.InsideElements] .= spl(H[mesh.InsideElements])

        # Check deviation from previous result
        div = mu0*maximum(abs.(H[mesh.InsideElements].-Hold[mesh.InsideElements]))
        println(att, " | mu0 |H(n)-H(n-1)| = ", div)

    end # Picard Iteration

    # Remove the lagrange multiplier
    pop!(u) # Remove last element of u

    println("Newton Raphson")

    # Newton-Raphson
    dmu::Vector{Float64} = zeros(mesh.nt)
    du::Vector{Float64} = zeros(mesh.nv+1)
    while div > maxDeviation && att < maxAtt # maxAtt

        att += 1
        Hold .= H

        # Update magnetic permeability
        mu[mesh.InsideElements] .= spl(H[mesh.InsideElements])
        
        # d/dH mu
        dmu[mesh.InsideElements] .= spl_dmu(H[mesh.InsideElements])

        # Check for nans
        
        if any(x -> !isfinite(x), mu) || any(x -> !isfinite(x), dmu)
            error("Nans/Infs in mu")
        end

        # Updated the compressed stiffness matrix
        n = 0
        for i in 1:4
            for j in 1:4
                n += 1
                for k in 1:mesh.nt
                    Acsc[(n-1)*mesh.nt + k] = Ak[n, k] * mu[k]
                end
            end
        end # Acsc

        # Update the global stiffness matrix
        A = sparse(rowIDs, colIDs, Acsc, mesh.nv, mesh.nv)

        # Tangential stiffness matrix
        At = tangentialStiffnessMatrix(mesh, Hfield, dmu)

        # Correction to the magnetic scalar potential
        du = [A+At Lag;Lag' 0]\[-RHS-A*u;0]
        
        # Norm of correction over the value
        println(att, " |du|/|u| = ", norm(du[1:mesh.nv])/norm(u))

        # Update the potential
        u .+= relax.*du[1:mesh.nv]

        # Magnetic field
        Hfield .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[:,k];

            Hx::Float64 = 0
            Hy::Float64 = 0
            Hz::Float64 = 0
            # Sum the contributions
            for nd in nds
                # obtain the element parameters
                _,b,c,d = abcd(mesh.p,nds,nd)

                Hx -= u[nd]*b
                Hy -= u[nd]*c
                Hz -= u[nd]*d
            end
            
            Hfield[1, k] = Hx
            Hfield[2, k] = Hy
            Hfield[3, k] = Hz
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(Hfield[:, k])
        end

        # Check deviation from previous result
        div = mu0*maximum(abs.(H[mesh.InsideElements].-Hold[mesh.InsideElements]))
        println(att, " | mu0 |H(n)-H(n-1)| = ", div)

    end # Newton iteration

    # Magnetic flux
    B::Vector{Float64} = mu.*H

    # Calculate the magnetostatic energy
    energy::Float64 = getEnergy(mesh, data, H, B)

    # Adjust the volume integral by the scale
    energy *= scale
    println("Energy (J): ", energy)

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
                        , lengthscale = 0.5
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    # Add a colorbar
    Colorbar(fig[1, 2], graph, label = "M (emu/g)")

    wait(display(fig))

end # end of main

meshSize  = 5.0
localSize = 0.1
showGmsh = false

main(meshSize, localSize, showGmsh)

