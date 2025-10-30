#=
    This showcases how to make many bodies, one inside the other
    while making a unified mesh for the entire volume.

    This example is for an iron sphere inside a gadolinium cube
    and with a bouding shell (air)
=#


include("../src/Femeko.jl")
include("../src/magneticProperties.jl")

# For plots
using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=true, saveMesh=false)
    #=
        Makes a model with cubes and spheres and refines the mesh on the spheres
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
            saveMesh  - Save mesh to a FEMCE compatible format 

    =#

    # vacuum magnetic permeability
    mu0 = pi*4e-7

    # Temperature
    T::Float64 = 293.0

    # Applied field | A/m
    Hext::Vector{Float64} = 1.0/mu0*[1.0,
                                     0.0,
                                     0.0]

    # Convergence criteria
    picardDeviation::Float64 = 1e-2
    maxDeviation::Float64 = 1e-6
    maxAtt::Int32 = 100
    relax::Float64 = 1.0 # Relaxation factor for N-R ]0, 1.0]

    # Data of magnetic materials
    materialProperties = Dict("Gd" => DATA(),
                              "Fe" => DATA())
    # Load Gd
    loadMaterial( materialProperties,
                       "Materials", # Folder with materials
                       "Gd_MFT",  # Data folder of target material
                       "Gd",        # Material name
                       7.9,
                       T)
    
    # Load Iron data
    materialProperties["Fe"].HofM = vec(readdlm("Materials/Pure_Iron_FEMM/H_Fe_extrap.dat"))  # A/m
    materialProperties["Fe"].B = vec(readdlm("Materials/Pure_Iron_FEMM/B_Fe_extrap.dat"))     # T

    # Get the permeability and its derivative
    materialPermeability(materialProperties["Fe"])

    # 3D Model
    gmsh.initialize()
    
    # Array of volume cell IDs
    cells = []
    cellLabels = []

    # Add a Fe sphere
    addSphere([0,0,0], 0.125, cells)
    push!(cellLabels, "Fe")

    # Add Gd cube
    addCuboid([0,0,0], [1.0, 1.0, 1.0], cells)
    push!(cellLabels, "Gd")

    # Unify model since there are cell intersections
    unifyModel(cells)

    # Create bounding shell
    box = addSphere([0,0,0], 5.0)
    push!(cellLabels, "Air")

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id = unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize, saveMesh)

    # Get element tags to then use GMSH 'get' functions
    t_tags, _ = gmsh.model.mesh.getElementsByType(4)

    # Store cell id of each element
    elementID::Vector{Int32} = zeros(mesh.nt)
    for k in 1:mesh.nt
        # element type , nodes of the element , dimension , id
        _, _, _, id = gmsh.model.mesh.getElement(t_tags[k])
        elementID[k] = id
    end
    
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
        nds = @view mesh.t[:, k]
        centroids[:, k] = mean(mesh.p[:, nds], 2)
    end

    # Check cell IDs
        # fig = Figure()
        # ax = Axis3(fig[1, 1], aspect = :data, title="Mesh and Cell Id")
        # scatterPlot = scatter!(ax, 
        #     centroids[1, :], # mesh.InsideElements
        #     centroids[2, :], # mesh.InsideElements
        #     centroids[3, :], # mesh.InsideElements
        #     color=elementID[:], # mesh.InsideElements
        #     colormap=:rainbow,  # :CMRmap :redsblues
        #     markersize=20)
        #     # markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements]))
        #     # color = mu0.*H[mesh.InsideElements], 
        # Colorbar(fig[1, 2], scatterPlot, label="Element ID") # Add a colorbar
        
        # # Display the figure (this will open an interactive window)
        # wait(display(fig)); return

    # Boundary conditions
    RHS::Vector{Float64} = BoundaryIntegral(mesh,mu0.*Hext,shell_id)

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
    div::Float64 = 1.0
    while div > picardDeviation && att < maxAtt 

        att += 1
        Hold .= H

        # Stiffness matrix
        A = stiffnessMatrix(mesh, mu)

        # Check condition number
        # println("Condition: ", cond(Matrix([A Lag; Lag' 0])))

        # Magnetic scalar potential
        u = [A Lag;Lag' 0]\[-RHS;0]

        # Magnetic field
        Hfield .= 0.0
        for k in 1:mesh.nt
            nds = mesh.t[:,k];

            # Sum the contributions
            for nd in nds
                # obtain the element parameters
                _,b,c,d = abcd(mesh.p,nds,nd)

                Hfield[1, k] -= u[nd]*b;
                Hfield[2, k] -= u[nd]*c;
                Hfield[3, k] -= u[nd]*d;
            end
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(Hfield[:, k])
        end

        # Update magnetic permeability
        for i in 1:length(cells)

            # Cell ID
            id = cells[i][2]

            # Get the data set of current cell ID
            key = cellLabels[id]

            if key == "Air"
                continue
            end

            spl = Spline1D(materialProperties[key].HofM,
                           materialProperties[key].mu
                           ) # ;bc="nearest") # nearest , extrapolate

            # Find all elements of current cell ID
            elements = findall(x -> x==id, elementID)
            
            # Interpolate the dataset for this elements
            mu[elements] .= spl(H[elements])

            idx = findall(findErr -> !isfinite(findErr), mu)
            if !isempty(idx)
                println(idx)
                error("Nans/Infs in mu")
            end
        end

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

        # Stiffness matrix
        A = stiffnessMatrix(mesh, mu)

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

        # Update magnetic permeability
        for i in 1:length(cells)

            # Cell ID
            id = cells[i][2]

            # Get the data set of current cell ID
            key = cellLabels[id]

            if key == "Air"
                continue
            end

            spl = Spline1D(materialProperties[key].HofM,
                           materialProperties[key].mu
                           ) # ;bc="nearest") # nearest , extrapolate

            # Find all elements of current cell ID
            elements = findall(x -> x==id, elementID)
            
            # Interpolate the dataset for this elements
            mu[elements] .= spl(H[elements])


            # d/dH mu

            # Get the data set of current cell ID
            key = cellLabels[id]

            if key == "Air"
                continue
            end

            spl = Spline1D(materialProperties[key].HofM,
                           materialProperties[key].dmu
                           ) # ;bc="extrapolate") # nearest , extrapolate

            dmu[elements] .= spl(H[elements])

            # Check for nans
            idx = findall(findErr -> !isfinite(findErr), mu)
            if !isempty(idx)
                println(idx)
                error("Nans/Infs in mu")
            end

            idx = findall(findErr -> !isfinite(findErr), dmu)
            if !isempty(idx)
                println(idx)
                error("Nans/Infs in dmu")
            end

        end # Data interpolation

        # Check deviation from previous result
        div = mu0*maximum(abs.(H[mesh.InsideElements].-Hold[mesh.InsideElements]))
        println(att, " | mu0 |H(n)-H(n-1)| = ", div)

    end # Newton iteration

    # Magnetic flux
    B::Vector{Float64} = mu.*H

    # Calculate the magnetostatic energy
    # >> not implemented for heterogeneous materials
    
    # Magnetization
    chi::Vector{Float64} = mu./mu0 .- 1.0

    M::Vector{Float64} = chi.*H

    Mfield::Matrix{Float64} = zeros(3, mesh.nt)
    Mfield[1, :] = chi.*Hfield[1, :]
    Mfield[2, :] = chi.*Hfield[2, :]
    Mfield[3, :] = chi.*Hfield[3, :]

    # Cell ID of Gd
    i = 2
    id = cells[i][2]

    # Find all elements of current cell ID
    elements = findall(x -> x==id, elementID)

    # Average magnetization of Gd
    M_avg::Float64 = 0.0
    volume::Float64 = 0.0
    for k in elements
        volume += mesh.VE[k]
        M_avg  += mesh.VE[k]*M[k]
    end
    M_avg /= volume
    println("\n<M> (emu/g) = ", M_avg/(7.9*1e3))


    # Plot result | Uncomment "using GLMakie"
    println("Generating plots...")

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)
    
    # graph = scatterPlot = scatter!(ax, 
    #     centroids[1,mesh.InsideElements],
    #     centroids[2,mesh.InsideElements],
    #     centroids[3,mesh.InsideElements], 
    #     color = M[mesh.InsideElements]./(7.9*1e3)
    #     , colormap=:turbo # :CMRmap :rainbow :CMRmap
    #     , markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements])
    #     # , markersize=20
    #     )

    ux::Vector{Float64} = Mfield[1, mesh.InsideElements]
    uy::Vector{Float64} = Mfield[2, mesh.InsideElements]
    uz::Vector{Float64} = Mfield[3, mesh.InsideElements]

    graph= arrows3d!(ax, centroids[1, mesh.InsideElements]
                , centroids[2, mesh.InsideElements]
                , centroids[3, mesh.InsideElements]
                , ux./maximum(M[mesh.InsideElements])
                , uy./maximum(M[mesh.InsideElements])
                , uz./maximum(M[mesh.InsideElements])
                , color = M[mesh.InsideElements]./(7.9*1e3)
                , lengthscale = 0.1
                , colormap = :turbo
                )

    Colorbar(fig[1, 2], graph, label="M (emu/g)")
    
    wait(display(fig))
    # save("fig.png",fig)

end

meshSize  = 2.5
localSize = 0.1
showGmsh = true

main(meshSize, localSize, showGmsh)

