#=
    This showcases how to make many bodies, one inside the other
    while making a unified mesh for the entire volume.

    This example is for an iron sphere inside a gadolinium cube
    and with a bouding shell (air)
=#


include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")
include("../src/magneticProperties.jl")

# For plots
using GLMakie

function main(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)
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
    T::Float64 = 293.0

    # Applied field
    Hext::Vector{Float64} = 1.2.*[1,0,0]./mu0    # A/m

    # Convergence criteria
    picardDeviation::Float64 = 1e-4
    maxDeviation::Float64 = 1e-6
    maxAtt::Int32 = 10
    relax::Float64 = 1.0 # Relaxation factor for N-R ]0, 1.0]

    # Data of magnetic materials
    materialProperties = Dict("Gd" => DATA(),
                              "Fe" => DATA())
    # Load Gd
    loadMaterial( materialProperties,
                       "Materials", # Folder with materials
                       "Gd_Ising",  # Data folder of target material
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

    # Import cad file
    # box = importCAD("../STEP_Models/cube.STEP", cells)
    # cellLabels = ["Gd"]
    # push!(cellLabels, "Air")

    # Add Gd cube
    addCuboid([0,0,0], [1.0, 1.0, 1.0], cells, true)
    cellLabels::Vector{String} = ["Gd"]

    # Add a Fe sphere
    addSphere([0,0,0], 0.125, cells, true)
    push!(cellLabels, "Fe")

    # Unify model since there are cell intersections
    unifyModel(cells)

    # Create bounding shell
    box = addSphere([0,0,0], 8.0)
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
        # gmsh.option.setNumber("Mesh.Clip", 1)
        # gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        # gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println(shell_id)
    println(cells)
    println(box)

    # Element centroids
    centroids::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[1,k] = sum(mesh.p[1,nds])/4
        centroids[2,k] = sum(mesh.p[2,nds])/4
        centroids[3,k] = sum(mesh.p[3,nds])/4
    end

    # Check cell IDs
        # fig = Figure()
        # ax = Axis3(fig[1, 1], aspect = :data, title="Mesh and Cell Id")
        # scatterPlot = scatter!(ax, 
        #     centroids[1,:], # mesh.InsideElements
        #     centroids[2,:], # mesh.InsideElements
        #     centroids[3,:], # mesh.InsideElements
        #     color=elementID[:], # mesh.InsideElements
        #     colormap=:rainbow,  # :CMRmap :viridis :redsblues
        #     markersize=20)
        #     # markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements]))
        #     # color = mu0.*H[mesh.InsideElements], 
        # Colorbar(fig[1, 2], scatterPlot, label="Element ID") # Add a colorbar
        
        # # Display the figure (this will open an interactive window)
        # wait(display(fig))

    # Boundary conditions
    RHS::Vector{Float64} = BoundaryIntegral(mesh,mu0.*Hext,shell_id)

    # Lagrange multiplier technique
    Lag::Vector{Float64} = lagrange(mesh)

    # Magnetic permeability
    mu::Vector{Float64} = zeros(mesh.nt) .+ mu0

    # FEM
    u::Vector{Float64} = zeros(mesh.nv+1)
    
    H_vectorField::Matrix{Float64} = zeros(mesh.nt,3)
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
        H_vectorField .= 0.0
        for k in 1:mesh.nt
            nds = mesh.t[:,k];

            # Sum the contributions
            for nd in nds
                # obtain the element parameters
                _,b,c,d = abcd(mesh.p,nds,nd)

                H_vectorField[k,1] -= u[nd]*b;
                H_vectorField[k,2] -= u[nd]*c;
                H_vectorField[k,3] -= u[nd]*d;
            end
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(H_vectorField[k,:])
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
        At = tangentialStiffnessMatrix(mesh, H_vectorField, dmu)

        # Correction to the magnetic scalar potential
        du = [A+At Lag;Lag' 0]\[-RHS-A*u;0]
        
        # Norm of correction over the value
        println(att, " |du|/|u| = ", norm(du[1:mesh.nv])/norm(u))

        # Update the potential
        u .+= relax.*du[1:mesh.nv]

        # Magnetic field
        H_vectorField .= 0.0
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
            
            H_vectorField[k,1] = Hx
            H_vectorField[k,2] = Hy
            H_vectorField[k,3] = Hz
        end

        # Magnetic field intensity
        for k in 1:mesh.nt
            H[k] = norm(H_vectorField[k,:])
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
    # energy::Float64 = getEnergy(mesh, materialProperties["Gd"], H, B)

    # Adjust the volume integral by the scale
    # energy *= scale
    # println("Energy (J): ", energy)
    
    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="With relax = "*string(relax))
    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        centroids[3,mesh.InsideElements], 
        color = mu0.*H[mesh.InsideElements], 
        colormap=:CMRmap, # :rainbow :CMRmap
        markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements]))
        # markersize=20)
    Colorbar(fig[1, 2], scatterPlot, label="|H|") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig))
    # save("relax_"*string(relax)*".png",fig)

end # end of main

meshSize  = 4.0
localSize = 0.1
showGmsh = true
saveMesh = false

main(meshSize,localSize,showGmsh,saveMesh)

