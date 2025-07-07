
include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

# For plots | Uncomment the plot section of "main()"
using GLMakie

# To read data from files
using DelimitedFiles

# Wrapper to Fortran dierckx | Interpolation functions
using Dierckx

mutable struct DATA
    # Magnetization data
    M # ::Matrix{Float64}

    # Magnetic field H
    HofM::Vector{Float64}
    
    # Temperature
    TofM::Vector{Float64}
    
    # Magnetic Flux
    B::Vector{Float64}

    # Density
    rho::Float64

    # Permeability mu = B/H
    mu::Vector{Float64}

    # d/dH mu (derivative of permeability)
    dmu::Vector{Float64}

    # Constructor
    DATA() = new()
end


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

    # Temperature
    T::Float64 = 293.0

    # Applied field
    Hext::Vector{Float64} = [1,0,0]./mu0    # A/m

    # Convergence criteria
    picardDeviation::Float64 = 1e-3
    maxDeviation::Float64 = 1e-6
    maxAtt::Int32 = 10

    # Data of magnetic materials
    folder::String = "Materials/"

    materialProperties = Dict("Gd" => DATA(),
                              "Fe" => DATA())


    # Load Gadolinium data
    materialProperties["Gd"].M = readdlm(folder*"Gd_MFT/M.dat")                 # emu/g
    materialProperties["Gd"].HofM = vec(readdlm(folder*"Gd_MFT/HofM.dat"))      # Oe
    materialProperties["Gd"].TofM = vec(readdlm(folder*"Gd_MFT/TofM.dat"))      # K
    materialProperties["Gd"].rho = 7.9 # g/cm3

    # Convert data units
    materialProperties["Gd"].M .*= materialProperties["Gd"].rho*1e3 # A/m
    materialProperties["Gd"].HofM .*= 1e-4/mu0                      # A/m

    # Interpolate data over the target temperature
    spl = Spline2D( materialProperties["Gd"].HofM,
                    materialProperties["Gd"].TofM,
                    materialProperties["Gd"].M)

    M::Vector{Float64} = zeros(length(materialProperties["Gd"].HofM))
    for i in 1:length(M)
        M[i] = spl(materialProperties["Gd"].HofM[i],T)
    end
    materialProperties["Gd"].M = M
    
    materialProperties["Gd"].B = mu0.*(materialProperties["Gd"].HofM .+
                                       materialProperties["Gd"].M)

    # Permeability
    materialProperties["Gd"].mu = materialProperties["Gd"].B./materialProperties["Gd"].HofM
    
    # Remove Inf
    idx = findall(x -> x==Inf, materialProperties["Gd"].mu)
    materialProperties["Gd"].mu[idx] .= 0.0
    materialProperties["Gd"].mu[idx] .= maximum(materialProperties["Gd"].mu)

    # d/dH mu
    dmu = gradient(materialProperties["Gd"].HofM, 
                   materialProperties["Gd"].B./materialProperties["Gd"].HofM)

    # Remove -Inf
    idx = findall(x -> x==-Inf, dmu)
    dmu[idx] .= 0.0
    dmu[idx] .= minimum(dmu)

    materialProperties["Gd"].dmu = dmu

    
    # Load Iron data
    materialProperties["Fe"].HofM = vec(readdlm(folder*"Pure_Iron_FEMM/H.dat"))  # A/m
    materialProperties["Fe"].B = vec(readdlm(folder*"Pure_Iron_FEMM/B.dat"))     # T

    # Permeability
    materialProperties["Fe"].mu = materialProperties["Fe"].B./materialProperties["Fe"].HofM
    
    # Remove Inf
    idx = findall(x -> x==Inf, materialProperties["Fe"].mu)
    materialProperties["Fe"].mu[idx] .= 0.0
    materialProperties["Fe"].mu[idx] .= maximum(materialProperties["Fe"].mu)

    # d/dH mu
    dmu = gradient(materialProperties["Fe"].HofM, 
                   materialProperties["Fe"].B./materialProperties["Fe"].HofM)

    # Remove -Inf
    idx = findall(x -> x==-Inf, dmu)
    dmu[idx] .= 0.0
    dmu[idx] .= minimum(dmu)

    materialProperties["Fe"].dmu = dmu

    # 3D Model
    gmsh.initialize()

    # Create an empty container
    # box = addCuboid([0,0,0],10.0*[1,1,1])
    box = addCuboid([0,0,0],8.25*[1,1,1])
    # box = addSphere([0,0,0],5*maximum(L))

    # Get how many surfaces compose the bounding shell
    temp = gmsh.model.getEntities(2)            # Get all surfaces of current model
    bounding_shell_n_surfaces = 1:length(temp)    # Get the number of surfaces in the bounding shell

    # List of cells inside the container
    cells = []
    cellLabels::Vector{String} = ["Air"]

    # Add Gadolinium
    addCuboid([0,0,0],[1.65,1.65,0.04],cells,true)
    push!(cellLabels,"Gd")

    # Add Iron
    # addCuboid([0,0,0.75],[1,1,0.5],cells,true)
    # push!(cellLabels,"Fe")


    # Fragment to make a unified geometry
    _, fragments = gmsh.model.occ.fragment([(3, box)], cells)
    gmsh.model.occ.synchronize()

    # Update container volume ID
    box = fragments[1][1][2]

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)
    
    # Get bounding shell surface id
    shell_id = gmsh.model.getAdjacencies(3, box)[2]

    # Must remove the surface Id of the interior surfaces
    shell_id = shell_id[bounding_shell_n_surfaces] # All other, are interior surfaces


    # Get element tags to then use GMSH 'get' functions
    t_tags, _ = gmsh.model.mesh.getElementsByType(4)

    # Store cell id of each element
    elementID::Vector{Int32} = zeros(mesh.nt)
    for k in 1:mesh.nt
        _, _, _, id = gmsh.model.mesh.getElement(t_tags[k])
        # element type , nodes of the element , dimension , id

        elementID[k] = id
    end
    

    if showGmsh
        # gmsh.option.setNumber("Mesh.Clip", 1)
        # gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        # gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

    # Element centroids
    centroids::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[1,k] = sum(mesh.p[1,nds])/4
        centroids[2,k] = sum(mesh.p[2,nds])/4
        centroids[3,k] = sum(mesh.p[3,nds])/4
    end

    # Boundary conditions
    RHS::Vector{Float64} = BoundaryIntegral(mesh,mu0.*Hext,shell_id)

    # Lagrange multiplier technique
    Lag::Vector{Float64} = lagrange(mesh)

    # Magnetic permeability
    mu::Vector{Float64} = zeros(mesh.nt) .+ mu0

    # FEM
    H_vectorField::Matrix{Float64} = zeros(mesh.nt,3)
    H::Vector{Float64} = zeros(mesh.nt)
    Hold::Vector{Float64} = zeros(mesh.nt)

    att::Int32 = 0
    div::Float64 = maxDeviation + 1.0
    while div > picardDeviation && att < maxAtt

        att += 1
        Hold .= H

        # Stiffness matrix
        A = stiffnessMatrix(mesh, mu)

        # Magnetic scalar potential
        u::Vector{Float64} = [A Lag;Lag' 0]\[-RHS;0]
        u = u[1:mesh.nv]

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
            spl = Spline1D(materialProperties[key].HofM,
                           materialProperties[key].mu)

            # Find all elements of current cell ID
            elements = findall(x -> x==id, elementID)
            
            # Interpolate the dataset for this elements
            mu[elements] .= spl(H[elements])
        end

        # Check deviation from previous result
        div = mu0*maximum(abs.(H.-Hold))
        println(att, " | mu0 |H(n)-H(n-1)| = ", div)

    end # Picard Iteration

    # Remove the lagrange multiplier
    pop!(u) # Remove last element of u

    # 



    
    
    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Magnetic field H")
    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        centroids[3,mesh.InsideElements], 
        color = mu0.*H[mesh.InsideElements], 
        colormap=:rainbow, 
        markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements]))

    Colorbar(fig[1, 2], scatterPlot, label="|H|") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig))
    
end # end of main

meshSize = 4.0
localSize = 0.1
showGmsh = false
saveMesh = false

main(meshSize,localSize,showGmsh,saveMesh)

