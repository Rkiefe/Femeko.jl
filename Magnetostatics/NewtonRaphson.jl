#=
    Newton Raphson (or Newton Galerkin) iteration method
    for magnetostatics    
=#


include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

# For plots
# using CairoMakie

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
    Hext::Vector{Float64} = [1.35,0,0]./mu0    # A/m

    # Convergence criteria
    picardDeviation::Float64 = 1e-4
    maxDeviation::Float64 = 1e-6
    maxAtt::Int32 = 2
    relax::Float64 = 0.001 # Relaxation factor for N-R ]0, 1.0]

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
    idx = findall(x -> !isfinite(x), materialProperties["Gd"].mu)
    materialProperties["Gd"].mu[idx] .= 0.0
    materialProperties["Gd"].mu[idx] .= maximum(materialProperties["Gd"].mu)

    # d/dH mu
    dmu = gradient(materialProperties["Gd"].HofM, 
                   materialProperties["Gd"].B./materialProperties["Gd"].HofM)

    # Remove -Inf
    idx = findall(x -> !isfinite(x), dmu)
    dmu[idx] .= 0.0
    dmu[idx] .= minimum(dmu)

    materialProperties["Gd"].dmu = dmu

    
    # Load Iron data
    materialProperties["Fe"].HofM = vec(readdlm(folder*"Pure_Iron_FEMM/H_Fe_extrap.dat"))  # A/m
    materialProperties["Fe"].B = vec(readdlm(folder*"Pure_Iron_FEMM/B_Fe_extrap.dat"))     # T

    # Permeability
    materialProperties["Fe"].mu = materialProperties["Fe"].B./materialProperties["Fe"].HofM
    
    # Remove Inf
    idx = findall(x -> !isfinite(x), materialProperties["Fe"].mu)
    materialProperties["Fe"].mu[idx] .= 0.0
    materialProperties["Fe"].mu[idx] .= maximum(materialProperties["Fe"].mu)

    # d/dH mu
    dmu = gradient(materialProperties["Fe"].HofM, 
                   materialProperties["Fe"].B./materialProperties["Fe"].HofM)

    # Remove -Inf
    idx = findall(x -> !isfinite(x), dmu)
    dmu[idx] .= 0.0
    dmu[idx] .= minimum(dmu)

    materialProperties["Fe"].dmu = dmu

    # 3D Model
    gmsh.initialize()

    # Dimensions of each plate
    spacing::Float64 = 1.25
    thick::Float64 = 0.5
    Leng::Float64 = 16.5

    L::Vector{Float64} = [Leng, Leng, thick]

    # Create an empty container
    # box = addCuboid([0,0,0],5*maximum(L)*[1,1,1])
    box = addSphere([0,0,0],5*maximum(L))

    # Get how many surfaces compose the bounding shell
    temp = gmsh.model.getEntities(2)            # Get all surfaces of current model
    bounding_shell_n_surfaces = 1:length(temp)    # Get the number of surfaces in the bounding shell

    # List of cells inside the container
    cells = []
    cellLabels::Vector{String} = ["Air"]


    # Add Gadolinium plates

    # 1
    y::Float64 = (spacing+2*thick)/2
    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # 2
    y += spacing+2*thick

    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # 3
    y += spacing+2*thick

    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # 4
    y += spacing+2*thick

    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # 5
    y = -(spacing+2*thick)/2

    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # 6
    y -= spacing+2*thick

    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # 7
    y -= spacing+2*thick

    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # 8
    y -= spacing+2*thick

    addCuboid([0,0,y], L, cells, true)
    push!(cellLabels,"Gd")

    # Add Iron
    thickIron::Float64 = 1.0
    addCuboid([-(Leng+thickIron)/2, 0.0 , 0],
              [thickIron, L[2], 2.1*y],
              cells,true)

    push!(cellLabels,"Fe")

    addCuboid([(Leng+thickIron)/2, 0.0 , 0],
              [thickIron, L[2], 2.1*y],
              cells,true)

    push!(cellLabels,"Fe")

    


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

    # return

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
    u::Vector{Float64} = zeros(mesh.nv+1)
    
    H_vectorField::Matrix{Float64} = zeros(mesh.nt,3)
    H::Vector{Float64} = zeros(mesh.nt)
    Hold::Vector{Float64} = zeros(mesh.nt)
    
    att::Int32 = 0
    div::Float64 = maxDeviation + 1.0
    while div > picardDeviation && att < 1 # maxAtt 

        att += 1
        Hold .= H

        # Stiffness matrix
        A = stiffnessMatrix(mesh, mu)

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
    
    # Average magnetic field of Gd
    H_avg::Float64 = 0.0
    volume::Float64 = 0.0
    for i in 1:length(cells)

        # Cell ID
        id = cells[i][2]

        # Get the data set of current cell ID
        key = cellLabels[id]

        if key == "Gd"

            # Find all elements of current cell ID
            elements = findall(x -> x==id, elementID)
        
            for k in elements
                volume += mesh.VE[k]
                H_avg += H[k]*mesh.VE[k]
            end

        end # Only add the Gd elements

    end # <H>
    H_avg /= volume
    println("mu_0 <H> = ", mu0*H_avg)

    # Open a file for writing (creates or overwrites)
    open("output.txt", "w") do file
        println(file, "relax = ", relax)
        println(file, "mu0 <H> = ", mu0*H_avg)
    end

    # Plot result | Uncomment "using GLMakie"
    # fig = Figure()
    # ax = Axis3(fig[1, 1], aspect = :data, title="With relax = "*string(relax))
    # scatterPlot = scatter!(ax, 
    #     centroids[1,mesh.InsideElements],
    #     centroids[2,mesh.InsideElements],
    #     centroids[3,mesh.InsideElements], 
    #     color = H[mesh.InsideElements], 
    #     colormap=:rainbow, 
    #     markersize=20) # 20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements])

    # Colorbar(fig[1, 2], scatterPlot, label="|H|") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    # wait(display(fig))
    # save("relax_"*string(relax)*".png",fig)

end # end of main

meshSize = 82.5
localSize = 1.0
showGmsh = false
saveMesh = false

main(meshSize,localSize,showGmsh,saveMesh)

