#=
    Newton Raphson (or Newton Galerkin) iteration method
    for magnetostatics    
=#


include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")
include("../src/magneticProperties.jl")

# For plots
using GLMakie

# To read data from files
using DelimitedFiles

# Wrapper to Fortran dierckx | Interpolation functions
using Dierckx

using IterativeSolvers



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
    Hext::Vector{Float64} = 1.2.*[1,0,0]./mu0    # A/m

    # Convergence criteria
    picardDeviation::Float64 = 1e-4
    maxDeviation::Float64 = 1e-6
    maxAtt::Int32 = 10
    relax::Float64 = 1.0 # Relaxation factor for N-R ]0, 1.0]

    # Data of magnetic materials
    folder::String = "Materials/"

    materialProperties = Dict("Gd" => DATA(),
                              "Fe" => DATA())
    # Load Gd
    loadMaterial( materialProperties,
                       "Materials",     # Folder with materials
                       "Gd_MFT",        # Data folder of target material
                       "Gd",            # Material name
                       7.9,
                       T)

    
    # Load Iron data
    materialProperties["Fe"].HofM = vec(readdlm(folder*"Pure_Iron_FEMM/H_Fe_extrap.dat"))  # A/m
    materialProperties["Fe"].B = vec(readdlm(folder*"Pure_Iron_FEMM/B_Fe_extrap.dat"))     # T

    # Get the permeability and its derivative
    materialPermeability(materialProperties, "Fe")



    # 3D Model
    gmsh.initialize()
    
    # Array of volume cell IDs
    cells = []

    # Import cad file
    # box = importCAD("../STEP_Models/cube.STEP", cells)
    # cellLabels = ["Gd"]
    # push!(cellLabels, "Air")


    # Add material
    addCuboid([0,0,0], [1.0, 1.0, 1.0], cells, true)
    cellLabels::Vector{String} = ["Gd"]

    # Create bounding shell
    box = addSphere([0,0,0], 5.0)
    push!(cellLabels, "Air")

    # Fragment to make a unified geometry
    fragments, _ = gmsh.model.occ.fragment(vcat(cells,[(3, box)]), [])
    gmsh.model.occ.synchronize()

    # Update cell ids
    cells = fragments[1:end-1]
    
    # Set the box to the last volume
    box = fragments[end][2]

    # Get bounding shell surface id
    shell_id = gmsh.model.getBoundary([(3, box)], false, false, false) # (dim, tag)
    shell_id = [s[2] for s in shell_id] # tag

    # Volume surface ids
    internal_surfaces = gmsh.model.getBoundary(cells, false, false, false) # (dim, tag)
    internal_surfaces = [s[2] for s in internal_surfaces] # tag

    shell_id = setdiff(shell_id, internal_surfaces) # Only the outer surfaces

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)

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

        # Check condition number
        # println("Condition: ", cond(Matrix([A Lag; Lag' 0])))

        # Magnetic scalar potential
        u = [A Lag;Lag' 0]\[-RHS;0]
        # u = minres([A Lag; Lag' 0], [-RHS; 0]
                    # ) # , verbose=true

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
        # du = minres([A+At Lag;Lag' 0], [-RHS-A*u;0])
        
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

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="With relax = "*string(relax))
    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        centroids[3,mesh.InsideElements], 
        color = mu0.*H[mesh.InsideElements], 
        colormap=:rainbow, 
        markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements]))
        # markersize=20)
    Colorbar(fig[1, 2], scatterPlot, label="|H|") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig))
    # save("relax_"*string(relax)*".png",fig)

end # end of main

meshSize  = 5.0
localSize = 0.1
showGmsh = true
saveMesh = false

main(meshSize,localSize,showGmsh,saveMesh)

