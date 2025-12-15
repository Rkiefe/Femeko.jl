#=
        Calculate the magnetostatic interaction between linear magnetic materials and a source field
    
    Using: 
        - linear Lagrange shape elements
        - The magnetostatic scalar potential
        - Magnetic permeability

    Why:
        - Simple to implement
        - Fast

    There is a more advanced Femeko implementation with Nedelec shape elements, which should be more stable
=#

include("../src/Femeko.jl")

using GLMakie

function main(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)
    #=
        Makes a model with cubes and spheres and refines the mesh on the spheres
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
            saveMesh  - Save mesh to a FEMCE compatible format 

    =#
    

    # Applied field
    mu0 = pi*4e-7      # vacuum magnetic permeability
    Hext::Vector{Float64} = [1.0, 0.0, 0.0]/mu0     # A/m

    # Dimensions
    L::Vector{Float64} = [20.0, 20.0, 20.0]
    
    # Relative magnetic permeability
    permeability::Float64 = mu0*(1 + 2) # mu0(1+chi)

    # Create a geometry
    gmsh.initialize()

    # List of cells inside the container
    cells = []

    # Import cad file
    # box = importCAD("../STEP_Models/cube.STEP", cells, true)

    # Or create your own geometry
    addCuboid([0,0,0], L, cells)
    box = addSphere([0,0,0], 5*maximum(L))

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id, _ = unifyModel(cells, box)
    
    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

    if showGmsh
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
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

    # FEM
    # Relative magnetic permeability 
    mu::Vector{Float64} = zeros(mesh.nt) .+ mu0
    mu[mesh.InsideElements] .= permeability

    # Boundary conditions
    RHS::Vector{Float64} = BoundaryIntegral(mesh, mu0*Hext, shell_id)

    # Lagrange multiplier technique
    Lag::Vector{Float64} = lagrange(mesh)

    # Stiffness matrix
    A = stiffnessMatrix(mesh, mu)

    # Magnetic scalar potential
    u::Vector{Float64} = [A Lag;Lag' 0]\[-RHS;0]
    u = u[1:mesh.nv]

    # Magnetic field
    H_vectorField::Matrix{Float64} = zeros(mesh.nt,3)
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
    H::Vector{Float64} = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = norm(H_vectorField[k,:])
    end

    B::Vector{Float64} = mu.*H

    # Magnetization
    chi::Float64 = permeability/mu0 - 1
    M_vectorField::Matrix{Float64} = chi.*H_vectorField[mesh.InsideElements,:]
    M::Vector{Float64} = chi.*H[mesh.InsideElements]

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Magnetic field H")
    scatterPlot = scatter!(ax, 
        centroids[1, mesh.InsideElements],
        centroids[2, mesh.InsideElements],
        centroids[3, mesh.InsideElements], 
        color = B[mesh.InsideElements], 
        colormap=:rainbow, 
        markersize=20) # 20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements])

    Colorbar(fig[1, 2], scatterPlot, label="|B| (T)") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    
    # save("H.png",fig)

end # end of main

meshSize = 30.0
localSize = 1.0
showGmsh = true

main(meshSize, localSize, showGmsh)

