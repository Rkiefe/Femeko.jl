#=
    Simulation:
        Simulates the magnetostatic interaction between a magnetic plate
        and a uniform external magnetic field
        The plate has a uniform, constant magnetic permeability
    
    Output:
        Expect a 3D tetrahedral mesh as a MESH() struct and one figure
        of the internal magnetic field of your geometry.
        
    
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

# For plots | Uncomment the plot section of "main()"
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
    Hext::Vector{Float64} = [1.35,0,0]     # T

    # Dimensions
    L::Vector{Float64} = [20.0, 20.0, 20.0]
    
    # Relative magnetic permeability
    permeability::Float64 = 2

    # Create a geometry
    gmsh.initialize()

    # List of cells inside the container
    cells = []

    # Import cad file
    # box = importCAD("../STEP_Models/cube.STEP", cells)

    # Or create your own geometry
    addCuboid([0,0,0], L, cells, true)
    box = addSphere([0,0,0], 5*maximum(L))

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

    # FEM
    # Relative magnetic permeability 
    mu::Vector{Float64} = ones(mesh.nt);
    mu[mesh.InsideElements] .= permeability

    # Boundary conditions
    RHS::Vector{Float64} = BoundaryIntegral(mesh,Hext,shell_id)

    # Lagrange multiplier technique
    Lag::Vector{Float64} = lagrange(mesh)

    # Stiffness matrix
    A = stiffnessMatrix(mesh,mu) 

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
    chi::Float64 = (permeability - 1)/mu0;
    M_vectorField::Matrix{Float64} = chi.*H_vectorField[mesh.InsideElements,:]
    M::Vector{Float64} = chi.*H[mesh.InsideElements]

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Magnetic field H")
    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        centroids[3,mesh.InsideElements], 
        color = B[mesh.InsideElements], 
        colormap=:rainbow, 
        markersize=20) # 20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements])

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    
    # save("H.png",fig)

end # end of main

meshSize = 30.0
localSize = 1.0
showGmsh = true
saveMesh = false

main(meshSize,localSize,showGmsh,saveMesh)

