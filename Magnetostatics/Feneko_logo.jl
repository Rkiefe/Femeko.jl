include("../gmsh_wrapper.jl")
include("../FEM.jl")

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
    Hext::Vector{Float64} = [1,0,0]     # T

    # Dimensions
    L::Vector{Float64} = [1.65,1.65,0.04]
    # L::Vector{Float64} = [1,1,1]
    
    # Relative magnetic permeability
    permeability::Float64 = 3

    # Create a geometry
    gmsh.initialize()

    # List of cells inside the container
    cells = []

    # Import cad file
    box, bounding_shell_n_surfaces = importCAD("../STEP_Models/Fennec_Fox.step",cells)

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

    # if showGmsh
        gmsh.fltk.run()
    # end
    gmsh.finalize()

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

    # FEM

    # Relative magnetic permeability 
    mu::Vector{Float64} = ones(mesh.nt);
    mu[mesh.InsideElements] .= permeability

    # Boundary conditions
    RHS = BoundaryIntegral(mesh,Hext,shell_id)

    # Lagrange multiplier technique
    Lag = lagrange(mesh)

    # Stiffness matrix
    A = stiffnessMatrix(mesh,mu) 
    
    # Extend the matrix for the Lagrange multiplier technique
    mat = [A Lag;Lag' 0]

    # Magnetic scalar potential
    u = mat\[-RHS;0]
    u = u[1:mesh.nv]

    #= Example of calling C++ to calculate the local stiffness matrix
        
        # Julia local stiffness matrix
        Ak::Matrix{Float64} = @time localStiffnessMatrix(mesh,mu)

        # C++ Local stiffness matrix
        t::Matrix{Int32} = mesh.t .- 1 # C++ index starts at 0 # C++ index starts at 0
        Ak::Matrix{Float64} = @time CstiffnessMatrix(mesh.p,t,mesh.VE,mu)
    
    =#

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

    # Magnetization
    chi::Vector{Float64} = mu .- 1;
    M_vectorField::Matrix{Float64} = zeros(mesh.nInside,3)
    M::Vector{Float64} = zeros(mesh.nInside)
    for ik in 1:mesh.nInside
        k = mesh.InsideElements[ik]
        
        M_vectorField[ik,1] = chi[k]*H_vectorField[k,1]
        M_vectorField[ik,2] = chi[k]*H_vectorField[k,2]
        M_vectorField[ik,3] = chi[k]*H_vectorField[k,3]

        M[ik] = chi[k]*H[k]
    end

    # Element centroids
    centroids::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        centroids[1,k] = sum(mesh.p[1,nds])/4
        centroids[2,k] = sum(mesh.p[2,nds])/4
        centroids[3,k] = sum(mesh.p[3,nds])/4
    end

    # Plot result | Uncomment "using GLMakie"
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Feneko.jl")
    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        centroids[3,mesh.InsideElements], 
        color = H[mesh.InsideElements], 
        colormap=:rainbow
        ) # markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements])

    Colorbar(fig[1, 2], scatterPlot, label="H field strength") # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    
    save("Feneko_logo.png",fig)

end # end of main

meshSize = 5*60
localSize = 5
showGmsh = true
saveMesh = false

main(meshSize,localSize,showGmsh,saveMesh)

