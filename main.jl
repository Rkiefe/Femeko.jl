#=
    Mesh generation in Julia using Gmsh
        . Can make tetrahedral meshes
        . Can extract surface triangles from each cell of the model
        . Can make local mesh refinements based on target cell

        . Can import STEP files and make volumes and surfaces out of them
        . Can make an automatic container for the STEP file
        . Can make local refinement, for the entire STEP file
=#

using Gmsh
include("gmsh_wrapper.jl")

# Linear basis function, f = a + bx + cy + dz
function abcd(p,nodes,nd)
    n1,n2,n3 = nodes[nodes .!= nd]
    x = @view p[1, [nd, n1, n2, n3]]
    y = @view p[2, [nd, n1, n2, n3]]
    z = @view p[3, [nd, n1, n2, n3]]

    M = [1.0 x[1] y[1] z[1];
         1.0 x[2] y[2] z[2];
         1.0 x[3] y[3] z[3];
         1.0 x[4] y[4] z[4]]
    r = M\[1;0;0;0]

    return r[1],r[2],r[3],r[4] # a b c d
end # Basis function coef.

function stiffnessMatrix(mesh,chi)
    A = zeros(mesh.nv,mesh.nv)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]

        b = zeros(4,1); c = zeros(4,1); d = zeros(4,1);
        for i in 1:length(nds)
            _,b[i],c[i],d[i] = abcd(mesh.p,nds,nds[i])
        end
        # A[nds,nds] += chi[k].*mesh.VE[k].*(b*b' + c*c' + d*d')
    
        # Calculate the local matrix contribution
        local_mat = chi[k] * mesh.VE[k] * (b*b' + c*c' + d*d')
        
        # Add to the global sparse matrix
        for i in 1:length(nds)
            for j in 1:length(nds)
                A[nds[i], nds[j]] += local_mat[i,j]
            end
        end
    end

    return A
end # Stiffness matrix

function main(meshSize=0,localSize=0,saveMesh=false)
    #=
        Makes a model with cubes and spheres and refines the mesh on the spheres
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
            saveMesh  - Save mesh to a FEMCE compatible format 

    =#
    
    # Create a geometry
    gmsh.initialize()

    # >> Model
    # Create an empty container
    box = addCuboid([0,0,0],[2,2,4])

    # List of cells inside the container
    cells = []

    # Add 1st sphere
    addSphere([0,0,-1],0.5,cells)

    # Add 2nd sphere
    addSphere([0,0,1],0.5,cells)

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
    shell_id = shell_id[1:6] # All other, are interior surfaces

    println("Number of elements ",size(mesh.t,2))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))

    gmsh.fltk.run()
    gmsh.finalize()

    # Build the stiffness matrix considering each cell to be a magnetic material
    # with linear susceptibility
    chi = zeros(mesh.nt,1);
    chi[mesh.InsideElements] .= 3 # Magnetic susceptibility of the cells inside the bounding shell
    
    # Stiffness matrix
    A = @inbounds stiffnessMatrix(mesh,chi)
    
end

function testCAD(meshSize=0,localSize=0,saveMesh=false)
    #=
        Imports a BREP, STEP or IGES file, makes a container automatically,
        sets the mesh size, sets each volume (excluding the container)
        to have local refinement and creates the 3D mesh
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
            saveMesh  - Save mesh to a FEMCE compatible format 
    =#
    gmsh.initialize()

    # List of cells inside the container
    cells = []

    # Import step file
    box = importCAD("STEP_Models/cube.step",cells)

    # Fragment to make a unified geometry
    gmsh.model.occ.fragment([(3, box)], cells)
    gmsh.model.occ.synchronize()

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)

    println("Number of elements ",size(mesh.t,2))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))


    gmsh.fltk.run()
    gmsh.finalize()
end

meshSize    =   10  # Target mesh element size
localSize   =   0   # Local target mesh element sie
saveMesh    = false # Save the mesh properties to .txt files

main(meshSize,localSize,saveMesh)    # Create your own model and mesh with local refinement
testCAD(20,1,false)                  # Import step file and make mesh with local refinement
