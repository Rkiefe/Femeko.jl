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

main(10,0.1,true)        # Create your own model and mesh with local refinement
testCAD(20,1,false)      # Import step file and make mesh with local refinement