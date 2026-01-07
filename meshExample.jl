#=
    An example of how you can add a step file, create a bounding shell
    and a mesh of the object using Femeko.jl
=#

include("src/Femeko.jl")

function userMade(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)
    #=
        Make your own model with cuboids and spheres
    =#
    gmsh.initialize()

    # Cuboid dimensions
    L::Vector{Float64} = [2.0, 2.0, 2.0]

    # List of volume cells
    cells = []

    # Add a cuboid
    addCuboid([0,0,0], L, cells) # position, dimensions, cell list, update cell list

    # Create a bounding shell
    box = addSphere([0,0,0], 5*maximum(L)) # Don't update the cell list and get the volume id

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id, box = unifyModel(cells, box)

    # Volume surface ids
    internal_surfaces = gmsh.model.getBoundary(cells, false, false, false) # (dim, tag)
    internal_surfaces = [s[2] for s in internal_surfaces] # tag

    # Generate Mesh
    meshExtend(0) # Extend or not the local refinement from the boundary to the volume
    mesh = Mesh(cells, meshSize, localSize, saveMesh)

    println("\nOuter shell ID: ", shell_id)
    println("Internal surfaces: ", internal_surfaces)
    println("Number of elements ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of nodes ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ne)
    println("")

    if showGmsh # Show GUI
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()

end

function meshCAD(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)
    #=
        Makes a model with cubes and spheres and refines the mesh on the spheres
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
            saveMesh  - Save mesh to a FEMCE compatible format 
    =#
    
    # Create a geometry
    gmsh.initialize()

    # List of cells inside the container
    cells = []

    # Import cad file
    box = importCAD("STEP_Models/Fennec_Fox.step", cells, true)

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id = unifyModel(cells, box)

    # Volume surface ids
    internal_surfaces = gmsh.model.getBoundary(cells, false, false, false) # (dim, tag)
    internal_surfaces = [s[2] for s in internal_surfaces] # tag

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)

    println("\nOuter shell ID: ", shell_id)
    println("Internal surfaces: ", internal_surfaces)
    println("Number of elements ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of nodes ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ne)
    println("")
    
    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

end # end of main

meshSize = 4
localSize = 0.1
showGmsh = true
saveMesh = false

userMade(meshSize, localSize, showGmsh, saveMesh)
meshCAD(0, 0, true, false)