#=
    An example of how you can add a step file, create a bounding shell
    and a mesh of the object using Femeko.jl
=#

include("src/gmsh_wrapper.jl")
using LinearAlgebra

using GLMakie

function userMade(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)
    #=
        Make your own model with cuboids and spheres
    =#
    gmsh.initialize()

    # Cuboid dimensions
    L::Vector{Float64} = [5.0, 1.0, 1.0]

    # >> Model

    # List of cells inside the container
    cells = []
    addCuboid([-3*L[1]/4,0,0], L, cells, true)

    # Create a bounding shell
    box = addSphere([0,0,0], 5*maximum(L), cells, true) # *[1,1,1]

    # Fragment to make a unified geometry
    gmsh.model.occ.fragment(cells, [])
    gmsh.model.occ.synchronize()
    
    # Get bounding shell surface id
    shell_id = gmsh.model.getBoundary([(3, box)], false, false, false)
    shell_id = [s[2] for s in shell_id]

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)

    gmsh.fltk.run()
    gmsh.finalize()

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

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
    box, bounding_shell_n_surfaces = importCAD("STEP_Models/Fennec_Fox.step",cells)

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

    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

end # end of main

meshSize = 0.0
localSize = 0.0
showGmsh = true
saveMesh = false

userMade(meshSize,localSize,showGmsh,saveMesh)
# meshCAD(0,0,true,false)

