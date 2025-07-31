#=
    An example of how you can add a step file, create a bounding shell
    and a mesh of the object using Femeko.jl
=#

include("src/gmsh_wrapper.jl")

function userMade(meshSize=0,localSize=0,showGmsh=true,saveMesh=false)
    #=
        Make your own model with cuboids and spheres
    =#
    gmsh.initialize()

    # Cuboid dimensions
    L::Vector{Float64} = [1.65, 1.65, 0.04]

    # 3D Model

    # List of cells inside the container
    cells = []

    # Add a cuboid
    addCuboid([0,0,0], L, cells, true) # position, dimensions, cell list, update cell list

    # Create a bounding shell
    box = addSphere([0,0,0],5*maximum(L)) # Don't update the cell list and get the volume id

    # Fragment to make a unified geometry
    fragments, _ = gmsh.model.occ.fragment(vcat(cells,[(3, box)]), [])
    gmsh.model.occ.synchronize()

    # Update cell ids
    cells = fragments[1:end-1]
    println(cells)
    
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
    mesh = Mesh(cells, meshSize, localSize, saveMesh)

    if showGmsh # Show GUI
        gmsh.fltk.run()
    end
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
    box = importCAD("STEP_Models/Fennec_Fox.step",cells)

    # Fragment to make a unified geometry
    fragments, _ = gmsh.model.occ.fragment(vcat(cells,[(3, box)]), [])
    gmsh.model.occ.synchronize()

    # Update cell ids
    cells = fragments[1:end-1]
    println(cells)
    
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

    println(shell_id)

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

meshSize = 4
localSize = 0.1
showGmsh = true
saveMesh = false

userMade(meshSize, localSize, showGmsh, saveMesh)
meshCAD(0, 0, true, false)

