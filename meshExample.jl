#=
    An example of how you can add a step file, create a bounding shell
    and a mesh of the object using Femeko.jl
=#

include("src/Femeko.jl")
using GLMakie

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
    extendLocalRefinement(0) # Extend or not the local refinement from the boundary to the volume
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

    # Example of Femeko 3D mesh interpolation
    r = rand(3) .- 0.5                          # 3D vector from -0.5 to 0.5
    P = L[1]*r[1], L[2]*r[2], L[3]*r[3]         # Random coordinate inside L
    k = findElement3D(mesh, P[1], P[2], P[3])   # Find the element k that contains the random coordinate P
    nds = mesh.t[:, k]                          # nodes of the element found

    # Define an example solution over the nodes
    T = zeros(mesh.nv)
    for nd in nds
        T[nd] = mesh.p[1, nd]^2
    end

    # Interpolate the solution over the coordinate P
    w = interp3Dmesh(mesh, r[1], r[2], r[3], T)
    
    # See the 3D interpolation
    println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1,1]
              # , aspect=:data
              )
    
    graph = scatter!(ax
                    , [mesh.p[1, nds]; r[1]]
                    , [mesh.p[2, nds]; r[2]]
                    , [mesh.p[3, nds]; r[3]]
                    , color = [T[nds]; w]
                    )
    
    Colorbar(fig[1, 2], graph)

    wait(display(fig))


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