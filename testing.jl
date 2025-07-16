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
    addCuboid([-2*L[1]/5,0,0], L, cells, true)

    # Create a bounding shell
    box = addSphere([0,0,0], 5*maximum(L)) # *[1,1,1]

    # Fragment to make a unified geometry
    fragments, _ = gmsh.model.occ.fragment(vcat(cells,[(3, box)]), [])
    gmsh.model.occ.synchronize()

    # Update cell ids
    cells = fragments[1:end-1]
    println(cells)
    
    # Set the box to the last volume
    box = fragments[end][2]
    println(box)

    # Get bounding shell surface id
    shell_id = gmsh.model.getBoundary([(3, box)], false, false, false)
    shell_id = [s[2] for s in shell_id]

    # Volume surface ids
    internal_surfaces = gmsh.model.getBoundary(cells, false, false, false)
    internal_surfaces = [s[2] for s in internal_surfaces]

    println(shell_id)
    println(internal_surfaces)    
    println(setdiff(shell_id, internal_surfaces))

    shell_id = setdiff(shell_id, internal_surfaces)

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)

    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))

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
    ax = Axis3(fig[1, 1], aspect = :data, title="Mesh")
    scatterPlot = scatter!(ax, 
        centroids[1,mesh.InsideElements],
        centroids[2,mesh.InsideElements],
        centroids[3,mesh.InsideElements], 
        color = :lightblue, 
        markersize=20)
        # markersize=20 .* mesh.VE[mesh.InsideElements]./maximum(mesh.VE[mesh.InsideElements]))
        # colormap=:rainbow, 

    # Colorbar(fig[1, 2], scatterPlot) # Add a colorbar
    
    # Display the figure (this will open an interactive window)
    wait(display(fig)) # This is required only if runing outside the repl
    

end

meshSize = 25.0
localSize = 1.0
showGmsh = true
saveMesh = false

userMade(meshSize,localSize,showGmsh,saveMesh)
# meshCAD(0,0,true,false)

