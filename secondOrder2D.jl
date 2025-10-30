#=
    An example on how to generate a 2nd order 2D triangular mesh
=#

include("src/Femeko.jl")

using GLMakie

function main(meshSize=0.0, showGmsh=false)
    
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.ElementOrder", 2)       # Set to quadratic
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 1)  # Dont conform at the boundary

    # Add object
    cells = []
    addDisk([0.0, 0.0], 0.5, cells)

    # Create a FEMEKO mesh data structure
    mesh = Mesh(cells, meshSize, 0.0, false, 2) # 2nd order triangular mesh
    
    println("")
    println("Number of elements ", mesh.nt)
    println("Number of nodes ", mesh.nv)
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("")

    if showGmsh
        gmsh.fltk.run()
    end
    gmsh.finalize()

    # Plot a single volume element
    nds = mesh.t[:, 1]
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    scatter!(ax, 
             mesh.p[1, nds], 
             mesh.p[2, nds], 
             mesh.p[3, nds])
    
    # The first 3 nodes are the nodes of the linear triangle
    scatter!(ax, 
             mesh.p[1, nds[1:3]], 
             mesh.p[2, nds[1:3]], 
             mesh.p[3, nds[1:3]])
    
    wait(display(fig))

end

main(0.1, true)