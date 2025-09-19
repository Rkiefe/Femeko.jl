#=
    An example on how to create a 2D model and generate a mesh
    with local refinement
=#

include("src/gmsh_wrapper.jl")
include("src/FEM.jl")

function main(meshSize=0.0, localSize=0.0, showGmsh=false)
	gmsh.initialize()

	cells = []

	# Add a 2D rectangle
	id = addRectangle([0,0,0], [2, 2], cells)

	# Add a container
	# box = addRectangle([0,0,0], [2, 4])
    box = addDisk([0,0,0], 4)

	# Combine the geometries
	gmsh.model.occ.fragment(vcat(cells,[(2,box)]), [])
	gmsh.model.occ.synchronize()

    # Generate mesh
	mesh::MESH = Mesh2D(cells, meshSize, localSize)

	# Run Gmsh GUI
    if showGmsh
	   gmsh.fltk.run()
    end
	gmsh.fltk.finalize()

end

main(4, 0.1, true)