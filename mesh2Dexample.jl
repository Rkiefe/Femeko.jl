#=
    An example on how to create a 2D model and generate a mesh
    with local refinement

    The model is a square with two holes inside and a bounding shell
    around the square

    Note: You set the holes last. 
    	. First add the square
    	. Then the bounding shell
    	. Unify the geometry so that the mesh process sees both geometries as 'one model'
		. Then add the holes
	Note: The holes only see the 'cells', this way you can't occCut the 'box'

=#

include("src/Femeko.jl")

function main(meshSize=0.0, localSize=0.0, showGmsh=false)
	gmsh.initialize()

	cells = []
	box = -1

	# Add a 2D rectangle
	id_rec = addRectangle([0,0,0], [4, 4], cells)

	# Add a bounding shell
	box = addDisk([0,0,0], 7.0)
	shell_id, box = unifyModel(cells, box)

	# Add a hole
    id = addDisk([-0.75,0,0], 0.2)
	cells = occCut(cells, [2, id])

    # Add another hole
    id = addDisk([0.75,0,0], 0.2)
    cells = occCut(cells, [2, id])

    # Generate mesh
	mesh::MESH = Mesh2D(cells, meshSize, localSize)

	println("\nBox cell ID: ", box)
	println("Inner cells: ", cells)
	println("Outer boundary curve ID: ", shell_id, "\n")

	# Run Gmsh GUI
    if showGmsh
	   gmsh.fltk.run()
    end
	gmsh.fltk.finalize()

end

meshSize = 2.0
localSize = 0.1
showGmsh = true

main(meshSize, localSize, showGmsh)