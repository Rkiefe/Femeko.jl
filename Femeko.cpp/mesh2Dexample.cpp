// Read the readme.md for instructions on how to setup and compile

#include "src/femeko.h"

void refineCell(std::vector<std::pair<int, int>> &cell,
				double localSize, double meshSize){
	
	// Get the boundary of each cell
	std::vector<std::pair<int, int>> cell_boundary;
	gmsh::model::getBoundary(cell, cell_boundary, false, false, false);

	// Create a distance field for local refinement
	int distance_field = gmsh::model::mesh::field::add("Distance");

	// if (cell_boundary[0].first < 2) { // 1 -> curves
	//         std::vector<int> curves;
	//         for (const auto& boundary : cell_boundary) {
	//             curves.push_back(boundary.second);
	//         }
	//         gmsh::model::mesh::field::setNumbers(distance_field, "CurvesList", curves);
	//     } else { // 2 -> faces
	//         std::vector<int> faces;
	//         for (const auto& boundary : cell_boundary) {
	//             faces.push_back(boundary.second);
	//         }
	//         gmsh::model::mesh::field::setNumbers(distance_field, "FacesList", faces);
	//     }
	    
	//     // Set number of sampling points over the surfaces
	//     gmsh::model::mesh::field::setNumber(distance_field, "Sampling", 500); // default is 100
	    
	//     // Enforce the distance field by a threshold field
	//     // You'll need to implement or include the setDistanceField function
	//     setDistanceField(distance_field, meshSize, localSize, 2 * meshSize);
}


int main()
{
	gmsh::initialize();

	double meshSize = 1.0; 	// Maximum mesh size
	double localSize = 0.1; // Local element size

	// Hold the label of each cell added
	std::vector<std::pair<int, int>> cells;

	{ // Add a rectangle
		std::vector<double> position = {0.0, 0.0};
		std::vector<double> dimensions = {2.0, 1.0};
		addRectangle(position, dimensions, cells);
	}

	// Add a disk as bounding shell
	int shellID;
	{ 
		std::vector<double> position = {0.0, 0.0};
		shellID = addDisk(position, 4.0);
	}

	{ // Fragment the geometry to make a conforming mesh
	    std::vector<std::pair<int, int> > outDimTags;
	    std::vector<std::vector<std::pair<int, int> > > outDimTagsMap;

	    // Fragment all objects in cells vector
	    gmsh::model::occ::fragment(cells, {{2, shellID}}, outDimTags, outDimTagsMap);
	    gmsh::model::occ::synchronize();
	    println("Fragmented geometry to create conforming mesh");

	    // Update the cells to the new IDs
	    cells = outDimTags;

	    // // Check the outDimTags
	    // for(auto dimTag : outDimTags){
	    // 	println(dimTag.first);
	    // 	println(dimTag.second);
	    // }
	}

	// Add local refinement
	refineCell(cells, localSize, meshSize);

	MESH2D mesh;
	Mesh2D(mesh, meshSize, localSize);

	println("Opening Gmsh GUI");
	gmsh::fltk::run();
	gmsh::finalize();
	return 0;
}