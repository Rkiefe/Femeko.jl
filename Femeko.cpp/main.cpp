// Read the readme.md for instructions on how to setup and compile

#include "src/femeko.h"

int main()
{
	gmsh::initialize();

	double meshSize = 0.1; 	// Maximum mesh size
	double localSize = 0.0; // Local element size

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
		shellID = addDisk(position, 0.25);
	}

	{ // Fragment the geometry to make a conforming mesh
	    std::vector<std::pair<int, int> > outDimTags;
	    std::vector<std::vector<std::pair<int, int> > > outDimTagsMap;

	    // Fragment all objects in cells vector
	    gmsh::model::occ::fragment(cells, {{2, shellID}}, outDimTags, outDimTagsMap);
	    gmsh::model::occ::synchronize();
	    println("Fragmented geometry to create conforming mesh");
	}

	MESH2D mesh;
	Mesh2D(mesh, meshSize, localSize);

	println("Opening Gmsh GUI");
	gmsh::fltk::run();
	gmsh::finalize();
	return 0;
}