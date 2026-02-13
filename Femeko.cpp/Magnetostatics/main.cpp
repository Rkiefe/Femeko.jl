/* 
	
	- compile with
	g++ main.cpp -o main.out -I ../gmsh-4.15.0-Linux64-sdk/include -L ../gmsh-4.15.0-Linux64-sdk/lib -l gmsh -Wl,-rpath,../gmsh-4.15.0-Linux64-sdk/lib
	
	- run with
	./main.out
*/

#include "../src/femeko.h"

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

	unifyModel(cells, shellID);

	println("Cells inside the bounding shell:");
	for(std::pair<int, int> cell : cells){
		println(cell.second);
	}
	
	// Add local refinement
	refineCell(cells, localSize, meshSize);

	// Create the mesh
	extendLocalRefinement(0.0);

	MESH2D mesh;
	Mesh2D(mesh, meshSize, localSize);

	// println("Opening Gmsh GUI");
	gmsh::fltk::run();
	gmsh::finalize();


	
	return 0;
}