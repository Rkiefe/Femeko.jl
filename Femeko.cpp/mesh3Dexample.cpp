/* !! Read the readme.md for instructions on how to setup and compile

	This creates a rectangle inside a disk and make a conforming mesh
	with local refinement on the boundaries if the rectangle

	- compile with
	g++ mesh3Dexample.cpp -o mesh3Dexample.out -I gmsh-4.15.0-Linux64-sdk/include -L gmsh-4.15.0-Linux64-sdk/lib -l gmsh -Wl,-rpath,gmsh-4.15.0-Linux64-sdk/lib
	
	- run with
	./mesh3Dexample.out
*/

#include "src/femeko.h"

int main()
{
	gmsh::initialize();

	double meshSize = 1.0; 	// Maximum mesh size
	double localSize = 0.1; // Local element size
	bool showGmsh = true; // Open gmsh GUI ?

	// Hold the label of each cell added
	std::vector<std::pair<int, int>> cells;

	{ // Add a prism
		std::vector<double> position = {0.0, 0.0, 0.0};
		std::vector<double> dimensions = {2.0, 1.0, 1.0};
		addCuboid(position, dimensions, cells);
	}

	// Add a disk as bounding shell
	int box;
	{ 
		std::vector<double> position = {0.0, 0.0, 0.0};
		box = addSphere(position, 4.0);
	}

	unifyModel(cells, box);

	// Show the cells that are inside the disk
	println("Cells inside the bounding shell:");
	for(std::pair<int, int> cell : cells){
		println(cell.second);
	}
	
	// Create the mesh
	// extendLocalRefinement(0.0);
	MESH mesh;
	Mesh(mesh, meshSize, localSize, cells);

	// Print some mesh properties
	print("\nNumber of elements: ");
	println(mesh.nt);
	
	print("Number of nodes: ");
	println(mesh.nv);

	// print("Number of boundary elements: ");
	// println(mesh.ns);

	// print("Number of elements in 'cells': ");
	// println(mesh.nInside);

	if(showGmsh){ gmsh::fltk::run(); }
	gmsh::finalize();

	return 0;
}