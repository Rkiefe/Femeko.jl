/*
	1st) get the gmsh sdk
		wget http://gmsh.info/bin/Linux/gmsh-4.15.0-Linux64-sdk.tgz
	
	2nd) unzip
		tar -xzf gmsh-4.15.0-Linux64-sdk.tgz
	
	3rd) COMPILE WITH
		g++ main.cpp -o main.out -I gmsh-4.15.0-Linux64-sdk/include -L gmsh-4.15.0-Linux64-sdk/lib -l gmsh -Wl,-rpath,gmsh-4.15.0-Linux64-sdk/lib 

*/
#include <iostream>
#include <vector>

#include "src/femeko.cpp"

int main()
{
	gmsh::initialize();

	// Hold the label of each cell added
	std::vector<int> cells;

	{ // Add a rectangle
		std::vector<double> position = {0.0, 0.0};
		std::vector<double> dimensions = {2.0, 1.0};
		addRectangle(position, dimensions, cells);
	}

	{ // Add a disk
		std::vector<double> position = {0.0, 0.0};
		addDisk(position, 0.25);
	}

	{ // Mesh
		double meshSize = 0.1;
		double localSize = 0.0;
		
		// Set maximum element size
		gmsh::option::setNumber("Mesh.MeshSizeMax", meshSize);

	    gmsh::model::mesh::generate(2);
	    gmsh::model::mesh::optimize();
	}
	println("\nFinished generating the mesh and optimizing it.");
	println("Extracting mesh info into Femeko format...");

	// Get the mesh nodes:
	std::vector<std::size_t> nodes;
	std::vector<double> coord, coordParam;
	gmsh::model::mesh::getNodes(nodes, coord, coordParam);

	int nv = nodes.size();

	println("Opening Gmsh GUI");
	gmsh::fltk::run();
	gmsh::finalize();
	return 0;
}