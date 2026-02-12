// Read the readme.md for instructions on how to setup and compile

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
	Eigen::MatrixXd p;
	{
		std::vector<std::size_t> nodes;
		std::vector<double> coord, coordParam;
		gmsh::model::mesh::getNodes(nodes, coord, coordParam);
		
		int nv = nodes.size();
		p = Eigen::Map<Eigen::MatrixXd>(coord.data(), 3, nv);
	}


	println("Opening Gmsh GUI");
	gmsh::fltk::run();
	gmsh::finalize();
	return 0;
}