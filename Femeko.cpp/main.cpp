// Read the readme.md for instructions on how to setup and compile

#include "src/femeko.cpp"

struct MESH2D
{
	int nv; // Number of nodes
	int nt; // Number of elements
	Eigen::MatrixXd p; // Mesh node coordinates | 3 by nv
	Eigen::MatrixXi t; // Mesh node connectivity| 3 by nt
};

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

	MESH2D mesh;

	// Get the mesh nodes:
	// Eigen::MatrixXd p;
	{
		std::vector<std::size_t> nodes;
		std::vector<double> coord, coordParam;
		gmsh::model::mesh::getNodes(nodes, coord, coordParam);
		
		int nv = nodes.size();
		mesh.p = Eigen::Map<Eigen::MatrixXd>(coord.data(), 3, nv);
	} // Get mesh node coordinates (x,y,z) by nv

	// Get mesh connectivity:
	// Eigen::MatrixXi t;
	{
		std::vector<std::size_t> elementTags, elementNodeTags;
		gmsh::model::mesh::getElementsByType(2, elementTags, elementNodeTags); // is type 2 for p1 triangles?

		int nt = elementTags.size();
		mesh.t = Eigen::MatrixXi::Zero(3, nt);
		int n = 0;
		for(int k = 0; k<nt; k++){ // For each element
			for(int i = 0; i<3; i++){ // For each node
				mesh.t(i, k) = elementNodeTags[n];
				n++;
			}
		}
	} // Mesh connectivity (3 by nt)
	

	println("Opening Gmsh GUI");
	gmsh::fltk::run();
	gmsh::finalize();
	return 0;
}