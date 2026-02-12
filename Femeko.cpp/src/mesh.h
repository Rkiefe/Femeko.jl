// Creates a mesh and outputs a struct with the data in a more familiar format

struct MESH2D
{
	int nv; // Number of nodes
	int nt; // Number of elements
	Eigen::MatrixXd p; // Mesh node coordinates | 3 by nv
	Eigen::MatrixXi t; // Mesh node connectivity| 3 by nt
};


void Mesh2D(MESH2D& mesh, // Populate this mesh struct
			  double meshSize, 
			  double localSize){

	{ // Mesh
		// Set maximum element size
		gmsh::option::setNumber("Mesh.MeshSizeMax", meshSize);

	    gmsh::model::mesh::generate(2);
	    gmsh::model::mesh::optimize();
	}
	println("\nFinished generating the mesh and optimizing it.");
	println("Extracting mesh info into Femeko format...");

	// Get the mesh nodes:
	{
		std::vector<std::size_t> nodes;
		std::vector<double> coord, coordParam;
		gmsh::model::mesh::getNodes(nodes, coord, coordParam);
		
		int nv = nodes.size();
		mesh.p = Eigen::Map<Eigen::MatrixXd>(coord.data(), 3, nv);
	} // Get mesh node coordinates (x,y,z) by nv

	// Get mesh connectivity:
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

	// return mesh;
}
