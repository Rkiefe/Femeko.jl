// Creates a mesh and outputs a struct with the data in a more familiar format

struct MESH2D
{
	int nv; // Number of nodes
	int nt; // Number of elements
	Eigen::MatrixXd p; // Mesh node coordinates | 3 by nv
	Eigen::MatrixXi t; // Mesh node connectivity| 3 by nt
};


// Set mesh extension to on/off (default off)
void extendLocalRefinement(double input){
	// 0 -> Dont extend mesh size, to prevent over-refinement inside entity
	// 1 -> Extend mesh refinement
	gmsh::option::setNumber("Mesh.MeshSizeFromPoints", input);
	gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", input);
	gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", input);
}

void setDistanceField(int distance_field, 
					  double meshSize, double localSize, 
					  double refineRange){

	// Create a threshold field that defines the refinement region
	int threshold_field = gmsh::model::mesh::field::add("Threshold");

	gmsh::model::mesh::field::setNumber(threshold_field, "InField", distance_field);
	gmsh::model::mesh::field::setNumber(threshold_field, "SizeMin", localSize);
	gmsh::model::mesh::field::setNumber(threshold_field, "SizeMax", meshSize);
	gmsh::model::mesh::field::setNumber(threshold_field, "DistMin", localSize);
	gmsh::model::mesh::field::setNumber(threshold_field, "DistMax", refineRange);

	// gmsh::model::mesh::field::setNumber(threshold_field, "Sigmoid", true);
	gmsh::model::mesh::field::setAsBackgroundMesh(threshold_field);
}

void refineCell(std::vector<std::pair<int, int>> &cell,
				double localSize, double meshSize){
	
	// Get the boundary of each cell
	std::vector<std::pair<int, int>> cell_boundary;
	gmsh::model::getBoundary(cell, cell_boundary, false, false, false);

	// Create a distance field for local refinement
	int distance_field = gmsh::model::mesh::field::add("Distance");

	if (cell_boundary[0].first < 2){ // 1 -> curves

		// Get the ID of each boundary
        std::vector<double> curves(cell_boundary.size()); // gmsh::...:setNumbers expects <double> for some reason
        for (int i = 0; i<cell_boundary.size(); i++){
            curves[i] = cell_boundary[i].second;
        }

        gmsh::model::mesh::field::setNumbers(distance_field, "CurvesList", curves);
    
    } else { // 2 -> faces
        
        // Get the ID of each boundary
        std::vector<double> faces(cell_boundary.size()); // gmsh::...:setNumbers expects <double> for some reason
        for (int i = 0; i<cell_boundary.size(); i++) {
            faces[i] = cell_boundary[i].second;
        }

        gmsh::model::mesh::field::setNumbers(distance_field, "FacesList", faces);
    }
    
    // Set number of sampling points over the surfaces
    gmsh::model::mesh::field::setNumber(distance_field, "Sampling", 500); // default is 100
    
    // Enforce the distance field by a threshold field
    setDistanceField(distance_field, meshSize, localSize, 2 * meshSize);

} // Local mesh refinement on target cells

// Generate 2D mesh
void Mesh2D(  MESH2D& mesh, // Populate this mesh struct
			  double meshSize, 
			  double localSize,
			  std::vector<std::pair<int, int>> &cells){

	// Set maximum element size
	gmsh::option::setNumber("Mesh.MeshSizeMax", meshSize);
	
	// Add local refinement
	if(cells.size()>0){
		refineCell(cells, localSize, meshSize);
	}
	
	// Generate the mesh
    gmsh::model::mesh::generate(2);
    gmsh::model::mesh::optimize();

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
