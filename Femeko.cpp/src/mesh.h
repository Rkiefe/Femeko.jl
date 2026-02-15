// Creates a mesh and outputs a struct with the data in a more familiar format

struct MESH
{
	int nv; // Number of nodes
	int nt; // Number of elements
	Eigen::MatrixXd p; // Mesh node coordinates | 3 by nv
	Eigen::MatrixXi t; // Mesh node connectivity| 3 by nt

	int ns; // Number of boundary elements
	Eigen::MatrixXi surfaceT; // Boundary elements connectivity (3 by ns)

	std::vector<int> InsideElements; // Elements tags for the volumes in 'cells'
	int nInside; // Number of elements in 'InsideElements'

	Eigen::MatrixXd normal; // Normal to boundary element

	std::vector<double> VE; // Volume of each element (triangle area in 2D)
	std::vector<double> AE; // Area of each element (line length in 2D)
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

// Area of surface triangle
double areaTriangle(Eigen::Ref<Eigen::MatrixXd> p, int i0, int i1, int i2){
	// Area by cross product 
	Eigen::Vector3d AB = p.col(i1) - p.col(i0);
	Eigen::Vector3d AC = p.col(i2) - p.col(i0);

	Eigen::Vector3d cross = AB.cross(AC);
	double area = 0.5 * cross.norm();

    return area;
} // Area of the 3D triangle

Eigen::Vector2d normalEdge(Eigen::Ref<Eigen::MatrixXd> p, 
						   Eigen::Vector3i nds){ // size 3 because an edge holds the Boundary ID

	double x = p(0, nds[1]) - p(0, nds[0]);
	double y = p(1, nds[1]) - p(1, nds[0]);

	double norm = std::sqrt(x*x + y*y);

	Eigen::Vector2d n = {-y/norm, x/norm}; // normal to edge
	return n;
}

// Volume of the tetrahedron devined by 4 points
double elementVolume(Eigen::Ref<Eigen::MatrixXd> p, 
					 Eigen::Vector4i nds){

	    // A = p[:, nds[0]]
	    // B = p[:, nds[1]]
	    // C = p[:, nds[2]]
	    // D = p[:, nds[3]]

	    // Compute vectors AB, AC, AD
	    std::vector<double> AB = {p(0, nds[1]) - p(0, nds[0]), 
	    					      p(1, nds[1]) - p(1, nds[0]),
	    					      p(2, nds[1]) - p(2, nds[0])}; // B - A
	    
	    std::vector<double> AC = {p(0, nds[2]) - p(0, nds[0]), 
	    					      p(1, nds[2]) - p(1, nds[0]),
	    					      p(2, nds[2]) - p(2, nds[0])};  // C - A
	    
		std::vector<double> AD = {p(0, nds[3]) - p(0, nds[0]), 
							      p(1, nds[3]) - p(1, nds[0]),
							      p(2, nds[3]) - p(2, nds[0])};  // D - A

	    // Compute AB ⋅ (AC × AD)
	    std::vector<double> cross_AC_AD = {AC[1]*AD[2] - AC[2]*AD[1],
	                   					   AC[2]*AD[0] - AC[0]*AD[2],
	                   					   AC[0]*AD[1] - AC[1]*AD[0]}; 

	    double aux = AB[0] * cross_AC_AD[0] + AB[1] * cross_AC_AD[1] + AB[2] * cross_AC_AD[2];

	    // Volume = (1/6) * |aux|
	    return std::abs(aux)/6.0;
} // Volume of tetrahedron

// Generate 2D mesh
void Mesh2D(  MESH& mesh, // Populate this mesh struct
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
		
		mesh.nv = nodes.size();
		mesh.p = Eigen::Map<Eigen::MatrixXd>(coord.data(), 3, mesh.nv);
	} // Get mesh node coordinates (x,y,z) by nv

	// Get mesh connectivity:
	{
		std::vector<std::size_t> elementTags, elementNodeTags;
		gmsh::model::mesh::getElementsByType(2, elementTags, elementNodeTags); // type 2 for p1 triangles

		mesh.nt = elementTags.size();
		mesh.t = Eigen::MatrixXi::Zero(3, mesh.nt);
		int n = 0;
		for(int k = 0; k<mesh.nt; k++){ // For each element
			for(int i = 0; i<3; i++){ // For each node
				mesh.t(i, k) = elementNodeTags[n]-1; // Gmsh node labels start at 1
				n++;
			}
		}

		// Get the element tags in 'cells'
		mesh.nInside = 0;
		mesh.InsideElements.reserve(mesh.nt); // At most this vector holds every tag

		if(cells.size()>0){ // Check which elements belong to 'cells' 

			for(int k = 0; k<mesh.nt; k++){
				
				size_t tag = elementTags[k]; // Input | element tag

		        int elementType; // output | type of element
		        std::vector<std::size_t> nodeTags; // Output | nodes of the element
		        int dim; // Output: geometric dimension of the element
		        int id; // Output: ID of the entity where element is classified
		        
		        // Get the cell ID of the current element 
		        gmsh::model::mesh::getElement(tag, elementType, nodeTags, dim, id);

		        // Sweep each cell ID and see if it matches
		        for(int cell = 0; cell<cells.size(); cell++){
					int cellID = cells[cell].second; // (dim, tag)
					if(cellID == id){ // Found a cell with the element tag
						mesh.InsideElements.push_back(k);
						mesh.nInside++;
					}
		        }
			
			} // Loop over the elements

		} // Check what element tags are in 'cells'
	
	} // Mesh connectivity (3 by nt)

	// Get the boundary connectivity (edges)
	{
		// Get the boundary element tags
		std::vector<std::size_t> elementTags, elementNodeTags;
		gmsh::model::mesh::getElementsByType(1, elementTags, elementNodeTags); // type 1 for p1 edges

		mesh.ns = elementTags.size();
		mesh.surfaceT = Eigen::MatrixXi::Zero(3, mesh.ns);
		for(int s = 0; s<mesh.ns; s++){
			size_t tag = elementTags[s]; // Input | element tag

			int elementType; // output | type of element
			std::vector<std::size_t> nodeTags; // Output | nodes of the element
			int dim; // Output: geometric dimension of the element
			int id; // Output: ID of the entity where element is classified

			// Get the cell ID of the current element 
			gmsh::model::mesh::getElement(tag, elementType, nodeTags, dim, id);

			// Store the nodes and the boundary ID in mesh.surfaceT
			mesh.surfaceT(0, s) = nodeTags[0]-1; // Gmsh labels starts at 1
			mesh.surfaceT(1, s) = nodeTags[1]-1; // ...
			mesh.surfaceT(2, s) = id;
		}
	}

	// Area of each element
	mesh.VE.assign(mesh.nt, 0.0); 
	for(int k = 0; k<mesh.nt; k++){
		mesh.VE[k] = areaTriangle(mesh.p, 
								  mesh.t(0, k), 
								  mesh.t(1, k),
								  mesh.t(2, k));
	}

	// Length and normal of each boundary element
	mesh.normal = Eigen::MatrixXd::Zero(2, mesh.ns);
	mesh.AE.assign(mesh.ns, 0.0);
	for(int s = 0; s<mesh.ns; s++){

		int i = mesh.surfaceT(0, s);
		int j = mesh.surfaceT(1, s);

		double dx = mesh.p(0, j) - mesh.p(0, i);
		double dy = mesh.p(1, j) - mesh.p(1, i);

		mesh.AE[s] = std::sqrt(dx*dx + dy*dy);
		mesh.normal.col(s) = normalEdge(mesh.p, mesh.surfaceT.col(s)); 
	}

} // Mesh2D()

// Generate 2D mesh
void Mesh( MESH& mesh, // Populate this mesh struct
		   double meshSize, 
		   double localSize,
		   std::vector<std::pair<int, int>> &cells){

	// Verify the dimension of the model
	if(cells.size()>0 && cells[0].first < 3){ // It is a 2D geometry
		Mesh2D(mesh, meshSize, localSize, cells);
		return;		
	}

	// Set maximum element size
	gmsh::option::setNumber("Mesh.MeshSizeMax", meshSize);
	
	// Add local refinement
	if(cells.size()>0){
		refineCell(cells, localSize, meshSize);
	}

	// Generate the mesh
    gmsh::model::mesh::generate(3);
    gmsh::model::mesh::optimize();

	println("\nFinished generating the mesh and optimizing it.");
	println("Extracting mesh info into Femeko format...");

	// Get the mesh nodes:
	{
		std::vector<std::size_t> nodes;
		std::vector<double> coord, coordParam;
		gmsh::model::mesh::getNodes(nodes, coord, coordParam);
		
		mesh.nv = nodes.size();
		mesh.p = Eigen::Map<Eigen::MatrixXd>(coord.data(), 3, mesh.nv);
	} // Get mesh node coordinates (x,y,z) by nv

	// Get mesh connectivity:
	{
		std::vector<std::size_t> elementTags, elementNodeTags;
		gmsh::model::mesh::getElementsByType(4, elementTags, elementNodeTags); // type 4 for p1 tetrahedrons

		mesh.nt = elementTags.size();
		mesh.t = Eigen::MatrixXi::Zero(4, mesh.nt);
		int n = 0;
		for(int k = 0; k<mesh.nt; k++){ // For each element
			for(int i = 0; i<4; i++){ // For each node
				mesh.t(i, k) = elementNodeTags[n]-1; // Gmsh node labels start at 1
				n++;
			}
		}

		// Get the element tags in 'cells'
		mesh.nInside = 0;
		mesh.InsideElements.reserve(mesh.nt); // At most this vector holds every tag

		if(cells.size()>0){ // Check which elements belong to 'cells' 

			for(int k = 0; k<mesh.nt; k++){
				
				size_t tag = elementTags[k]; // Input | element tag

		        int elementType; // output | type of element
		        std::vector<std::size_t> nodeTags; // Output | nodes of the element
		        int dim; // Output: geometric dimension of the element
		        int id; // Output: ID of the entity where element is classified
		        
		        // Get the cell ID of the current element 
		        gmsh::model::mesh::getElement(tag, elementType, nodeTags, dim, id);

		        // Sweep each cell ID and see if it matches
		        for(int cell = 0; cell<cells.size(); cell++){
					int cellID = cells[cell].second; // (dim, tag)
					if(cellID == id){ // Found a cell with the element tag
						mesh.InsideElements.push_back(k);
						mesh.nInside++;
					}
		        }
			
			} // Loop over the elements

		} // Check what element tags are in 'cells'
	
	} // Mesh connectivity (3 by nt)

	// Get the boundary connectivity (triangles)
	{
		// Get the boundary element tags
		std::vector<std::size_t> elementTags, elementNodeTags;
		gmsh::model::mesh::getElementsByType(2, elementTags, elementNodeTags); // type 1 for p1 edges

		mesh.ns = elementTags.size();
		mesh.surfaceT = Eigen::MatrixXi::Zero(4, mesh.ns); // 3 nodes + boundary ID
		for(int s = 0; s<mesh.ns; s++){
			size_t tag = elementTags[s]; // Input | element tag

			int elementType; // output | type of element
			std::vector<std::size_t> nodeTags; // Output | nodes of the element
			int dim; // Output: geometric dimension of the element
			int id; // Output: ID of the entity where element is classified

			// Get the cell ID of the current element 
			gmsh::model::mesh::getElement(tag, elementType, nodeTags, dim, id);

			// Store the nodes and the boundary ID in mesh.surfaceT
			mesh.surfaceT(0, s) = nodeTags[0]-1; // Gmsh labels starts at 1
			mesh.surfaceT(1, s) = nodeTags[1]-1; // ...
			mesh.surfaceT(2, s) = nodeTags[2]-1; // ...
			mesh.surfaceT(3, s) = id;
		}
	}

	mesh.VE.assign(mesh.nt, 0.0); 
	for(int k = 0; k<mesh.nt; k++){
	    mesh.VE[k] = elementVolume(mesh.p, mesh.t.col(k));
	} // Volume of each element

	// Area and normal of each boundary element
	mesh.normal = Eigen::MatrixXd::Zero(2, mesh.ns);
	mesh.AE.assign(mesh.ns, 0.0);
	for(int s = 0; s<mesh.ns; s++){

		mesh.AE[s] = areaTriangle(mesh.p, 
								  mesh.surfaceT(0, s), 
								  mesh.surfaceT(1, s),
								  mesh.surfaceT(2, s));

		// mesh.normal.col(s) = normalEdge(mesh.p, mesh.surfaceT.col(s)); 
	}

} // Mesh2D()