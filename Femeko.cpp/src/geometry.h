template<typename T>
void println(T input){
	std::cout << input << std::endl;
}

template<typename T>
void print(T input){
	std::cout << input;
}

// Add a center-based rectangle
int addRectangle(std::vector<double> &position, 
				 std::vector<double> &dimensions,
				 std::vector<std::pair<int, int>> &cells){

	double x = position[0] - dimensions[0]/2.0;
	double y = position[1] - dimensions[1]/2.0;

	int id = gmsh::model::occ::addRectangle(x, y, 0.0, dimensions[0], dimensions[1]);
    gmsh::model::occ::synchronize(); // Sync kernel before exiting

    cells.push_back({2, id}); // Store the new cell (dim, tag)
    return id;
}

// Add a center-based rectangle. Does not upodate a array of cell IDs
int addRectangle(std::vector<double> &position, 
				 std::vector<double> &dimensions){

	double x = position[0] - dimensions[0]/2.0;
	double y = position[1] - dimensions[1]/2.0;

	int id = gmsh::model::occ::addRectangle(x, y, 0.0, dimensions[0], dimensions[1]);
    gmsh::model::occ::synchronize(); // Sync kernel before exiting
    return id;
}

// Add center-based disk and include the new cell ID
int addDisk(std::vector<double> &position, double radius,
			std::vector<std::pair<int, int>> &cells){

	int id = gmsh::model::occ::addDisk(position[0], position[1], 0.0,
	                               	   radius, radius);

	gmsh::model::occ::synchronize(); // Sync kernel before exiting
    cells.push_back({2, id}); // Store the new cell (dim, tag)
	return id;
}

// Add center-based disk. Does not upodate a array of cell IDs
int addDisk(std::vector<double> &position, double radius){
	int id = gmsh::model::occ::addDisk(position[0], position[1], 0.0,
	                               	   radius, radius);
	gmsh::model::occ::synchronize(); // Sync kernel before exiting
	return id;
}

// Add a cuboid based on its center and include the new cell ID
int addCuboid(std::vector<double> &position, 
			  std::vector<double> &dimensions,
			  std::vector<std::pair<int, int>> &cells){

	double x = position[0] - dimensions[0]/2.0;
	double y = position[1] - dimensions[1]/2.0;
	double z = position[2] - dimensions[2]/2.0;

	int id = gmsh::model::occ::addBox(x, y, z, 
									 dimensions[0], 
									 dimensions[1], 
									 dimensions[2]);

    cells.push_back({3, id}); // Store the new cell (dim, tag)
	gmsh::model::occ::synchronize(); // Sync kernel before exiting
    return id;
} // add a cuboid and update the 'cells'

// Add a cuboid based on its center. Does not update the array of cells
int addCuboid(std::vector<double> &position, 
			  std::vector<double> &dimensions){

	double x = position[0] - dimensions[0]/2.0;
	double y = position[1] - dimensions[1]/2.0;
	double z = position[2] - dimensions[2]/2.0;

	int id = gmsh::model::occ::addBox(x, y, z, 
									 dimensions[0], 
									 dimensions[1], 
									 dimensions[2]);

	gmsh::model::occ::synchronize(); // Sync kernel before exiting
    return id;
}

// Add a sphere and upodate the array of cells
int addSphere(std::vector<double> &position, double radius,
			  std::vector<std::pair<int, int>> &cells){

	int id = gmsh::model::occ::addSphere(position[0], 
										 position[1],
										 position[2],
										 radius);

	cells.push_back({3, id});
	gmsh::model::occ::synchronize(); // Sync kernel before exiting
	return id;
}

// Add a sphere. Does not upodate the array of cells
int addSphere(std::vector<double> &position, double radius){

	int id = gmsh::model::occ::addSphere(position[0], 
										 position[1],
										 position[2],
										 radius);

	gmsh::model::occ::synchronize(); // Sync kernel before exiting
	return id;
}


// Fragment the geometry to make a conforming mesh
void unifyModel(std::vector<std::pair<int, int>> &cells,
				int box)
{ 

	int dim = cells[0].first; // Get dimension of the model

    // Unify model by fragment()
    std::vector<std::pair<int, int> > outDimTags;
    std::vector<std::vector<std::pair<int, int> > > tagsMap;

    gmsh::model::occ::fragment(cells, {{dim, box}}, outDimTags, tagsMap);
    gmsh::model::occ::synchronize();
    println("Fragmented geometry to create conforming mesh");

    // Update the cells to the new IDs
    cells = tagsMap[0];
}
