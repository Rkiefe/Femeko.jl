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

    cells.push_back({2, id}); // Add the ID to the cell (dim, tag)
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
    cells.push_back({2, id}); // Add the ID to the cell (dim, tag)
	return id;
}

// Add center-based disk. Does not upodate a array of cell IDs
int addDisk(std::vector<double> &position, double radius){
	int id = gmsh::model::occ::addDisk(position[0], position[1], 0.0,
	                               	   radius, radius);
	gmsh::model::occ::synchronize(); // Sync kernel before exiting
	return id;
}
