// Read the readme.md for instructions on how to setup and compile

#include "src/femeko.h"

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
	MESH2D mesh;
	Mesh2D(mesh, meshSize, localSize);

	// println("Opening Gmsh GUI");
	gmsh::fltk::run();
	gmsh::finalize();
	return 0;
}