/* 
	
	- compile with
	g++ main.cpp -o main.out -I ../gmsh-4.15.0-Linux64-sdk/include -L ../gmsh-4.15.0-Linux64-sdk/lib -l gmsh -Wl,-rpath,../gmsh-4.15.0-Linux64-sdk/lib
	
	- run with
	./main.out
*/

#include "../src/femeko.h"

Eigen::Vector3d linearBasis(Eigen::Ref<Eigen::MatrixXd> p,
							Eigen::Vector3i nodes,
							int nd){



	int nd1 = -1;
	int nd2 = -1;
	for(int i : nodes){
		if(i != nd){
			if(nd1<0){
				nd1 = i;
			}else if(nd2<0){
				nd2 = i;
			}
		}
	} // Find the nodes that are not 'nd'

	// Build the matrix equation for the a b c d
	Eigen::Matrix3d M(3, 3);
	M.row(0) << 1.0, p(0, nd),  p(1, nd);
	M.row(1) << 1.0, p(0, nd1), p(1, nd1);
	M.row(2) << 1.0, p(0, nd2), p(1, nd2);

	Eigen::Vector3d b(1.0, 0.0, 0.0);
	Eigen::Vector3d r = M.colPivHouseholderQr().solve(b);

	return r;
}


int main()
{
	gmsh::initialize();

	double meshSize = 1.0; 	// Maximum mesh size
	double localSize = 0.1; // Local element size
	bool showGmsh = false; // Open gmsh GUI ?

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

	// Show the cells that are inside the disk
	println("Cells inside the bounding shell:");
	for(std::pair<int, int> cell : cells){
		println(cell.second);
	}
	
	// Create the mesh
	extendLocalRefinement(0.0);
	MESH2D mesh;
	Mesh2D(mesh, meshSize, localSize, cells);

	// Print some mesh properties
	print("\nNumber of elements: ");
	println(mesh.nt);
	print("Number of nodes: ");
	println(mesh.nv);
	print("Number of elements in 'cells': ");
	println(mesh.nInside);

	if(showGmsh){ gmsh::fltk::run(); }
	gmsh::finalize();

	// int k = 0;
	// Eigen::Vector3d abc = linearBasis(mesh.p, mesh.t.col(k), mesh.t(0, k));
	
	// Get the linear basis function for each node of each element
	Eigen::MatrixXd P1a(3, mesh.nt);
	Eigen::MatrixXd P1b(3, mesh.nt);
	Eigen::MatrixXd P1c(3, mesh.nt);

	for(int k = 0; k<mesh.nt; k++){ // Each element
		for(int i = 0; i<3; i++){ // Each node of current element
			Eigen::Vector3d abc = linearBasis(mesh.p, mesh.t.col(k), mesh.t(i, k));
			
			P1a(i, k) = abc(0);
			P1b(i, k) = abc(1);
			P1c(i, k) = abc(2);
		}
	}

	// Build the stiffness matrix
	Eigen::MatrixXd Ak(9, mesh.nt); // Local stiffness matrix
	for(int k = 0; k<mesh.nt; k++){
		Eigen::Matrix3d aux = mesh.VE[k]*(P1b.col(k) * P1b.col(k).transpose() + 
										  P1c.col(k) * P1c.col(k).transpose());

		Ak.col(k) = Eigen::Map<Eigen::VectorXd>(aux.data(), 9);
	}



	return 0;
}