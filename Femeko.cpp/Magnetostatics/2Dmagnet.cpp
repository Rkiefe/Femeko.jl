/* 
	Calculates the demagnetizing field of a given magnetization field M
	
	- compile with
	g++ 2Dmagnet.cpp -o 2Dmagnet.out -I ../gmsh-4.15.0-Linux64-sdk/include -L ../gmsh-4.15.0-Linux64-sdk/lib -l gmsh -Wl,-rpath,../gmsh-4.15.0-Linux64-sdk/lib
	
	- run with
	./main.out
*/

#include "../src/femeko.h"
#include <fstream> // Save results to .txt

int main()
{
	gmsh::initialize();

	// Mesh settings
	double meshSize = 1.0; 	// Maximum mesh size
	double localSize = 0.1; // Local element size
	bool showGmsh = false; // Open gmsh GUI ?

	// Applied field
	double mu0 = pi*4e-7; // Vacuum magnetic permeability
	Eigen::Vector2d Bext = {1.0, 0.0}; // Not used in this simulation

	// Linear solver
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;

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

	// Centroids of each element
	Eigen::MatrixXd centers(2, mesh.nt);
	for(int k = 0; k<mesh.nt; k++){
		centers.col(k) << 0.0, 0.0;
		for(int i = 0; i<3; i++){
			int nd = mesh.t(i, k);
			centers(0, k) += mesh.p(0, nd)/3.0;
			centers(1, k) += mesh.p(1, nd)/3.0;
		} // Average position of the nodes of the element
	}

	// Magnetic permeability
	std::vector<double> mu(mesh.nt, 1.0);

	// Magnetization field
	Eigen::MatrixXd M = Eigen::MatrixXd::Zero(2, mesh.nt);
	for(int k : mesh.InsideElements){
		M(0, k) = 1.0;
	}

	print("\nBuilding the P1 coefficients for each node of each element... ");	
	// Get the linear basis function for each node of each element
	Eigen::MatrixXd P1a(3, mesh.nt);
	Eigen::MatrixXd P1b(3, mesh.nt);
	Eigen::MatrixXd P1c(3, mesh.nt);

	for(int k = 0; k<mesh.nt; k++){ // Each element
		for(int i = 0; i<3; i++){ // Each node of current element
			Eigen::Vector3d abc = linearBasis2D(mesh.p, mesh.t.col(k), mesh.t(i, k));
			
			P1a(i, k) = abc(0);
			P1b(i, k) = abc(1);
			P1c(i, k) = abc(2);
		}
	}
	println("Done.");

	print("Building the element-wise stiffness matrix... ");
	// Local stiffness matrix (dense)
	Eigen::MatrixXd Ak = localStiffnessMatrix2D(mesh, P1b, P1c);
	println("Done.");

	// Global stiffness matrix (sparse)
	Eigen::SparseMatrix<double> A = stiffnessMatrix2D(mesh, Ak, mu);

	// Build the load vector
	Eigen::VectorXd RHS = Eigen::VectorXd::Zero(mesh.nv);
	for(int k = 0; k<mesh.nt; k++){
		for(int i = 0; i<3; i++){
			int nd = mesh.t(i, k);
			RHS[nd] += ( M(0, k)*P1b(i, k) + M(1, k)*P1c(i, k) ) * mesh.VE[k];
		}
	}

	// Solve for the magnetic scalar potential
	print("Calculating the magnetostatic potential by the Conjugate Gradient method... ");
	Eigen::VectorXd u = solver.compute(A) .solve(RHS);
	println("Done.");

	// Calculate the magnetic field from the scalar potential
	Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(2, mesh.nt);
	for(int k = 0; k<mesh.nt; k++){

		// Sum contributions
		for(int i = 0; i<3; i++){
			int nd = mesh.t(i, k);

			Hd(0, k) -= u[nd]*P1b(i, k);
			Hd(1, k) -= u[nd]*P1c(i, k);
		}
	}

	print("Saving to .txt files... ");
	{ // Magnetization
		std::ofstream file("M.txt");
		if (file.is_open()) {
		    file << M << std::endl;
		    file.close();
		}
	}
	{ // Mesh centroids
		std::ofstream file("centers.txt");
		if (file.is_open()) {
		    file << centers << std::endl;
		    file.close();
		}
	}
	{ // Magnetic field
		std::ofstream file("Hd.txt");
		if (file.is_open()) {
		    file << Hd << std::endl;
		    file.close();
		}
	}
	println("Done.");

	return 0;
}