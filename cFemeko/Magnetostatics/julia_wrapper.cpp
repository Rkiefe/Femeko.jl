// Wrapper to C++ magnetostatic simulation to be called in julia

#include "../../src/FEMc.cpp" // Include Femeko C++ FEM functions

// Magnetostatic simulation of a paramagnet
Eigen::VectorXd magnetostatics(
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> t, 
	Eigen::Ref<Eigen::MatrixXi> surfaceT, 
	Eigen::Ref<Eigen::MatrixXd> normal, 
	double* VE,
	double* mu,
	std::vector<double>& F,
	int shell_id)
{

	int nv = p.cols();
	int nt = t.cols();

	// Run the boundary integral
	Eigen::VectorXd RHS = BoundaryIntegral(p, surfaceT, normal, F, shell_id);
	
	// Run the lagrange multiplier technique
	Eigen::VectorXd Lag = lagrange(t, VE, nv, nt);

	// Global sparse stiffness matrix
	Eigen::SparseMatrix<double> A = stiffnessMatrix(p, t, VE, mu);

	// Create the extended matrix
	int n = A.rows();
	int newSize = n + 1;
	std::vector<Eigen::Triplet<double>> triplets;
	
	// Reserve space: non-zeros from A + 2 * non-zeros in Lag + bottom-right zero
	triplets.reserve(A.nonZeros() + 2*Lag.size() + 1);

	// Add original matrix A
	for (int k = 0; k < A.outerSize(); ++k) {
	    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
	        triplets.emplace_back(it.row(), it.col(), it.value());
	    }
	}

	// Add Lag (i, n) and Lag' (n, i)
	for (int i = 0; i < n; ++i) {
	    triplets.emplace_back(i, n, Lag(i));
	    triplets.emplace_back(n, i, Lag(i));
	}

	// Explicitly set bottom-right element to 0
	triplets.emplace_back(n, n, 0.0);

	// Build the sparse matrix
	Eigen::SparseMatrix<double> mat(newSize, newSize);
	mat.setFromTriplets(triplets.begin(), triplets.end());
	mat.makeCompressed();  // Optimize storage

	// Extend the RHS | [-RHS; 0]
    Eigen::VectorXd RHS_ext(n + 1);
    RHS_ext.head(n) = -RHS;
    RHS_ext(n) = 0.0;

    // Solve using SparseLU
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(mat);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }

    // Magnetostatic potential
    Eigen::VectorXd u = solver.solve(RHS_ext);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed");
    }

    return u;
}

extern "C"{
	
	/*
		Takes Julia input arrays and converts to eigen vectors
		and std vectors
	*/ 
	void cMagnetoStatics(double* u_in,
						double* p_in,
						int* t_in,
						int* surfaceT_in,
						double* normal_in,
						int nv,
						int nt,
						int ne,
						double* VE,
						double* mu,
						double* F_in,
						int shell_id)
	{

		// Node coordinates
		Eigen::Map<Eigen::MatrixXd> p(p_in,3,nv);
		
		// Element connectivity
		Eigen::Map<Eigen::MatrixXi> t(t_in,4,nt);

		// Surface element node connectivity
		Eigen::Map<Eigen::MatrixXi> surfaceT(surfaceT_in,4,ne);
		
		// Surface normals
		Eigen::Map<Eigen::MatrixXd> normal(normal_in,3,ne);

		// Source field
		std::vector<double> F(F_in, F_in + 3);  
		// Example of source field {1.0,2.0,3.0}


		// Run
		Eigen::VectorXd u = magnetostatics(
			 p, 
			 t, 
			 surfaceT, 
			 normal,
			 VE,
			 mu, 
			 F,
			 shell_id);

		// Update the input u
		for (int i = 0; i < nv; i++)
		{
			// std::cout << u(i) << std::endl;
			u_in[i] = u(i);
		}

	} // Wrapper to C++ magnetic field simulation

} // extern C