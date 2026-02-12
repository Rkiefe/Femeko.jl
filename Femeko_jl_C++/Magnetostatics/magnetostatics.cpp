// Wrapper to C++ magnetostatic simulation to be called in julia
#include "../../src/FEMc.cpp"

// Magnetostatic simulation of a paramagnet
Eigen::VectorXd magnetostatics(
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> t, 
	Eigen::Ref<Eigen::MatrixXi> surfaceT, 
	Eigen::Ref<Eigen::MatrixXd> normal, 
	double* VE,
	double* mu,
	double* Hext,
	int shell_id)
{

	int nv = p.cols();
	int nt = t.cols();

	// Run the boundary integral
	Eigen::VectorXd RHS = BoundaryIntegral(p, surfaceT, normal, Hext, shell_id);
	
	// Global sparse stiffness matrix
	Eigen::SparseMatrix<double> A = stiffnessMatrix(p, t, VE, mu);

    // Solve using conjugate gradient
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Conjugate gradient failed to compute stiffness matrix");
    }

    // Scalar potential
    Eigen::VectorXd u = solver.solve(-RHS);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Conjugate gradient method failed to solve the linear system");
    }

    return u;
}

extern "C"{
	
	/*
		Takes Julia input arrays and converts to eigen vectors
		and std vectors
	*/ 
	void scalarPotential(double* u_in, 		// scalar potential on each node
						 double* p_in, 		// Node coordinates 3 by nv
						 int* t_in,			// Element connectivity 4 by nt
						 int* surfaceT_in,	// Surface connectivity 3 by ns
						 double* normal_in, // Surface normals 3 by ns
						 int nv, // Number of nodes
						 int nt, // Number of elements
						 int ns, // Number of surface elements
						 double* VE, // Volume of each element
						 double* mu, // Permeability of each element
						 double* Hext, // External field
						 int shell_id) // ID of the bounding shell
	{

		// Node coordinates
		Eigen::Map<Eigen::MatrixXd> p(p_in, 3, nv);
		
		// Element connectivity
		Eigen::Map<Eigen::MatrixXi> t(t_in, 4, nt);

		// Surface element node connectivity
		Eigen::Map<Eigen::MatrixXi> surfaceT(surfaceT_in, 4, ns);
		
		// Surface normals
		Eigen::Map<Eigen::MatrixXd> normal(normal_in, 3, ns);

		// Source field
		// std::vector<double> F(Hext, Hext + 3);   // Map the input vector to a std::vector with pointer arithmetic (pointer + size)

		// Run
		Eigen::VectorXd u = magnetostatics(
			 p, 
			 t, 
			 surfaceT, 
			 normal,
			 VE,
			 mu, 
			 Hext,
			 shell_id);

		// Update the input u
		for (int i = 0; i < nv; i++)
		{
			u_in[i] = u(i);
		}

	} // Wrapper to C++ magnetic field simulation

} // extern C