#include "LL.cpp"

extern "C"{

	// Sends raw Julia pointers to a C++ object
	void LandauLifshitz(double* p_ptr,  // Node coordinates 3 by nv
						int* t_ptr, 	// Node connectivity 3 by nt
						double* VE, 	// Volume of each element
						int* InsideElements_ptr, // Array of element labels in magnetic volume
						int* InsideNodes_ptr, 	 // Array of node labels in magnetic volume
						int nv, int nt,
						int nInside, int nInsideNodes,
						double* M_ptr) // Magnetization field
	{

		// Map the mesh info to Eigen
		Eigen::Map<Eigen::MatrixXd> p(p_ptr, 3, nv);
		Eigen::Map<Eigen::MatrixXi> t(t_ptr, 4, nt);
		Eigen::Map<Eigen::VectorXi> InsideElements(InsideElements_ptr, nInside);
		Eigen::Map<Eigen::VectorXi> InsideNodes(InsideNodes_ptr, nInsideNodes);

		// Map the magnetization vector field
		Eigen::Map<Eigen::MatrixXd> M(M_ptr, 3, nv);

		// Create the micromagnetics solver object
		LL solver(p, t, InsideElements, InsideNodes, VE, M);
		solver.run(); // Run the simulation

	} // Sends raw Julia pointers to a C++ object

} // Julia to C++ wrapper 
