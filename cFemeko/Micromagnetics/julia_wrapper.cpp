// Wrapper to C++ micromagnetics simulation to be called in julia

#include "../src/BEMc.cpp" // Include Femeko C++ FEM functions

extern "C"{
	
	// Get the demagnetizing field on the mesh nodes
	void demag(double* Hd_in,
			   double* p_in,
			   int* t_in,
			   int* surfaceT_in,
			   double* normal_in,
			   double* areaT,
			   int nv,
			   int nt,
			   int ne,
			   double* VE,
			   double* Vn,
			   double* m_in)
	{

		// Node coordinates
		Eigen::Map<Eigen::MatrixXd> p(p_in,3,nv);
		
		// Element connectivity
		Eigen::Map<Eigen::MatrixXi> t(t_in,4,nt);

		// Surface element node connectivity
		Eigen::Map<Eigen::MatrixXi> surfaceT(surfaceT_in,4,ne);
		
		// Surface normals
		Eigen::Map<Eigen::MatrixXd> normal(normal_in,3,ne);

		// Magnetization
		Eigen::Map<Eigen::MatrixXd> m(m_in,3,nv);

		// Map the input magnetostatic field
		Eigen::Map<Eigen::MatrixXd> Hd_out(Hd_in,3,nv);

		// Run
		BEMdmag(
				  Hd_out,
				  p, 
				  t, 
				  surfaceT, 
				  normal,
				  areaT,
				  VE,
				  Vn,
				  m);

	} // Wrapper to C++ demag field

} // extern C