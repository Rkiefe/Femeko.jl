// Wrapper to the C++ code, to be called in Julia

#include "FEMc.cpp"

extern "C"{

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