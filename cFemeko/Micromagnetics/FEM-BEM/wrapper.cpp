// Wrapper to C++ micromagnetics simulation to be called in julia

// #include "LandauLifshitz.cpp"
// #include "SteepestDescent.cpp"
#include "Micromagnetics.cpp"

extern "C"{
	
	// Run the landau lifshitz solver
	void LL(double* m_in,
			double* p_in,
		    int* t_in,
		    int* surfaceT_in,
		    double* normal_in,
		    double* AE, double* VE,
		    int nv, int nt, int ne,
		    double* m_out,
		    int maxAtt)
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

		// Map the input to the proper structure
		Eigen::Map<Eigen::Vector3d> m_result(m_out);


		// Create a microamgnetics object
		Micromagnetics micro(p, t, surfaceT, normal, AE, VE);
		
		// Set the properties
		micro.Ms 	 = 800e3;
		micro.Aexc 	 = 13e-12;
		micro.Aan 	 = 0.0;
		// micro.uan
		micro.Hext(0) = 0.1/micro.mu0;
		micro.maxAtt = maxAtt;

		// Run the solver
		micro.LandauLifshitz(m, p, t, surfaceT, normal, AE, VE);

		// Update the output
		for(int i = 0; i<3; i++)
		{
			m_out[i] = micro.M_avg(i,micro.att-1);
		}

	} // Wrapper to C++ Landau Lifshitz solver

} // extern C