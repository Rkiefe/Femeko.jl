// Wrapper to C++ micromagnetics simulation to be called in julia

#include "../src/BEMc.cpp" // Include Femeko C++ FEM functions

// Calculate the magnetic scalar potential
Eigen::VectorXd BEMdmag(
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> t, 
	Eigen::Ref<Eigen::MatrixXi> surfaceT, 
	Eigen::Ref<Eigen::MatrixXd> normal,
	double* areaT,
	double* VE,
	Eigen::Ref<Eigen::MatrixXd> m)
{

	int nv = p.cols();
	int nt = t.cols();
	int ne = surfaceT.cols();

	// Make the FEM-BEM matrices

	// Stiffness matrix
	Eigen::MatrixXd A = denseStiffnessMatrix(p,t,VE); 		// nv by nv
	Eigen::MatrixXd B = Bmatrix(p,surfaceT,areaT); 			// nv by ne
	Eigen::MatrixXd C = Cmatrix(p,surfaceT,normal,areaT); 	// ne by nv
	Eigen::MatrixXd D = Dmatrix(p,surfaceT,areaT); 			// ne by ne

	// Extend the matrix
	// LHS = [-A B; C D]
	Eigen::MatrixXd LHS(nv+ne,nv+ne);
	LHS.block(0,0,nv,nv) 		= -A;
	LHS.block(0, nv, nv, ne) 	= B;
	LHS.block(nv, 0, ne, nv) 	= C;
	LHS.block(nv, nv, ne, ne) 	= D;

	// Righ hand side of FEM-BEM:
	Eigen::VectorXd RHS = Eigen::VectorXd::Zero(nv+ne);
	for (int s = 0; s<ne; s++)
	{
		// Average magnetization on the surface element
		Eigen::Vector3d mavg = (  m.col(surfaceT(0,s)) 
								+ m.col(surfaceT(1,s)) 
								+ m.col(surfaceT(2,s))
								)/3.0;
		
		double integral = areaT[s]/3.0 * mavg.dot(normal.col(s));

		// Update RHS
		for(int i = 0; i<3; i++){
			RHS(surfaceT(i,s)) -= integral; 
		}
	} // RHS of FEM-BEM magnetostatic potential linear eq.

	// Magnetostatic scalar potential | Solve LHS * u = RHS
    Eigen::VectorXd u = LHS.partialPivLu().solve(RHS);

    // Check solution
    if (LHS.partialPivLu().info() != Eigen::Success) {
        // std::cerr << "Solving failed!" << std::endl;
        throw std::runtime_error("Linear solver failed for FEM-BEM magnetostatic potential");
    }

	// Eigen::VectorXd u = Eigen::VectorXd::Zero(nv);

	return u;
} // Magnetostatic scalar potential

extern "C"{
	
	/*
		Takes Julia input arrays and converts to eigen vectors
		and std vectors
	*/ 
	void demag(double* u_in,
			   double* p_in,
			   int* t_in,
			   int* surfaceT_in,
			   double* normal_in,
			   double* areaT,
			   int nv,
			   int nt,
			   int ne,
			   double* VE,
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

		// Run
		Eigen::VectorXd u = BEMdmag(
			 p, 
			 t, 
			 surfaceT, 
			 normal,
			 areaT,
			 VE,
			 m);

		// Update the input u
		for (int i = 0; i < nv; i++)
		{
			// std::cout << u(i) << std::endl;
			u_in[i] = u(i);
		}

	} // Wrapper to C++ magnetic field simulation

} // extern C