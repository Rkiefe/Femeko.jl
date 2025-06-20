// Wrapper to C++ micromagnetics simulation to be called in julia

#include "../src/BEMc.cpp" // Include Femeko C++ FEM functions

// Calculate the magnetic scalar potential
void BEMdmag(
	Eigen::Ref<Eigen::MatrixXd> Hd,
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> t, 
	Eigen::Ref<Eigen::MatrixXi> surfaceT, 
	Eigen::Ref<Eigen::MatrixXd> normal,
	double* areaT,
	double* VE,
	double* Vn,
	Eigen::Ref<Eigen::MatrixXd> m)
{

	/*		Inputs
		Hd 		 -> Input and output  | Magnetostatic field, 3 by nv
		p  		 -> Nodes coordiantes | 3 by nv
		t  		 -> Volume element node connectivity | 4 by nt
		surfaceT -> Surface element node connectivity | 4 by ne (3 nodes + boundary label)
		normal 	 -> Normal to surface element | 3 by ne
		areaT 	 -> Area of surface element
		VE 		 -> Volume of volume element
		Vn 		 -> Total volume of elements with node i | Vector of length nv
	*/ 

	int nv = p.cols(); 			// Number of mesh nodes
	int nt = t.cols(); 			// Number of volume elements
	int ne = surfaceT.cols();   // Number of surface elements

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

    // Calculate the Magnetostatic field
    Eigen::MatrixXd Hdk = Eigen::MatrixXd::Zero(3,nt);
    for (int k = 0; k<nt; k++)
    {
    	for(int i = 0; i<4; i++)
    	{
	    	Eigen::Vector4d r = abcd(p,t.col(k),t(i,k));
	    	// phi = a + bx + cy + dz
    		Hdk(0,k) -= u(t(i,k))*r(1);
    		Hdk(1,k) -= u(t(i,k))*r(2);
    		Hdk(2,k) -= u(t(i,k))*r(3);
    	}
    } // Magnetostatic field on each element of the mesh

    // Update the input magnetostatic field on the nodes 
    // by the average contributions
    
    for(int k = 0; k<nt; k++)
    {
    	// Add to each node, the element contribution
    	for(int i = 0; i<4; i++)
    	{
    		int nd = t(i,k);
    		Hd.col(nd) += VE[k]*Hdk.col(k);
    	}

    } // Element contribution to the magnetic field on the nodes

    // Divide by the total volume of each node
    for(int i = 0; i<nv; i++)
    {
    	Hd.col(i) /= Vn[i];
    }

} // Magnetostatic scalar potential

extern "C"{
	
	/*
		Takes Julia input arrays and converts to eigen vectors
		and std vectors
	*/ 
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

		// // Update the input Hd
		// for (int i = 0; i < nv; i++)
		// {
		// 	Hd_out.col(i) = Hd.col(i);
		// }

	} // Wrapper to C++ magnetic field simulation

} // extern C