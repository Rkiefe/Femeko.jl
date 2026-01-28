#include "LL.h"

// Contructor - Create the micromagnetics solver
LL::LL(Eigen::Ref<Eigen::MatrixXd> p_input, 
	   Eigen::Ref<Eigen::MatrixXi> t_input,
	   Eigen::Ref<Eigen::VectorXi> InsideElements_input,
	   Eigen::Ref<Eigen::VectorXi> InsideNodes_input,
	   double* VE_ptr,
	   Eigen::Ref<Eigen::MatrixXd> M_input)
	 : p(p_input), t(t_input), 
	   InsideElements(InsideElements_input), 
	   InsideNodes(InsideNodes_input), 
	   VE(VE_ptr), M(M_input)
{
	/*
		1) Stores locally the input mesh data
		2) Calculates the FEM basis function of each node and element
		3) Calculates the sparse stiffness matrix
		4) Pre-computes the conjugate gradient solver for that stiffness matrix
	*/

	// Get the linear basis function of each element and node
	std::cout << "Calculating the linear basis function over each element and node" << std::endl;
	linearBasis();

	// Stiffness matrix
	std::cout << "Building the sparse stiffness matrix" << std::endl;
	Eigen::SparseMatrix<double> A = stiffnessMatrix(p, t, VE);

	// Process the stiffness matrix with the linear solver
	std::cout << "Applying the linear solver pre-computation" << std::endl;
	CG.compute(A);

	// Compute the exchange stiffness matrix
	std::cout << "Calculating the exchange stiffness matrix" << std::endl;
	exchangeStiffness();

	// Run the micromagnetics solver
	run();

} // LL constructor

// Run the micromagnetic simulation
void LL::run(){

	// Get each magnetic field contribution
	magnetostaticField();
	exchangeField();
	anisotropyField();

	// Update the effective field
	updateEffectiveField(); // Updated 'H'

	// Prepare a copy of the magnetization field of current iteration
	Eigen::MatrixXd Mold = M; // Copy of M

	// Preapre a copy of H
	Eigen::MatrixXd Hold = H; // Copy of H

	// Norm of the magnetization on each node
	Eigen::VectorXd Mnorm = Eigen::VectorXd::Zero(p.cols());
	for(int i = 0; i<InsideNodes.size(); i++){
		int nd = InsideNodes(i);
		Mnorm(nd) = M.col(nd).norm();
	}

} // Run the micromagnetic solver

// FEM linear basis function of each node and element
void LL::linearBasis(){
	a = Eigen::MatrixXd::Zero(4, t.cols());
	b = Eigen::MatrixXd::Zero(4, t.cols());
	c = Eigen::MatrixXd::Zero(4, t.cols());
	d = Eigen::MatrixXd::Zero(4, t.cols());
	volumes = Eigen::VectorXd::Zero(p.cols());

	for(int k = 0; k<t.cols(); k++){
		for(int i = 0; i<4; i++){
			
			int nd = t(i, k); // Global mesh node label

			Eigen::Vector4d r = abcd(p, t.col(k), nd);
			a(i, k) = r(0);
			b(i, k) = r(1);
			c(i, k) = r(2);
			d(i, k) = r(3);

			volumes(nd) += VE[k];
		}
	} // Loop over the elements

} // FEM linear basis function of each node and element

// Get the stiffness matrix of the exchange field
void LL::exchangeStiffness(){

	// Create the local stiffness matrix
	Eigen::MatrixXd Ak = Eigen::MatrixXd::Zero(16, t.cols());
	for(int ik = 0; ik<InsideElements.size(); ik++){
		int k = InsideElements(ik);
		Eigen::Matrix4d aux = VE[k]*(b.col(k)*b.col(k).transpose() + c.col(k)*c.col(k).transpose() + d.col(k)*d.col(k).transpose());
		Ak.col(k) = Eigen::Map<Eigen::VectorXd>(aux.data(), 16);
	} // Local exchange stiffness matrix

	// Temporary storage for triplets (row, col, value)
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(16 * t.cols());

	// Update the global matrix
	for(int ik = 0; ik<InsideElements.size(); ik++){
		int k = InsideElements(ik);
		int n = -1;
		for (int i = 0; i < 4; ++i) {
		    for (int j = 0; j < 4; ++j) {
		        n++;
		        
		        // Add contribution to global matrix
		        triplets.emplace_back(t(i, k), t(j, k), Ak(n, k));
		    }
		}
	} // Loop over the elements

	// Build global stiffness matrix from triplets
	AEXC = Eigen::SparseMatrix<double>(p.cols(), p.cols());
	AEXC.setFromTriplets(triplets.begin(), triplets.end());
	// AEXC.makeCompressed();

} // Get the stiffness matrix of the exchange field

// Get the demagnetizing field
void LL::magnetostaticField(){
	
	int nv = p.cols();
	int nt = t.cols();
	
	Eigen::VectorXd RHS = Eigen::VectorXd::Zero(nv); // Load vector
	for(int ik = 0; ik<InsideElements.size(); ik++){ // sweep the elements in the magnetic volume
		
		int k = InsideElements(ik); // Global element label

		// Mean magnetization over each element node
		double mx, my, mz;
		mx = my = mz = 0.0;
		for(int i = 0; i<4; i++){
			int nd = t(i, k); // Global node label
			
			mx += M(0, nd);
			my += M(1, nd);
			mz += M(2, nd);
		}
		mx /= 4; my /= 4; mz /= 4;
		
		// Update load vector
		for(int i = 0; i<4; i++)
		{
			RHS(t(i, k)) += VE[k]*(b(i, k)*mx + c(i, k)*my + d(i, k)*mz);
		}
	} // Boundary conditions for the demagnetizing field

	// Magnetostatic scalar potential
	Eigen::VectorXd u = CG.solve(RHS); // CG.info() output 0 if CG succeeds

	// Calculate the demagnetizing field from the scalar potential
	Eigen::MatrixXd Hdk = Eigen::MatrixXd::Zero(3, nt);
	for(int k = 0; k<nt; k++){
		for(int i = 0; i<4; i++){
			
			int nd = t(i, k); // Global node label

			Hdk(0, k) -= u(nd)*b(i, k);
			Hdk(1, k) -= u(nd)*c(i, k);
			Hdk(2, k) -= u(nd)*d(i, k);
		}
	} // Demag field over the elements

	// Map the element-wise demag field to the mesh nodes
	Hd = Eigen::MatrixXd::Zero(3, nv);
	for(int k = 0; k<nt; k++){
		for(int i = 0; i<4; i++){
			int nd = t(i, k);
			Hd(0, nd) += Hdk(0, k)*VE[k];
			Hd(1, nd) += Hdk(1, k)*VE[k];
			Hd(2, nd) += Hdk(2, k)*VE[k];
		}
	}

	for(int nd = 0; nd<nv; nd++){
		Hd(0, nd) /= volumes(nd);
		Hd(1, nd) /= volumes(nd);
		Hd(2, nd) /= volumes(nd);
	} // Average over the nodes

} // Demagnetizing field from scalar potential

// Exchange Field (3 by nv) 
void LL::exchangeField(){
	double mu0 = pi*4e-7; // Vacuum magnetic permeability
	Hexc = Eigen::MatrixXd::Zero(3, p.cols());
	
	for(int i = 0; i<3; i++){
		double coef = -2*mu0*Aexc/(0.25 * Ms*Ms *scale*scale);
	    
	    Eigen::VectorXd temp = AEXC * M.row(i).transpose();
	    temp = temp.array() / volumes.array();
		Hexc.row(i) = (coef * temp).transpose();
	}
} // Exchange field over the mesh nodes

// Anisotropy field (3 by nv)
void LL::anisotropyField(){
	double mu0 = pi*4e-7; // Vacuum magnetic permeability
	Han = Eigen::MatrixXd::Zero(3, p.cols());
	for(int i = 0; i<InsideNodes.size(); i++){

		int nd = InsideNodes(i); // Global node label

		double Msqrd = M(0, nd)*M(0, nd) + M(1, nd)*M(1, nd) + M(2, nd)*M(2, nd); // |M|^2
	    double temp = mu0 * 2.0*Aan/Msqrd;
	    temp *= M(0, nd)*uan(0) + M(1, nd)*uan(1) + M(2, nd)*uan(2); // M dot easy-axis
	    
	    Han.col(nd) = temp * uan; // Update the field on the mesh node 'nd'
	} // Update the field on each mesh node
} // Anisotropy field

void LL::updateEffectiveField(){
	
	// Reset to all zeros
	H = Eigen::MatrixXd::Zero(3, p.cols());
	
	H.colwise() += Hext; // Add the external field
	H += Hd; 	// Add the demagnetizing field
	H += Hexc; 	// Add the exchange field
	H += Han;	// Add the anisotropy field

}

// Destructor
LL::~LL(){std::cout << "Micromagnetics solver leaving scope \n";}
