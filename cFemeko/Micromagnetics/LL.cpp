#include "LL.h"

// Contructor - Create the micromagnetics solver
LL::LL(Eigen::Ref<Eigen::MatrixXd> p_input, 
	   Eigen::Ref<Eigen::MatrixXi> t_input,
	   Eigen::Ref<Eigen::VectorXi> InsideElements_input,
	   Eigen::Ref<Eigen::VectorXi> InsideNodes_input,
	   double* VE_ptr,
	   Eigen::Map<Eigen::MatrixXd> M_input)
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

	std::cout << "Found an out of bounds access when the mesh has local refinement. Do not run at this stage." << std::endl;

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
	std::cout << "Running the magnetic field contributions" << std::endl;
	
	std::cout << "Running the demag field" << std::endl;
	magnetostaticField();

	std::cout << "Running the exchange field" << std::endl;
	exchangeField();

	std::cout << "Running the anisotropy field" << std::endl;
	anisotropyField();

	// Update the effective field
	std::cout << "Updating the total field" << std::endl;
	updateEffectiveField(); // Updated 'H' matrix

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

	std::cout << "Running Landau-Lifshitz simulation" << std::endl;

	int frame = 0;
	double torque = 1e2; // Max of |dM/dt|

	// while(frame < 1){ // For testing
	while(torque > maxTorque && frame < maxSteps){

		torque = 0.0; // Reset maximum torque

		// Evolve the magnetization of each node
		for(int i = 0; i<InsideNodes.size(); i++){
			int nd = InsideNodes(i); // Global node label
	
			// New magnetization
			Eigen::Vector3d M2 = step(M.col(nd), Mold.col(nd), H.col(nd), Hold.col(nd));

			Mold.col(nd) = M.col(nd); // Update the old magnetization value
			M.col(nd) = M2; // Update the new magnetization value

			// Store the old magnetic field H
			Hold.col(nd) = H.col(nd);

			// Update the magnetizaiton norm
			Mnorm(nd) = M.col(nd).norm();

			// Check the torque term dM/dt
			Eigen::Vector3d M_nd, H_nd;
			M_nd = M.col(nd); H_nd = H.col(nd);

			Eigen::Vector3d dM_dt = M_nd.cross(H_nd) + alfa* M_nd.cross( M_nd.cross(H_nd) );
			torque = std::max(torque, dM_dt.norm()); // Store the largest torque value of the mesh

			// Store the average Mx,y,z of the mesh
			M_time.col(frame) += M_nd;

		} // Update M

		// Average Mx,y,z of current time frame
		M_time.col(frame) /= InsideNodes.size(); 

		// Get each magnetic field contributions for the new magnetization field
		magnetostaticField();
		exchangeField();
		anisotropyField();

		updateEffectiveField(); // Update the effective field

		if(verbose){
			std::cout << frame << "/" << maxSteps << " |dM/dt| = " << torque << std::endl;	
		} 

		frame++;	  // Update iteration step

	} // Time step loop

	std::cout << "Simulation finished. Saving..." << std::endl;

	std::ofstream file("M_time.txt");
    if (file.is_open()) {
        file << M_time << std::endl;
        file.close();
    } // Save the M(time)

    std::ofstream file2("Mfield.txt");
    if (file2.is_open()) {
        file2 << M << std::endl;
        file2.close();
    } // Save the M

    std::ofstream file3("Hfield.txt");
    if (file3.is_open()) {
        file3 << H << std::endl;
        file3.close();
    } // Save the H

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

	std::cout << "Building the local exchange stiffness matrix" << std::endl;

	// Create the local stiffness matrix
	Eigen::MatrixXd Ak = Eigen::MatrixXd::Zero(16, t.cols());
	for(int ik = 0; ik<InsideElements.size(); ik++){
		int k = InsideElements(ik);

		std::cout << k << " | " << t.cols() << std::endl;


		Eigen::Matrix4d aux = VE[k]*(b.col(k)*b.col(k).transpose() + c.col(k)*c.col(k).transpose() + d.col(k)*d.col(k).transpose());
		Ak.col(k) = Eigen::Map<Eigen::VectorXd>(aux.data(), 16);
	} // Local exchange stiffness matrix

	std::cout << "Building the triplets" << std::endl;

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

	std::cout << "Assembling the sparse exchange stiffness matrix" << std::endl;

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

	    double temp = mu0 * 2.0*Aan/M.col(nd).norm();
	    temp *= M.col(nd).dot(uan); // M dot easy-axis
	    
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

// Update the direction of the magnetization
Eigen::Vector3d LL::step(Eigen::Vector3d M,
					 Eigen::Vector3d Mold,
					 Eigen::Vector3d H,
					 Eigen::Vector3d Hold){

	double d = timeStep/2.0;

	// M(n+1/2)
	Eigen::Vector3d M12 = 3/2*M - 1/2*Mold;

	// M(n+1)
	Eigen::Vector3d M2 = {0.0, 0.0, 0.0};

	// 1) Initial guess of the new magnetic field
	Eigen::Vector3d H12 = 3/2 *H - 0.5 *Hold;

	// Htild from Oriano et al 2008
	Eigen::Vector3d Htild = alfa* M12.cross(H12) + H12;

	// Repeat until M(n+1) doesn't change
	Eigen::Vector3d aux = M; // Copy of M
	double err = 1.0;
	int att = 0;
	while(err > 1e-6 && att < 1000){
	    
	    att += 1;

	    // 2) M (n+1) from M(n) and H(n+1/2)
	    Eigen::Matrix3d mat;
	    mat << 1, d*Htild(2), -d*Htild(1),
	           -d*Htild(2), 1, d*Htild(0),
	           d*Htild(1), -d*Htild(0), 1;
	    
	    Eigen::Vector3d RHS = M - d* M.cross(Htild);

	    M2 = mat.colPivHouseholderQr().solve(RHS); // M(n+1)

	    // 3) M (n + 1/2)
	    M12 = 0.5*(M + M2);

	    // 4) Calculate H~ (n+1/2) from M(n+1/2)
	    if(precession){
	    	Htild = alfa* M12.cross(H12) + H12;	
	    } else{
	    	Htild =	alfa* M12.cross(H12);
	    }
	    
	    // Max difference between new M(n+1) and old M(n+1)
	    err = (M2-aux).norm();
	    // std::cout << "M(n+1) iteration error: " << err << std::endl;

	    // Update M(n+1) from last iteration
	    aux = M2;

	} // Get new magnetization value

	return M2;
}

// Destructor
LL::~LL(){std::cout << "Micromagnetics solver leaving scope \n";}



