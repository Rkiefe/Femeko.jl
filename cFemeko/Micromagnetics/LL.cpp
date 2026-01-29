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

	// Get the linear basis function of each element and node
	std::cout << "Calculating the linear basis function over each element and node" << std::endl;
	linearBasis();

	// Stiffness matrix
	std::cout << "Building the sparse stiffness matrix" << std::endl;
	A = stiffnessMatrix(p, t, VE);

	// Process the stiffness matrix with the linear solver
	std::cout << "Applying the linear solver pre-computation" << std::endl;
	CG.compute(A);

	// Compute the exchange stiffness matrix
	std::cout << "Building the exchange stiffness matrix" << std::endl;
	exchangeStiffness();

	// Pre-set the magnetic fields to all zeros
	H = Eigen::MatrixXd::Zero(3, p.cols());
	Hd = Eigen::MatrixXd::Zero(3, p.cols());
	Han = Eigen::MatrixXd::Zero(3, p.cols());
	Hexc = Eigen::MatrixXd::Zero(3, p.cols());

	// Pre-set the output M(t)
	M_time = Eigen::MatrixXd(3, maxSteps);

	// Run the micromagnetics solver
	run();

} // LL constructor

void LL::run(){

	// Update the effective field H
	effectiveField();

	// Store a copy of the input magnetization and field
	Eigen::MatrixXd Mold = M;
	Eigen::MatrixXd Hold = H;
	
	int frame = 0;
	double torque = 1.0;

	while (frame < maxSteps){

		// Update the effective field
		effectiveField();

		// Update the magnetization direction
		for(int i = 0; i<InsideNodes.size(); i++){

			int nd = InsideNodes(i); // Global node label

			// Store the average magnetization
			M_time.col(frame) += M.col(nd)* 1.0/p.cols();

			// Update magnetization
			Eigen::Vector3d M2 = step(M.col(nd), Mold.col(nd), 
									  H.col(nd), Hold.col(nd));

			Mold.col(nd) = M.col(nd);
			M.col(nd) = M2;

			// Check the torque |dM/dt|
			Eigen::Vector3d dM_dt = getTorque(M.col(nd), H.col(nd));
			torque = std::max(torque, dM_dt.norm());
		}

		// Store the old magnetic field
		Hold = H;

		// Show the maximum torque value
		if(verbose)
		{
			std::cout << frame << "/" << maxSteps;
			std::cout << " |dM/dt| = " << torque << std::endl;
		}

		frame++;
	
	} // Time step

	// Save outputs to file
	std::ofstream file("M_time.txt");
	if (file.is_open()) {
	    file << M_time << std::endl;
	    file.close();
	}

}

void LL::linearBasis(){
	int nv = p.cols();
	int nt = t.cols();

	a = Eigen::MatrixXd::Zero(4, nt);
	b = Eigen::MatrixXd::Zero(4, nt);
	c = Eigen::MatrixXd::Zero(4, nt);
	d = Eigen::MatrixXd::Zero(4, nt);

	volumes = Eigen::VectorXd::Zero(nv);

	for(int k = 0; k<nt; k++){
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
}

// Exchange field stiffness matrix
void LL::exchangeStiffness(){

	// Local exchange stiffness matrix
	Eigen::MatrixXd Ak = Eigen::MatrixXd::Zero(16, t.cols());
	for(int ik = 0; ik<InsideElements.size(); ik++){

		int k = InsideElements(ik); // Global element label

		// Local matrix
		Eigen::Matrix4d aux = VE[k]*(b.col(k) * b.col(k).transpose() + c.col(k) * c.col(k).transpose() + d.col(k) * d.col(k).transpose());
		Ak.col(k) = Eigen::Map<Eigen::VectorXd>(aux.data(), 16);
	
	} // Local exchange stiffness matrix
	
	// Temporary storage for triplets (row, col, value)
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(16 * t.cols());

	// Update the global matrix
	for(int k = 0; k<t.cols(); k++){
	    int n = -1;
	    for (int i = 0; i < 4; ++i) {
	        for (int j = 0; j < 4; ++j) {
	            n++;
	            
	            // Add contribution to global matrix
	            triplets.emplace_back(t(i,k), t(j,k), Ak(n, k));
	        }
	    }
	} // Loop over the elements

	// Build global stiffness matrix from triplets
	Eigen::SparseMatrix<double> AEXC(p.cols(), p.cols());
	AEXC.setFromTriplets(triplets.begin(), triplets.end());

} // Exchange stiffness matrix


// Update the direction of the magnetization
Eigen::Vector3d LL::step(Eigen::Vector3d M,
					 Eigen::Vector3d Mold,
					 Eigen::Vector3d H,
					 Eigen::Vector3d Hold){

	double d = timeStep/2.0;

	// M(n+1/2)
	Eigen::Vector3d M12 = 1.5*M - 0.5*Mold;

	// M(n+1)
	Eigen::Vector3d M2 = {0.0, 0.0, 0.0};

	// 1) Initial guess of the new magnetic field
	Eigen::Vector3d H12 = 1.5*H - 0.5*Hold;

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
	    mat << 1.0, d*Htild(2), -d*Htild(1),
	           -d*Htild(2), 1.0, d*Htild(0),
	           d*Htild(1), -d*Htild(0), 1.0;
	    
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

void LL::effectiveField(){

	// Reset H to zeros
	H = Eigen::MatrixXd::Zero(3, p.cols());
	
	// Compute each component of H
	magnetostaticField();	// Compute the demag field
	exchangeField(); 		// Compute the exchange field
	anisotropyField(); 		// Compute the anisotropy field

	// Sum all contributions on each node
	for(int i = 0; i<p.cols(); i++){
		H.col(i) += Hext;
		H.col(i) += Hd.col(i);
		H.col(i) += Han.col(i);
		H.col(i) += Hexc.col(i);
	}
}

void LL::magnetostaticField(){
// Calculates the scalar potential u
// then the magnetic field, constant on each element
// and then maps the field to the nodes by the box method 

	// Load vector
	Eigen::VectorXd RHS = Eigen::VectorXd::Zero(p.cols());
	for(int ik = 0; ik<InsideElements.size(); ik++){

		int k = InsideElements(ik); // Global element label
		
		// Average magnetization over the element nodes
		Eigen::Vector3d m_avg = {0.0, 0.0, 0.0};
		for(int i = 0; i<4; i++){
			int nd = t(i, k); // Global node label
			m_avg(0) += M(0, nd)/4.0; // <Mx>
			m_avg(1) += M(1, nd)/4.0; // <My>
			m_avg(2) += M(2, nd)/4.0; // <Mz>
		} // <M>_k

		// Update the load vector
		for(int i = 0; i<4; i++){
			int nd = t(i, k); // Global node label
			RHS(nd) += VE[k]*( b(i, k)*m_avg(0) + c(i, k)*m_avg(1) + d(i, k)*m_avg(2) );
		}

	} // Loop over the elements in the magnetic domain  
	
	// Solve the linear system
	Eigen::VectorXd u = CG.solve(RHS); // Scalar potential

	// Demagnetizing field constant on each element
	Eigen::MatrixXd Hdk = Eigen::MatrixXd::Zero(3, t.cols());
	for(int k = 0; k<t.cols(); k++){
		for(int i = 0; i<4; i++){
			Hdk(0, k) -+ b(i, k)*u(t(i, k));
			Hdk(1, k) -+ c(i, k)*u(t(i, k));
			Hdk(2, k) -+ d(i, k)*u(t(i, k));
		}
	} // Hd constant on each element

	// Map the magnetic field to the nodes
	Hd = Eigen::MatrixXd::Zero(3, p.cols());
	for(int k = 0; k<t.cols(); k++){
		for(int i = 0; i<4; i++){
			int nd = t(i, k);
			Hd(0, nd) += VE[k]*Hdk(0, k);
			Hd(1, nd) += VE[k]*Hdk(1, k);
			Hd(2, nd) += VE[k]*Hdk(2, k);
		}
	} // Multiply Hd_k by the volume of k
	for(int nd = 0; nd < p.cols(); nd++){
		Hd(0, nd) /= volumes(nd);
		Hd(1, nd) /= volumes(nd);
		Hd(2, nd) /= volumes(nd);
	} // Divide by the total volume of 'k' containing 'nd' (box method)

} // Demagnetizing field on the nodes of the mesh

void LL::exchangeField(){

	// Hexc = Eigen::MatrixXd::Zero(3, p.cols()); // Don't need to reset because it is over written
	double coef = -2.0*mu0/(0.25 * Ms*Ms * scale*scale);
	for(int i = 0; i<3; i++){ // x,y,z 
		
		Eigen::VectorXd temp =  coef * Aexc*M.row(i);
		
		for(int nd = 0; nd<p.cols(); nd++){
			Hexc(i, nd) = temp(nd)/volumes(nd); // Overwrite the exchange field
		} // Loop over each node

	} // Loop over x,y,z 

} // Exchange field

void LL::anisotropyField(){
	
	Han = Eigen::MatrixXd::Zero(3, p.cols()); // Don't need to reset because it is overwritten
	double coef = mu0*2.0*Aan;

	for(int i = 0; i<InsideNodes.size(); i++){
		int nd = InsideNodes(i);
		double Msqrd = M.col(nd).dot(M.col(nd));
		Han.col(nd) = coef/Msqrd * (M.col(nd).dot(uan)) * uan;
	}

} // Anisotropy field

Eigen::Vector3d LL::getTorque(Eigen::Vector3d m, 
							  Eigen::Vector3d h){
	
	// Check the torque term dM/dt
	Eigen::Vector3d dM_dt = m.cross(h) + alfa* m.cross( m.cross(h) );
	return dM_dt;
}

LL::~LL(){
	std::cout << "Micromagnetics solver leaving scope" << std::endl;
}

int main(){
	std::cout << "Hello world!" << std::endl;
}