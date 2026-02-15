/*
		Logic for 2D FEM
	
	Requirements:
		Eigen
		mesh.h
*/

// Linear basis function in global coordinates
Eigen::Vector3d linearBasis2D(Eigen::Ref<Eigen::MatrixXd> p,
							Eigen::Vector3i nodes,
							int nd){

// Phi of node 'nd' = a + b*x + c*y
// where Phi at node 'nd' = 1.0 ; and 0.0 every other node

	// Find the nodes that are node the target node 'nd'
	int nd1 = -1;
	int nd2 = -1;
	for(int i : nodes){
		if(i != nd){
			if(nd1<0){
				nd1 = i;
			}else if(nd2<0){
				nd2 = i;
			}
		}
	} // Find the nodes that are not 'nd'

	// Build the matrix equation for the a b c d
	Eigen::Matrix3d M(3, 3);
	M.row(0) << 1.0, p(0, nd),  p(1, nd);
	M.row(1) << 1.0, p(0, nd1), p(1, nd1);
	M.row(2) << 1.0, p(0, nd2), p(1, nd2);

	Eigen::Vector3d b(1.0, 0.0, 0.0);
	Eigen::Vector3d r = M.colPivHouseholderQr().solve(b);

	return r;
} // Linear basis function in global coordinates

// Build the local stiffness matrix
// Requires the pre-computed linear basis coefficients P1b and P1c (does not need P1a)
Eigen::MatrixXd localStiffnessMatrix2D(MESH &mesh,
									   Eigen::Ref<Eigen::MatrixXd> P1b,
									   Eigen::Ref<Eigen::MatrixXd> P1c){
	
	Eigen::MatrixXd Ak(9, mesh.nt); // Local stiffness matrix
	for(int k = 0; k<mesh.nt; k++){
		Eigen::Matrix3d aux = mesh.VE[k]*(P1b.col(k) * P1b.col(k).transpose() + 
										  P1c.col(k) * P1c.col(k).transpose());

		Ak.col(k) = Eigen::Map<Eigen::VectorXd>(aux.data(), 9);
	}

	return Ak;
} // Local stiffness matrix from pre-computed P1 coef.

// Build the sparse stiffness matrix 
// Requires the pre-computed local stiffness matrix
Eigen::SparseMatrix<double> stiffnessMatrix2D(MESH &mesh,
											  Eigen::Ref<Eigen::MatrixXd> Ak,
											  std::vector<double> &mu){

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(9 * mesh.nt); // At most, the matrix is fully populated

	// Update the global matrix
	for(int k = 0; k<mesh.nt; k++){
	    int n = 0;
	    for (int i = 0; i < 3; i++) {
	        for (int j = 0; j < 3; j++) {
	            
	            // Add contribution to global matrix
	            triplets.emplace_back(mesh.t(i, k), mesh.t(j, k), Ak(n, k)*mu[k]);
	            n++;
	        }
	    }
	} // Loop over the elements

	// Build global stiffness matrix from triplets
	Eigen::SparseMatrix<double> A(mesh.nv, mesh.nv);
	A.setFromTriplets(triplets.begin(), triplets.end());

	return A;
} // sparse stiffness matrix from the local dense matrix
