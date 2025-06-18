// Wrapper to the C++ code, to be called in Julia

#include "FEMc.cpp"

extern "C"{
	// Function called from Julia
	void abcd_julia(double* pJL, int* nodes, int nd, int nv, double* result){
		// Convert the pointer array to Eigen matrix
		Eigen::Map<Eigen::MatrixXd> p(pJL,3,nv);
		Eigen::Map<Eigen::Vector4i> nds(nodes);
		
		// Run C++ abcd
		Eigen::Vector4d r = abcd(p,nds,nd);

		// Update the output
		result[0] = r[0];
		result[1] = r[1];
		result[2] = r[2];
		result[3] = r[3];
	}

	void stiffnessMatrix(double* AkJL, double* pJL, int* tJL, int nv, int nt, double* VE, double* mu){
		// Convert the pointer array to Eigen matrix
		Eigen::Map<Eigen::MatrixXd> Ak(AkJL,16,nt);
		Eigen::Map<Eigen::MatrixXd> p(pJL,3,nv);
		Eigen::Map<Eigen::MatrixXi> t(tJL,4,nt);

		// Run C++ local stiffness matrix
		localStiffnessMatrix(Ak,p,t,VE,mu);
	}

} // extern C