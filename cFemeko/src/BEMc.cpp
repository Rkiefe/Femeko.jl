#include "FEMc.cpp"

/* 	   FEM-BEM approach
	| A B | |u| = -|RHS|
	| C D | |p|	   | 0 |
*/

// A | Stiffness matrix using dense matrix for BEM
Eigen::MatrixXd denseStiffnessMatrix(
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> t, 
	double* VE,
	double* mu)
{

	// Same as FEM stiffnessMatrix but using a dense version of the matrix
	// This matrix in the cope of FEM-BEM should only be called once 

	int nv = p.cols();
	int nt = t.cols();

	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(nv,nv);

	// Populate the global stiffness matrix directly
	for(int k = 0; k<t.cols(); k++){
		Eigen::Vector4d b,c,d;
		for(int i = 0; i<4; i++){
			Eigen::Vector4d r = abcd(p,t.col(k),t(i,k));
			b(i) = r[1];
			c(i) = r[2];
			d(i) = r[3]; 
		}
		Eigen::Matrix4d aux = VE[k]*mu[k]*(b*b.transpose() + c*c.transpose() + d*d.transpose());
		
		// Update the global stiffness matrix
		for(int i = 0; i<4; i++){
			for(int j = 0; j<4; j++){
				A(t(i,k),t(j,k)) += aux(i,j);
			}
		} // updated the global stiffness matrix
	} // End of loop over the mesh elements 

	return A;
} // Dense stiffness matrix for BEM

// B | dont know how to name this one






int main(int argc, char const *argv[])
{
	/* code */
	return 0;
}