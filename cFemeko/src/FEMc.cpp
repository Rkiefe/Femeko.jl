/*
	C++ version of Femeko.jl to speed up calculations
	and avoid the garbage collector
*/

#include <iostream>	
#include <vector>
#include <omp.h>
// #include <algorithm> // For std::find

// Eigen for linear algebra
#include "../../extern/eigen/Eigen/Dense"
#include "../../extern/eigen/Eigen/Sparse"

// FEM basis function
Eigen::Vector4d abcd(Eigen::Ref<Eigen::MatrixXd> p,Eigen::Vector4i nodes,int nd){
	
	int nds[3] = {0,0,0};		// All other nodes of the element
	{ // Get nodes different than nd
	int n = 0;
	for(int i = 0; i<4; i++){
		if(nodes(i)!=nd){
			nds[n] = nodes(i);
			n++;
		}
	}
	} // Get nodes different than nd

	// Build the matrix equation for the a b c d
	Eigen::Matrix4d M = Eigen::Matrix4d::Zero();
	M.row(0) << 1,p(0,nd),p(1,nd),p(2,nd);
	M.row(1) << 1,p(0,nds[0]),p(1,nds[0]),p(2,nds[0]);
	M.row(2) << 1,p(0,nds[1]),p(1,nds[1]),p(2,nds[1]);
	M.row(3) << 1,p(0,nds[2]),p(1,nds[2]),p(2,nds[2]);

	Eigen::Vector4d b(1.0,0,0,0);
	Eigen::Vector4d r = M.colPivHouseholderQr().solve(b);

	return r;
} // FEM basis function

// Area of surface triangle
double areaTriangle(Eigen::Ref<Eigen::MatrixXd> p, int i0, int i1, int i2){
	// Area by cross product 
	Eigen::Vector3d AB = p.col(i1) - p.col(i0);
	Eigen::Vector3d AC = p.col(i2) - p.col(i0);

	Eigen::Vector3d cross = AB.cross(AC);
	double area = 0.5 * cross.norm();

    return area;
} // Area of the 3D triangle

// Lagrange multiplier technique | Volume integral of basis function
Eigen::VectorXd lagrange(
	Eigen::Ref<Eigen::MatrixXi> t, 
	double* VE, 
	int nv, int nt)
{
	Eigen::VectorXd C = Eigen::VectorXd::Zero(nv);
	for(int k = 0; k<nt; k++){
		Eigen::Vector4i nds = t.col(k); // Nodes of element k

		// Update C
		for(int i = 0; i<4; i++){
			C(nds(i)) += VE[k]/4.0;
		}
	} // End of loop over the mesh elements

	return C;
} // Lagrange multiplier integral of linear basis function

// Dense, element-wise stiffness matrix
Eigen::MatrixXd localStiffnessMatrix(
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> t, 
	double* VE, double* mu)
{
	
	// Create the local stiffness matrix
	Eigen::MatrixXd Ak = Eigen::MatrixXd::Zero(16,t.cols());

	#pragma omp parallel for
	for(int k = 0; k<t.cols(); k++){
		Eigen::Vector4d b,c,d;
		for(int i = 0; i<4; i++){
			Eigen::Vector4d r = abcd(p,t.col(k),t(i,k));
			b(i) = r[1];
			c(i) = r[2];
			d(i) = r[3]; 
		}
		Eigen::Matrix4d aux = VE[k]*mu[k]*(b*b.transpose() + c*c.transpose() + d*d.transpose());
		Ak.col(k) = Eigen::Map<Eigen::VectorXd>(aux.data(), 16);
	}
	return Ak;

} // Local stiffness matrix

// Gloval stiffness matrix
Eigen::SparseMatrix<double> stiffnessMatrix(
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> t, 
	double* VE, double* mu)
{
	// First calculate the local stiffness matrix
	Eigen::MatrixXd Ak = localStiffnessMatrix(p, t, VE, mu);
	
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
	Eigen::SparseMatrix<double> A(p.cols(), p.cols());
    A.setFromTriplets(triplets.begin(), triplets.end());
    // A.makeCompressed();
    
    return A;
} // Sparse, global stiffness matrix

// Boundary integral | Vector field boundary conditions
Eigen::VectorXd BoundaryIntegral(
	Eigen::Ref<Eigen::MatrixXd> p, 
	Eigen::Ref<Eigen::MatrixXi> surfaceT, 
	Eigen::Ref<Eigen::MatrixXd> normal, 
	Eigen::Ref<Eigen::Vector3d> F, 
	std::vector<int>& shell_id)
{

	int nv = p.cols(); 			// Number of mesh nodes
	int ne = surfaceT.cols(); 	// Number of surface elements

	Eigen::VectorXd RHS = Eigen::VectorXd::Zero(nv);

	for (int s = 0; s < ne; s++)
	{
		// Only integrate over the outer shell (shell_id)
		if (!( surfaceT(3,s) == shell_id[0] || surfaceT(3,s) == shell_id[1] || surfaceT(3,s) == shell_id[2] )) {
		    continue;
		}
		
		// Area of surface triangle
		double areaT = areaTriangle(p, surfaceT(0,s), surfaceT(1,s), surfaceT(2,s));
		
		for(int i = 0; i<3; i++){
			RHS(surfaceT(i,s)) += (normal(0,s)*F(0) + normal(1,s)*F(1) + normal(2,s)*F(2))*areaT/3; 
		} // Update RHS

	} // End of loop of surface elements

	return RHS;
} // Boundary integral (vector field)



// For testing
int main(int argc, char const *argv[])
{

	// Mesh data
	int nt = 2; // Number of elements
	int nv = 5; // Number of nodes
	int ne = 6; // Number of surface elements

	// Node coordinates
	Eigen::MatrixXd p = Eigen::MatrixXd::Zero(3,nv);

	p.row(0) << 0,1,0,0,0;
	p.row(1) << 0,0,1,0,-1;
	p.row(2) << 0,0,0,1,0;
	
	// Element connectivity
	Eigen::MatrixXi t = Eigen::MatrixXi::Zero(4,nt);
	t.col(0) << 0, 1, 2, 3;
	t.col(1) << 0, 1, 3, 4;

	// Volume of each element
	std::vector<double> VE = {1.0/6.0, 1.0/6.0};

	// Surface element node connectivity
	Eigen::MatrixXi surfaceT = Eigen::MatrixXi::Zero(4,ne);
	
	surfaceT.col(0) << 0,1,2,0;
	surfaceT.col(1) << 0,2,3,0;
	surfaceT.col(2) << 1,2,3,0;
	surfaceT.col(3) << 4,1,3,0;
	surfaceT.col(4) << 4,1,0,0;
	surfaceT.col(5) << 4,0,3,0;

	// Surface normals
	Eigen::MatrixXd normal = Eigen::MatrixXd::Zero(3,ne);
	Eigen::Vector3d F(1.0,2.0,3.0);
	std::vector<int> shell_id = {0,1,2};

	// Define constants
	std::vector<double> mu = {1.5, 1.5};

	// Run the boundary integral
	Eigen::VectorXd RHS = BoundaryIntegral(p, surfaceT, normal, F, shell_id);
	
	// Run the lagrange multiplier technique
	Eigen::VectorXd Lag = lagrange(t, VE.data(), nv, nt);

	// Global sparse stiffness matrix
	Eigen::SparseMatrix<double> A = stiffnessMatrix(p, t, VE.data(), mu.data());

	// Create the extended matrix
	Eigen::SparseMatrix<double> mat(A.rows() + 1, A.cols() + 1);

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(A.nonZeros() + 2 * Lag.nonZeros());

	// Add A's elements
	for (int k = 0; k < A.outerSize(); ++k) {
	    for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
	        triplets.emplace_back(it.row(), it.col(), it.value());
	    }
	}

	// Add Lag's elements (right column)
	for (int i = 0; i < Lag.size(); ++i) {
	    if (Lag[i] != 0) {
	        triplets.emplace_back(i, A.cols(), Lag[i]);
	    }
	}

	// Add Lag's elements transposed (bottom row)
	for (int i = 0; i < Lag.size(); ++i) {
	    if (Lag[i] != 0) {
	        triplets.emplace_back(A.rows(), i, Lag[i]);
	    }
	}

	mat.setFromTriplets(triplets.begin(), triplets.end());
	mat.makeCompressed();

	std::cout << "Working" << std::endl;
	return 0;
}




