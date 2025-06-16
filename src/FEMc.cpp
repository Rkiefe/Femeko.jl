/*
	C++ version of Femeko.jl to speed up calculations
*/

#include <iostream>	
#include <vector>
#include <omp.h>
// #include <algorithm> // For std::find

// Eigen for linear algebra
#include "../extern/eigen/Eigen/Dense"

// FEM basis function
Eigen::Vector4d abcd(Eigen::MatrixXd p,Eigen::Vector4i nodes,int nd){
	
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

// Dense, element-wise stiffness matrix
void localStiffnessMatrix(Eigen::Ref<Eigen::MatrixXd> Ak, Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXi> t, double* VE, double* mu){
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
} // Local stiffness matrix


/*
	-- In progress --
*/
// // Boundary integral | Vector field boundary conditions
// Eigen::VectorXd BoundaryIntegral(Eigen::Ref<Eigen::MatrixXd> p, Eigen::Ref<Eigen::MatrixXd> surfaceT, Eigen::Ref<Eigen::MatrixXd> normal, Eigen::Ref<Eigen::Vector3d> F, std::vector<int>& shell_id){

// 	int nv = p.cols(); 			// Number of mesh nodes
// 	int ne = surfaceT.cols(); 	// Number of surface elements

// 	Eigen::VectorXd RHS = Eigen::VectorXd::Zero(nv);

// 	for (int s = 0; s < ne; s++)
// 	{
// 		// Only integrate over the outer shell (shell_id)
// 		if (!( surfaceT(3,s) == shell_id[0] || surfaceT(3,s) == shell_id[1] || surfaceT(3,s) == shell_id[2] )) {
// 		    continue;
// 		}
		
// 		// Area of surface triangle
// 		double areaT = areaTriangle(p, surfaceT(0,s), surfaceT(1,s), surfaceT(2,s));
// 		std::cout << areaT << std::endl;
		
// 		return RHS; 
// 	}

// 	return RHS;
// }

// int main(int argc, char const *argv[])
// {

// 	int nv = 4;
// 	int ne = 4;
// 	Eigen::MatrixXd p = Eigen::MatrixXd::Zero(3,nv);

// 	p.row(0) << 0,1,0,0;
// 	p.row(1) << 0,0,1,0;
// 	p.row(2) << 0,0,0,1;

// 	Eigen::MatrixXd surfaceT = Eigen::MatrixXd::Zero(4,ne);
	
// 	surfaceT.col(0) << 0,1,2,0;
// 	surfaceT.col(1) << 0,1,3,0;
// 	surfaceT.col(2) << 0,2,3,0;
// 	surfaceT.col(3) << 1,2,3,0;

// 	Eigen::MatrixXd normal = Eigen::MatrixXd::Zero(3,ne);
// 	Eigen::Vector3d F(1.0,2.0,3.0);
// 	std::vector<int> shell_id = {0,1,2};

// 	Eigen::VectorXd RHS = BoundaryIntegral(p, surfaceT, normal, F, shell_id);

// 	std::cout << "Working" << std::endl;
// 	return 0;
// }


// Wrapper to the C++ code, to be called in Julia
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

