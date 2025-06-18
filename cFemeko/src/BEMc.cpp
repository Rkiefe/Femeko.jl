#include "FEMc.cpp"

/* 	   FEM-BEM approach
	| A B | |u| = -|RHS|
	| C D | |p|	   | 0 |
*/

/* 
	The integral of the BEM matrices is handled with a simple quadrature
	on the midpoints of the element edges. The triangle is subdivided
 	into 4 smaller triangles and the integral is evaluated on each midpoint.

 	It cannot be the center of the original element as it holds a singularity
 	thanks to the Greens function
*/

// Get midpoints of input triangle
Eigen::Matrix3d midPoints(Eigen::Ref<Eigen::Matrix3d> r)
{
	Eigen::Matrix3d t(3,3);
	t.col(0) = 0.5*(r.col(0) + r.col(1));
	t.col(1) = 0.5*(r.col(1) + r.col(2));
	t.col(2) = 0.5*(r.col(0) + r.col(2));

	return t;	
} // Triangle midpoints

// Subdivide the element into 4 and get the midpoints of the edges
Eigen::MatrixXd subtriangle(Eigen::Ref<Eigen::MatrixXd> p)
{
	Eigen::MatrixXd r(3,15);

	// Thre first three coordinates remain
	r.block<3,3>(0,0) = p;

	// Midpoints of large triangle
	r.block<3,3>(0,3) = midPoints(r.block<3,3>(0,0));

	// -- Mid points of each 4 smaller triangles --
	// 1st small triangle
	Eigen::Matrix3d temp(3,3);
	temp.col(0) = r.col(0);
	temp.col(1) = r.col(3);
	temp.col(2) = r.col(5);
	r.block<3,3>(0,6) = midPoints(temp);

	// 2nd small triangle
	temp.col(0) = r.col(3);
	temp.col(1) = r.col(1);
	temp.col(2) = r.col(4);
	r.block<3,3>(0,9) = midPoints(temp);

	// 3nd small triangle
	temp.col(0) = r.col(4);
	temp.col(1) = r.col(2);
	temp.col(2) = r.col(5);
	r.block<3,3>(0,12) = midPoints(temp);
	
	// (the 4th center one is complete by doing the other 3) 

	return r;
} // Get element quadrature nodes 


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
Eigen::MatrixXd Bmatrix(
	Eigen::Ref<Eigen::MatrixXd> p,
	Eigen::Ref<Eigen::MatrixXi> surfaceT)
{
	int nv = p.cols(); 			// Number of mesh nodes
	int ne = surfaceT.cols(); 	// Number of surface elements
	
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nv,ne);
	
	// Update the matrix
	for(int s = 0; s<ne; s++){
		
		// Area of the surface element 
		double areaT = areaTriangle(p, surfaceT(0,s), surfaceT(1,s), surfaceT(2,s));

		// For each node of the surface triangle...
		for(int i = 0; i<3; i++){
			// int nd = surfaceT(i,s);
			B(surfaceT(i,s),s) += areaT/3; 
		}

	} // End of loop over surface elements

	return B;
}

// C | also dont know how to name this function
Eigen::MatrixXd Cmatrix(
	Eigen::Ref<Eigen::MatrixXd> p,
	Eigen::Ref<Eigen::MatrixXi> surfaceT)
{

	// Pre compute the values of the linear basis function over the quadrature
	Eigen::MatrixXd phi(3,15);
	phi.row(0) << 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.75, 0.5, 0.75, 0.25, 0.0, 0.25, 0.0, 0.25, 0.25;
	phi.row(1) << 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.25, 0.25, 0.0, 0.75, 0.75, 0.5, 0.25, 0.0, 0.25;
	phi.row(2) << 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.25, 0.25, 0.0, 0.25, 0.25, 0.75, 0.75, 0.5;
	
	int nv = p.cols();
	int ne = surfaceT.cols();

	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(ne,nv);
	for(int m = 0; m < ne; m++){

		Eigen::Vector3d xm;

		// Update C over the nodes of m
		for(int i = 0; i<3; i++){
			int nd = surfaceT(i,m);
			C(m,nd) = 1.0/6.0;
		}


		for(int s = 0; s < ne; s++){

		} // Loop over surface elements | s

	} // Loop over surface elements | m


	return C;
}



int main(int argc, char const *argv[])
{
	/* code */
	return 0;
}