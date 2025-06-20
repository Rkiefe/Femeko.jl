
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

#include "FEMc.cpp"

// Get midpoints of input triangle
Eigen::Matrix3d midPoints(Eigen::Ref<const Eigen::Matrix3d> r)
{
	Eigen::Matrix3d t(3,3);
	t.col(0) = 0.5*(r.col(0) + r.col(1));
	t.col(1) = 0.5*(r.col(1) + r.col(2));
	t.col(2) = 0.5*(r.col(0) + r.col(2));

	return t;	
} // Triangle midpoints

// Subdivide the element into 4 and get the midpoints of the edges
Eigen::MatrixXd subtriangle(Eigen::Ref<const Eigen::Matrix3d> p)
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
	double* VE)
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
		Eigen::Matrix4d aux = VE[k]*(b*b.transpose() + c*c.transpose() + d*d.transpose());
		
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
	Eigen::Ref<Eigen::MatrixXi> surfaceT,
	double* AE)
{
	int nv = p.cols(); 			// Number of mesh nodes
	int ne = surfaceT.cols(); 	// Number of surface elements
	
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nv,ne);
	
	// Update the matrix
	for(int s = 0; s<ne; s++){
		
		// Area of the surface element 
		// double areaT = areaTriangle(p, surfaceT(0,s), surfaceT(1,s), surfaceT(2,s));
		double areaT = AE[s];

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
	Eigen::Ref<Eigen::MatrixXi> surfaceT,
	Eigen::Ref<Eigen::MatrixXd> normal,
	double* areaT)
{

	// Pre compute the values of the linear basis function over the quadrature
	Eigen::MatrixXd phi(3,15);
	phi.row(0) << 1.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.75, 0.5, 0.75, 0.25, 0.0, 0.25, 0.0, 0.25, 0.25;
	phi.row(1) << 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.25, 0.25, 0.0, 0.75, 0.75, 0.5, 0.25, 0.0, 0.25;
	phi.row(2) << 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.25, 0.25, 0.0, 0.25, 0.25, 0.75, 0.75, 0.5;
	
	int nv = p.cols(); 			// Number of mesh nodes
	int ne = surfaceT.cols(); 	// Number of surface elements

	double one_six = 1.0/6.0; // Avoid repeats
	double one_over_4pi = 1/(4*pi);

	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(ne,nv);
	for(int m = 0; m < ne; m++){

		// Center of surface element m
		Eigen::Vector3d xm(0.0,0.0,0.0);

		// Update C over the nodes of m
		for(int i = 0; i<3; i++){
			int nd = surfaceT(i,m);
			C(m,nd) = one_six;

			// Calculate the centroid of the surface element m
			xm += p.col(nd);
		}
		xm /= 3.0;

		// Now add the boundary integral
		for(int s = 0; s < ne; s++){

			// Get the element node coordinates xyz into a matrix
			Eigen::Matrix3d triangle_coord;
			triangle_coord.col(0) = p.col(surfaceT(0,s));
			triangle_coord.col(1) = p.col(surfaceT(1,s)); 
			triangle_coord.col(2) = p.col(surfaceT(2,s)); 

			// Get the nodes of the quadrature for this element
			Eigen::MatrixXd r = subtriangle(triangle_coord);
			
			// Temporary vector to hold the integral value
			Eigen::Vector3d integral(0.0,0.0,0.0);
			
			// distance between centroid xm and the quadrature point
			Eigen::Vector3d distance(0.0,0.0,0.0);
			
			// Evaluate the integral over the quadrature points
			for(int quad = 0; quad<r.cols(); quad++){ // 15 quadrature nodes
				Eigen::Vector3d y = r.col(quad);
				distance = y - xm;

				integral += (distance.dot(normal.col(s))/pow(distance.norm(),3)) * phi.col(quad);
			}
			// Update C
			C(m,surfaceT(0,s)) += (one_over_4pi*areaT[s]/r.cols()) * integral(0);
			C(m,surfaceT(1,s)) += (one_over_4pi*areaT[s]/r.cols()) * integral(1);
			C(m,surfaceT(2,s)) += (one_over_4pi*areaT[s]/r.cols()) * integral(2);

		} // Loop over surface elements | s
	} // Loop over surface elements | m

	return C;
} // C matrix of BEM 

// D matrix | i give up naming this matrices (its the last one :D)
Eigen::MatrixXd Dmatrix(
	Eigen::Ref<Eigen::MatrixXd> p,
	Eigen::Ref<Eigen::MatrixXi> surfaceT,
	double* areaT)
{
	int ne = surfaceT.cols(); // Number of surface elements
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(ne,ne);

	double one_over_4pi = 1/(4*pi);

	for (int m = 0; m<ne; m++)
	{
		// Centroid of surface element m
		Eigen::Vector3d xm(0.0,0.0,0.0);

		// Go over each node of the surface element to get the average position | centroid
		for(int i = 0; i<3; i++){
			xm += p.col(surfaceT(i,m));
		}
		xm = xm/3.0;

		for (int n = 0; n<ne; n++)
		{

			// Get the element node coordinates xyz into a matrix
			Eigen::Matrix3d triangle_coord;
			triangle_coord.col(0) = p.col(surfaceT(0,n));
			triangle_coord.col(1) = p.col(surfaceT(1,n)); 
			triangle_coord.col(2) = p.col(surfaceT(2,n)); 

			// Get the nodes of the quadrature for this element n
			Eigen::MatrixXd r = subtriangle(triangle_coord);

			// Go over the quadrature points
			double integral = 0.0;
			for(int quad = 0; quad<r.cols(); quad++){
				Eigen::Vector3d y = r.col(quad);
				double distance = (y-xm).norm();

				integral += 1/distance;
			}

			// Update D
			D(m,n) += one_over_4pi * integral * areaT[n]/r.cols();

		} // End of loop over surface elements n
	} // End of loop over surface elements m

	return D;
} // D matrix | BEM

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

