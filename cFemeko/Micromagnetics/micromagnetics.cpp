#include "../../src/FEMc.cpp"

class LL
{
public:

	// Mesh data
	Eigen::MatrixXd p;
	Eigen::MatrixXi t;
	Eigen::VectorXi InsideElements;
	Eigen::VectorXi InsideNodes;
	double* VE;

	// FEM lienar basis function F = a + bx + cy + dz
	Eigen::MatrixXd a, b, c, d;

	Eigen::SparseMatrix<double> AEXC; // Exchange field stiffness matrix

	// Material properties
	Eigen::MatrixXd M; // M field, 3 by nv
	double Aexc = 0.0;
	double Aan = 0.0;
	Eigen::Vector3d uan = {1.0, 0.0, 0.0};  // Easy axis

	Eigen::Vector3d Hext = {0.0, 0.0, 0.0}; // Applied field (Tesla)

	// Solver properties
	int maxSteps = 100000; 	 // Max number of time steps
	double timeStep = 0.01;  // Time step in giro seconds
	double maxTorque = 1e-5; // Max torque |dM/dt|

	double alfa = 1.0;
	bool precession = true;  // Consider precession ? 

	// 		-- Outputs --
	Eigen::MatrixXd M_time = Eigen::MatrixXd::Zero(3, maxSteps); // <M>x,y,z over time

	// Conjugate Gradient solver for the magnetostatic scalar potential
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> CG; 
	
	Eigen::MatrixXd Hd; // Demagnetizing field, 3 by nv
	Eigen::MatrixXd Hexc; // Exchange field, 3 by nv




	// Create the micromagnetics solver
	LL(Eigen::Ref<Eigen::MatrixXd> p_input, 
	   Eigen::Ref<Eigen::MatrixXi> t_input,
	   Eigen::Ref<Eigen::VectorXi> InsideElements_input,
	   Eigen::Ref<Eigen::VectorXi> InsideNodes_input,
	   double* VE_ptr,
	   Eigen::Ref<Eigen::MatrixXd> M_input)
	: p(p_input), 
	  t(t_input), 
	  InsideElements(InsideElements_input), 
	  InsideNodes(InsideNodes_input), 
	  VE(VE_ptr), 
	  M(M_input)
	{
		/*
			0) Runs a 'setup()' that updates the solver properties
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
		exchangeStiffness();

		// Run the micromagnetics solver
		run();

	} // LL constructor

	// Destructor
	~LL(){std::cout << "Micromagnetics solver leaving scope \n";}

	// Run the micromagnetics solver
	void run(){

		// Update the demagnetizing field
		magnetostaticField();
	
	} // Run the micromagnetic solver
	

	// Get the linear basis function over each node and element
	void linearBasis(){
		a = Eigen::MatrixXd::Zero(4, t.cols());
		b = Eigen::MatrixXd::Zero(4, t.cols());
		c = Eigen::MatrixXd::Zero(4, t.cols());
		d = Eigen::MatrixXd::Zero(4, t.cols());
		for(int k = 0; k<t.cols(); k++){

			for(int i = 0; i<4; i++){
				Eigen::Vector4d r = abcd(p, t.col(k), t(i, k));
				a(i, k) = r(0);
				b(i, k) = r(1);
				c(i, k) = r(2);
				d(i, k) = r(3);
			}
		} // Loop over the elements

	} // FEM linear basis function of each node and element

	// Exchange field stiffness matrix
	void exchangeStiffness(){

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
		Eigen::SparseMatrix<double> AEXC(p.cols(), p.cols());
		AEXC.setFromTriplets(triplets.begin(), triplets.end());
		// AEXC.makeCompressed();

	}

	// Get the demagnetizing field
	void magnetostaticField(){

		// Update the boundary conditions for the solver
		std::cout << "Applying the boundary conditions" << std::endl;
		
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
		std::cout << "Solving the scalar potential" << std::endl;
		Eigen::VectorXd u = CG.solve(RHS); // CG.info() output 0 if CG succeeds

		// Calculate the demagnetizing field from the scalar potential
		std::cout << "Updating the demagnetizing field" << std::endl;
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
		std::vector<double> volumes(nv, 0.0); 
		for(int k = 0; k<nt; k++){
			for(int i = 0; i<4; i++){
				int nd = t(i, k);
				volumes[nd] += VE[k];
				Hd(0, nd) += Hdk(0, k)*VE[k];
				Hd(1, nd) += Hdk(1, k)*VE[k];
				Hd(2, nd) += Hdk(2, k)*VE[k];
			}
		}

		for(int nd = 0; nd<nv; nd++){
			Hd(0, nd) /= volumes[nd];
			Hd(1, nd) /= volumes[nd];
			Hd(2, nd) /= volumes[nd];
		} // Average over the nodes

	} // Demagnetizing field from scalar potential

	void exchangeField(){
// 	function exchangeField(self::LL, M::Matrix{Float64}, AEXC, Volumes::Vector{Float64})
// 		# AEXC -> Stiffness matrix of the exchange field
// 		# Aexc -> Exchange constant
// 		mu0::Float64 = pi*4e-7 # Vaccum magnetic permeability
// 		Hexc::Matrix{Float64} = zeros(3, self.mesh.nv)
// 		for i in 1:3
// 		    Hexc[i, :] = -2*mu0*self.Aexc* AEXC*M[i, :]./(0.25*self.Ms^2 *self.scale^2 *Volumes)
// 		end 

// 		return Hexc
// 	end # Exchange field (T)
	}


};

extern "C"{

	// Sends raw Julia pointers to a C++ object
	void LandauLifshitz(double* p_ptr,  // Node coordinates 3 by nv
						int* t_ptr, 	// Node connectivity 3 by nt
						double* VE, 	// Volume of each element
						int* InsideElements_ptr, // Array of element labels in magnetic volume
						int* InsideNodes_ptr, 	 // Array of node labels in magnetic volume
						int nv, int nt,
						int nInside, int nInsideNodes,
						double* M_ptr) // Magnetization field
	{

		// Map the mesh info to Eigen
		Eigen::Map<Eigen::MatrixXd> p(p_ptr, 3, nv);
		Eigen::Map<Eigen::MatrixXi> t(t_ptr, 4, nv);
		Eigen::Map<Eigen::VectorXi> InsideElements(InsideElements_ptr, nInside);
		Eigen::Map<Eigen::VectorXi> InsideNodes(InsideNodes_ptr, nInsideNodes);

		// Map the magnetization vector field
		Eigen::Map<Eigen::MatrixXd> M(M_ptr, 3, nv);


		// Create the micromagnetics solver object
		LL solver(p, t, InsideElements, InsideNodes, VE, M);

	} // Sends raw Julia pointers to a C++ object

} // Julia to C++ wrapper 

// int main()
// {
// 	std::cout << "Hello world! \n";
// }