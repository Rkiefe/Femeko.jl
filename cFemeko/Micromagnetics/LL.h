#include "../../src/FEMc.cpp"
#include <fstream> // Save matrices to .txt files

class LL
{
public:

	// Mesh properties
	Eigen::MatrixXd p; // Node coordinates, 3 by nv
	Eigen::MatrixXi t; // Element connectivity 3 by nt
	Eigen::VectorXi InsideElements; // Array of element labels in magnetic volume
	Eigen::VectorXi InsideNodes;    // Array of node labels in magnetic volume
	double* VE; 					// Array of element volumes (size nt)
	Eigen::VectorXd volumes; // Volume of elements surrounding each node (size = nv)
	Eigen::MatrixXd a, b, c, d; // FEM basis coefficients

	// Linear solver for the Magnetic scalar potential
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> CG; // Conjugate Gradient

	// Stiffness matrix of magnetostatic field
	Eigen::SparseMatrix<double> A;

	// Material properties
	double mu0 = pi*4e-7;	// Vacuum magnetic permeability
	double Ms = mu0*860e3;	// Saturation magnetization (Tesla)
	double Aexc = 12e-13; 		// Exchange energy (J/m)
	double scale = 1e-9; 	// Scale of the material
	double Aan = 0.0; 		// Anisotropy energy (J/m3)
	Eigen::Vector3d uan = {1.0, 0.0, 0.0};  // Easy axis
	
	Eigen::Vector3d Hext = {0.0, mu0*50e3, 0.0};  // External field (tesla)

 	// Solver properties
	int maxSteps = 7035; 	 // Max number of time steps
	double timeStep = 0.01;  // Time step in giro seconds
	double maxTorque = 1e-5; // Max torque |dM/dt|
	double alfa = 0.1/Ms;	 // Damping coef.
	bool precession = true;  // Consider precession ? 
	bool verbose = true; 	 // Show torque values of each simulation step

	// Output
	Eigen::Map<Eigen::MatrixXd> M;
	Eigen::MatrixXd M_time;

	Eigen::MatrixXd H;
	Eigen::MatrixXd Hd;
	Eigen::MatrixXd Han;
	Eigen::MatrixXd Hexc;

	// Create the micromagnetics solver
	LL(Eigen::Ref<Eigen::MatrixXd> p_input, 
	   Eigen::Ref<Eigen::MatrixXi> t_input,
	   Eigen::Ref<Eigen::VectorXi> InsideElements_input,
	   Eigen::Ref<Eigen::VectorXi> InsideNodes_input,
	   double* VE_ptr,
	   Eigen::Map<Eigen::MatrixXd> M_input);

	// Compute the FEM basis coefficients
	void linearBasis();

	// Compute the total magnetic field H
	void effectiveField();

	// Compute the demag field
	void magnetostaticField();

	// Compute the exchange field
	void exchangeField();

	// Compute the anisotropy field
	void anisotropyField();

	// Update the magnetization direction
	Eigen::Vector3d step(Eigen::Vector3d M,
					 		 Eigen::Vector3d Mold,
					 		 Eigen::Vector3d H,
					 		 Eigen::Vector3d Hold);

	Eigen::Vector3d getTorque(Eigen::Vector3d m, 
							  Eigen::Vector3d h);

	// Run the micromagnetic simulation
	void run();

	// Destructor
	~LL();
	
};