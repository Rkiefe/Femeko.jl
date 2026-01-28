#pragma once
#include "../../src/FEMc.cpp"

class LL
{
public:

	// Mesh data
	Eigen::MatrixXd p; // Node coordinates, 3 by nv
	Eigen::MatrixXi t; // Element connectivity 3 by nt
	Eigen::VectorXi InsideElements; // Array of element labels in magnetic volume
	Eigen::VectorXi InsideNodes;    // Array of node labels in magnetic volume
	double* VE; 					// Array of element volumes (size nt)
	Eigen::VectorXd volumes; // Volume of elements surrounding each node (size = nv)

	Eigen::MatrixXd a, b, c, d; 		// FEM linear basis function F = a + bx + cy + dz, 4 by nt
	Eigen::SparseMatrix<double> AEXC; 	// Exchange field stiffness matrix
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> CG; // Conjugate Gradient solver for the magnetostatic scalar potential

	// Material properties
	double Ms = 1.0;		// Saturation magnetization (Tesla)
	double Aexc = 0.0; 		// Exchange energy (J/m)
	double scale = 1e-9; 	// Scale of the material
	double Aan = 0.0; 		// Anisotropy energy (J/m3)
	Eigen::Vector3d uan = {1.0, 0.0, 0.0};  // Easy axis

 	// Solver properties
	int maxSteps = 100000; 	 // Max number of time steps
	double timeStep = 0.01;  // Time step in giro seconds
	double maxTorque = 1e-5; // Max torque |dM/dt|
	double alfa = 1.0;
	bool precession = true;  // Consider precession ? 

	// 		-- Outputs --
	Eigen::MatrixXd M_time = Eigen::MatrixXd::Zero(3, maxSteps); // <M>x,y,z over time

	Eigen::Vector3d Hext = {0.0, 0.0, 0.0}; // Applied field (Tesla)
	Eigen::MatrixXd Hd; 	// Demagnetizing field, 3 by nv
	Eigen::MatrixXd Hexc; 	// Exchange field, 3 by nv
	Eigen::MatrixXd Han; 	// Anisotropy field, 3 by nv
	Eigen::MatrixXd M; 		// M field, 3 by nv
	Eigen::MatrixXd H; 		// Effective field, 3 by nv (sum of all H fields)

	// Create the micromagnetics solver
	LL(Eigen::Ref<Eigen::MatrixXd> p_input, 
	   Eigen::Ref<Eigen::MatrixXi> t_input,
	   Eigen::Ref<Eigen::VectorXi> InsideElements_input,
	   Eigen::Ref<Eigen::VectorXi> InsideNodes_input,
	   double* VE_ptr,
	   Eigen::Ref<Eigen::MatrixXd> M_input);

	// Destructor
	~LL();

	// Run the micromagnetics solver
	void run();
	
	// Get the linear basis function over each node and element
	void linearBasis();

	// Get the Exchange field stiffness matrix
	void exchangeStiffness();

	// Get the demagnetizing field on each node
	void magnetostaticField();

	// Get the Exchange Field (3 by nv) 
	void exchangeField();

	// Get the Anisotropy field (3 by nv)
	void anisotropyField();

	// Creates a zero(3, nv) array and adds each magnetic field
	void updateEffectiveField();

	// Update the magnetization direction
	Eigen::Vector3d step(Eigen::Vector3d M,
						 Eigen::Vector3d Mold,
						 Eigen::Vector3d H,
						 Eigen::Vector3d Hold);
};
