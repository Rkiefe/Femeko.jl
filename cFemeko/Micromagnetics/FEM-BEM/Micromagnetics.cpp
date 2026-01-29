/*
	This is a micromagnetics solver in the form using OOP
	This way, the solver calculates the FEM-BEM matrices once
	and it is easier to call them without having to pass pointers
	around. The variable name immediately corresponds to the matrix

*/

#include <iostream>
#include <memory>

// Include Finite Element package
#include "../../../src/BEMc.cpp"


class Micromagnetics{
private:
	Eigen::MatrixXd A; // Stiffness matrix
	Eigen::PartialPivLU<Eigen::MatrixXd> lu; // FEM_BEM factorized matrix

	std::vector<double> Vn;
	std::vector<double> Lagrange;

public:
	// double pi = 3.14159265358979311600;
	double mu0 = pi*4e-7; 				// Vacuum magnetic permeability
	double giro = 2.210173e5 /mu0; 		// Gyromagnetic ratio (rad T-1 s-1)
	double precession = 1.0; 			// Include the precession term ? 1.0 : 0.0
	double damp = 1.0; 					// Damping parameter (dimensionless [0,1])

	double scl = 1e-9; 	// Scale of model (nm)

	// Simulation settings
	double Ms, Aexc, Aan;
	
	Eigen::Vector3d uan; 	// Easy axis direction
	Eigen::Vector3d Hext;  // Applied field

	double totalTime = 1.0/0.0; // inf | Not the best way to define it, but for now it works
	double maxTorque = 1e-14;   // <|dm/dt|>
	int maxAtt = 5000;

	double dt = 0.0117;			// giro seconds | Time step

	// Simulation output
	Eigen::MatrixXd M_avg;
	Eigen::MatrixXd Heff;
	int att = 0; 			// Store the last iteration of the solver (+1)


	// Constructor
	Micromagnetics(
		Eigen::Ref<Eigen::MatrixXd> p,
		Eigen::Ref<Eigen::MatrixXi> t,
		Eigen::Ref<Eigen::MatrixXi> surfaceT,
		Eigen::Ref<Eigen::MatrixXd> normal,
		double* AE,
		double* VE)
		
	{
		uan << 1.0, 0.0, 0.0;
		Hext << 0.0, 0.0, 0.0;

		// Calculate the FEM-BEM matrices
		FEM_BEM_Matrices(p, t, surfaceT, normal, AE, VE);

		// Now for the integral of the basis function and the
		// total volume of each node
		this->Vn.assign(p.cols(),0.0);
		this->Lagrange.assign(p.cols(),0.0);
		for(int k = 0; k<t.cols(); k++)
		{
			for(int i = 0; i<4; i++)
			{
				Vn[t(i,k)] 		 += VE[k];
				Lagrange[t(i,k)] += VE[k]/4;
			}
		}

	} // Constructor


	// ----- Methods -----
	void FEM_BEM_Matrices(
		Eigen::Ref<Eigen::MatrixXd> p,
		Eigen::Ref<Eigen::MatrixXi> t,
		Eigen::Ref<Eigen::MatrixXi> surfaceT,
		Eigen::Ref<Eigen::MatrixXd> normal,
		double* AE,
		double* VE)
	{
		int nv = p.cols(); 			// Number of mesh nodes
		int nt = t.cols(); 			// Number of volume elements
		int ne = surfaceT.cols(); 	// Number of surface elements

		// Make the FEM-BEM matrices
		this->A = denseStiffnessMatrix(p,t,VE); 		// nv by nv
		Eigen::MatrixXd B = Bmatrix(p,surfaceT,AE); 			// nv by ne
		Eigen::MatrixXd C = Cmatrix(p,surfaceT,normal,AE); 	// ne by nv
		Eigen::MatrixXd D = Dmatrix(p,surfaceT,AE); 			// ne by ne

		// Extend the matrix
		// LHS = [-A B; C D]
		Eigen::MatrixXd LHS(nv+ne,nv+ne);
		LHS.block(0,0,nv,nv) 		= -A;
		LHS.block(0, nv, nv, ne) 	= B;
		LHS.block(nv, 0, ne, nv) 	= C;
		LHS.block(nv, nv, ne, ne) 	= D;

		// Factorize LHS once
		this->lu = LHS.partialPivLu(); // Compute and store the factorized LHS
		// Maybe replace partialPivLu() with colPivHouseholderQr() , slower but more accurate

	} // FEM BEM matrices and node volume / integral


	// Landau Lifshitz time step
	void timeStep(
		Eigen::Ref<Eigen::Vector3d> m,
		Eigen::Ref<Eigen::Vector3d> H,
		Eigen::Ref<Eigen::Vector3d> Hold,
		Eigen::Ref<Eigen::Vector3d> Hef,
		double dt,
		double giro,
		double damp,
		double precession)
	{
		double d = dt*giro/2.0;
		
		// m (n+1)
		Eigen::Vector3d mNew, m12, H12;

		// 1) Initial guess of the new magnetic field
		H12 = 1.5 *H - 0.5 *Hold;

		// Find new m(n+1) until it converges
		Eigen::Vector3d oldGuess = m;
		Eigen::PartialPivLU<Eigen::Matrix3d> luSolver;

		double err = 1.0;
		int trial = 0;
		Eigen::Matrix3d mat;
		while(err > 1e-6 && trial < 1000)
		{
			trial++;

			// 2) m (n+1) from m(n) and H(n+1/2)
			mat << 1.0,      d*H12(2),  -d*H12(1),
		          	-d*H12(2),  1.0,     d*H12(0),
		           	d*H12(1), -d*H12(0),  1.0;

		    luSolver.compute(mat);

		    mNew = luSolver.solve(m - d * m.cross(H12));

		    // mNew = mat.lu().solve(m - d * m.cross(H12));  // Best for 3x3 systems

		    // 3) m (n+1/2)
		    m12 = 0.5*(m + mNew);

		    // 4) Get H(n+1/2) from m (n+1/2)
		    H12 = precession*Hef + damp*m12.cross(Hef);

		    // Max difference between mNew and oldGuess
		    for(int i = 0; i<3; i++)
		    {
		    	err = std::max(err,std::abs(mNew(i) - oldGuess(i)));
		    }

		    // Update oldGuess
		    oldGuess = mNew;
		}

		// Update input m to the new m
		m = mNew;
	} // Landau Lifshitz time step

	// Steepest descent | new magnetization 
	void nextM(
		Eigen::Ref<Eigen::Vector3d> m,
		Eigen::Ref<Eigen::Vector3d> Heff,
		double dt)
	{
		double d = dt/2.0;
		
		Eigen::Vector3d H12 = m.cross(Heff);

		// m (n+1)
		Eigen::Matrix3d mat;
		mat << 1.0,      d*H12(2),  -d*H12(1),
	          	-d*H12(2),  1.0,     d*H12(0),
	           	d*H12(1), -d*H12(0),  1.0;

		Eigen::PartialPivLU<Eigen::Matrix3d> luSolver;
	    luSolver.compute(mat);

	    // Check for failure
	    if (luSolver.info() != Eigen::Success) {
	        throw std::runtime_error("nextM | compute(mat) failed");
	    }

	    // Update input m
	    Eigen::Vector3d mNew = luSolver.solve(m - d * m.cross(H12));

	    // Check for failure
	    if (luSolver.info() != Eigen::Success) {
	        throw std::runtime_error("nextM | luSolver failed");
	    }

	    // Update m only if solver worked
	    m = mNew;
	} // Steepest descent time step

	// Micromagnetics with fixed time step
	void LandauLifshitz(
		Eigen::Ref<Eigen::MatrixXd> m,
		Eigen::Ref<Eigen::MatrixXd> p,
		Eigen::Ref<Eigen::MatrixXi> t,
		Eigen::Ref<Eigen::MatrixXi> surfaceT,
		Eigen::Ref<Eigen::MatrixXd> normal,
		double* AE,
		double* VE)
	{
		// Process the mesh input data
		int nv = p.cols();
		int nt = t.cols();
		int ne = surfaceT.cols();
		
		// Prepare output
		this->M_avg = Eigen::MatrixXd::Zero(3,maxAtt);

		// Calculate initial magnetic field
		Heff = Eigen::MatrixXd::Zero(3,nv);
		for(int i = 0; i<nv; i++)
		{
			Heff.col(i) += mu0 *Hext;
		}

		// Demagnetizing field
		Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(3,nv);
		BEMdmag(Hd, p, t, surfaceT, normal, AE, VE, Vn.data(), lu, m);
		Hd *= mu0*Ms;

		// Exchange field
		Eigen::MatrixXd Hexc(3,nv);
		for(int i = 0; i<3; i++)
		{
			// Update the entire xyz row
			Hexc.row(i) = -2*Aexc * (A*(m.row(i).transpose())).transpose();
			
			// Update each node to be at the correct units
			for(int nd = 0; nd<nv; nd++)
			{
				Hexc(i,nd) /= Ms*pow(scl,2)*Lagrange[nd];
			} 
		}

		// Anisotropy field
		Eigen::MatrixXd Han(3,nv);
		for(int i = 0; i<nv; i++)
		{
			Han.col(i) = 2*Aan/Ms *uan.dot(m.col(i)) * uan;
		}

		// Effective field
		Heff += Hd + Hexc + Han;

		// H = Heff + damp M cross Heff
		Eigen::MatrixXd H(3,nv);
		for(int i = 0; i<nv; i++)
		{
			// Eigen::Vector3d aux = m.col(i).head<3>().cross(Heff.col(i).head<3>());
			H.col(i) = precession * Heff.col(i) + damp * m.col(i).head<3>().cross(Heff.col(i).head<3>());
		}

		// Time iteration
		Eigen::MatrixXd Hold = H;

		double time = 0.0;
		this->att = 0;
		while (time < totalTime && att < this->maxAtt)
		{
			// New magnetization
			#pragma omp parallel for
			for(int i = 0; i<nv; i++)
			{
				timeStep(m.col(i).head<3>(),
						 H.col(i),
						 Hold.col(i),
						 Heff.col(i),
						 dt,
						 1.0,
						 damp,
						 precession); // Considering time was normalized by giro
			

			} // New magnetization

			// -- New magnetic field --
			Hold = H; // Store the old one

			// Reset effective field
			for(int i = 0; i<nv; i++)
			{
				Heff.col(i) = mu0 *Hext;
			}

			// Demagnetizing field
			Hd.setZero();
			BEMdmag(Hd, p, t, surfaceT, normal, AE, VE, Vn.data(), lu, m);
			Hd *= mu0*Ms;

			// Exchange field
			for(int i = 0; i<3; i++)
			{
				// Update the entire xyz row
				Hexc.row(i) = -2*Aexc * (A*(m.row(i).transpose())).transpose();
				
				// Update each node to be at the correct units
				for(int nd = 0; nd<nv; nd++)
				{
					Hexc(i,nd) /= Ms*pow(scl,2)*Lagrange[nd];
				}
			}

			// Anisotropy field
			for(int i = 0; i<nv; i++)
			{
				Han.col(i) = 2*Aan/Ms *uan.dot(m.col(i)) * uan;
			}

			// Add all field contributions
			Heff += Hd + Hexc + Han;

			// New H = Heff + damp M cross Heff
			for(int i = 0; i<nv; i++)
			{
				// Eigen::Vector3d aux = m.col(i).head<3>().cross(Heff.col(i).head<3>());
				H.col(i) = precession * Heff.col(i) + damp * m.col(i).head<3>().cross(Heff.col(i).head<3>());
			}

			// Average magnetization
			this->M_avg.col(att) = m.rowwise().mean();

			// <|dm/dt|> , "Torque"
			double dtau = 0.0;
			for(int i = 0; i<nv; i++){
				dtau += (m.col(i).head<3>().cross(Heff.col(i).head<3>())).dot(m.col(i).head<3>().cross(Heff.col(i).head<3>()));
			}


			att++;
			time += dt;

			// Log results
			if(!(att%100))
			{
				std::cout << att << std::endl;
			}

			// return att;

		} // Time iteration

	} // Micromagnetics with fixed time step

	// Micromagnetics energy minimization by steepest descent
	void SteepestDescent(
		Eigen::Ref<Eigen::MatrixXd> m,
		Eigen::Ref<Eigen::MatrixXd> p,
		Eigen::Ref<Eigen::MatrixXi> t,
		Eigen::Ref<Eigen::MatrixXi> surfaceT,
		Eigen::Ref<Eigen::MatrixXd> normal,
		double* AE,
		double* VE)
	{

		// Process the mesh input data
		int nv = p.cols();
		int nt = t.cols();
		int ne = surfaceT.cols();

		// Prepare output
		this->M_avg = Eigen::MatrixXd::Zero(3,maxAtt);

		// Calculate initial magnetic field
		this->Heff = Eigen::MatrixXd::Zero(3,nv);
		for(int i = 0; i<nv; i++)
		{
			Heff.col(i) += mu0 *Hext;
		}

		// Demagnetizing field
		Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(3,nv);
		BEMdmag(Hd, p, t, surfaceT, normal, AE, VE, Vn.data(), lu, m);
		Hd *= mu0*Ms;

		// Exchange field
		Eigen::MatrixXd Hexc(3,nv);
		for(int i = 0; i<3; i++)
		{
			// Update the entire xyz row
			Hexc.row(i) = -2*Aexc * (A*(m.row(i).transpose())).transpose();
			
			// Update each node to be at the correct units
			for(int nd = 0; nd<nv; nd++)
			{
				Hexc(i,nd) /= Ms*pow(scl,2)*Lagrange[nd];
			} 
		}

		// Anisotropy field
		Eigen::MatrixXd Han(3,nv);
		for(int i = 0; i<nv; i++)
		{
			Han.col(i) = 2*Aan/Ms *uan.dot(m.col(i)) * uan;
		}

		// Effective field
		Heff += Hd + Hexc + Han;

		// H = M cross Heff
		Eigen::MatrixXd H(3,nv);
		for(int i = 0; i<nv; i++)
		{
			// Eigen::Vector3d aux = m.col(i).head<3>().cross(Heff.col(i).head<3>());
			H.col(i) = m.col(i).head<3>().cross(Heff.col(i).head<3>());
		}

		// New magnetization
		Eigen::MatrixXd mOld = m; // Store the old magnetization

		#pragma omp parallel for
		for(int i = 0; i<nv; i++)
		{
			timeStep(m.col(i).head<3>(),
					 H.col(i),
					 H.col(i),
					 Heff.col(i),
					 dt,
					 1.0,
					 1.0,
					 0.0); // Considering time was normalized by giro
		} // New magnetization

		// Store initial guess of magnetic field
		Eigen::MatrixXd Hold = H;
		Eigen::MatrixXd HeffOld = Heff;

		// <magnetization>
		

		// Time iteration
		this->att = 0; 
		double dtau = 2*maxTorque;
		while ( dtau > maxTorque && att < maxAtt)
		{

			// -- New magnetic field --

			// Reset effective field to Hext
			for(int i = 0; i<nv; i++)
			{
				Heff.col(i) = mu0 *Hext;
			}

			// Demagnetizing field
			Hd.setZero();
			BEMdmag(Hd, p, t, surfaceT, normal, AE, VE, Vn.data(), lu, m);
			Hd *= mu0*Ms;

			// Exchange field
			for(int i = 0; i<3; i++)
			{
				// Update the entire xyz row
				Hexc.row(i) = -2*Aexc * (A*(m.row(i).transpose())).transpose();
				
				// Update each node to be at the correct units
				for(int nd = 0; nd<nv; nd++)
				{
					Hexc(i,nd) /= Ms*pow(scl,2)*Lagrange[nd];
				}
			}

			// Anisotropy field
			for(int i = 0; i<nv; i++)
			{
				Han.col(i) = 2*Aan/Ms *uan.dot(m.col(i)) * uan;
			}

			// Add all field contributions
			Heff += Hd + Hexc + Han;

			// New H = M cross Heff
			for(int i = 0; i<nv; i++)
			{
				// Eigen::Vector3d aux = m.col(i).head<3>().cross(Heff.col(i).head<3>());
				H.col(i) = m.col(i).head<3>().cross(Heff.col(i).head<3>());
			}

			// Calculate steepest descent step size
			double snN, snD, snD2;
			snN = 0.0;
			snD = 0.0;
			snD2 = 0.0;
			Eigen::Vector3d sn, gn2, gn1;
			for(int i = 0; i<nv; i++)
			{
				sn = m.col(i) - mOld.col(i);
				
				// m cross (m cross Heff)
				gn2 = m.col(i).head<3>().cross( m.col(i).head<3>().cross(Heff.col(i).head<3>()) );
		
				// m old cross (m old cross Heff old)
				gn1 = mOld.col(i).head<3>().cross( mOld.col(i).head<3>().cross(HeffOld.col(i).head<3>()) );

				snN += sn.dot(sn);
				snD += sn.dot(gn2-gn1);
				snD2 += (gn2-gn1).dot(gn2-gn1);
			}

			double tau1 = snN/snD;
			double tau2 = snD/snD2;

			dt = att%2 > 0 ? tau1 : tau2;

			// New magnetization
			#pragma omp parallel for
			for(int i = 0; i<nv; i++)
			{
				// Update mOld to store the latest m before iteration
				mOld.col(i) = m.col(i);
				
				// Update m
				nextM(m.col(i).head<3>(),
					  Heff.col(i),
					  dt); // Considering time was normalized by giro

			} // New magnetization

			// Update old magnetic fields
			Hold = H;
			HeffOld = Heff;

			// Average magnetization
			M_avg.col(att) = m.rowwise().mean();

			// <|dm/dt|> , "Torque"
			dtau = 0.0;
			for(int i = 0; i<nv; i++){
				dtau += (m.col(i).head<3>().cross(Heff.col(i).head<3>())).dot(m.col(i).head<3>().cross(Heff.col(i).head<3>()));
			}
			dtau /= nv;

			// Update tic
			att++;
			
			if(!(att%100)) // Log results every 100 iterations
			{
				std::cout << dtau << std::endl;
			}

		} // Time iteration

	} // Steepest descent solver


};

// For testing only
// int main(int argc, char const *argv[])
// {
// 	Eigen::MatrixXd m, p, normal, M_avg;
// 	Eigen::MatrixXi t, surfaceT;
// 	double* AE = (double*) malloc(3*sizeof(double));
// 	double* VE = (double*) malloc(3*sizeof(double));

// 	Micromagnetics micro(p,t,surfaceT,normal,AE,VE);

// 	free(AE);
// 	free(VE);
// 	return 0;
// }