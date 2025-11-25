#=
	Landau-Lifshizt equation solver
=#

include("../src/Femeko.jl")
include("../src/magneticProperties.jl")

using IterativeSolvers

mutable struct LL
	# Methods
	run::Function 		# Run the LL for t = 0 to totalTime

	# Time iteration settings
	timeStep::Float64 	# Time step of the solver (normalized by the giromagnetic ratio)
	totalTime::Float64 	# Total time of the LL simulation
	nSteps::Int   		# The number of expected time steps
	
	maxTorque::Float64  # Tolerance, solver stops when |dM/dt| < maxTorque
	# maxSteps::Int 		# Maximum allowed number of time steps on the minimization

	# Landau-Lifshitz-Bernoulli coefficients
	alfa::Float64 # Damping term
	precession::Bool # Consider the precession term ?

	Hext::Vector{Float64} 	# External field (Tesla)
	Ms::Float64  # Saturation magnetization (Tesla)

	Aan::Float64 # Anisotropy constant (J/m3)
	uan::Vector  # Easy axis

	Aexc::Float64 # Exchange constant (J/m)

	scale::Float64 # scale of the geometry | nm: 1e-9 m

	# Mesh data
	mesh::MESH

	# Constructor
	LL() = new(run, 0.1, 0.0, 10_000, 1e-5, 1.0, true)
	# LL() = new(run) # No default parameters, only the method run()
end

function run(  self::LL
			 , M::Matrix{Float64}
			 , A 					# Global stiffness matrix
			 , AEXC 				# Stiffness matrix for the exchange field
			 , b::Matrix{Float64}	# Basis function coefs. _,b,c,d
			 , c::Matrix{Float64}
			 , d::Matrix{Float64}
			 , Volumes::Vector{Float64} # Volume of elements of each node
			 , showTorque::Bool = false  # Print max|dM/dt| of each simulation step
			 )

	# Set the number of iteration steps based on the parameters
	if self.totalTime > self.timeStep
		self.nSteps = floor(self.totalTime/self.timeStep) + 1
	
	# else # -> The solver runs until dM/dt = 0 or reaches the max number of steps
		# Use the default self.nSteps as the max number of steps
	end

	# Initial magnetostatic field | Nodes
	Hdfield = magnetostaticField(self.mesh, M, A, b, c, d, Volumes)

	# Anisotropy field
	Han::Matrix{Float64} = anisotropyField(self, M)

	# Exchange field
	Hexc::Matrix{Float64} = exchangeField(self, M, AEXC, Volumes)

	# Prepare simulation
	Mold::Matrix{Float64} = deepcopy(M)
	H::Matrix{Float64} = zeros(3, self.mesh.nv)
	Hold::Matrix{Float64} = zeros(3, self.mesh.nv)

	Mnorm::Vector{Float64} = zeros(self.mesh.nv)
	for nd in self.mesh.InsideNodes
		Mnorm[nd] = norm(M[:, nd])
	end

	# <M> (in time)
	Mx::Vector{Float64} = zeros(self.nSteps)
	My::Vector{Float64} = zeros(self.nSteps)
	Mz::Vector{Float64} = zeros(self.nSteps)

	# Run simulation
	println("Running LL simulation")
	frame::Int32 = 0
	torque::Float64 = Inf # maximum value of |dM/dt|

	while (self.totalTime > self.timeStep && frame < self.nSteps) || (self.totalTime < self.timeStep && torque > self.maxTorque && frame < self.nSteps)

		frame += 1
		torque = 0.0 # reset maximum torque

		# Evolve the magnetization on each node
		for nd in self.mesh.InsideNodes

			# Get the total magnetic field on current node
			H[:, nd] = self.Hext + Hdfield[:, nd] + Han[:, nd]

			# Update magnetization
			M2= step(  M[:, nd]    	# M(n)
				 	 , Mold[:, nd] 	# M(n-1)
				 	 , H[:, nd]    	# H(n)
				 	 , Hold[:, nd] 	# H(n-1)
		         	 , self.timeStep# Time step
		         	 , self.alfa 	# LL damping
		         	 , self.precession) # Consider Precession

			Mold[:, nd] .= M[:, nd]
			M[:, nd] .= M2

			# Store the old magnetic field
			Hold[:, nd] .= H[:, nd]

			# Update the magnetizaiton norm
			Mnorm[nd] = norm(M[:, nd])

			# Check the torque term dM/dt
			dM_dt =   cross(M[:, nd], H[:, nd]) 
					+ self.alfa*cross(M[:, nd], cross(M[:, nd], H[:, nd])) 

			torque = max(torque, norm(dM_dt)) # Store the largest torque value
		
		end # Update M

		# Get the new magnetostatic field from the new magnetization
		Hdfield = magnetostaticField(self.mesh, M, A, b, c, d, Volumes)

		# Anisotropy field with the new M
		Han = anisotropyField(self, M)

		# Exchange field
		Hexc = exchangeField(self, M, AEXC, Volumes)

		# Average x,y,z components in Tesla
		Mx[frame] = mean(M[1, self.mesh.InsideNodes])
		My[frame] = mean(M[2, self.mesh.InsideNodes])
		Mz[frame] = mean(M[3, self.mesh.InsideNodes])

		# Print the largest torque value of each iteration step
		if showTorque
			println("|dM/dt| = ", torque)
		end

	end # Runs the LL simulation
	println("Simulation finished")
	
	# Cut up to the last frame
	Mx = Mx[1:frame]
	My = My[1:frame]
	Mz = Mz[1:frame]

	self.nSteps = frame

	return M, Mnorm, Mx, My, Mz, H
end

# Demagnetizing field in units of M
function magnetostaticField(mesh::MESH, M::Matrix{Float64}
							, A # Sparse matrix
							, b::Matrix{Float64}
							, c::Matrix{Float64}
							, d::Matrix{Float64}
							, Volumes::Vector{Float64})

	# Load vector 
	RHS::Vector{Float64} = zeros(mesh.nv)
	for k in mesh.InsideElements
		nds = @view mesh.t[:, k]

		# Mean magnetization on the element
		f::Vector{Float64} = mean(M[:, nds], 2)
		
		# Update load vector
		for i in 1:4
			RHS[nds[i]] += mesh.VE[k]*dot([b[i, k], c[i, k], d[i, k]], f)
		end

	end # Load vector

	# Magnetostatic potential
	u = cg(A, RHS) # Conjugate gradient solver

	# Demagnetizing field
	Hdfield::Matrix{Float64} = zeros(3, mesh.nt)
	for k in 1:mesh.nt
		nds = @view mesh.t[:, k]
		for i in 1:4
			Hdfield[1, k] -= b[i, k]*u[nds[i]]
			Hdfield[2, k] -= c[i, k]*u[nds[i]]
			Hdfield[3, k] -= d[i, k]*u[nds[i]]
		end
	end

	# Magnetic field on the nodes
	Hdfield_nodes::Matrix{Float64} = zeros(3, mesh.nv)
	
	for k in mesh.InsideElements
		nds = @view mesh.t[:, k]
		Hdfield_nodes[:, nds] .+= mesh.VE[k]*Hdfield[:, k]
	end

	Hdfield_nodes[1, :] ./= Volumes
	Hdfield_nodes[2, :] ./= Volumes 
	Hdfield_nodes[3, :] ./= Volumes 

	return Hdfield_nodes
end # Demagnetizing field

# Anisotropy field (tesla)
function anisotropyField( self::LL, M::Matrix{Float64})
	
	mu0::Float64 = pi*4e-7 # Vaccum magnetic permeability
	Han::Matrix{Float64} = zeros(3, self.mesh.nv)
	for i in self.mesh.InsideNodes
		Msqrd::Float64 = dot(M[:, i], M[:, i])
	    Han[:, i] = mu0*2*self.Aan/Msqrd *dot(M[:, i], self.uan) .*self.uan
	end
	return Han
end # Anisotropy field (T)

# Exchange field (tesla)
function exchangeField(self::LL, M::Matrix{Float64}, AEXC, Volumes::Vector{Float64})
	
	mu0::Float64 = pi*4e-7 # Vaccum magnetic permeability
	Hexc::Matrix{Float64} = -2*self.Aexc.* (AEXC*M')'

	for i in 1:3
		Hexc[i, :] ./= (0.25*self.Ms^2 *self.scale^2) .* Volumes
	end

	return Hexc
end # Exchange field (T)

# Find new magnetization after time iteration
function step(  M::Vector{Float64},   	# M(n)
			  	Mold::Vector{Float64},	# M(n-1)
			  	H::Vector{Float64},   	# H(n)
			  	Hold::Vector{Float64},	# H(n-1)
              	dt::Float64, 			# Time step
              	alfa::Float64=1.0, 		# LL damping
              	precession::Bool=true 	# Consider the precession term
              )
    #=
        Repeats the search of a new magnetization until the solution is stable
    =#

    d::Float64 = dt/2

    Mnorm = norm(M)
    Mnorm_old = norm(Mold)

    # M(n+1/2)
    M12::Vector{Float64} = 3/2*M - 1/2*Mold

    # M(n+1)
    M2::Vector{Float64} = [0.0, 0.0, 0.0]

    # 1) Initial guess of the new magnetic field
    H12::Vector{Float64} = 3/2 *H - 0.5 *Hold

    # Htild from Oriano 2008
    Htild::Vector{Float64} = alfa.*cross(M12, H12) + H12

    # Repeat until M(n+1) doesn't change
    aux::Vector{Float64} = deepcopy(M)
    err::Float64 = Inf
    att::Int32 = 0
    while err > 1e-6 && att < 1_000
        att += 1

        # 2) M (n+1) from M(n) and H(n+1/2)
        mat::Matrix{Float64} = [1 d*Htild[3] -d*Htild[2];
                                -d*Htild[3] 1 d*Htild[1];
                                d*Htild[2] -d*Htild[1] 1]
        
        M2 = mat\( M - d*cross(M, Htild) ) # M(n+1)

        # 3) M (n + 1/2)
        M12 = 0.5*(M + M2)

        # 4) Calculate H~ (n+1/2) from M(n+1/2)
        Htild = precession ?  alfa.*cross(M12, H12) + H12 : alfa.*cross(M12, H12)

        # Max difference between new M(n+1) and old M(n+1)
        err = maximum(abs.(M2-aux))
        # println("M(n+1) iteration error: ", err)

        # Update M(n+1) from last iteration
        aux .= M2

    end # Get new magnetization value

    return M2
end # Find new magnetization after time iteration
