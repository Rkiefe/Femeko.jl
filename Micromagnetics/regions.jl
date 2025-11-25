#=
	Replicates the mumax3
=#

include("LL.jl") 	# Include Landau-Lifshitz solver
using GLMakie 		# Include Makie for plots

function exchangeField( M::Matrix{Float64}, 
						Ms::Float64, 
						Aexc::Float64,  # Exchange constant
						AEXC, 			# Exchange stiffness matrix
						Volumes::Vector{Float64},
						scale::Float64)

	mu0::Float64 = pi*4e-7 # Vaccum magnetic permeability
	Hexc::Matrix{Float64} = -2*Aexc.* (AEXC*M')'

	for i in 1:3
		Hexc[i, :] ./= (0.25*Ms^2 *scale^2) .* Volumes
	end
	return Hexc
end

function main(meshSize::Float64=0.0, localSize::Float64=0.0, showGmsh::Bool=true)

	mu0::Float64 = pi*4e-7 			# Vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
	
	# Update LL parameters
	Hext = [0.0, 0.0, 0.0] 		# Applied field (T)
	Ms::Float64 = mu0 * 800e3	# Mag. saturation (T)
	scale = 1e-9				# scale of the geometry | nm: 1e-9 m

	Aexc = 12e-13	# Exchange   (J/m)

	Aan = [1e5, 2e5] # Anisotropy constant J/m3 of each region
	uan = [[0,1,0], [1,0,0]] 	# Easy axis
	
	timeStep = 0.05		# Time step (normalized by the gyromagnetic ratio)
	totalTime = 10.0 	# Stop when time > total time (normalized by the gyromagnetic ratio)
	# maxTorque = 1e-5   # If 'totalTime' is not provided, it minimizes the Energy/M state
	
	alfa = 1.0 	# damping

	showTorque::Bool = true # Print the dM/dt term on each time step
	
	# Create a 3D Model
	gmsh.initialize()
	cells = []
	addCylinder([0,0,0], [0,0,4], 256.0, cells)
	box = addSphere([0,0,0], 2500.0) # Add a bounding shell

	unifyModel(cells, box)

	# Create a mesh
	mesh = Mesh(cells, meshSize, localSize)

	# Print number of elements and nodes
	println("\nNumber of elements: ", mesh.nt)
	println("Number of internal elements: ", mesh.nInside)
	println("Number of internal nodes: ", mesh.nInsideNodes, "\n")

	if showGmsh
		gmsh.option.setNumber("Mesh.Clip", 1)
		gmsh.option.setNumber("Mesh.VolumeFaces", 1)
		gmsh.option.setNumber("General.ClipWholeElements", 1)
		gmsh.fltk.run()
	end
	gmsh.finalize()

	left = []
	right = []
	for nd in mesh.InsideNodes
		# x of the node is on the left
		if mesh.p[1, nd] < 0.0
			push!(left, nd)
		else # node is on the right
			push!(right, nd)
		end
	end

	# Show left and right regions
	# fig = Figure()
	# ax = Axis(fig[1, 1]
	# 			# , aspect = DataAspect()
	# 		  )
	# scatter!(ax, 
	#          mesh.p[1, left], 
	#          mesh.p[2, left], 
	#          mesh.p[3, left], label="left")

	# scatter!(ax, 
	#          mesh.p[1, right], 
	#          mesh.p[2, right], 
	#          mesh.p[3, right], label="right")

    # axislegend() # position=:rb
    # wait(display(fig)); return

	# FEM
	print("Calculating the Lagrange shape elements... ")
	b::Matrix{Float64} = zeros(4, mesh.nt)
	c::Matrix{Float64} = zeros(4, mesh.nt)
	d::Matrix{Float64} = zeros(4, mesh.nt)
	for k in 1:mesh.nt
		nds = @view mesh.t[:, k]
		for i in 1:4
			_, b[i, k], c[i, k], d[i, k] = abcd(mesh.p, nds, nds[i])
		end
	end
	println("Done.")

	# Stiffness matrix | Exchange field 
    AEXC = spzeros(mesh.nv, mesh.nv)

	Ak::Matrix{Float64} = zeros(4*4, mesh.nt)
	for k in mesh.InsideElements
	    nds = @view mesh.t[:,k]
	    Ak[:, k] = vec( mesh.VE[k]*( b[:, k]*b[:, k]' 
										+ c[:, k]*c[:, k]'
										+ d[:, k]*d[:, k]' ) )
	end
    
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            AEXC += sparse(mesh.t[i,:], mesh.t[j,:], Ak[n,:], mesh.nv, mesh.nv)
        end
    end # Stiffness matrix for Exchange Field

	# Node volumes
	Volumes::Vector{Float64} = zeros(mesh.nv)
	for k in mesh.InsideElements
		nds = @view mesh.t[:, k]
		Volumes[nds] .+= mesh.VE[k]
	end

	# Global stiffness matrix
	A = stiffnessMatrix(mesh)

	# Initial magnetization state
	M::Matrix{Float64} = zeros(3, mesh.nv)
	
	# Left side
	direction = [-1, 1, 0]/norm([-1, 1, 0])
	for nd in left
		M[:, nd] = Ms*direction
	end

	# Right side
	direction = [1, 1, 0]/norm([1, 1, 0])
	for nd in right
		M[:, nd] = Ms*direction
	end

	# Run
	nSteps::Int = 0
	if totalTime > timeStep
		nSteps = floor(totalTime/timeStep) + 1
	
	else # -> The solver runs until dM/dt = 0 or reaches the max number of steps
		# Use the default self.nSteps as the max number of steps
		nSteps = 10_000
	end

	# Initial magnetostatic field | Nodes
	Hdfield::Matrix{Float64} = magnetostaticField(mesh, M, A, b, c, d, Volumes)

	# Anisotropy field
	Han::Matrix{Float64} = zeros(3, mesh.nv)
	for i in left
		Msqrd::Float64 = dot(M[:, i], M[:, i])
	    Han[:, i] = mu0*2*Aan[1]/Msqrd *dot(M[:, i], uan[1]) .* uan[1]
	end
	for i in right
		Msqrd::Float64 = dot(M[:, i], M[:, i])
	    Han[:, i] = mu0*2*Aan[2]/Msqrd *dot(M[:, i], uan[2]) .* uan[2]
	end

	# Exchange field
	Hexc::Matrix{Float64} = exchangeField(M, Ms, Aexc, AEXC, Volumes, scale)

	# Prepare simulation
	Mold::Matrix{Float64} = deepcopy(M)
	H::Matrix{Float64} = zeros(3, mesh.nv)
	Hold::Matrix{Float64} = zeros(3, mesh.nv)

	Mnorm::Vector{Float64} = zeros(mesh.nv)
	for nd in mesh.InsideNodes
		Mnorm[nd] = norm(M[:, nd])
	end

	# <M> (in time)
	Mx::Vector{Float64} = zeros(nSteps)
	My::Vector{Float64} = zeros(nSteps)
	Mz::Vector{Float64} = zeros(nSteps)
	
	# Run simulation
	println("Running LL simulation")
	frame::Int32 = 0
	torque::Float64 = Inf # maximum value of |dM/dt|

	while (totalTime > timeStep && frame < nSteps) || (totalTime < timeStep && torque > maxTorque && frame < nSteps)

		frame += 1
		torque = 0.0 # reset maximum torque

		# Evolve the magnetization on each node
		for nd in mesh.InsideNodes

			# Get the total magnetic field on current node
			H[:, nd] = Hext + Hdfield[:, nd] + Han[:, nd]

			# Update magnetization
			M2= step(  M[:, nd]    	# M(n)
				 	 , Mold[:, nd] 	# M(n-1)
				 	 , H[:, nd]    	# H(n)
				 	 , Hold[:, nd] 	# H(n-1)
		         	 , timeStep# Time step
		         	 , alfa 	# LL damping
		         	 , true) # Consider Precession

			Mold[:, nd] .= M[:, nd]
			M[:, nd] .= M2

			# Store the old magnetic field
			Hold[:, nd] .= H[:, nd]

			# Update the magnetizaiton norm
			Mnorm[nd] = norm(M[:, nd])

			# Check the torque term dM/dt
			dM_dt =   cross(M[:, nd], H[:, nd]) 
					+ alfa*cross(M[:, nd], cross(M[:, nd], H[:, nd])) 

			torque = max(torque, norm(dM_dt)) # Store the largest torque value
		
		end # Update M

		# Get the new magnetostatic field from the new magnetization
		Hdfield = magnetostaticField(mesh, M, A, b, c, d, Volumes)

		# Anisotropy field
		Han = zeros(3, mesh.nv)
		for i in left
			Msqrd::Float64 = dot(M[:, i], M[:, i])
		    Han[:, i] = mu0*2*Aan[1]/Msqrd *dot(M[:, i], uan[1]) .* uan[1]
		end
		for i in right
			Msqrd::Float64 = dot(M[:, i], M[:, i])
		    Han[:, i] = mu0*2*Aan[2]/Msqrd *dot(M[:, i], uan[2]) .* uan[2]
		end

		# Exchange field
		Hexc = exchangeField(M, Ms, Aexc, AEXC, Volumes, scale)

		# Average x,y,z components in Tesla
		Mx[frame] = mean(M[1, mesh.InsideNodes])
		My[frame] = mean(M[2, mesh.InsideNodes])
		Mz[frame] = mean(M[3, mesh.InsideNodes])

		# Print the largest torque value of each iteration step
		if showTorque
			println("|dM/dt| = ", torque)
		end

	end # Runs the LL simulation
	println("Simulation finished")
	
	# Plot the M(t)
	fig = Figure()
	ax = Axis3(fig[1,1]
				, aspect = :data
				)
	
	graph = arrows3d!(  ax
	                    , mesh.p[1, mesh.InsideNodes]
	                    , mesh.p[2, mesh.InsideNodes]
	                    , mesh.p[3, mesh.InsideNodes]
	                    , M[1, mesh.InsideNodes]
	                    , M[2, mesh.InsideNodes]
	                    , M[3, mesh.InsideNodes]
	                    , color = Mnorm[mesh.InsideNodes]./(mu0*1e3)
	                    , lengthscale = 50.0
	                    , colormap = :rainbow,  # :CMRmap :viridis :redsblues :turbo :rainbow
	                )
	
	# Colorbar(fig[1, 2], graph, label = "M (kA/m)"
	#          #, vertical = false
	#          )

    # save("M_time_permalloy.png", fig)
	wait(display(fig))

end

meshSize = 500.0
localSize = 75.0
showGmsh = false
main(meshSize, localSize, showGmsh)