#=
	See Femeko.jl/Micromagnetics/permalloy.jl
	This passes the solver logic to a C++ implementation
	
	/!\ /!\ Compile micromagnetics.cpp with 
		g++ -O3 -fPIC -shared -o micromagnetics.so micromagnetics.cpp
	
	or add multithreading with
		g++ -O3 -fPIC -shared -fopenmp -o micromagnetics.so micromagnetics.cpp
=#

include("../../src/Femeko.jl") 	# Include Landau-Lifshitz solver
using IterativeSolvers

using GLMakie 		# Include Makie for plots
using DelimitedFiles

meshSize = 200.0
localSize = 20.0
showGmsh = false

function main( meshSize::Float64=0.0
			 , localSize::Float64=0.0
			 , showGmsh::Bool=true)

	mu0::Float64 = pi*4e-7 			# Vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
	
	# Update LL parameters in the LL.h file
		Hext = [0.0, mu0* 50e3, 0.0] # Applied field (T)
		Ms = mu0 * 860e3	# Mag. saturation (T)
		# scale = 1e-9		# scale of the geometry | nm: 1e-9 m
		# Aexc = 12e-13	# Exchange   (J/m)
		# Aan = 0.0 		# Anisotropy constant J/m3
		# uan = [1,0,0] 	# Easy axis
		timeStep = 0.01 	# Time step (normalized by the gyromagnetic ratio)
		# nSteps = 7035
		# totalTime = 70.35 	# Stop when time > total time (normalized by the gyromagnetic ratio)
		# alfa = 0.1/Ms 	# damping
		
	# Create a 3D Model
	gmsh.initialize()
	cells = []
	addCuboid([0,0,0], [100, 100, 5], cells)
	box = addSphere([0,0,0], 500.0) # Add a bounding shell

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

	# Element centroids
	centroids::Matrix{Float64} = zeros(3, mesh.nt) # Element centroids
	for k in 1:mesh.nt
	    nds = mesh.t[:, k]
	    centroids[:, k] = mean(mesh.p[:, nds], 2)
	end

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

	# Node volumes
	Volumes::Vector{Float64} = zeros(mesh.nv)
	for k in mesh.InsideElements
		nds = @view mesh.t[:, k]
		Volumes[nds] .+= mesh.VE[k]
	end

	# Initial magnetization
	M::Matrix{Float64} = zeros(3, mesh.nv)
	M[1, mesh.InsideNodes] .= Ms

	# -- Solve the demag field with julia --

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

	# Global stiffness matrix
	A = stiffnessMatrix(mesh)

	# Magnetostatic potential
	u = cg(A, RHS) # Conjugate gradient solver

	# Demagnetizing field
	Hd::Matrix{Float64} = zeros(3, mesh.nt)
	for k in 1:mesh.nt
		nds = @view mesh.t[:, k]
		for i in 1:4
			Hd[1, k] -= b[i, k]*u[nds[i]]
			Hd[2, k] -= c[i, k]*u[nds[i]]
			Hd[3, k] -= d[i, k]*u[nds[i]]
		end
	end

	# Map the demag field to the mesh nodes
	Hd_nodes::Matrix{Float64} = zeros(3, mesh.nv)
	
	for k in mesh.InsideElements
		nds = @view mesh.t[:, k]
		Hd_nodes[:, nds] .+= mesh.VE[k]*Hd[:, k]
	end

	Hd_nodes[1, :] ./= Volumes
	Hd_nodes[2, :] ./= Volumes 
	Hd_nodes[3, :] ./= Volumes 
	
	H = zeros(mesh.nv)
	for nd in 1:mesh.nv
		H[nd] = norm(Hd_nodes[:, nd])
	end

	# Plot the results
	println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data, title="Julia")

    graph = arrows3d!(  ax
                        , mesh.p[1, mesh.InsideNodes]
                        , mesh.p[2, mesh.InsideNodes]
                        , mesh.p[3, mesh.InsideNodes]
                        , Hd_nodes[1, mesh.InsideNodes]
                        , Hd_nodes[2, mesh.InsideNodes]
                        , Hd_nodes[3, mesh.InsideNodes]
                        , color = H[mesh.InsideNodes]
                        , lengthscale = 100.0
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    # Add a colorbar
    Colorbar(fig[2, 1], graph, vertical=false)
    display(fig)

	# -- Solve the demag field with C++ --

	# Convert Julia 1 indexing to C++ 0 indexing
	t::Matrix{Int32} = mesh.t .- 1
	InsideElements::Vector{Int32} = mesh.InsideElements .- 1
	InsideNodes::Vector{Int32} = mesh.InsideNodes .- 1

	# Call C++ solver
	@ccall "./micromagnetics.so".LandauLifshitz(
	      mesh.p::Ptr{Float64}      # Node coordinates, 3 by nv
	    , t::Ptr{Int32}             # Node connectivity, 4 by nt
	    , mesh.VE::Ptr{Float64}  	# Volume of each element
	    , InsideElements::Ptr{Int32}
	    , InsideNodes::Ptr{Int32}
	    , mesh.nv::Int32         # Number of nodes
	    , mesh.nt::Int32         # Number of elements
	    , mesh.nInside::Int32
	    , mesh.nInsideNodes::Int32
	    , M::Ptr{Float64}       # 
	)::Cvoid

	# Add C++ results
	Hd_cpp = readdlm("Hd.txt") # T
	H_cpp = zeros(mesh.nv)
	for nd in 1:mesh.nv
		H_cpp[nd] = norm(Hd_cpp[:, nd])
	end

	println("Updating plot with C++ results...")
    ax = Axis3(fig[1, 2], aspect = :data, title="C++ with Eigen")

    graph = arrows3d!(  ax
                        , mesh.p[1, mesh.InsideNodes]
                        , mesh.p[2, mesh.InsideNodes]
                        , mesh.p[3, mesh.InsideNodes]
                        , Hd_cpp[1, mesh.InsideNodes]
                        , Hd_cpp[2, mesh.InsideNodes]
                        , Hd_cpp[3, mesh.InsideNodes]
                        , color = H_cpp[mesh.InsideNodes]
                        , lengthscale = 100.0
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                      )

    # Add a colorbar
    Colorbar(fig[2, 2], graph, vertical=false)
    wait(fig.scene)

end

main(meshSize, localSize, showGmsh)