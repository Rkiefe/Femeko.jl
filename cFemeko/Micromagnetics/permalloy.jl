#=
	See Femeko.jl/Micromagnetics/permalloy.jl
	This passes the solver logic to a C++ implementation
=#

include("../../src/Femeko.jl") 	# Include Landau-Lifshitz solver
using GLMakie 		# Include Makie for plots

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
		nSteps = 7035
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

	# Initial magnetization
	M::Matrix{Float64} = zeros(3, mesh.nv)
	M[1, mesh.InsideNodes] .= Ms

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


end

main(meshSize, localSize, showGmsh)