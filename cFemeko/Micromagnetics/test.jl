include("../../src/Femeko.jl")

function main(meshSize::Float64 = 0.0
			, localSize::Float64 = 0.0
			, showGmsh::Bool = false)
	
	gmsh.initialize()	

	# Create 3D model
	cells = []
	addCuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], cells)
	box = addSphere([0.0, 0.0, 0.0], 5.0)

	_, box = unifyModel(cells, box) # unify volume

	# Generate mesh
	mesh = Mesh(cells, meshSize, localSize)
	
	# View mesh in Gmsh gui
	showGmsh ? gmsh.fltk.run() : nothing
	gmsh.finalize()

	# Normalized 3D magnetization field
	M = zeros(3, mesh.nv)
	M[:, mesh.InsideNodes] .= rand(3, mesh.nInsideNodes) .- 0.5
	for i in mesh.InsideNodes
		M[:, i] ./= norm(M[:, i])
	end

	# Shift to 0 index for C++
	t::Matrix{Int32} = mesh.t .- 1
	InsideElements::Vector{Int32} = mesh.InsideElements .- 1
	InsideNodes::Vector{Int32} = mesh.InsideNodes .- 1

	# Send mesh data to the C++ micromagnetic solver
	@ccall "./micromagnetics.so".LandauLifshitz(
		  mesh.p::Ptr{Float64}
	    , t::Ptr{Int32}
	    , mesh.VE::Ptr{Float64}   
	    , InsideElements::Ptr{Int32}
	    , InsideNodes::Ptr{Int32}
	    , mesh.nv::Int32, mesh.nt::Int32
	    , mesh.nInside::Int32, mesh.nInsideNodes::Int32
    	, M::Ptr{Float64}
    	)::Cvoid

end

main(0.0, 0.0, false)

# g++ -O3 -fPIC -shared -fopenmp -o micromagnetics.so micromagnetics.cpp
# g++ -fPIC -shared -fopenmp -o micromagnetics.so micromagnetics.cpp