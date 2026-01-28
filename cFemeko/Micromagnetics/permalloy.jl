#= 
	Same as the Femeko.jl permalloy.jl, but with C++.
	Make sure to update the LL.h header file with the correct values and compile with
		g++ -O3 -fPIC -shared -fopenmp -o micromagnetics.so micromagnetics.cpp
=#

include("../../src/Femeko.jl")
using GLMakie
using DelimitedFiles

meshSize = 200.0
localSize = 20.0
showGmsh = false

function main(meshSize::Float64 = 0.0
			, localSize::Float64 = 0.0
			, showGmsh::Bool = false)
	
	gmsh.initialize()	

	# Create 3D model
	cells = []
	addCuboid([0.0, 0.0, 0.0], [100.0, 100.0, 5.0], cells)
	box = addSphere([0.0, 0.0, 0.0], 500.0)

	_, box = unifyModel(cells, box) # unify volume

	# Generate mesh
	mesh = Mesh(cells, meshSize, localSize)
	
	# View mesh in Gmsh gui
	if showGmsh
		gmsh.option.setNumber("Mesh.Clip", 1)
		gmsh.option.setNumber("Mesh.VolumeFaces", 1)
		gmsh.option.setNumber("General.ClipWholeElements", 1)
		gmsh.fltk.run()
	end
	gmsh.finalize()

	mu0::Float64 = pi*4e-7 			# Vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
	timeStep = 0.01
	nSteps = 70350

	# Normalized 3D magnetization field
	M = zeros(3, mesh.nv)
	M[1, mesh.InsideNodes] .= 1.0807
	# M[:, mesh.InsideNodes] .= rand(3, mesh.nInsideNodes) .- 0.5
	# for i in mesh.InsideNodes
	# 	M[:, i] ./= norm(M[:, i])
	# end

	# Shift to 0 index for C++
	t::Matrix{Int32} = mesh.t .- 1
	InsideElements::Vector{Int32} = mesh.InsideElements .- 1
	InsideNodes::Vector{Int32} = mesh.InsideNodes .- 1

	# Output plot
	println("Generating plots...")
	fig = Figure()
	ax = Axis3(fig[1, 1], aspect = :data, title="Initial M")
	graph = arrows3d!(  ax
	                    , mesh.p[1, mesh.InsideNodes]
	                    , mesh.p[2, mesh.InsideNodes]
	                    , mesh.p[3, mesh.InsideNodes]
	                    , M[1, mesh.InsideNodes]
	                    , M[2, mesh.InsideNodes]
	                    , M[3, mesh.InsideNodes]
	                    , color = M[3, mesh.InsideNodes]
	                    , lengthscale = 10.0
	                    , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
	                  )

	# Colorbar(fig[1, 2], graph, label = "M (emu/g)") # Add a colorbar
	# display(fig)
	# wait(display(fig))
	# display(GLMakie.Screen(), fig)


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


	# Plot new magnetization
	println("Adding new M to the plot")
	ax = Axis3(fig[1, 2], aspect = :data, title="Final M")
	graph = arrows3d!(  ax
	                    , mesh.p[1, mesh.InsideNodes]
	                    , mesh.p[2, mesh.InsideNodes]
	                    , mesh.p[3, mesh.InsideNodes]
	                    , M[1, mesh.InsideNodes]
	                    , M[2, mesh.InsideNodes]
	                    , M[3, mesh.InsideNodes]
	                    , color = M[3, mesh.InsideNodes]
	                    , lengthscale = 10.0
	                    , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
	                  )
	# Colorbar(fig[1, 2], graph, label = "M (emu/g)") # Add a colorbar	
	display(GLMakie.Screen(), fig)
	# display(fig)



	# Read the output
	Mxyz = readdlm("M_time.txt")

	println("Making a new window with M(t)")
	fig = Figure()
	ax = Axis(fig[1,1]) #, aspect =  DataAspect()
	scatter!(ax, (timeStep*1e9/giro)*(0:nSteps-1), Mxyz[1, :], label="Mx")
	scatter!(ax, (timeStep*1e9/giro)*(0:nSteps-1), Mxyz[2, :], label="My")
	scatter!(ax, (timeStep*1e9/giro)*(0:nSteps-1), Mxyz[3, :], label="Mz")
    axislegend(position=:rb)
	
	# display(GLMakie.Screen(), fig)
	wait(display(fig))
	# wait(fig.scene)

end

main(meshSize, localSize, showGmsh)


