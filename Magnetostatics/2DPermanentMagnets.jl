#=
	2D magnetostatic simulation of two permanent magnets intearcting
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)
	gmsh.initialize()
	cells = []

	# Vaccum magnetic permeability
	mu0::Float64 = pi*4e-7

	# Magnets
	M0::Vector{Float64} = [0.0, 1.0] 	# Magnetic field intensity of the magnet | Tesla

	# Add two magnets
	magnetCells::Vector{Int32} = [0, 0] # Store the cell id of each magnet

	magnetCells[1] = addRectangle([0.0,  1.0], [1.0, 1.0], cells)
	magnetCells[2] = addRectangle([0.0, -1.0], [1.0, 1.0], cells)

	# Add a bounding shell
	box = addDisk([0.0, 0.0], 5.0)
	shell_id = unifyModel(cells, box)

    # Generate mesh
	mesh::MESH = Mesh2D(cells, meshSize, localSize)

	# Get element tags to then use GMSH 'get' functions
	t_tags, _ = gmsh.model.mesh.getElementsByType(2) # 3 node first order triangle

    # Store cell id of each element
    elementID::Vector{Int32} = zeros(mesh.nt)
    for k in 1:mesh.nt
        # element type , nodes of the element , dimension , id
        _, _, _, id = gmsh.model.mesh.getElement(t_tags[k])
        elementID[k] = id
    end

	println("\nBox cell ID: ", box)
	println("Inner cells: ", cells)
	println("Outer boundary curve ID: ", shell_id, "\n")

	# Run Gmsh GUI
    if showGmsh
	   gmsh.fltk.run()
    end
	gmsh.fltk.finalize()

	# Element centroid
	centroids::Matrix{Float64} = zeros(2, mesh.nt)
	for k in 1:mesh.nt
		nds = @view mesh.t[:, k]
		centroids[:, k] = mean(mesh.p[1:2, nds], 2)
	end

	# Lagrange multiplier technique
    Lag::Vector{Float64} = zeros(mesh.nv)
    for k = 1:mesh.nt
        nds = @view mesh.t[:,k]
        Lag[nds] .+= mesh.VE[k]/3
    end

	# Stiffness matrix
	A = stiffnessMatrix2D(mesh)

	# Load vector from magnets
	RHS = zeros(mesh.nv)
	for k in 1:mesh.nt

		# Only integrate over the magnet volume
		if !(elementID[k] in magnetCells)
			continue
		end

		nds = @view mesh.t[:, k]
		for i in 1:3
			_, b, c = abc(mesh.p, nds, nds[i])

			RHS[nds[i]] += mesh.VE[k]*( M0[1]*b + M0[2]*c )
		end
	end # Load vector

	# Magnetic scalar potential
	u = [A Lag;Lag' 0]\[RHS;0]

	# Magnetic field
	Hfield::Matrix{Float64} = zeros(2, mesh.nt)
	for k in 1:mesh.nt
	    nds = @view mesh.t[:,k];

	    # Sum the contributions
	    for nd in nds
	        _, b, c = abc(mesh.p, nds, nd)

	        Hfield[1, k] -= u[nd]*b;
	        Hfield[2, k] -= u[nd]*c;
	    end
	end

	# Magnetic field intensity
	H::Vector{Float64} = zeros(mesh.nt)
	for k in 1:mesh.nt
		H[k] = norm(Hfield[:, k])
	end

	# Plot result
	println("Generating plots...")

	fig = Figure()
	ax = Axis(fig[1,1], aspect = DataAspect())

	graph = arrows2d!(	ax,
						centroids[1,:],
						centroids[2,:],
						Hfield[1,:]./maximum(H),
						Hfield[2,:]./maximum(H)
						, color = H./mu0
						, lengthscale = 0.1
						, colormap = :turbo
					)

	Colorbar(fig[1, 2], graph, label = "H (A/m)")

	wait(display(fig))

end

meshSize = 2.0
localSize = 0.1
showGmsh = false

main(meshSize, localSize, showGmsh)