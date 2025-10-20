#=
    	Internal magnetic field using FEM/BEM for non-linear magnetic materials

    This implementation is based on https://doi.org/10.1016/j.jmmm.2012.01.016
    but adapted for non-linear magnetic materials and any magnetic field souce (non conservative)        	

	The non-linearity is handled by the Fixed-Point iteration method

	This example considers different magnetic bodies, interacting at a distance, without meshing
	the space between them - showcasing the advantage of FEM-BEM
=#


include("../../src/gmsh_wrapper.jl")
include("../../src/BEM.jl")
include("../../src/magneticProperties.jl")

using GLMakie

# Goes to each volume cell and updates the permeability on the 
# mesh elements of that cell
function updatePermeability(mu::Vector{Float64},
							H::Vector{Float64}, 
							materialProperties, 
							cellLabels, 
							elementID)

	# Update magnetic permeability
	for (id, key) in enumerate(cellLabels) # cell number, string with the data label

	    spl = Spline1D(materialProperties[key].HofM,
	                   materialProperties[key].mu
	                   ) # ;bc="nearest") # nearest , extrapolate

	    # Find all elements of current cell ID
	    elements = findall(x -> x==id, elementID)
	    
	    # Interpolate the dataset for this elements
	    mu[elements] .= spl(H[elements])

	    idx = findall(findErr -> !isfinite(findErr), mu)
	    if !isempty(idx)
	        println(idx)
	        error("Nans/Infs in mu")
	    end
	end

end # Update element-wise permeability


function main(meshSize=0.0, showGmsh=false)
	gmsh.initialize()

	# Applied field
	mu0::Float64 = pi*4e-7 			# Vaccum magnetic permeability
	Hext::Vector{Float64} = 1.2.*[	1.0, 
							 		0.0, 
							 		0.0]/mu0
	T::Float64 = 293.0 # 400.0 # K

	# Numerical method settings
	maxAtt::Int32 = 100
	maxDeviation::Float64 = 1e-7

	# Create geometry
	spacing::Float64 = 0.1
	cells = []

	L::Vector{Float64} = [1.0, 1.0, 0.2]
	addCuboid([0.0, 0.0, 0.0], L, cells, true)
	cellLabels = ["Gd"]

	L2::Vector{Float64} = [1.0, 0.1, L[3]]
	addCuboid([0.0, spacing + (L[2] + L2[2])/2, 0.0], L2, cells, true)
	push!(cellLabels, "Fe")

	addCuboid([0.0, -spacing - (L[2] + L2[2])/2, 0.0], L2, cells, true)
	push!(cellLabels, "Fe")

	unifyModel(cells) # In case of intersections, you must unify the model

	# Generate mesh
	mesh::MESH = Mesh([], meshSize, 0.0)

	# Get element tags to then use GMSH 'get' functions
	t_tags, _ = gmsh.model.mesh.getElementsByType(4)

	# Store cell id of each element
	elementID::Vector{Int32} = zeros(mesh.nt)
	for k in 1:mesh.nt
	    # element type , nodes of the element , dimension , id
	    _, _, _, id = gmsh.model.mesh.getElement(t_tags[k])
	    elementID[k] = id
	end

	println("\nNumber of elements ",size(mesh.t,2))
    println("Number of surface elements ",size(mesh.surfaceT,2))

	# Run Gmsh GUI
    if showGmsh
		# gmsh.option.setNumber("Mesh.Clip", 1)
		# gmsh.option.setNumber("General.ClipWholeElements", 1)
		# gmsh.option.setNumber("Mesh.VolumeFaces", 1)
	   	gmsh.fltk.run()
    end
	gmsh.fltk.finalize()

	# Element centroids
	centroids::Matrix{Float64} = zeros(3, mesh.nt)
	for k in 1:mesh.nt
		nds = @view mesh.t[1:4, k]
		centroids[:, k] = mean(mesh.p[1:3, nds], 2)
	end

	# Check cell IDs
        # fig = Figure()
        # ax = Axis3(fig[1, 1], aspect = :data, title="Mesh and Cell Id")
        # scatterPlot = scatter!(ax, 
        #     centroids[1,:], # mesh.InsideElements
        #     centroids[2,:], # mesh.InsideElements
        #     centroids[3,:], # mesh.InsideElements
        #     color=elementID[:], # mesh.InsideElements
        #     colormap=:rainbow,  # :CMRmap :viridis :redsblues
        #     markersize=20)
        # Colorbar(fig[1, 2], scatterPlot, label="Element ID") # Add a colorbar
        
        # # Display the figure (this will open an interactive window)
        # wait(display(fig)); return

	# Getting the volume element from the surface element
	surface2element = surface2volume(mesh)

	# Load Data
	density::Float64 = 7.9 # g/cm3
						
	# Data of magnetic materials
	materialProperties = Dict("Gd" => DATA(),
	                          "Fe" => DATA())
	# Load Gd
	loadMaterial( materialProperties,
                  "../Materials", # Folder with materials
	              "Gd_MFT",  	# Data folder of target material
	              "Gd",        # Material name
	              density,
	              T)

	# Load Iron data
	materialProperties["Fe"].HofM = vec(readdlm("../Materials/Pure_Iron_FEMM/H_Fe_extrap.dat"))  # A/m
	materialProperties["Fe"].B = vec(readdlm("../Materials/Pure_Iron_FEMM/B_Fe_extrap.dat"))     # T

	# Get the permeability and its derivative
	materialPermeability(materialProperties["Fe"])
	
	# BEM matrices
	@time begin 
		C = Cmatrix(mesh)
		D = Dmatrix(mesh)
		B::Matrix{Float64} = zeros(mesh.nv, mesh.ne)
		for s in 1:mesh.ne
		    nds = mesh.surfaceT[1:3, s]
		    B[nds,s] .-= mu0*mesh.AE[s]/3
		end
	end

	# Magnetic permeability
	mu::Vector{Float64} = zeros(mesh.nt) .+ mu0

	# Find the magnetic field
	Hfield::Matrix{Float64} = zeros(3, mesh.nt)
	H::Vector{Float64} = zeros(mesh.nt)
	Hold::Vector{Float64} = zeros(mesh.nt)

	div::Float64 = 1.0
	att::Int32 = 0
	@time while att < maxAtt && div > maxDeviation

		att += 1
		Hold .= H

		# Stiffness matrix
		A = denseStiffnessMatrix(mesh, mu)  # ij
		LHS::Matrix{Float64} = [A B; C D]; # Final FEM-BEM matrix

		# # !! From Bruckner 2012, the load is a surface integral
		# RHS::Vector{Float64} = zeros(mesh.nv + mesh.ne)
		# for s in 1:mesh.ne
		#     nds = @view mesh.surfaceT[1:3, s]
		#     k = surface2element[s]
		#     RHS[nds] .+= (mu[k]-mu0)*dot(mesh.normal[:, s], Hext)*mesh.AE[s]/3
		# end

		# !! For non-linear materials, the RHS is a volume integral
		RHS::Vector{Float64} = zeros(mesh.nv + mesh.ne)
		for k in 1:mesh.nt
		    nds = @view mesh.t[1:4, k]
		    for i = 1:4
		    	_, bi, ci, di = abcd(mesh.p, nds, nds[i])
		    	RHS[nds[i]] += mesh.VE[k]* (mu[k]-mu0)*(Hext[1]*bi + Hext[2]*ci + Hext[3]*di)
		    end
		end

		# Magnetic scalar potential
		u = LHS\RHS

		# Magnetic vector field
		Hfield .= 0.0
		for k in 1:mesh.nt
		    nds = @view mesh.t[:, k]

		    Hfield[:, k] = Hext

		    for j in 1:4
		        _, b, c, d = abcd(mesh.p, nds, nds[j])

		        Hfield[1,k] -= u[nds[j]]*b;
		        Hfield[2,k] -= u[nds[j]]*c;
		        Hfield[3,k] -= u[nds[j]]*d;
		    end
		end

		# Magnetic field intensity
		H .= 0.0
		for k in 1:mesh.nt
			H[k] = norm(Hfield[:, k])
		end

		# Update magnetic permeability
		updatePermeability( mu,
							H, 
							materialProperties, 
							cellLabels, 
							elementID)

		# Check interpolation
		idx = findall(findErr -> !isfinite(findErr), mu)
		if !isempty(idx)
		    println(idx)
		    error("Nans or Infs in mu")
		end

		# Check deviation from previous result
		div = mu0*maximum(abs.(H-Hold))
		println(att, " | mu0 |H(n)-H(n-1)| = ", div)
	end

	# Magnetization
	chi::Vector{Float64} = mu./mu0 .- 1.0

	M::Vector{Float64} = chi.*H

	Mfield::Matrix{Float64} = zeros(3, mesh.nt)
	Mfield[1, :] = chi.*Hfield[1, :]
	Mfield[2, :] = chi.*Hfield[2, :]
	Mfield[3, :] = chi.*Hfield[3, :]

	# Plot vector field
	println("Generating plots...")
	fig = Figure()
	ax = Axis3(fig[1, 1], aspect = :data, title="Spacing: "*string(10*spacing)*" mm")

	ux::Vector{Float64} = Mfield[1,:]
	uy::Vector{Float64} = Mfield[2,:]
	uz::Vector{Float64} = Mfield[3,:]

	for (id, key) in enumerate(cellLabels) # cell number, string with the data label

	    # Find all elements of current cell ID
	    elements = findall(x -> x==id, elementID)
	    
	    # Normalize the magnetization by the max of the cell
	    ux[elements] ./= maximum(M[elements])
	    uy[elements] ./= maximum(M[elements])
	    uz[elements] ./= maximum(M[elements])

	end

	# Plot only the magnetization of a given cell ID
	ID::Int32 = 1
	elements = findall(x -> x==ID, elementID)

	graph = arrows3d!(	  centroids[1, elements]
				  		, centroids[2, elements]
				  		, centroids[3, elements]
				  		, ux[elements]
				  		, uy[elements]
				  		, uz[elements]
                  		, color = M[elements]./(density*1e3)
                  		, lengthscale = 0.1
                  		, colormap = :turbo
                  		# , colorrange = (53.0, 63.8)
                  		)

	Colorbar(    fig[1, 2], graph
			   , label = "M (emu/g)"
	         # , vertical = false
	         )

	wait(display(fig))
end

main(0.0, true)