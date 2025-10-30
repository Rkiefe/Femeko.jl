#=
    	Internal magnetic field using FEM/BEM for non-linear magnetic materials
    
    This implementation is based on https://doi.org/10.1016/j.jmmm.2012.01.016
    but adapted for non-linear magnetic materials and any magnetic field souce (non conservative)        	

	The non-linearity is handled by the Fixed-Point iteration method
=#

include("../../src/Femeko.jl")
include("../../src/BEM.jl")
include("../../src/magneticProperties.jl")

using GLMakie

function main(meshSize=0.0, showGmsh=false)
	gmsh.initialize()

	# Applied field | A/m
	mu0::Float64 = pi*4e-7 			# Vaccum magnetic permeability

	Hext::Vector{Float64} = [1.0, 
							 0.0, 
							 0.0]*1.2/mu0

	# Temperature of the magnetic material | K
	T::Float64 = 293.0

	# Numerical method settings
	maxAtt::Int32 = 100
	maxDeviation::Float64 = 1e-7

	# Create geometry
	addCuboid([0,0,0], [1.0, 1.0, 1.0]) 	
	# addSphere([0,0,0], 0.5)

	# Generate mesh
	mesh::MESH = Mesh([], meshSize, 0.0)

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

	# Each surface element has a corresponding volume element. Not needed anymore
	# surface2element = surface2volume(mesh) # maps the surface element ID to a volume element ID

	# Load Data
	density::Float64 = 7.9 # g/cm3
						
	data = DATA()
	loadMaterial( data,
               	  "../Materials", # Folder with materials
               	  "Gd_MFT",    # Data folder of target material
               	  "Gd",        # Material name
               	  density,
               	  T)

	spl = Spline1D(data.HofM, data.mu)
	
	# Magnetic relative permeability
	mu::Vector{Float64} = zeros(mesh.nt) .+ mu0

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
		# Hfield = zeros(3, mesh.nt) .+ Hext
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
		mu = spl(H)

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

	# Average magnetization
	M_avg::Float64 = 0.0
	volume::Float64 = 0.0
	for k in 1:mesh.nt
		volume 	+= mesh.VE[k]
		M_avg 	+= mesh.VE[k]*M[k]
	end
	M_avg /= volume

	println("<M> (emu/g) = ", M_avg/(density*1e3))

	# Plot vector field
	println("Generating plots...")
	fig = Figure()
	ax = Axis3(fig[1, 1], aspect = :data)

	u::Vector{Float64} = Mfield[1,:]./M
	v::Vector{Float64} = Mfield[2,:]./M
	w::Vector{Float64} = Mfield[3,:]./M

	graph = arrows3d!(	  centroids[1,:]
				  		, centroids[2,:]
				  		, centroids[3,:]
				  		, u
				  		, v
				  		, w
                  		, color = M./(density*1e3)
                  		, lengthscale = 0.1
                  		, colormap = :turbo
                  		)

	Colorbar(fig[1, 2], graph, label = "M (emu/g)"
	         #, vertical = false
	         )

	wait(display(fig))
end

main(0.0, false)