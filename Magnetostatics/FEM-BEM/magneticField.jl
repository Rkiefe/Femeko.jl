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
# using IterativeSolvers

meshSize = 1.0
localSize = 0.0
showGmsh = false

function main(meshSize=0.0, localSize=0.0, showGmsh=false, verbose=true)
	gmsh.initialize()

	# Applied field | A/m
	mu0::Float64 = pi*4e-7 			# Vaccum magnetic permeability

	Hext::Vector{Float64} = [1.2, 
							 0.0, 
							 0.0]/mu0

	# Temperature of the magnetic material | K
	T::Float64 = 293.0

	# Numerical method settings
	maxAtt::Int32 = 100
	maxDeviation::Float64 = 1e-7

	# Create geometry
	cells = []
	addCuboid([0,0,0], [1.0, 1.0, 1.0], cells) 	
	# addSphere([0,0,0], 0.5, cells)

	# Generate mesh
	localSize > 0.0 ? extendLocalRefinement() : nothing # keep local refinement on the boundary
	mesh::MESH = Mesh(cells, meshSize, localSize)

	println("\nNumber of elements ", mesh.nt)
    println("Number of surface elements ", mesh.ns)

	# Run Gmsh GUI
    if showGmsh
		gmsh.option.setNumber("Mesh.Clip", 1)
		gmsh.option.setNumber("General.ClipWholeElements", 1)
		gmsh.option.setNumber("Mesh.VolumeFaces", 1)
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
	println("Building the BEM matrices and the FEM-BEM coupling")
	@time begin 
		C = Cmatrix(mesh)
		D = Dmatrix(mesh)
		B::Matrix{Float64} = zeros(mesh.nv, mesh.ns)
		for s in 1:mesh.ns
		    nds = mesh.surfaceT[1:3, s]
		    B[nds,s] .-= mu0*mesh.AE[s]/3
		end
	end

	println("Building the element-wise stiffness matrix")
	Ak = @time localStiffnessMatrix(mesh)

	# Find the magnetic field
	Hfield::Matrix{Float64} = zeros(3, mesh.nt)
	H::Vector{Float64} = zeros(mesh.nt)
	Hold::Vector{Float64} = zeros(mesh.nt)

	div::Float64 = 1.0
	att::Int32 = 0
	println("Running the Fixed-Point iteration method")
	@time while att < maxAtt && div > maxDeviation

		att += 1
		Hold .= H

		# Update the stiffness matrix
	    A::Matrix{Float64} = zeros(mesh.nv, mesh.nv)
	    for k in 1:mesh.nt
	    	nds = @view mesh.t[:, k]
	    	n = 0
	    	for i in 1:4
	    		for j in 1:4
	    			n += 1
	    			A[nds[i], nds[j]] += Ak[n, k]*mu[k]
	    		end
    		end  
	    end

		LHS = [A B; C D]; # Final FEM-BEM matrix

		# # !! From Bruckner 2012, the load is a surface integral
		# RHS::Vector{Float64} = zeros(mesh.nv + mesh.ns)
		# for s in 1:mesh.ns
		#     nds = @view mesh.surfaceT[1:3, s]
		#     k = surface2element[s]
		#     RHS[nds] .+= (mu[k]-mu0)*dot(mesh.normal[:, s], Hext)*mesh.AE[s]/3
		# end

		# !! For non-linear materials, the RHS is a volume integral
		RHS::Vector{Float64} = zeros(mesh.nv + mesh.ns)
		for k in 1:mesh.nt
		    nds = @view mesh.t[1:4, k]
		    for i = 1:4
		    	_, bi, ci, di = abcd(mesh.p, nds, nds[i])
		    	RHS[nds[i]] += mesh.VE[k]* (mu[k]-mu0)*(Hext[1]*bi + Hext[2]*ci + Hext[3]*di)
		    end
		end

		# Magnetic scalar potential
		u = LHS\RHS
		residue = norm(RHS - LHS*u)

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
        verbose ? println(att, " | mu0 |H(n)-H(n-1)| = ", div, " , |y-Ax| = ", residue) : nothing

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

	graph = arrows3d!(	  centroids[1,:]
				  		, centroids[2,:]
				  		, centroids[3,:]
				  		, Mfield[1,:]
				  		, Mfield[2,:]
				  		, Mfield[3,:]
                  		, color = M./(density*1e3)
                  		, lengthscale = 0.1/maximum(M)
                  		, colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                  		)

	Colorbar(fig[1, 2], graph, label="M (emu/g)")
	wait(display(fig))
end

main(meshSize, localSize, showGmsh)