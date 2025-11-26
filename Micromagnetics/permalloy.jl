#=
	Solves the Landau-Lifshitz equation for a permalloy, replicating the results from
	OOMMF reported in this article from Oriano et al. 2008
	https://doi.org/10.1109/TMAG.2008.2001666

	The solver is based on Oriano et al. 2008, but with a lot of modifications
=#

include("LL.jl") 	# Include Landau-Lifshitz solver
using GLMakie 		# Include Makie for plots

function main(meshSize::Float64=0.0, localSize::Float64=0.0, showGmsh::Bool=true)

	mu0::Float64 = pi*4e-7 			# Vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
	
	# A struct holding all the info needed for the LL solver
	ll = LL()

	# Update LL parameters
	ll.Hext = [0.0, mu0* 50e3, 0.0] # Applied field (T)
	ll.Ms = mu0 * 860e3	# Mag. saturation (T)
	ll.scale = 1e-9		# scale of the geometry | nm: 1e-9 m

	ll.Aexc = 12e-13	# Exchange   (J/m)

	ll.Aan = 0.0 		# Anisotropy constant J/m3
	ll.uan = [1,0,0] 	# Easy axis
	
	ll.timeStep = 0.01 	# Time step (normalized by the gyromagnetic ratio)
	ll.totalTime = 70.35 	# Stop when time > total time (normalized by the gyromagnetic ratio)

	ll.alfa = 0.1 	# damping
	
	# Create a 3D Model
	gmsh.initialize()
	cells = []
	addCuboid([0,0,0], [100, 100, 5], cells)
	box = addSphere([0,0,0], 500.0) # Add a bounding shell

	unifyModel(cells, box)

	# Create a mesh
	ll.mesh = Mesh(cells, meshSize, localSize)

	# Print number of elements and nodes
	println("\nNumber of elements: ", ll.mesh.nt)
	println("Number of internal elements: ", ll.mesh.nInside)
	println("Number of internal nodes: ", ll.mesh.nInsideNodes, "\n")

	if showGmsh
		gmsh.option.setNumber("Mesh.Clip", 1)
		gmsh.option.setNumber("Mesh.VolumeFaces", 1)
		gmsh.option.setNumber("General.ClipWholeElements", 1)
		gmsh.fltk.run()
	end
	gmsh.finalize()

	# FEM
	print("Calculating the Lagrange shape elements... ")
	b::Matrix{Float64} = zeros(4, ll.mesh.nt)
	c::Matrix{Float64} = zeros(4, ll.mesh.nt)
	d::Matrix{Float64} = zeros(4, ll.mesh.nt)
	for k in 1:ll.mesh.nt
		nds = @view ll.mesh.t[:, k]
		for i in 1:4
			_, b[i, k], c[i, k], d[i, k] = abcd(ll.mesh.p, nds, nds[i])
		end
	end
	println("Done.")

	# Stiffness matrix | Exchange field 
    AEXC = spzeros(ll.mesh.nv, ll.mesh.nv)

	Ak::Matrix{Float64} = zeros(4*4, ll.mesh.nt)
	for k in ll.mesh.InsideElements
	    nds = @view ll.mesh.t[:,k]
	    Ak[:, k] = vec( ll.mesh.VE[k]*( b[:, k]*b[:, k]' 
										+ c[:, k]*c[:, k]'
										+ d[:, k]*d[:, k]' ) )
	end
    
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            AEXC += sparse(ll.mesh.t[i,:], ll.mesh.t[j,:], Ak[n,:], ll.mesh.nv, ll.mesh.nv)
        end
    end # Stiffness matrix for Exchange Field

	# Node volumes
	Volumes::Vector{Float64} = zeros(ll.mesh.nv)
	for k in ll.mesh.InsideElements
		nds = @view ll.mesh.t[:, k]
		Volumes[nds] .+= ll.mesh.VE[k]
	end

	# Global stiffness matrix
	A = stiffnessMatrix(ll.mesh)

	# Initial magnetization state
	M::Matrix{Float64} = zeros(3, ll.mesh.nv)

	# Uniform initial magnetization
	M[1, ll.mesh.InsideNodes] .= ll.Ms

	# Run
	M, Mnorm, Mx, My, Mz, _ = ll.run(ll, M, A, AEXC, b, c, d, Volumes, true)

	# Plot the M(t)
	fig = Figure()
	ax = Axis(fig[1,1], title="<M> (emu/g)", xlabel="Time (ns)", ylabel="M (kA/m)")
	scatter!(ax, (ll.timeStep*1e9/giro)*(0:ll.nSteps-1), Mx./(mu0*1e3), label="Mx")
	scatter!(ax, (ll.timeStep*1e9/giro)*(0:ll.nSteps-1), My./(mu0*1e3), label="My")
	scatter!(ax, (ll.timeStep*1e9/giro)*(0:ll.nSteps-1), Mz./(mu0*1e3), label="Mz")
    axislegend(position=:rb)

    save("M_time_permalloy.png", fig)
	wait(display(fig))

end

meshSize = 200.0
localSize = 20.0
showGmsh = false
main(meshSize, localSize, showGmsh)