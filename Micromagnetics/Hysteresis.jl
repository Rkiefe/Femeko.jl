#=
    Solves the Landau-Lifshitz equation for a permalloy, replicating the results from
    OOMMF reported in this article from Oriano et al. 2008
    https://doi.org/10.1109/TMAG.2008.2001666

    The solver is based on Oriano et al. 2008, but with a lot of modifications
=#

include("LL.jl")    # Include Landau-Lifshitz solver
using GLMakie       # Include Makie for plots

function main(meshSize::Float64=0.0, localSize::Float64=0.0, showGmsh::Bool=true)

    mu0::Float64 = pi*4e-7          # Vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
    
    # A struct holding all the info needed for the LL solver
    ll = LL()

    # Update LL parameters
    ll.Hext = [0.0, 0.0, 0.0] # Applied field (T)
    ll.Ms = mu0 * 800e3 # Mag. saturation (T)
    ll.scale = 1e-9     # scale of the geometry | nm: 1e-9 m

    ll.Aexc = 13e-12    # Exchange   (J/m)

    ll.Aan = 0.0        # Anisotropy constant J/m3
    ll.uan = [1,0,0]    # Easy axis
    
    ll.timeStep = 0.1  # Time step (normalized by the gyromagnetic ratio)
    # ll.totalTime = 70.35    # Stop when time > total time (normalized by the gyromagnetic ratio)
    ll.maxTorque = 1e-5   # If 'totalTime' is not provided, it minimizes the Energy/M state

    # ll.alfa = 0.1   # damping
    
    # Create a 3D Model
    gmsh.initialize()
    cells = []
    addCuboid([0,0,0], [512, 128, 30], cells)
    box = addSphere([0,0,0], 2500.0) # Add a bounding shell

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

    # Random initial magnetization state
    M::Matrix{Float64} = zeros(3, ll.mesh.nv)
    theta::Vector{Float64} = 2*pi.*rand(ll.mesh.nv)
    phi::Vector{Float64} = pi.*rand(ll.mesh.nv)
    for i in ll.mesh.InsideNodes
        M[1, i] = ll.Ms*sin(phi[i])*cos(theta[i])
        M[2, i] = ll.Ms*sin(phi[i])*sin(theta[i])
        M[3, i] = ll.Ms*cos(phi[i])
    end

    # Run
    Hspan = [range(0.0, 0.1, 101); range(0.1, -0.1, 201); range(-0.1, 0.1, 201)]
    M_H::Vector{Float64} = zeros(length(Hspan))
    for (i, h) in enumerate(Hspan)
        
        # Set the applied field
        ll.Hext[1] = h

        # Run
        M, _, Mx, _, _, _ = ll.run(ll, M, A, AEXC, b, c, d, Volumes, false)
        M_H[i] = Mx[end] # Save the last value
    
    end # Hysteresis loop
    
    # Plot the M(t)
    fig = Figure()
    ax = Axis(fig[1,1], title="<M> (emu/g)", xlabel="Time (s)", ylabel="M (kA/m)")
    scatter!(ax, Hspan, M_H)

    # save("M_H.png", fig)
    wait(display(fig))
end

meshSize = 2000.0
localSize = 20.0
showGmsh = true
# main(meshSize, localSize, showGmsh)