#=
    Uses the mixed Finite-Element / Boundary Element Method
    
    Solves the time dependent Landau-Lifshitz equation
    Considering a permalloy like in the Fig. 2 of this article
        https://doi.org/10.1109/TMAG.2008.2001666
=#

include("LandauLifshitz.jl")
using GLMakie

function main(meshSize=0.0)
    gmsh.initialize()

    # Constants
    mu0::Float64 = pi*4e-7          # vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
    dt::Float64 = 0.67e-13           # Time step (s)
    totalTime::Float64 = 0.4        # Total time of spin dynamics simulation (ns)
    damp::Float64 = 0.1             # Damping parameter (dimensionless [0,1])
    precession::Float64 = 1.0       # Include precession or not (0 or 1)

    # Dimension of the magnetic material 
    L::Vector{Float64} = [100,100,5] # [512,128,30]
    scl::Float64 = 1e-9                 # scale of the geometry | (m -> nm)
    
    # Conditions
    Ms::Float64   = 860e3               # Magnetic saturation (A/m)
    Aexc::Float64 = 13e-12              # Exchange   (J/m)
    Aan::Float64  = 0                   # Anisotropy (J/m3)
    uan::Vector{Float64}  = [1,0,0]     # easy axis direction
    Hap::Vector{Float64}  = [0,50e3,0] # A/m

    # Convergence criteria | Only used when totalTime != Inf
    maxTorque::Float64 = 0              # Maximum difference between current and previous <M>
    maxAtt::Int32 = 15_000              # Maximum number of iterations in the solver
    
    # -- Create a geometry --

    # Magnetic body
    # addSphere([0,0,0],50)
    addCuboid([0,0,0],L)

    # Generate Mesh
    mesh = Mesh([], meshSize, 0.0, false)

    # Finalize Gmsh and show mesh properties
    # gmsh.fltk.run()
    gmsh.finalize()

    println("Number of elements ", mesh.nt)
    println("Number of nodes ", mesh.nv)
    println("Number of surface elements ", mesh.ns)
    
    # Volume of elements of each mesh node | Needed for the demagnetizing field
    Vn::Vector{Float64} = zeros(mesh.nv)

    # Integral of basis function over the domain | Needed for the exchange field
    nodeVolume::Vector{Float64} = zeros(mesh.nv)
    
    for k in 1:mesh.nt
        Vn[mesh.t[:,k]]         .+= mesh.VE[k]
        nodeVolume[mesh.t[:,k]] .+= mesh.VE[k]/4
    end

    # FEM/BEM matrices
    A = denseStiffnessMatrix(mesh)  # ij
    B = Bmatrix(mesh)        # in
    C = Cmatrix(mesh)        # mj
    D = Dmatrix(mesh)        # mn

    LHS::Matrix{Float64} = [A B; C D]; # Final BEM matrix

    # Initial magnetization field
    m::Matrix{Float64} = zeros(3,mesh.nv)
    m[1,:] .= 1

    # Landau Lifshitz
    m, Heff, M_avg, E_time, torque_time = LandauLifshitz(mesh, m, Ms,
                                                        Hap, Aexc, Aan,
                                                        uan, scl, damp, giro,
                                                        A, LHS, Vn, nodeVolume,
                                                        dt, precession, maxTorque,
                                                        maxAtt, totalTime)

    time::Vector{Float64} = 1e9*dt .* (1:size(M_avg,2))

    fig = Figure()
    ax = Axis(  fig[1,1],
                xlabel = "Time (ns)", 
                ylabel = "<M> (kA/m)",
                title = "Micromagnetic simulation",
                yticks = range(-1500,1500,5))

    scatter!(ax,time,Ms/1000 .*M_avg[1,:], label = "M_x")
    scatter!(ax,time,Ms/1000 .*M_avg[2,:], label = "M_y")
    scatter!(ax,time,Ms/1000 .*M_avg[3,:], label = "M_z")
    axislegend()

    # save("M_time_Sphere.png",fig)
    wait(display(fig))

    # save("M_time_permalloy.png",fig)
end
main()