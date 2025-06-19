# For plots
# using CairoMakie

include("../../src/gmsh_wrapper.jl")
include("SteepestDescent.jl")

function main()
    meshSize::Float64 = 0

    # Constants
    mu0::Float64 = pi*4e-7          # vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
    dt::Float64 = 1e-12             # Time step (s)
    totalTime::Float64 = Inf        # Total time of spin dynamics simulation (ns)
    damp::Float64 = 1.0             # Damping parameter (dimensionless [0,1])
    precession::Float64 = 0.0       # Include precession or not (0 or 1)

    # Dimension of the magnetic material 
    L::Vector{Float64} = [512,128,30]
    scl::Float64 = 1e-9                 # scale of the geometry | (m -> nm)
    
    # Conditions
    Ms::Float64   = 800e3               # Magnetic saturation (A/m)
    Aexc::Float64 = 13e-12              # Exchange   (J/m)
    Aan::Float64  = 0                   # Anisotropy (J/m3)
    uan::Vector{Float64}  = [1,0,0]     # easy axis direction
    Hap::Vector{Float64}  = [0.1/mu0,0,0]     # A/m

    Hext::Vector{Float64} = vcat(0:1e-3:0.1,0.1:-1e-3:-0.1,-0.1:1e-3:0.1)./mu0
    Hap[1] = Hext[1]

    # Convergence criteria | Only used when totalTime != Inf
    maxTorque::Float64 = 1e-14          # Maximum difference between current and previous <M>
    maxAtt::Int32 = 5_000               # Maximum number of iterations in the solver
    
    # -- Create a geometry --
    gmsh.initialize()

    # Magnetic body
    # addSphere([0,0,0],50)
    addCuboid([0,0,0],L)

    # Generate Mesh
    mesh = Mesh([],meshSize,0,false)

    # Finalize Gmsh and show mesh properties
    # gmsh.fltk.run()
    gmsh.finalize()
        
    # -----------------------

    println("Number of elements ",mesh.nt)
    println("Number of nodes ",mesh.nv)
    println("Number of surface elements ",mesh.ne)
    println("Number of Inside elements ",mesh.nInside)
    println("Number of Inside nodes ",mesh.nInsideNodes)
    # viewMesh(mesh)
    # return
    
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

    LHS::Matrix{Float64} = [-A B; C D]; # Final BEM matrix

    # Initial magnetization field
    m::Matrix{Float64} = zeros(3,mesh.nv)
    m[1,:] .= 1
    theta::Vector{Float64} = 2*pi.*rand(mesh.nv)
    phi::Vector{Float64} = pi.*rand(mesh.nv)
    for i in 1:mesh.nv
        m[:,i] = [sin(theta[i])*cos(phi[i]), sin(theta[i])*sin(phi[i]), cos(theta[i])]
        m[:,i] = m[:,i]./norm(m[:,i])
    end

    Heff::Matrix{Float64} = Matrix{Float64}(undef,0,0)

    # Landau Lifshitz
    m, Heff, M_avg, E_time, torque_time = LandauLifshitz(mesh, m, Ms,
                                                        Hap, Aexc, Aan,
                                                        uan, scl, damp, giro,
                                                        A, LHS, Vn, nodeVolume,
                                                        dt, precession, maxTorque,
                                                        maxAtt, totalTime)

    # M_H::Matrix{Float64} = zeros(3,length(Hext))
    # for iH in 1:length(Hext)
        # Hap[1] = Hext[iH]
        m, Heff, M_avg = SteepestDescent(mesh, m, Ms, Heff,
                                        Hap, Aexc, Aan,
                                        uan, scl,
                                        A, LHS, Vn, nodeVolume,
                                        maxTorque, maxAtt)
        
        # M_H[:,iH] = M_avg[:,end]    
    # end

    # fig = Figure()
    # ax = Axis(  fig[1,1],
    #             xlabel = "B applied",
    #             ylabel = "<M>")

    # scatter!(ax, mu0.*Hext, M_H[1,:])

    # save("M_H_"*string(mesh.nv)*".png",fig)
    # display(fig)
end

main()
