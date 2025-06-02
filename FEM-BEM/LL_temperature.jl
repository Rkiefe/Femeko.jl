#=
    Uses the mixed Finite-Element / Boundary Element Method

    Includes a thermal field to the Landau-Lifshitz equation
=#

# For plots
using GLMakie

include("../gmsh_wrapper.jl")
include("LandauLifshitz.jl")

function main()
    meshSize::Float64 = 0

    # Constants
    mu0::Float64 = pi*4e-7          # vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
    dt::Float64 = 0.67e-13          # Time step (s)
    totalTime::Float64 = Inf        # Total time of spin dynamics simulation (ns)
    damp::Float64 = 0.1             # Damping parameter (dimensionless [0,1])
    precession::Float64 = 1.0       # Include precession or not (0 or 1)

    # Temperature
    T::Float64  = 0.0 # K

    # Dimension of the magnetic material 
    L::Vector{Float64} = [100,100,5] # [512,128,30]
    scl::Float64 = 1e-9                 # scale of the geometry | (m -> nm)
    
    # Conditions
    Ms::Float64   = 860e3               # Magnetic saturation (A/m)
    Aexc::Float64 = 13e-12              # Exchange   (J/m)
    Aan::Float64  = 0                   # Anisotropy (J/m3)
    uan::Vector{Float64}  = [1,0,0]     # easy axis direction
    Hap::Vector{Float64}  = [0,50e3,0]  # A/m

    # Convergence criteria | Only used when totalTime != Inf
    maxTorque::Float64 = 1e-10              # Maximum difference between current and previous <M>
    maxAtt::Int32 = 5_000              # Maximum number of iterations in the solver
    
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

    println("Number of nodes ",mesh.nv)
    println("Number of surface elements ",mesh.ne)
    println("Number of elements ",mesh.nt)
    println("Number of Inside elements ",mesh.nInside)
    println("Number of Inside nodes ",mesh.nInsideNodes)
    # viewMesh(mesh)
    # return

    # Pre-calculate the area of each surface triangle
    areaT::Vector{Float64} = zeros(mesh.ne)
    for s in 1:mesh.ne
        nds = mesh.surfaceT[1:3,s]
        areaT[s] = areaTriangle(mesh.p[1,nds],mesh.p[2,nds],mesh.p[3,nds])
    end

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
    B = Bmatrix(mesh, areaT)        # in
    C = Cmatrix(mesh, areaT)        # mj
    D = Dmatrix(mesh, areaT)        # mn

    LHS::Matrix{Float64} = [-A B; C D]; # Final BEM matrix

    # Initial magnetization field
    m::Matrix{Float64} = zeros(3,mesh.nv)
    # m[1,:] .= 1
    begin # Random initial magnetization
        theta::Vector{Float64} = 2*pi.*rand(mesh.nv)
        phi::Vector{Float64} = pi.*rand(mesh.nv)
        for i in 1:mesh.nv
            m[:,i] = [sin(theta[i])*cos(phi[i]), sin(theta[i])*sin(phi[i]), cos(theta[i])]
            m[:,i] = m[:,i]./norm(m[:,i])
        end
    end # Random initial magnetization

    # M(T)
    mOld::Matrix{Float64} = deepcopy(m)
    Tspan::Vector{Float64} = range(0.0,0.15,50)
    M_T::Matrix{Float64} = zeros(3,length(Tspan))
    for iT in 1:length(Tspan)
        T = Tspan[iT]

        # Reset m
        m = deepcopy(mOld)

        # Landau Lifshitz
        _, _, M_avg = LandauLifshitz(mesh, m, Ms,
                                    Hap, Aexc, Aan,
                                    uan, scl, damp, giro,
                                    A, LHS, Vn, nodeVolume, areaT,
                                    dt, precession, maxTorque,
                                    maxAtt, totalTime, T)

        M_T[:,iT] = M_avg[:,end]
    end

    # time::Vector{Float64} = (1e9*dt) * (1:size(M_avg,2))

    fig = Figure()
    ax = Axis(  fig[1,1],
                xlabel = "T (K)", 
                ylabel = "<M> (kA/m)",
                title = "Thermal noise")

    scatter!(ax, Tspan, M_T[1,:], label = "M_x")
    scatter!(ax, Tspan, M_T[2,:], label = "M_y")
    scatter!(ax, Tspan, M_T[3,:], label = "M_z")
    axislegend()

    save("M_T.png",fig)
    wait(display(fig))
end

main()