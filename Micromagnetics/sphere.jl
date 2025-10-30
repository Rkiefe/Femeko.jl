#=
    Solves the Landau-Lifshitz for a high stress-test example compared
    against an analytic solution
        A sphere without an exchange field and no damping
        -> a sinusoidal behavior of the <M> over time.

    The result is heavily dependent on the mesh quality
    and the bounding shell size and the time step
=#

include("../src/Femeko.jl")
include("LandauLifshitz.jl")

# For plots
using GLMakie

function main(showGmsh=true)
    meshSize::Float64 = 2500
    localSize::Float64 = 1

    # Constants
    mu0::Float64 = pi*4e-7          # vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
    dt::Float64 = 1e-12             # Time step (s)
    totalTime::Float64 = 0.1        # Total time of spin dynamics simulation (ns)
    damp::Float64 = 0.0             # Damping parameter (dimensionless [0,1])
    precession::Float64 = 1.0       # Include precession or not (0 or 1)

    # Dimension of the magnetic material 
    L::Vector{Float64} = [512,128,30]   # (rectangle)
    scl::Float64 = 1e-9                 # scale of the geometry | (m -> nm)
    
    # Conditions
    Ms::Float64   = 1400e3              # Magnetic saturation (A/m)
    Aexc::Float64 = 0                   # Exchange   (J/m)
    Aan::Float64  = 500e3               # Anisotropy (J/m3)
    uan::Vector{Float64}  = [1,0,0]     # easy axis direction
    Hap::Vector{Float64}  = [0,400e3,0] # A/m

    # Convergence criteria | Only used when totalTime != Inf
    maxTorque::Float64 = 1e-6           # Maximum difference between current and previous <M>
    maxAtt::Int32 = 15_000              # Maximum number of iterations in the solver
    
    # Create a geometry
    gmsh.initialize()
    cells = [] # Array of volume cell IDs

    # Model
    addSphere([0,0,0], 50.0, cells)
    box = addSphere([0,0,0], 5.0*50.0) # Create bounding shell

    # Unify the volumes for a single geometry and get the bounding shell
    shell_id = unifyModel(cells, box)

    # Generate Mesh
    mesh = Mesh(cells, meshSize, localSize, false)
    mesh.shell_id = shell_id

    println("Number of elements ", mesh.nt)
    println("Number of nodes ", mesh.nv)
    println("Number of surface elements ", mesh.ne)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Bounding shell: ", mesh.shell_id)
    
    # Finalize Gmsh and print mesh properties
    if showGmsh
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.finalize()
    # viewMesh(mesh)
    
    # Volume of elements of each mesh node | Needed for the demagnetizing field
    Vn::Vector{Float64} = zeros(mesh.nv)

    # Integral of basis function over the domain | Needed for the exchange field
    nodeVolume::Vector{Float64} = zeros(mesh.nv)
    
    for ik in 1:mesh.nInside
        k = mesh.InsideElements[ik]
        Vn[mesh.t[:,k]]         .+= mesh.VE[k]
        nodeVolume[mesh.t[:,k]] .+= mesh.VE[k]/4
    end
    Vn = Vn[mesh.InsideNodes]
    nodeVolume = nodeVolume[mesh.InsideNodes]

    # Magnetization field
    m::Matrix{Float64} = zeros(3,mesh.nv)
    m[1,mesh.InsideNodes] .= 1

    # Dirichlet boundary condition
    fixed::Vector{Int32} = findNodes(mesh,"face",mesh.shell_id)
    free::Vector{Int32} = setdiff(1:mesh.nv,fixed)

    # Stiffness matrix | Demagnetizing field
    AD = stiffnessMatrix(mesh)

    # Stiffness matrix | Exchange field
    A = spzeros(mesh.nv,mesh.nv)
    Ak::Matrix{Float64} = zeros(4*4,mesh.nt)
    b::Vector{Float64} = zeros(4)
    c::Vector{Float64} = zeros(4)
    d::Vector{Float64} = zeros(4)
    aux::Matrix{Float64} = zeros(4,4)
    
    for ik in 1:mesh.nInside
        k = mesh.InsideElements[ik]
        for i in 1:4
            _,b[i],c[i],d[i] = abcd(mesh.p,mesh.t[:,k],mesh.t[i,k])
        end
        aux = mesh.VE[k]*(b*b' + c*c' + d*d')
        Ak[:,k] = aux[:] # vec(aux)
    end

    # Update sparse global matrix
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            A += sparse(Int.(mesh.t[i,:]),Int.(mesh.t[j,:]),Ak[n,:],mesh.nv,mesh.nv)
        end
    end
    
    # Remove all exterior nodes
    A = A[mesh.InsideNodes,mesh.InsideNodes]

    # Landau Lifshitz
    m, Heff::Matrix{Float64},
    M_avg::Matrix{Float64},
    E_time::Vector{Float64}, 
    torque_time::Vector{Float64} = @time LandauLifshitz(  
                                                    mesh, m, Ms,
                                                    Hap, Aexc, Aan,
                                                    uan, scl, damp, giro,
                                                    A, AD, Vn, nodeVolume,
                                                    free, fixed,
                                                    dt, precession, maxTorque,
                                                    maxAtt, totalTime
                                                 )

    time::Vector{Float64} = 1e9*dt .* (1:size(M_avg,2))

    fig = Figure()
    ax = Axis(  fig[1,1],
                xlabel = "Time (ns)", 
                ylabel = "<M> (kA/m)",
                title = "Micromagnetic simulation",
                yticks = range(-1500,1500,5))

    scatter!(ax,time,Ms/1000 .*M_avg[1,:], label = "M_x")
    # scatter!(ax,time,Ms/1000 .*M_avg[2,:], label = "M_y")
    # scatter!(ax,time,Ms/1000 .*M_avg[3,:], label = "M_z")
    axislegend()

    # save("M_time_Sphere.png",fig)
    wait(display(fig))
end

main()