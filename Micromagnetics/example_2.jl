#=
    Minimizes the energy by steepest descent and by the 
    Landau-Lifshitz equation.

    This example is the same as 'example_1' but without time dynamics.
    The system is minimized directly towards a stable magnetization and low energy
    resulting in the same <M> at the end, as expected.
=#


# For plots
using GLMakie

include("../src/gmsh_wrapper.jl")
include("SteepestDescent.jl")

function main()
    meshSize::Float64 = 500
    localSize::Float64 = 5

    # Constants
    mu0::Float64 = pi*4e-7          # vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
    dt::Float64 = 0.1e-12             # Time step (s)
    totalTime::Float64 = 0.1        # Total time of spin dynamics simulation (ns)
    damp::Float64 = 0.0             # Damping parameter (dimensionless [0,1])
    precession::Float64 = 1.0       # Include precession or not

    # Dimension of the magnetic material 
    L::Vector{Float64} = [100,100,5]   # (rectangle)
    scl::Float64 = 1e-9                 # scale of the geometry | (m -> nm)
    
    # Conditions
    Ms::Float64   = 860e3              # Magnetic saturation (A/m)
    Aexc::Float64 = 13e-12                 # Exchange   (J/m)
    Aan::Float64  = 0.0               # Anisotropy (J/m3)
    uan::Vector{Float64}  = [1,0,0]     # easy axis direction
    Hap::Vector{Float64}  = [0,50e3,0] # A/m

    # Convergence criteria | Only used when totalTime != Inf
    maxTorque::Float64 = 1e-6           # Maximum difference between current and previous <M>
    maxAtt::Int32 = 1_000              # Maximum number of iterations in the solver
    
    # -- Create a geometry --
        gmsh.initialize()

        # >> Model
        # Create an empty container
        container = addSphere([0,0,0],5*maximum(L))

        cells = [] # List of cells inside the container

        # Get how many surfaces compose the bounding shell
        temp = gmsh.model.getEntities(2)                # Get all surfaces of current model
        bounding_shell_n_surfaces = 1:length(temp)      # Get the number of surfaces in the bounding shell

        # Add another object inside the container
        addCuboid([0,0,0],L,cells,true)

        # Fragment to make a unified geometry
        _, fragments = gmsh.model.occ.fragment([(3, container)], cells)
        gmsh.model.occ.synchronize()

        # Update container volume ID
        container = fragments[1][1][2]

        # Generate Mesh
        mesh = Mesh(cells,meshSize,localSize,false)
        
        # Get bounding shell surface id
        mesh.shell_id = gmsh.model.getAdjacencies(3, container)[2]

        # Must remove the surface Id of the interior surfaces
        mesh.shell_id = mesh.shell_id[bounding_shell_n_surfaces] # All other, are interior surfaces

        # Finalize Gmsh and show mesh properties
        # gmsh.fltk.run()
        gmsh.finalize()
    
    # -----------------------

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    # println("Bounding shell: ",mesh.shell_id)
    # return

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
    # begin
    #     theta::Vector{Float64} = 2*pi.*rand(mesh.nInsideNodes)
    #     phi::Vector{Float64} = pi.*rand(mesh.nInsideNodes)
    
    #     for i in 1:mesh.nInsideNodes
    #         nd = mesh.InsideNodes[i]
    #         m[:,nd] =   [
    #                         sin(theta[i])*cos(phi[i])
    #                         sin(theta[i])*sin(phi[i])
    #                         cos(theta[i])
    #                     ]

    #         m[:,nd] = m[:,nd]./norm(m[:,nd])
    #     end
    # end # Random initial magnetization

    # Dirichlet boundary condition
    fixed::Vector{Int32} = findNodes(mesh,"face",mesh.shell_id)
    free::Vector{Int32} = setdiff(1:mesh.nv,fixed)

    # Stiffness matrix | Demagnetizing field
    AD = stiffnessMatrix(mesh)

    # Stiffness matrix | Exchange field
    A = spzeros(mesh.nv,mesh.nv)

    begin
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

    end # End of Stiffness Matrix of Exchange field

    # Landau Lifshitz
    # m, Heff::Matrix{Float64},
    # M_avg::Matrix{Float64},
    # E_time::Vector{Float64}, 
    # torque_time::Vector{Float64} = LandauLifshitz(  
    #                                                 mesh, m, Ms,
    #                                                 Hap, Aexc, Aan,
    #                                                 uan, scl, damp, giro,
    #                                                 A, AD, Vn, nodeVolume,
    #                                                 free, fixed,
    #                                                 dt, precession, maxTorque,
    #                                                 maxAtt, totalTime
    #                                              )

    Heff::Matrix{Float64} = Matrix{Float64}(undef,0,0)

    m, Heff,
    M_avg::Matrix{Float64},
    E_time::Vector{Float64}, 
    torque_time::Vector{Float64} = @time SteepestDescent(
                                                    mesh, m, Ms, Heff,
                                                    Hap, Aexc, Aan,
                                                    uan, scl,
                                                    A, AD, Vn, nodeVolume,
                                                    free, fixed,
                                                    maxTorque, maxAtt
                                                  )

    step::Vector{Int32} = 1:size(M_avg,2)

    fig = Figure()
    ax = Axis(  fig[1,1],
                xlabel = "Iteration step", 
                ylabel = "Energy density",
             )
    scatter!(ax,step,E_time)

    ax = Axis(  fig[1,2],
                xlabel = "Iteration step", 
                ylabel = "Log of <|dm/dt|>",
             )
    scatter!(ax,step,log.(10,torque_time))

    ax = Axis(  fig[1,3],
                xlabel = "Iteration step", 
                ylabel = "<M> (kA/m)",
             )
    scatter!(ax,step,Ms/1000 .*M_avg[1,:], label = "M_x")
    scatter!(ax,step,Ms/1000 .*M_avg[2,:], label = "M_y")
    scatter!(ax,step,Ms/1000 .*M_avg[3,:], label = "M_z")
    axislegend()

    # save("M_time_Sphere.png",fig)
    wait(display(fig))
end

main()