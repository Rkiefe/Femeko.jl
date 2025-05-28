#=
    Uses the mixed Finite-Element / Boundary Element Method

    Solves the Landau-Lifshitz for a high stress-test example compared
    against an analytic solution
        A sphere without an exchange field and no damping
        -> a sinusoidal behavior of the <M> over time.
=#

# For plots
using GLMakie

include("../gmsh_wrapper.jl")
include("../FEMjl.jl")
# include("LandauLifshitz.jl")

function demagnetizingField(mesh::MESH,m::Matrix{Float64},areaT::Vector{Float64},LHS::Matrix{Float64})
    
    RHS::Vector{Float64} = zeros(mesh.nv + mesh.ne)
    for s in 1:mesh.ne
        nds = @view mesh.surfaceT[1:3,s]
        RHS[nds] .+= dot(mesh.normal[:,s],mean(m[:,nds],2))*areaT[s]/3
    end

    # Magnetic scalar potential
    u = LHS\-RHS

    # Demag field | Elements
    Hdk::Matrix{Float64} = zeros(3,mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:,k] # all nodes of that element

        # Sum the contributions
        for j in 1:4
            _, b::Float64, c::Float64, d::Float64 = abcd(mesh.p,nds,nds[j])

            Hdk[1,k] = Hdk[1,k] - u[nds[j]]*b;
            Hdk[2,k] = Hdk[2,k] - u[nds[j]]*c;
            Hdk[3,k] = Hdk[3,k] - u[nds[j]]*d;
        end
    end

    # Demag field | Nodes
    Hd::Matrix{Float64} = zeros(3,mesh.nv)
    for k in 1:mesh.nt
        nds = mesh.t[:,k]
        Hd[:,nds] .+= mesh.VE[k].*Hdk[:,k]
    end

    return Hd
end


function Bmatrix(mesh::MESH, areaT::Vector{Float64})
    B::Matrix{Float64} = zeros(mesh.nv,mesh.ne)
    for s in 1:mesh.ne
        nds = mesh.surfaceT[1:3,s]
        B[nds,s] .+= areaT[s]/3
    end
    return B
end # BEM matrix Min

function Cmatrix(mesh::MESH, areaT::Vector{Float64})
    # Basis function value on the quadrature points
    phi = [1.0000         0         0
                 0    1.0000         0
                 0         0    1.0000
            0.5000    0.5000         0
                 0    0.5000    0.5000
            0.5000         0    0.5000
            0.7500    0.2500         0
            0.5000    0.2500    0.2500
            0.7500         0    0.2500
            0.2500    0.7500         0
                 0    0.7500    0.2500
            0.2500    0.5000    0.2500
                 0    0.2500    0.7500
            0.2500         0    0.7500
            0.2500    0.2500    0.5000]';

    # Mmj
    C::Matrix{Float64} = zeros(mesh.ne,mesh.nv)
    for m in 1:mesh.ne
        nds = mesh.surfaceT[1:3,m]

        C[m,nds] .= 1/6
        xm::Vector{Float64} = mean(mesh.p[:,nds],2)

        # Now add the boundary integral
        for s in 1:mesh.ne
            nds_j = mesh.surfaceT[1:3,s]

            # Nodes of the quadrature
            p::Matrix{Float64} = subtriangle(mesh.p[1:3,nds_j])
            
            aux::Vector{Float64} = [0,0,0]
            for quad in 1:size(p,2)
                y::Vector{Float64} = p[:,quad]
                r::Vector{Float64} = y-xm
                
                aux .+= dot(r,mesh.normal[:,s])/(norm(r)^3) .* phi[:,quad]
            end

            C[m,nds_j] .+= 1/(4*pi) * areaT[s]/size(p,2) .* aux
        end
    end

    return C
end # BEM matrix Mmj

function Dmatrix(mesh::MESH, areaT::Vector{Float64})

    # Mmn
    D::Matrix{Float64} = zeros(mesh.ne,mesh.ne)
    for m in 1:mesh.ne

        nds = mesh.surfaceT[1:3,m]                       # Nodes of the surface triangle m
        xm::Vector{Float64} = mean(mesh.p[1:3,nds],2)    # Center of edge

        for n in 1:mesh.ne
            nds = mesh.surfaceT[1:3,n]                   # Nodes of the surface triangle n

            # Quadrature coordinates
            p = subtriangle(mesh.p[1:3,nds])
            
            aux::Float64 = 0;
            for quad = 1:size(p,2)
                y::Vector{Float64} = p[:,quad]
                r::Float64 = norm(xm-y)

                aux += 1/r
            end

            D[m,n] += 1/(4*pi) * aux * areaT[n]/size(p,2)
        end
    end
    
    return D
end # BEM matrix Mmn

function surface2element(mesh::MESH,nds_j)
    k::Int32 = 0
    for ik in 1:mesh.nt
        nds = mesh.t[:,ik]
        inside::Vector{Bool} = Vector{Bool}(undef,3)
        for i in 1:length(nds_j)
            inside[i] =  nds_j[i] in nds ? true : false
        end
        if all(inside)
            k = ik
            println("Found the element that contains the surface triangle")
            break
        end
    end

    return k
end

function midpoints(r::Matrix{Float64})
    # All 3 mid points of the triangle edges
    t::Matrix{Float64} = 0.5.*(r[:,[1,2,1]] + r[:,[2,3,3]])
    return t
end

function subtriangle(p::Matrix{Float64})
    #=
        input: p[1:3,s_nds] -> xyz coordinates of surface triangle nodes
         s_nds -> nodes of surface triangle (3)
        
        Subdivides a triangle into 4 and calculates the midpoints
        of each edge of the small triangles

        output: 15 xyz coordinates 
    =#
    r::Matrix{Float64} = zeros(3,15) # 15 nodes

    # Large triangle | First 3 nodes
    r[:,1:3] .= p

    # Midpoints of large triangle
    r[:,4:6] = midpoints(r[:,1:3])

    # Mid points of each 4 small triangles
    r[:,7:9] = midpoints(r[:,[1,4,6]])
    r[:,10:12] = midpoints(r[:,[4,2,5]])
    r[:,13:15] = midpoints(r[:,[5,3,6]])

    return r
end # 4 sub-triangles


function main()
    meshSize::Float64 = 0

    # Constants
    mu0::Float64 = pi*4e-7          # vacuum magnetic permeability
    giro::Float64 = 2.210173e5 /mu0 # Gyromagnetic ratio (rad T-1 s-1)
    dt::Float64 = 1e-12             # Time step (s)
    totalTime::Float64 = 0.1        # Total time of spin dynamics simulation (ns)
    damp::Float64 = 0.0             # Damping parameter (dimensionless [0,1])
    precession::Float64 = 1.0       # Include precession or not (0 or 1)

    # Dimension of the magnetic material 
    L::Vector{Float64} = [512,128,30]
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

    println("Number of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))
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
    m[1,:] .= 1
    # begin # Random initial magnetization
       #  theta::Vector{Float64} = 2*pi.*rand(mesh.nv)
       #  phi::Vector{Float64} = pi.*rand(mesh.nv)
       #  for i in 1:mesh.nv
       #      m[:,i] = [sin(theta[i])*cos(phi[i]), sin(theta[i])*sin(phi[i]), cos(theta[i])]
       #      m[:,i] = m[:,i]./norm(m[:,i])
       #  end
    # end # Random initial magnetization

    # Hd::Matrix{Float64} = demagnetizingField(mesh,m,areaT,LHS)
    return
    
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