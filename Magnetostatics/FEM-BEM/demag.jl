#=
    Gets the demagnetizing field of a magnetization field using 
    the mixed-element FEM-BEM coupling
    
    This implementation is based on https://doi.org/10.1016/j.jmmm.2012.01.016
    but adapted to take magnetization fields (M) instead of a source field and a permeability

    The units of the demagnetizing field are the same as the units of the input M
=#

include("../../src/Femeko.jl")
include("../../src/BEM.jl")

using GLMakie

meshSize = 1.0
localSize = 0.0
showGmsh = false

function main(meshSize=0.0, localSize=0.0, showGmsh=false, verbose=true)
    gmsh.initialize()

    # Applied field | A/m
    mu0::Float64 = pi*4e-7          # Vaccum magnetic permeability
    Ms::Float64 = 1.0 # Tesla

    # Numerical method settings
    maxAtt::Int32 = 100
    maxDeviation::Float64 = 1e-7

    # Create geometry
    cells = []
    addCuboid([0,0,0], [1.0, 1.0, 0.1], cells)  
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

    # Initial magnetization state
    M::Matrix{Float64} = zeros(3, mesh.nv)
    M[1, :] .= Ms

    # BEM matrices
    println("Building the BEM matrices and the FEM-BEM coupling")
    @time begin 
        B = Bmatrix(mesh)
        C = Cmatrix(mesh)
        D = Dmatrix(mesh)
    end

    println("Building the global stiffness matrix")
    A::Matrix{Float64} = zeros(mesh.nv, mesh.nv)
    b::Matrix{Float64} = zeros(4, mesh.nt)
    c::Matrix{Float64} = zeros(4, mesh.nt)
    d::Matrix{Float64} = zeros(4, mesh.nt)
    @time for k in 1:mesh.nt

        nds = @view mesh.t[:, k]
        for i in 1:4
            _,
            b[i, k],
            c[i, k],
            d[i, k] = abcd(mesh.p, nds, nds[i]) 

        end

        for i in 1:4
            for j in i:4
                A[nds[i], nds[j]] += (b[i, k]*b[j, k] +
                                      c[i, k]*c[j, k] +
                                      d[i, k]*d[j, k]) * mesh.VE[k] 
                
                A[nds[j], nds[i]] = A[nds[i], nds[j]]
            end
        end
    end
    
    # Full FEM-BEM coupled matrix
    LHS = [A -B; C D]

    # Load vector 
    RHS::Vector{Float64} = zeros(mesh.nv + mesh.ns)
    for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k]
        for i = 1:4
            RHS[nds[i]] +=  (M[1, nds[i]]*b[i, k] + 
                             M[2, nds[i]]*c[i, k] + 
                             M[3, nds[i]]*d[i, k]) * mesh.VE[k]
        end
    end

    # Magnetostatic potential
    u = LHS\RHS

    # Demagnetizing field
    Hdfield::Matrix{Float64} = zeros(3, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]
        for i in 1:4
            Hdfield[1, k] -= b[i, k]*u[nds[i]]
            Hdfield[2, k] -= c[i, k]*u[nds[i]]
            Hdfield[3, k] -= d[i, k]*u[nds[i]]
        end
    end

    H::Vector{Float64} = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = norm(Hdfield[:, k])
    end


    # Plot vector field
    println("Generating plots...")
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = :data)

    graph = arrows3d!(    centroids[1,:]
                        , centroids[2,:]
                        , centroids[3,:]
                        , Hdfield[1,:]
                        , Hdfield[2,:]
                        , Hdfield[3,:]
                        , color = H
                        , colormap = :CMRmap  # :CMRmap :viridis :redsblues :turbo :rainbow
                        , lengthscale = 0.5
                        )

    Colorbar(fig[1, 2], graph, label="Hd (T)")
    wait(display(fig))
end

main(meshSize, localSize, showGmsh)
