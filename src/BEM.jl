#=
    Includes the BEM matrices and quadrature logic for magnetostatic simulations
    with the mixed FEM-BEM formulation
=#

include("FEM.jl")

function denseStiffnessMatrix(mesh::MESH, F::Vector{Float64}=ones(mesh.nt))
    #= 
        This is is the same as the FEM stiffness matrix, but using the dense Matrix
        data type instead of spzeros.
    =#
    A::Matrix{Float64} = zeros(mesh.nv,mesh.nv)

    b::Vector{Float64} = zeros(4)
    c::Vector{Float64} = zeros(4)
    d::Vector{Float64} = zeros(4)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        for i in 1:4
            _,b[i],c[i],d[i] = abcd(mesh.p,nds,nds[i])
        end
        A[nds,nds] .+= mesh.VE[k]*(b*b' + c*c' + d*d')*F[k]
    end

    return A
end

function Bmatrix(mesh::MESH)
    B::Matrix{Float64} = zeros(mesh.nv, mesh.ns)
    for s in 1:mesh.ns
        nds = mesh.surfaceT[1:3, s]
        B[nds,s] .+= mesh.AE[s]/3
    end
    return B
end # BEM matrix Min

function Cmatrix(mesh::MESH)
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
    C::Matrix{Float64} = zeros(mesh.ns, mesh.nv)
    for m in 1:mesh.ns
        nds = mesh.surfaceT[1:3,m]

        C[m,nds] .= 1/6
        xm::Vector{Float64} = mean(mesh.p[:,nds], 2)

        # Now add the boundary integral
        for s in 1:mesh.ns
            nds_j = mesh.surfaceT[1:3,s]

            # Nodes of the quadrature
            p::Matrix{Float64} = subtriangle(mesh.p[1:3,nds_j])
            
            aux::Vector{Float64} = [0,0,0]
            for quad in 1:size(p,2)
                y::Vector{Float64} = p[:,quad]
                r::Vector{Float64} = y-xm
                
                aux .+= dot(r,mesh.normal[:,s])/(norm(r)^3) .* phi[:,quad]
            end

            C[m,nds_j] .+= 1/(4*pi) * mesh.AE[s]/size(p,2) .* aux
        end
    end

    return C
end # BEM matrix Mmj

function Dmatrix(mesh::MESH)

    # Mmn
    D::Matrix{Float64} = zeros(mesh.ns,mesh.ns)
    for m in 1:mesh.ns

        nds = mesh.surfaceT[1:3,m]                       # Nodes of the surface triangle m
        xm::Vector{Float64} = mean(mesh.p[1:3,nds],2)    # Center of edge

        for n in 1:mesh.ns
            nds = mesh.surfaceT[1:3,n]                   # Nodes of the surface triangle n

            # Quadrature coordinates
            p = subtriangle(mesh.p[1:3,nds])
            
            aux::Float64 = 0;
            for quad = 1:size(p,2)
                y::Vector{Float64} = p[:,quad]
                r::Float64 = norm(xm-y)

                aux += 1/r
            end

            D[m,n] += 1/(4*pi) * aux * mesh.AE[n]/size(p,2)
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
