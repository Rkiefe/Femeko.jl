using LinearAlgebra, SparseArrays

# FEM linear basis function
function abcd(p::Matrix{Float64},nodes::AbstractVector,nd::Int32)
    n1,n2,n3 = nodes[nodes .!= nd]
    x = @view p[1, [nd, n1, n2, n3]]
    y = @view p[2, [nd, n1, n2, n3]]
    z = @view p[3, [nd, n1, n2, n3]]

    r::Vector{Float64} =[1.0 x[1] y[1] z[1];
                         1.0 x[2] y[2] z[2];
                         1.0 x[3] y[3] z[3];
                         1.0 x[4] y[4] z[4]]\[1;0;0;0]

    return r[1],r[2],r[3],r[4] # a b c d
end # Basis function coef.

# 2D linear basis function
function abc(p::Matrix{Float64}, nds::AbstractVector, nd::Int32)

    nd1::Int32 = 0
    nd2::Int32 = 0
    for i in nds
        if i != nd 
            if nd1==0
                nd1 = i
            else nd2==0
                nd2 = i
            end
        end
    end

    M::Matrix{Float64} = [1 p[1,nd]  p[2,nd];
                          1 p[1,nd1] p[2,nd1];
                          1 p[1,nd2] p[2,nd2]]
    S::Vector{Float64} = M\[1,0,0]

    return S[1], S[2], S[3] # f = a + bx + cy
end # 2D Linear basis function


# 2D quadratic basis function
function quadraticBasis2D(p::Matrix{Float64}, nodes::AbstractVector,nd::Int32)
    # f = S[1] + S[2]*x + S[3]*y + S[4]*x^2 + S[5]*x*y + S[6]*y^2

    nds = nodes[nodes .!= nd]

    x = @view p[1, [nd, nds[1], nds[2], nds[3], nds[4], nds[5]] ]
    y = @view p[2, [nd, nds[1], nds[2], nds[3], nds[4], nds[5]] ]
    z = @view p[3, [nd, nds[1], nds[2], nds[3], nds[4], nds[5]] ]

    S::Vector{Float64} = [1.0 x[1] y[1] x[1]^2  x[1]*y[1] y[1]^2;
                          1.0 x[2] y[2] x[2]^2  x[2]*y[2] y[2]^2;
                          1.0 x[3] y[3] x[3]^2  x[3]*y[3] y[3]^2;
                          1.0 x[4] y[4] x[4]^2  x[4]*y[4] y[4]^2;
                          1.0 x[5] y[5] x[5]^2  x[5]*y[5] y[5]^2;
                          1.0 x[6] y[6] x[6]^2  x[6]*y[6] y[6]^2]\[1;0;0;0;0;0]

    return S
end


# Sparse, global stiffness matrix
function stiffnessMatrix(mesh::MESH, f::Vector{Float64}=ones(mesh.nt))
    A = spzeros(mesh.nv,mesh.nv)

    # Local stiffness matrix
    Ak::Matrix{Float64} = localStiffnessMatrix(mesh,f)

    # Update sparse global matrix
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            A += sparse(mesh.t[i,:],mesh.t[j,:],Ak[n,:],mesh.nv,mesh.nv)
        end
    end

    return A
end # Sparse, global stiffness matrix

# Local stiffnessmatrix in 100% Julia
function localStiffnessMatrix(mesh::MESH,f::Vector{Float64})
    Ak::Matrix{Float64} = zeros(4*4,mesh.nt)
    b::Vector{Float64} = zeros(4)
    c::Vector{Float64} = zeros(4)
    d::Vector{Float64} = zeros(4)
    aux::Matrix{Float64} = zeros(4,4)
    
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        for i in 1:4
            _,b[i],c[i],d[i] = abcd(mesh.p,nds,nds[i])
        end
        aux = mesh.VE[k]*f[k]*(b*b' + c*c' + d*d')
        Ak[:,k] = aux[:] # vec(aux)
    end

    return Ak
end # Local stiffnessmatrix in 100% Julia

# Local 2D quadratic stiffness matrix
function quadraticLocalStiffnessMatrix2D(mesh::MESH, mu::Vector{Float64})
    # Local stiffness matrix
    Ak::Matrix{Float64} = zeros(36, mesh.nt) # 6 x 6
    temp::Matrix{Float64} = zeros(6,6)

    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]

        temp .= 0
        for i in 1:length(nds)
            Si = quadraticBasis2D(mesh.p, nds, nds[i])
            
            for j in i:length(nds)
                Sj = quadraticBasis2D(mesh.p, nds, nds[j])

                # 6 node quadrature
                aux::Float64 = 0.0
                for n in 1:6
                    dxi::Float64 = Si[2] + 2*Si[4]*mesh.p[1,nds[n]] + Si[5]*mesh.p[2,nds[n]]
                    dyi::Float64 = Si[3] + Si[5]*mesh.p[1,nds[n]] + 2*Si[6]*mesh.p[2,nds[n]] 
                    dxj::Float64 = Sj[2] + 2*Sj[4]*mesh.p[1,nds[n]] + Sj[5]*mesh.p[2,nds[n]]
                    dyj::Float64 = Sj[3] + Sj[5]*mesh.p[1,nds[n]] + 2*Sj[6]*mesh.p[2,nds[n]]
                    aux += dxi*dxj + dyi*dyj
                end 
                aux /= 6

                temp[i,j] = mu[k]*aux*mesh.VE[k]
                temp[j,i] = temp[i,j] # It is symmetric
            end
        end

        Ak[:,k] = vec(temp)
    end # Local stiffness matrix

    return Ak
end

# Tangential stiffness matrix for Newton Iteration
function tangentialStiffnessMatrix(mesh::MESH, H_vec::Matrix{Float64}, dmu::Vector{Float64})

    At = spzeros(mesh.nv,mesh.nv)
    Ak::Matrix{Float64} = zeros(16, mesh.nt)
    b::Vector{Float64} = zeros(4)
    c::Vector{Float64} = zeros(4)
    d::Vector{Float64} = zeros(4)
    aux::Matrix{Float64} = zeros(4,4)

    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        for i in 1:4
            _,b[i],c[i],d[i] = abcd(mesh.p,nds,nds[i])
        end

        H::Float64 = norm(H_vec[k,:])
        
        for i in 1:4
            for j in i:4
                aux[i,j] = mesh.VE[k]*dmu[k]/H *
                           (H_vec[k,1]*b[i] + H_vec[k,2]*c[i] + H_vec[k,3]*d[i]) *
                           (H_vec[k,1]*b[j] + H_vec[k,2]*c[j] + H_vec[k,3]*d[j])
            
                aux[j,i] = aux[i,j]
            end
        end

        Ak[:,k] = aux[:] # vec(aux)
    end

    # Update sparse global matrix
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            At += sparse(mesh.t[i,:],mesh.t[j,:],Ak[n,:],mesh.nv,mesh.nv)
        end
    end

    return At
end # Newton iteration matrix

# Lagrange multiplier technique
function lagrange(mesh::MESH)
    C::Vector{Float64} = zeros(mesh.nv)
    for k in 1:mesh.nt
        nds::AbstractVector{Int32} = @view mesh.t[:,k];   # Nodes of that element
        C[nds] .+= mesh.VE[k]/4
    end
    return C
end # Lagrange multiplier technique

# Boundary condition
function BoundaryIntegral(mesh::MESH,F::Vector{Float64},shell_id)
    RHS::Vector{Float64} = zeros(mesh.nv);
    for s in 1:mesh.ne

        # Only integrate over the outer shell
        if !(mesh.surfaceT[4,s] in shell_id)
            continue
        end
        
        RHS[mesh.surfaceT[1:3,s]] .+= dot(mesh.normal[:,s],F)*mesh.AE[s]/3;
    end

    return RHS
end # Boundary conditions

# Mean function
function mean(arr::Vector,dimension=1)
    m::Real = 0
    for x in arr
        m += x
    end

    return m/length(arr)
end # End of mean for vectors

function mean(arr::Matrix,dimension=1)
    if dimension == 1 
        d = size(arr,1)
        m = 0 .*arr[1,:]
        for i in 1:d
            m .+= arr[i,:]
        end

    else
        d = size(arr,2)
        m = 0 .*arr[:,1]
        for i in 1:d
            m .+= arr[:,i]
        end
    end

    return m./d
end # End of mean for 2D matrices

# Linear interp function
function interp1(x::Vector, y::Vector, xq::Real)
    if minimum(x) > xq || maximum(x) < xq
        error("Interp1 | xq is out of bounds of x")
        return
    end

    yq::Real = NaN
    if xq == x[1]
        yq = y[1]
        return yq
    elseif xq == x[end]
        yq = y[end]
        return y[end]
    end

    for i in 2:length(x)
        if x[i] > xq && x[i-1] < xq
            yq = y[i-1] + (y[i]-y[i-1])/(x[i]-x[i-1]) *(xq-x[i-1])
            return yq
        end
    end
end # Interp a single value in a dataset

# Linear interp over an array
function interp1(x::Vector, y::Vector, xq::Vector)
    yq::Vector = zeros(length(xq))
    for i in 1:length(xq)
        yq[i] = interp1(x, y, xq[i])
    end
    return yq
end # Interp over the entire input array

# 1D Gradient function
function gradient(x::Vector{Float64},y::Vector{Float64})
    #=
        Forward and backward finite difference in the limits of the dataset ;
        Central finite difference on everything else
    =#

    if length(x) != length(y)
        error("Input x and y must have the same length")
    end
    
    if length(x) < 2
        error("Input data must have more than 1 element")
    end

    size::Int32 = length(x)
    z::Vector{Float64} = zeros(size)


    z[1] = (y[2] - y[1])/(x[2] - x[1])
    z[size] = (y[size] - y[size-1])/(x[size] - x[size-1])

    for i in 2:size-1
        z[i] = (y[i+1] - y[i-1])/(x[i+1] - x[i-1])
    end

    return z
end