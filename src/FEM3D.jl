# Finite element method logic for 3D tetrahedral meshes


# 3D linear basis function
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


# 3D quadratic basis function
function quadraticBasis(mesh::MESH, nds::AbstractVector, nd::Integer)

    # f(x,y,z) = S[1] + S[2]*x + S[3]*y + S[4]*z 
    #          + S[5]*x^2 + S[6]*x*y + S[7]*x*z
    #          + S[8]*y^2 + S[9]*y*z + S[10]*z^2

    nodes = nds[nds .!= nd] # All other nodes of the element

    # x y z coordinates on the order 'nd' then the other 'nodes'
    xo = @view mesh.p[1, [nd; nodes]]
    yo = @view mesh.p[2, [nd; nodes]]
    zo = @view mesh.p[3, [nd; nodes]]

    S = [ones(10) xo yo zo xo.^2 xo.*yo xo.*zo yo.^2 yo.*zo zo.^2]\ [1.0; zeros(9)]

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
function localStiffnessMatrix(mesh::MESH,f::Vector{Float64}=ones(mesh.nt))
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

        H::Float64 = norm(H_vec[:, k])
        
        for i in 1:4
            for j in i:4
                aux[i,j] = mesh.VE[k]*dmu[k]/H *
                           (H_vec[1, k]*b[i] + H_vec[2, k]*c[i] + H_vec[3, k]*d[i]) *
                           (H_vec[1, k]*b[j] + H_vec[2, k]*c[j] + H_vec[3, k]*d[j])
            
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

# Local stiffness matrix with Nedelec shape elements
function nedelecLocalStiffness(mesh::MESH)

    # Local stiffness matrix
    Ak::Matrix{Float64} = zeros(36, mesh.nt) # 6x6 edges per volume element
    for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k] # Nodes of the linear volume element

        # Hat shape element for each of the 4 nodes
        hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
        for i in 1:4
            hat[1, i], 
            hat[2, i], 
            hat[3, i], 
            hat[4, i] = abcd(mesh.p, nds, nds[i])
        end

        curlN::Matrix{Float64} = zeros(3, 6)
        for ie in 1:6
            # Global node labels of the edge
            ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

            # Length of edge
            edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])
            
            # Curl of Nedelec shape element
            curlN[:, ie] = 2.0*edgeLength*cross(hat[2:4, ndi], hat[2:4, ndj])
        end

        n = 0
        # Update the stiffness matrix
        for ie in 1:6 # For each edge of the tetrahedron
            for je in 1:6
                n += 1
                Ak[n, k] = mesh.VE[k]*dot(curlN[:, ie], curlN[:, je])
            end
        end     # Loop over the edges of the element
        
    end # Loop over the volume element labels

    return Ak
end # Local stiffness matrix with Nedelec shape elements

# Global tangential stiffness matrix with Nedelec shape elements
function nedelecTangentialStiffness(mesh::MESH, 
                                    dnu::Vector{Float64},    # d/dB nu
                                    Bfield::Matrix{Float64}, # Vector field B
                                    B::Vector{Float64}       # Norm of B
                                    )
    
    # Global sparse tangential stiffness matrix
    At = spzeros(mesh.ne, mesh.ne)

    # Local tangential stiffness matrix
    Alocal::Matrix{Float64} = zeros(6, 6)    # For a single element
    Ak::Matrix{Float64} = zeros(36, mesh.nt) # For all elements
    for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k] # Nodes of the linear volume element

        # Hat shape element for each of the 4 nodes
        hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
        for i in 1:4
            hat[1, i], 
            hat[2, i], 
            hat[3, i], 
            hat[4, i] = abcd(mesh.p, nds, nds[i])
        end

        curlN::Matrix{Float64} = zeros(3, 6)
        for ie in 1:6
            # Global node labels of the edge
            ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

            # Length of edge
            edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])
            
            # Curl of Nedelec shape element
            curlN[:, ie] = 2.0*edgeLength*cross(hat[2:4, ndi], hat[2:4, ndj])
        end

        # Update local tangential stiffness matrix
        for i in 1:6
            for j in i:6 # Force symmetry in the tangential stiffness matrix

                Alocal[i, j] = mesh.VE[k] * dot(curlN[:, i], Bfield[:, k])  *
                                            dot(curlN[:, j], Bfield[:, k])  *
                                            dnu[k]/B[k]

                Alocal[j, i] = Alocal[i, j] # At should be symmetric
            end
        end

        # Update local tangential stiffness matrix
        n = 0
        for i in 1:6
            for j in 1:6
                n += 1
                Ak[n, k] = Alocal[i, j]
            end
        end

    end # Loop over the volume element labels

    # Update global sparse stiffness matrix
    n = 0
    for i in 1:6
        edge1 = mesh.edge2localMap[mesh.t[4+i, :]]

        for j in 1:6
            edge2 = mesh.edge2localMap[mesh.t[4+j, :]]

            n += 1
            At += sparse(edge1, edge2, Ak[n, :], mesh.ne, mesh.ne)
        end
    end

    return At
end # Tangential stiffness matrix

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
    for s in 1:mesh.ns

        # Only integrate over the outer shell
        if !(mesh.surfaceT[4,s] in shell_id)
            continue
        end
        
        RHS[mesh.surfaceT[1:3,s]] .+= dot(mesh.normal[:,s],F)*mesh.AE[s]/3;
    end

    return RHS
end # Boundary conditions

# Updates the Newton-Raphson iteration with a line search
# where the full step size is reduced to minimize the new residue
function lineSearch(u, du, A, RHS, minStep::Float64=1e-4)
    # Update the solution with a weight: u_new = u + alf*du
    alf = 1.0 # Initial weight (full step)

    # Residual before update of u
    residual_old = norm(RHS - A*u) 

    # New solution trial
    u_trial = u + du
    residue = norm(RHS - A*u_trial)
    
    # Reduce the step size until the new solution is more accurate than the old solution
    while residue > residual_old && alf > minStep
        alf *= 0.9 # Reduce the size of the weight by 10 %
        
        # New trial solution
        u_trial = u + alf*du
        residue = norm(RHS - A*u_trial)
    end
    
    # Update the solution
    if alf < minStep
        # Don't update u otherwise it will increase the residue
        println("Warning: N-R could not decrease the residue further")
    else
        u .= u_trial
    end

    return residue
end # Adapt the N-R step size based on the residual


function quadraticLocalStiffnessMatrix(mesh::MESH
                                       , S # Quadratic basis functions for every node and element
                                       , mu::Vector{Float64} = ones(mesh.nt) # Viscosity
                                      )
# S -> quadratic basis function. 10 by 10 by nt
# S[:, 1, 2] is the coef. of the basis function for the node 1 on
# the second tetrahedron

    Ak::Matrix{Float64} = zeros(100, mesh.nt)
    temp::Matrix{Float64} = zeros(10, 10)     # Local element wise stiffness matrix

    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        for i in 1:10
            Si = @view S[:, i, k]
            
            for j in i:10
                Sj = @view S[:, j, k]

                # 10 node quadrature
                aux::Float64 = 0.0
                for n in 1:10

                    x::Float64 = mesh.p[1, nds[n]] 
                    y::Float64 = mesh.p[2, nds[n]] 
                    z::Float64 = mesh.p[3, nds[n]]
                    
                    dxi::Float64 = Si[2] + 2*Si[5] *x + Si[6]*y + Si[7]*z
                    dyi::Float64 = Si[3] + 2*Si[8] *y + Si[6]*x + Si[9]*z
                    dzi::Float64 = Si[4] + 2*Si[10]*z + Si[7]*x + Si[9]*y

                    dxj::Float64 = Sj[2] + 2*Sj[5] *x + Sj[6]*y + Sj[7]*z
                    dyj::Float64 = Sj[3] + 2*Sj[8] *y + Sj[6]*x + Sj[9]*z
                    dzj::Float64 = Sj[4] + 2*Sj[10]*z + Sj[7]*x + Sj[9]*y

                    aux += dxi*dxj + dyi*dyj + dzi*dzj
                end 
                aux /= 10

                temp[i,j] = mesh.VE[k]*mu[k]*aux
                temp[j,i] = temp[i,j] # Stiffness matrix is symmetric
            end
        end # Local element wise stiffness matrix

        # Update the local stiffness matrix
        Ak[:, k] = vec(temp)

    end # Local stiffness matrix

    return Ak

end # Quadratic order, element-wise stiffness matrix

# Local element-wise divergence matrix
function localDivergenceMatrix(mesh::MESH, S)
    
    # Considering a combination of linear and quadratic Lagrange elements

    # S -> quadratic basis function. 10 by 10 by nt
    # S[:, 1, 2] is the coef. of the basis function for the node 1 on
    # the second tetrahedron

    # Local divergence matrix
    B1k::Matrix{Float64} = zeros(40, mesh.nt) # 4 by 10 by nt
    B2k::Matrix{Float64} = zeros(40, mesh.nt) # ...
    B3k::Matrix{Float64} = zeros(40, mesh.nt) # ...

    # Local Divergence matrix
    B1temp::Matrix{Float64} = zeros(4, 10) # Element wise matrix
    B2temp::Matrix{Float64} = zeros(4, 10) # ...
    B3temp::Matrix{Float64} = zeros(4, 10) # ...

    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]

        for i in 1:4
            a, b, c, d = abcd(mesh.p, nds[1:4], nds[i]) # Linear basis function
            
            for j in 1:10
                Sj = @view S[:, j, k]
                
                # 10 Node quadrature
                b1::Float64 = 0.0
                b2::Float64 = 0.0
                b3::Float64 = 0.0
                for n in 1:10

                    x::Float64 = mesh.p[1, nds[n]]
                    y::Float64 = mesh.p[2, nds[n]]
                    z::Float64 = mesh.p[3, nds[n]]

                    b1 -= (a + b*x + c*y + d*z)*                    # Linear
                          (Sj[2] + 2*Sj[5] *x + Sj[6]*y + Sj[7]*z)  # Quadratic
                    
                    b2 -= (a + b*x + c*y + d*z)*                    # Linear
                          (Sj[3] + 2*Sj[8] *y + Sj[6]*x + Sj[9]*z)  # Quadratic

                    b3 -= (a + b*x + c*y + d*z)*                    # Linear
                          (Sj[4] + 2*Sj[10]*z + Sj[7]*x + Sj[9]*y)  # Quadratic


                end # 10 node quadrature (quadratic nodes)

                # Element wise divergence matrix
                B1temp[i, j] = mesh.VE[k]*b1/10
                B2temp[i, j] = mesh.VE[k]*b2/10
                B3temp[i, j] = mesh.VE[k]*b3/10

            end # Quadratic nodes loop
        end # Linear nodes loop

        # Update local divergence matrix
        B1k[:, k] = vec(B1temp')
        B2k[:, k] = vec(B2temp')
        B3k[:, k] = vec(B3temp')
    
    end # Loop over the elements

    return B1k, B2k, B3k

end # Divergence matrix

function massMatrix(mesh::MESH)
    # Local mass matrix
    Mlocal::Matrix{Float64} = 1/20 * [2 1 1 1;
                                      1 2 1 1;
                                      1 1 2 1;
                                      1 1 1 2]; # 3D element-wise mass matrix

    Mk::Matrix{Float64} = zeros(16, mesh.nt)
    for k in 1:mesh.nt
        Mk[:,k] = mesh.VE[k]*Mlocal[:];
    end

    # Create the sparse global matrix
    M = spzeros(mesh.nv,mesh.nv)
    n = 0
    for i in 1:4
        for j in 1:4
            n += 1
            M += sparse(  mesh.t[i, :]
                        , mesh.t[j, :]
                        , Mk[n, :]
                        , mesh.nv, mesh.nv)
        end
    end

    return M
end