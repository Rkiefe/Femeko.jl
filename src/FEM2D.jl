# Finite element method logic for 2D triangular meshes


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

# 2D Global stiffness matrix
function stiffnessMatrix2D(mesh::MESH, mu::Vector{Float64}=ones(mesh.nt))
    # Global sparse stiffness matrix
    A = spzeros(mesh.nv,mesh.nv)

    # Local stiffness matrix
    Ak::Matrix{Float64} = zeros(9, mesh.nt)
    b::Vector{Float64} = [0,0,0]
    c::Vector{Float64} = [0,0,0]
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        for i in 1:3
            _, b[i], c[i] = abc(mesh.p, nds, nds[i])
        end

        aux = mesh.VE[k]*mu[k]*(b*b' + c*c')
        Ak[:,k] = aux[:]
    end

    # Update sparse global matrix
    n = 0
    for i in 1:3
        for j in 1:3
            n += 1
            A += sparse(mesh.t[i,:],mesh.t[j,:],Ak[n,:],mesh.nv,mesh.nv)
        end
    end

    return A
end # 2D Global stiffness matrix


# Local 2D quadratic stiffness matrix
function quadraticLocalStiffnessMatrix2D(mesh::MESH)
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

                temp[i,j] = aux*mesh.VE[k]
                temp[j,i] = temp[i,j] # It is symmetric
            end
        end

        Ak[:,k] = vec(temp)
    end # Local stiffness matrix

    return Ak
end

function quadraticStiffnessMatrix2D(  mesh::MESH
                                    , mu::Vector{Float64} = ones(mesh.nt))

    # Local stiffness matrix
    Ak::Matrix{Float64} = quadraticLocalStiffnessMatrix2D(mesh)

    # Build the sparse stiffness matrix
    A = spzeros(mesh.nv, mesh.nv)    
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            A += sparse(  mesh.t[i, :]
                        , mesh.t[j, :]
                        , Ak[n, :] .* mu
                        , mesh.nv, mesh.nv)
        end
    end

    return A
end

function massMatrix2D(mesh::MESH)
    Mlocal::Matrix{Float64} = 1/12 *[2 1 1;
                                     1 2 1;
                                     1 1 2]

    M = spzeros(mesh.nv,mesh.nv)
    Mk::Matrix{Float64} = zeros(9, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        Mk[:,k] = mesh.VE[k]*Mlocal[:];
    end

    # Update sparse global matrix
    n = 0
    for i in 1:3
        for j in 1:3
            n += 1
            M += sparse(mesh.t[i,:],mesh.t[j,:],Mk[n,:],mesh.nv,mesh.nv)
        end
    end

    return M
end

function quadraticMassMatrix2D(mesh::MESH)
    # 6-point Gaussian quadrature for triangles (exact precision)
    # Coordinates and weights for reference triangle (0,0), (1,0), (0,1)
    weights::Vector{Float64}, 
    points::Matrix{Float64} = GaussQuadrature2D(4)

    Mlocal = zeros(36, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]
        vertices = mesh.p[:, nds[1:3]]  # Triangle vertices
        
        # Jacobian for affine transformation from reference to physical triangle
        J11 = vertices[1, 2] - vertices[1, 1]
        J12 = vertices[1, 3] - vertices[1, 1]
        J21 = vertices[2, 2] - vertices[2, 1]
        J22 = vertices[2, 3] - vertices[2, 1]
        
        # detJ = abs(J11 * J22 - J12 * J21) # = 2x the area of the triangle. Note: The reference triangle area = 0.5

        # Precompute basis coefficients for all 6 nodes
        S::Matrix{Float64} = zeros(6 ,6)  
        for i in 1:6
            S[:, i] = quadraticBasis2D(mesh.p, nds, nds[i])
        end
        
        # Initialize local mass matrix
        Mk = zeros(6, 6)
        
        # Loop over quadrature points
        for q in 1:length(weights)
            xi, eta = points[q, 1:2]
            
            # Transform from reference to physical coordinates
            x = vertices[1, 1] + J11 * xi + J12 * eta
            y = vertices[2, 1] + J21 * xi + J22 * eta
            
            # Evaluate all basis functions at this quadrature point
            phi = zeros(6)
            for i in 1:6
                phi[i] = S[1, i] + S[2, i]*x + S[3, i]*y + S[4, i]*x^2 + S[5, i]*x*y + S[6, i]*y^2
            end
            
            # Accumulate to mass matrix
            w = weights[q] * mesh.VE[k] # detJ/2
            for i in 1:6
                for j in i:6  # Only compute upper triangle
                    Mk[i, j] += w * phi[i] * phi[j]
                    Mk[j, i] = Mk[i, j] # Fill lower triangle (symmetric)
                end
            end
        end
        
        Mlocal[:, k] = vec(Mk)
    end
    
    # Assemble global sparse mass matrix
    M = spzeros(mesh.nv, mesh.nv)
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            M += sparse(mesh.t[i, :], mesh.t[j, :], Mlocal[n, :], mesh.nv, mesh.nv)
        end
    end
    
    return M
end

function quadraticConvectionMatrix2D(mesh::MESH, u::Matrix{Float64})

    # Quadratic shape function times its gradient -> Order 3 integrand
    weights::Vector{Float64}, 
    points::Matrix{Float64} = GaussQuadrature2D(3)

    # Local convection matrix
    Clocal = zeros(36, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:, k]
        vertices = mesh.p[:, nds[1:3]]  # Triangle vertices
        
        # Jacobian for affine transformation from reference to physical triangle
        J11 = vertices[1, 2] - vertices[1, 1]
        J12 = vertices[1, 3] - vertices[1, 1]
        J21 = vertices[2, 2] - vertices[2, 1]
        J22 = vertices[2, 3] - vertices[2, 1]
        
        # detJ = abs(J11 * J22 - J12 * J21)  # = 2x the area of the triangle. Note: The reference triangle area = 0.5

        # Precompute basis coefficients for all 6 nodes
        S::Matrix{Float64} = zeros(6 ,6)  
        for i in 1:6
            S[:, i] = quadraticBasis2D(mesh.p, nds, nds[i])
        end
        
        # Initialize element-wise matrix
        Ck = zeros(6, 6)

        # Loop over quadrature points
        for q in 1:length(weights)
            xi, eta = points[q, 1:2]
            
            # Transform from reference to physical coordinates
            x = vertices[1, 1] + J11 * xi + J12 * eta
            y = vertices[2, 1] + J21 * xi + J22 * eta
            
            # Evaluate on the current quadrature point the
                #  2nd order Lagrange shape function
                #  The velocity field
                #  The gradient of the 2nd order Lagrange shape function
            
            phi::Vector{Float64} = zeros(6)         # Shape function 
            ux::Float64 = 0.0; uy::Float64 = 0.0    # Velocity field
            gradPhi::Matrix{Float64} = zeros(2, 6)  # Gradient of shape function
            
            for i in 1:6 # All nodes of the element
            
                # Basis function on the quadrature point
                phi[i] = S[1, i] + S[2, i]*x + S[3, i]*y + S[4, i]*x^2 + S[5, i]*x*y + S[6, i]*y^2
                
                # Velocity on the quadrature point
                ux += u[1, nds[i]]*phi[i]
                uy += u[2, nds[i]]*phi[i]

                # Gradient of basis function on quadrature point
                gradPhi[1, i] = S[2, i] + 2*S[4, i]*x + S[5, i]*y
                gradPhi[2, i] = S[3, i] + 2*S[6, i]*y + S[5, i]*x
            end

            # Accumulate to local matrix
            w = weights[q] * mesh.VE[k] # detJ/2
            for i in 1:6
                grad_i = gradPhi[1, i] * ux + gradPhi[2, i] * uy
                for j in 1:6
                    Ck[i, j] += w * grad_i * phi[j]
                end
            end
        end # Loop over quadrature points
        
        Clocal[:, k] = vec(Ck)
    end # Loop over the element labels
    
    # Sparse global convection matrix
    C = spzeros(mesh.nv, mesh.nv)
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            C += sparse(  mesh.t[i, :]
                        , mesh.t[j, :]
                        , Clocal[n, :]
                        , mesh.nv, mesh.nv)
        end
    end

    return C
end

# 2D Sparse divergence matrix
function divergenceMatrix2D(mesh::MESH, vertexID::Vector{Int32}, nVertices::Int32)
    
    # Global sparse divergence matrix
    B1 = spzeros(nVertices, mesh.nv) # Vertices x Nodes
    B2 = spzeros(nVertices, mesh.nv) # Vertices x Nodes
    
    # Local dense divergence matrix
    B1k::Matrix{Float64}, 
    B2k::Matrix{Float64} = localDivergenceMatrix2D(mesh::MESH)

    # Update sparse global matrix
    n = 0
    for i in 1:3
        for j in 1:6
            n += 1
            B1 += sparse(vertexID[mesh.t[i,:]],      # Convert original node ID to sorted node ID
                         vertexID[mesh.t[j,:]],      # Convert original node ID to sorted node ID
                         B1k[n,:], nVertices, mesh.nv)

            B2 += sparse(vertexID[mesh.t[i,:]],      # Convert original node ID to sorted node ID
                         vertexID[mesh.t[j,:]],      # Convert original node ID to sorted node ID
                         B2k[n,:], nVertices, mesh.nv)
        end
    end

    return B1, B2
end # Sparse divergence matrix

# 2D Local dense divergence matrix
function localDivergenceMatrix2D(mesh::MESH)

    # Local divergence matrix
    B1k::Matrix{Float64} = zeros(18, mesh.nt)
    B2k::Matrix{Float64} = zeros(18, mesh.nt)

    B1_local::Matrix{Float64} = zeros(3,6) # Temporary for assembly of the local matrix
    B2_local::Matrix{Float64} = zeros(3,6) # ...

    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]

        for i in 1:3
            a::Float64, b::Float64, c::Float64 = abc(mesh.p, nds[1:3], nds[i])
            for j in 1:6
                S = quadraticBasis2D(mesh.p, nds, nds[j])
                
                # 6 Node quadrature
                b1::Float64 = 0.0
                b2::Float64 = 0.0
                for n in 1:6
                    b1 -= (a + b*mesh.p[1,nds[n]] + c*mesh.p[2,nds[n]])*            # Linear
                          (S[2] + 2*S[4]*mesh.p[1,nds[n]] + S[5]*mesh.p[2,nds[n]])  # Quadratic
                    
                    b2 -= (a + b*mesh.p[1,nds[n]] + c*mesh.p[2,nds[n]])*            # Linear
                          (S[3] + S[5]*mesh.p[1,nds[n]] + 2*S[6]*mesh.p[2,nds[n]])  # Quadratic
                end # 6 node quadrature (quadratic nodes)

                B1_local[i,j] = mesh.VE[k]*b1/6
                B2_local[i,j] = mesh.VE[k]*b2/6
            
            end # Quadratic nodes loop
        end # Linear nodes loop

        # Update local divergence matrix
        B1k[:,k] = vec(B1_local')
        B2k[:,k] = vec(B2_local')

    end # Element loop

    return B1k, B2k
end # 2D Local dense divergence matrix

# 1D quadratic mass matrix, used to calculate the surface integral of over 2D triangles
function quadraticMassMatrix1D(mesh::MESH, h::Vector{Float64} = ones(mesh.ns))


    # mesh.surfaceT[1:3, s] >> node1, node2, midpoint
    Mk::Matrix{Float64} = 1.0/30.0 *[4  2 -1;
                                     2  16 2;
                                     -1 2  4]

    # Local mass matrix
    M = spzeros(mesh.nv, mesh.nv)
    for s in 1:mesh.ns # For each edge (3 nodes)
        
        nds = @view mesh.surfaceT[:, s]

        for i in 1:3
            for j in i:3 # Symmetry
                M[nds[i], nds[j]] += mesh.AE[s]*Mk[i,j]*h[s]
                M[nds[j], nds[i]] = M[nds[i], nds[j]]  # Mass matrix is symmetric 
            end
        end

    end

    return M
end

# All 3 mid points of the triangle edges
function midpoints(r::Matrix{Float64})
    t::Matrix{Float64} = 0.5.*(r[:,[1,2,1]] + r[:,[2,3,3]])
    return t
end

# Subdivides a triangle into 4 smaller triangles
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

function GaussQuadrature2D(precision::Integer)
# From Larson 2013 FEM book "Theory implementation and applications"
# on 2D reference element quadrature points and weights

# weights -> vector of weight values by point
# points -> matrix n by 2 with the x,y coordinates of each quadrature point

weights::Vector{Float64} = zeros(precision)
points::Matrix{Float64} = zeros(precision, 2)

    if precision == 1

        weights = [1.0]
        points = [1.0/3.0 1.0/3.0]
    
    elseif precision == 2
    
        weights = [1.0/3.0, 1.0/3.0, 1.0/3.0]
        points =  [1.0/6.0 1.0/6.0;
                   2.0/3.0 1.0/6.0;
                   1.0/6.0 2.0/3.0]
    
    elseif precision == 3
    
        weights = [-27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0]
        points = [1.0/3.0 1.0/3.0;
                  0.2     0.2;
                  0.6     0.2;
                  0.2     0.6]
    
    elseif precision == 4

        weights = [ 0.223381589678011
                    0.223381589678011
                    0.223381589678011
                    0.109951743655322
                    0.109951743655322
                    0.109951743655322]

        points = [0.445948490915965 0.445948490915965;
                  0.445948490915965 0.108103018168070;
                  0.108103018168070 0.445948490915965;
                  0.091576213509771 0.091576213509771;
                  0.091576213509771 0.816847572980459;
                  0.816847572980459 0.091576213509771]

    else
        @error "Requested precision for 2D Gaussian quadrature is not available. Limit is 4"
        return nothing
    end

    return weights, points 

end


function GaussQuadrature1D(order::Integer)
# From wikipedia | https://en.wikipedia.org/wiki/Gaussian_quadrature

weights::Vector{Float64} = zeros(order)
points::Vector{Float64} = zeros(order)    

    # Returns (points, weights) for 1D Gaussian quadrature on [-1, 1]
    if order == 1

         points = [0.0] 
         weights = [2.0]

    elseif order == 2

         points = [-1/sqrt(3), 1/sqrt(3)]
         weights = [1.0, 1.0]

    elseif order == 3

         points = [-sqrt(3/5), 0.0, sqrt(3/5)]
         weights = [5/9, 8/9, 5/9]

    elseif order == 4

        points = [
                  -sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0)),
                  -sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)),
                   sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0)),
                   sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0))
                 ]

        weights = [
                    18.0-sqrt(30)/36.0, 
                    18.0+sqrt(30)/36.0, 
                    18.0+sqrt(30)/36.0, 
                    18.0-sqrt(30)/36.0, 
                  ] 

    else
        error("Requested precision for 1D Gaussian quadrature is not available. Limit is 4")
        return nothing
    end
    
    return points, weights
end