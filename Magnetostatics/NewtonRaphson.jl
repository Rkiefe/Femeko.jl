#=
    Newton-Raphson iteration method for magnetostatic simulations of nonlinear
    magnetic properties
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

using GLMakie

function main(mesh::MESH, data)

    # Magnetic permeability
    mu0::Float64 = pi*4e-7
    mu::Vector{Float64} = zeros(mesh.nt) .+ mu0;

    # Boundary conditions
    RHS::Vector{Float64} = BoundaryIntegral(mesh,Hext,shell_id)

    # Lagrange multiplier technique
    Lag::Vector{Float64} = lagrange(mesh)

    # -- Picard iteration method --
    H_vec::Matrix{Float64} = zeros(3,mesh.nt)
    H::Vector{Float64} = zeros(mesh.nt)
    H_old::Vector{Float64} = zeros(mesh.nt)

    err::Float64 = maxDeviation
    att::Int32 = 0
    while err > maxDeviation && att < maxAtt
        att += 1

        # Store the last iteration magnetic field
        H_old .= H
        
        # Set magnetic permeability
        mu[mesh.InsideElements] .= interp1(data.H,data.B,H[mesh.InsideElements])./H[mesh.InsideElements]

        # Stiffness matrix
        A = stiffnessMatrix(mesh,mu)

        # Magnetic scalar potential
        u::Vector{Float64} = [A Lag;Lag' 0]\[-RHS;0]
        u = u[1:mesh.nv]
        
        # New magnetic field 
        H_vec .= zeros(3,mesh.nt) # Reset to 0
        for k in 1:mesh.nt
            nds = mesh.t[:,k] # all nodes of that element

            # Sum the contributions
            for ind = 1:numel(nds)
                nd = nds[ind]

                # obtain the element parameters
                _,b,c,d = abcd(mesh.p,nds,nd)

                H_vec[1,k] -= u[nd]*b;
                H_vec[2,k] -= u[nd]*c;
                H_vec[3,k] -= u[nd]*d;
            end
        end

        # |H| and check the deviation to the previous solution
        err = 0.0
        for k in 1:mesh.nt
            H[k] = norm(H_vec[:,k])
            err = max(err, mu0*(H[k] - H_old[k])) 
        end

    end # Picard iteration

    # -- Newton Raphson --

    # Define d_dH mu
    dmu_dH::Vector{Float64} = gradient(data.H, data.B./data.H)
    # Missing!! | Finding Inf and NaN in dmu_dH

    # this code is not complete. Not tests have been done and 
    # the Newton raphson iteration is missing

end
