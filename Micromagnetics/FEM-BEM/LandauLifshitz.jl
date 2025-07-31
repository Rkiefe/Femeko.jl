#=
    3D Finite element implementation of the Landau-Lifshitz equation with a semi-implciti
    time step, based on https://doi.org/10.1109/TMAG.2008.2001666 (Oriano 2008)
=#

# Include FEM functions
include("../../src/BEM.jl")

# For the thermal field
using Random
import Distributions as Dist

# Demagntizing field with FEM/BEM
function BEMdmag(mesh::MESH,m::Matrix{Float64},LHS::Matrix{Float64})
    
    RHS::Vector{Float64} = zeros(mesh.nv + mesh.ne)
    for s in 1:mesh.ne
        nds = @view mesh.surfaceT[1:3,s]
        RHS[nds] .+= dot(mesh.normal[:,s],mean(m[:,nds],2))*mesh.AE[s]/3
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
end # Demagntizing field with FEM/BEM

# Find new magnetization after time iteration
function timeStep(m::AbstractVector{Float64},H::AbstractVector{Float64},Hold::AbstractVector{Float64},
                  Heff::AbstractVector{Float64},
                  dt::Float64,
                  giro::Float64,damp::Float64=1.0,
                  precession::Float64=1.0)
    #=
        Repeats the search of a new magnetization until the solution is stable
    =#

    d::Float64 = dt*giro/2

    # m (n+1)
    m2::Vector{Float64} = zeros(3)

    # 1) Initial guess of the new magnetic field
    H12::Vector{Float64} = 3/2 *H - 0.5 *Hold

    # Find new m(n+1) until it doesn't change
    aux::Vector{Float64} = zeros(3)
    aux .= m

    err::Float64 = 1.0
    att::Int32 = 0
    mat::Matrix{Float64} = zeros(3,3)
    while err > 1e-6
        att += 1

        # 2) m (n+1) from m (n) and H (n+1/2)
        mat[1,1] = 1; mat[1,2] = d*H12[3]; mat[1,3] = -d*H12[2] 
        mat[2,1] = -d*H12[3]; mat[2,2] = 1; mat[2,3] = d*H12[1]
        mat[3,1] = d*H12[2];  mat[3,2] = -d*H12[1];  mat[3,3] = 1;

        m2 = mat\(m - d*cross(m,H12))

        # 3) m (n + 1/2)
        m12::Vector{Float64} = 0.5*(m + m2)

        # 4) Calculate H (n+1/2) from m (n+1/2)
        H12 = precession.*Heff + damp*cross(m12,Heff)

        # Max difference between m2 and previous m2
        err = maximum(abs.(m2-aux))
        # println(err)

        aux .= m2
        if att > 1_000
            println("Time step did not converge in ",att," steps")
            break
        end
    end

    # Update m directly and exit
    m .= m2

end # Find new magnetization after time iteration

function LandauLifshitz(mesh::MESH, m::Matrix{Float64}, Ms::Float64,
                        Hap::Vector{Float64}, Aexc::Float64, Aan::Float64,
                        uan::Vector{Float64}, scl::Float64, damp::Float64, giro::Float64,
                        A::Matrix{Float64}, LHS::Matrix{Float64}, Vn::Vector{Float64}, nodeVolume::Vector{Float64},
                        dt::Float64, precession::Float64, maxTorque::Float64,
                        maxAtt::Int32, totalTime::Float64=Inf, T::Float64=0.0)

    mu0::Float64 = pi*4e-7

    # Thermal constant
    Kb::Float64 = 1.380649e-23   # J/K, Boltzmann constant
    Cth::Float64 = sqrt(2*damp*Kb*T/(mu0*Ms*giro*dt*scl^3))
    rng = Dist.Normal(0.0,1.0) # Normal distribution for the thermal noise

    # -- Initial Magnetic Field --
    
    # Applied field
    Hext::Matrix{Float64} = zeros(3,mesh.nv) .+ mu0.*Hap

    # Demagnetizing field
    Hd::Matrix{Float64} = BEMdmag(mesh,m,LHS)

    # Exchange field
    Hexc::Matrix{Float64} = -2*Aexc.* (A * m')'

    # Thermal field
    Hth::Matrix{Float64} = zeros(3, mesh.nv)
    rand!(rng,Hth)
    Hth .*= Cth
    
    # Correct units of Demag, Exchange and Thermal fields
    for i in 1:3
        @view(Hd[i,:])     .*= mu0*Ms./Vn
        @view(Hexc[i,:])   ./= Ms*scl^2 .*nodeVolume
        @view(Hth[i,:])    ./= sqrt.(Vn)
    end

    # Anisotropy field
    Han::Matrix{Float64} = zeros(3,mesh.nv)
    for i in 1:mesh.nv
        Han[:,i] = 2*Aan/Ms *dot(m[:,i],uan) .*uan
    end

    # Effective field
    Heff::Matrix{Float64} = Hext + Hd + Hexc + Han + Hth
    H::Matrix{Float64} = zeros(3,mesh.nv)
    for i in 1:mesh.nv
        H[:,i] = precession*Heff[:,i] + damp*cross(m[:,i],Heff[:,i])
    end

    # -- Energy density --
    E::Float64 = 0.0        # Total
    Eext::Float64 = 0.0     # External field
    Ed::Float64 = 0.0       # Magnetostatic
    Eexc::Float64 = 0.0     # Exchange
    Ean::Float64 = 0.0      # Anisotropy

    for i in 1:mesh.nv
        Eext    -= mu0*Ms*dot(m[:,i],Hext[:,i])
        Ed      -= 0.5*mu0*Ms*dot(m[:,i],Hd[:,i])
        Eexc    -= 0.5*mu0*Ms*dot(m[:,i],Hexc[:,i])
        Ean     -= 0.5*mu0*Ms*dot(m[:,i],Han[:,i])
    end
    E = Eext + Ed + Eexc + Ean

    # Time iteration
    t::Float64 = 0.0 # Initial time (s)

    Hold::Matrix{Float64} = deepcopy(H)

    M_avg::Matrix{Float64} = zeros(3,maxAtt)
    E_time::Vector{Float64} = zeros(maxAtt)
    torque_time::Vector{Float64} = zeros(maxAtt)

    att::Int32 = 0
    while 1e9*t < totalTime && att < maxAtt
        t += dt
        att += 1

        # New magnetization
        for i in 1:mesh.nv
            timeStep(@view(m[:,i]), @view(H[:,i]), @view(Hold[:,i]),
                  @view(Heff[:,i]),
                  dt, giro, damp,
                  precession)
        end


        # -- New magnetic field --
        copyto!(Hold,H)    # Store the old magnetic field

        # Applied field | Is constant so don't update
        # Hext = zeros(3,mesh.nv) .+ mu0.*Hap

        # Demagnetizing field
        Hd = BEMdmag(mesh,m,LHS)

        # Exchange field
        Hexc = -2*Aexc.* (A * m')'

        # Thermal field
        rand!(rng,Hth)
        Hth .*= Cth

        # Correct units of Demag, Exchange and Thermal fields
        for i in 1:3
            @view(Hd[i,:])     .*= mu0*Ms./Vn
            @view(Hexc[i,:])   ./= Ms*scl^2 .*nodeVolume
            @view(Hth[i,:])    ./= sqrt.(Vn)
        end

        # Anisotropy field
        for i in 1:mesh.nv
            Han[:,i] = 2*Aan/Ms *dot(m[:,i],uan) .*uan
        end

        # Effective field
        Heff = Hext + Hd + Hexc + Han + Hth
        for i in 1:mesh.nv
            H[:,i] = precession*Heff[:,i] + damp*cross(m[:,i],Heff[:,i])
        end

        # -- Energy density --
        Eext = 0.0      # External field
        Ed = 0.0        # Magnetostatic
        Eexc = 0.0      # Exchange
        Ean = 0.0       # Anisotropy
        for i in 1:mesh.nv
            Eext    -= mu0*Ms*dot(m[:,i],Hext[:,i])
            Ed      -= 0.5*mu0*Ms*dot(m[:,i],Hd[:,i])
            Eexc    -= 0.5*mu0*Ms*dot(m[:,i],Hexc[:,i])
            Ean     -= 0.5*mu0*Ms*dot(m[:,i],Han[:,i])
        end
        E = Eext + Ed + Eexc + Ean
        E_time[att] = E

        # <|dm/dt|> , "Torque"
        dtau::Float64 = 0.0
        for i in 1:mesh.nv
            dtau += norm(cross(m[:,i],Heff[:,i]))
        end
        dtau /= mesh.nv
        torque_time[att] = dtau

        # Average magnetization
        # M_avg[:,att] = mean(m,2)
        M_avg[:,att] .= sum(m, dims=2) ./ size(m, 2)

        # Check if <|dm/dt|> is less than maxTorque 
        if dtau < maxTorque
            println("Converged, returning from LandauLifshitz()")
            break
        end

        # Print <|dm/dt|> every n iteration
        if mod(att,100) < 1
            println("t (ns) = ", 1e9*t, " , <|dm/dt|> = ", dtau)
        end

    end # End of time iteration

    # Remove excess zeros
    M_avg = M_avg[:,1:att]
    E_time = E_time[1:att]
    torque_time = torque_time[1:att]

    return m, Heff, M_avg, E_time, torque_time

end