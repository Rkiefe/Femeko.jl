#=
    3D Finite element implementation of the Landau-Lifshitz equation with a semi-implicit
    time step, based on https://doi.org/10.1109/TMAG.2008.2001666 (Oriano 2008)

    The energy is minimized by a steepest descent algorithm from 

    https://doi.org/10.1063/1.4862839
                and
    https://doi.org/10.1063/1.4896360

=#
include("LandauLifshitz.jl")

# Next step in magnetization by steepest descent
function nextM(M::Vector{Float64},Heff::Vector{Float64},dt::Float64)

    # Semi-implicit time step
    d::Float64 = dt/2;
    h12::Vector{Float64} = cross(M,Heff)

    mat::Matrix{Float64} = [1 d*h12[3] -d*h12[2];
                           -d*h12[3] 1 d*h12[1];
                           d*h12[2] -d*h12[1] 1]

    M .= mat\(M-d*cross(M,h12))
end # New magnetization with steepest descent


function SteepestDescent(mesh::MESH, m::Matrix{Float64}, Ms::Float64, Heff::Matrix{Float64},
                        Hap::Vector{Float64}, Aexc::Float64, Aan::Float64,
                        uan::Vector{Float64}, scl::Float64,
                        A::Matrix{Float64}, LHS::Matrix{Float64}, Vn::Vector{Float64}, nodeVolume::Vector{Float64},
                        maxTorque::Float64, maxAtt::Int32)

    mu0::Float64 = pi*4e-7

    # Applied field

    # -- Initial Magnetic Field --
    Hext::Matrix{Float64} = zeros(3,mesh.nv) .+ mu0.*Hap
    Hd::Matrix{Float64} = zeros(3,mesh.nv)
    Hexc::Matrix{Float64} = zeros(3,mesh.nv)
    Han::Matrix{Float64} = zeros(3,mesh.nv)

    if isempty(Heff) # Get the initial effective field
        # Demagnetizing field
        Hd = BEMdmag(mesh,m,LHS)

        # Exchange field
        Hexc = -2*Aexc.* (A * m')'

        # Correct units of Demag and Exchange fields
        for i in 1:3
            Hd[i,:]     .*= mu0*Ms./Vn
            Hexc[i,:]   ./= Ms*scl^2 .*nodeVolume
        end

        # Anisotropy field
        Han = zeros(3,mesh.nv)
        for i in 1:mesh.nv
            Han[:,i] = 2*Aan/Ms *dot(m[:,i],uan) .*uan
        end

        # Effective field
        Heff = Hext + Hd + Hexc + Han

    # else
        # Use input effective field
    end

    H::Matrix{Float64} = zeros(3,mesh.nv)
    for i in 1:mesh.nv
        H[:,i] = cross(m[:,i],Heff[:,i])
    end

    # New magnetization
    mOld::Matrix{Float64} = deepcopy(m)
    for i in 1:mesh.nv
        m[:,i] = timeStep(m[:,i], H[:,i], H[:,i], Heff[:,i],
                          0.03, 1.0, 1.0, 0.0)
    end

    HeffOld::Matrix{Float64} = zeros(3,mesh.nv)
    Hold::Matrix{Float64}    = zeros(3,mesh.nv) 

    # -- Energy density --
    E::Float64 = 0.0        # Total
    Eext::Float64 = 0.0     # External field
    Ed::Float64 = 0.0       # Magnetostatic
    Eexc::Float64 = 0.0     # Exchange
    Ean::Float64 = 0.0      # Anisotropy


    # Energy minimization by steepest descent
    M_avg::Matrix{Float64} = zeros(3,maxAtt)
    E_time::Vector{Float64} = zeros(maxAtt)
    torque_time::Vector{Float64} = zeros(maxAtt)

    dtau::Float64 = 2*maxTorque
    att::Int32 = 0
    while dtau > maxTorque && att < maxAtt
        att += 1

        HeffOld = deepcopy(Heff)
        Hold = deepcopy(H)

        # -- New magnetic field --
        
        # Demagnetizing field
        Hd = BEMdmag(mesh,m,LHS)

        # Exchange field
        Hexc = -2*Aexc.* (A * m')'

        # Correct units of Demag and Exchange fields
        for i in 1:3
            Hd[i,:]     .*= mu0*Ms./Vn
            Hexc[i,:]   ./= Ms*scl^2 .*nodeVolume
        end

        # Anisotropy field
        Han = zeros(3,mesh.nv)
        for i in 1:mesh.nv
            Han[:,i] = 2*Aan/Ms *dot(m[:,i],uan) .*uan
        end

        # Effective field
        Heff = Hext + Hd + Hexc + Han

        # H = m cross Heff
        for i in 1:mesh.nv
            H[:,i] = cross(m[:,i],Heff[:,i])
        end

        # -- Energy density --
        E = 0.0     # Total
        Eext = 0.0  # External field
        Ed = 0.0    # Magnetostatic
        Eexc = 0.0  # Exchange
        Ean = 0.0   # Anisotropy
        for i in 1:mesh.nv
            Eext    -= mu0*Ms*dot(m[:,i],Hext[:,i])
            Ed      -= 0.5*mu0*Ms*dot(m[:,i],Hd[:,i])
            Eexc    -= 0.5*mu0*Ms*dot(m[:,i],Hexc[:,i])
            Ean     -= 0.5*mu0*Ms*dot(m[:,i],Han[:,i])
        end
        E = Eext + Ed + Eexc + Ean
        E_time[att] = Eext + Ed + Eexc + Ean

        snN::Float64 = 0.0
        snD::Float64  = 0.0
        snD2::Float64 = 0.0
        sn::Vector{Float64} = zeros(3)
        gn1::Vector{Float64} = zeros(3)
        gn2::Vector{Float64} = zeros(3)

        for i in 1:mesh.nv
            sn = m[:,i] - mOld[:,i]
            
            gn2 = cross(m[:,i],cross(m[:,i],Heff[:,i]))
            gn1 = cross(mOld[:,i],cross(mOld[:,i],HeffOld[:,i]))

            snN  += dot(sn,sn)
            snD  += dot(sn,gn2-gn1)
            snD2 += dot(gn2-gn1,gn2-gn1)
        end

        tau1::Float64 = snN/snD
        tau2::Float64 = snD/snD2

        dt::Float64 = mod(att,2) > 0 ? tau1 : tau2

        # New magnetization
        for i in 1:mesh.nv
            mOld[:,i] = m[:,i]
            m[:,i] = nextM(m[:,i],Heff[:,i],dt)
        end

        # Average magnetization
        M_avg[:,att] = mean(m,2)

        # <|dm/dt|> , "Torque"
        dtau = 0.0
        for i in 1:mesh.nv
            dtau += norm(cross(m[:,i],Heff[:,i]))
        end
        dtau /= mesh.nv
        torque_time[att] = dtau

        # Print <|dm/dt|> every n iteration
        if mod(att,100) < 1
            println(att, ": <|dm/dt|> = ", dtau)
        end

    end # End of energy minimization

    # Remove excess zeros
    M_avg = M_avg[:,1:att]
    E_time = E_time[1:att]
    torque_time = torque_time[1:att]

    return m, Heff, M_avg, E_time, torque_time

end


