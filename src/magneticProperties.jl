#=
    Handles magnetic materials and their datasets
    Requires FEM.jl for the gradient and trapz functions
=#
# To read data from files
using DelimitedFiles

# Wrapper to Fortran dierckx | Interpolation functions
using Dierckx

mutable struct DATA
    # Magnetization data
    M # A/m

    # Magnetic field H
    HofM::Vector{Float64} # A/m
    
    # Temperature
    TofM::Vector{Float64} # K
    
    # Magnetic Flux
    B::Vector{Float64}  # T

    # Density
    rho::Float64 # g/cm3

    # Permeability mu = B/H
    mu::Vector{Float64}

    # d/dH mu (derivative of permeability)
    dmu::Vector{Float64}

    # Constructor
    DATA() = new()
end

function loadMaterial( materialProperties,      # Dict or DATA
                       folder::String,          # Folder with materials
                       data::String,            # Data folder of target material
                       key::String,             # Material name
                       density::Float64,        # Density (g/cm3)
                       T::Float64)              # Temperature (K)
    
    # vacuum magnetic permeability
    mu0 = pi*4e-7

    # If the user just wants the DATA struct, and not update a hashtable:
    if typeof(materialProperties) == DATA
        # Update material properties directly

        # Load material properties
        materialProperties.M = readdlm(folder*"/"*data*"/M.dat")                 # emu/g
        materialProperties.HofM = vec(readdlm(folder*"/"*data*"/HofM.dat"))      # Oe
        materialProperties.TofM = vec(readdlm(folder*"/"*data*"/TofM.dat"))      # K
        materialProperties.rho = density # g/cm3

        # Convert data units
        materialProperties.M .*= materialProperties.rho*1e3 # A/m
        materialProperties.HofM .*= 1e-4/mu0                      # A/m

        # Interpolate data over the target temperature
        spl = Spline2D( materialProperties.HofM,
                        materialProperties.TofM,
                        materialProperties.M)

        M = zeros(length(materialProperties.HofM))::Vector{Float64}
        for i in 1:length(M)
            M[i] = spl(materialProperties.HofM[i], T)
        end
        materialProperties.M = M
        
        # Magnetic flux density
        materialProperties.B = mu0.*(materialProperties.HofM .+
                                     materialProperties.M)

        # Get the permeability and its derivative
        materialPermeability(materialProperties)

        return
    end

    # Else, update the hashtable

    # Load material properties
    materialProperties[key].M = readdlm(folder*"/"*data*"/M.dat")                 # emu/g
    materialProperties[key].HofM = vec(readdlm(folder*"/"*data*"/HofM.dat"))      # Oe
    materialProperties[key].TofM = vec(readdlm(folder*"/"*data*"/TofM.dat"))      # K
    materialProperties[key].rho = density # g/cm3

    # Convert data units
    materialProperties[key].M .*= materialProperties[key].rho*1e3  # A/m
    materialProperties[key].HofM .*= 1e-4/mu0                      # A/m

    # Interpolate data over the target temperature
    spl = Spline2D( materialProperties[key].HofM,
                    materialProperties[key].TofM,
                    materialProperties[key].M)

    M = zeros(length(materialProperties[key].HofM))::Vector{Float64}
    for i in 1:length(M)
        M[i] = spl(materialProperties[key].HofM[i], T)
    end
    materialProperties[key].M = M
    
    # Magnetic flux density
    materialProperties[key].B = mu0.*(materialProperties[key].HofM .+
                                       materialProperties[key].M)

    # Get the permeability and its derivative
    materialPermeability(materialProperties[key])

end

# Magnetic permeability from dataset
function materialPermeability(data::DATA)

    # Sets the permeability of the material from the dataset
    # And the derivate of the permeability for the Newton Rapshon methods
    # And cleans NaN and Inf

    # Permeability
    data.mu = data.B./data.HofM
    
    # Remove Inf
    idx = findall(x -> !isfinite(x), data.mu)
    data.mu[idx] .= 0.0
    data.mu[idx] .= maximum(data.mu)

    # d/dH mu
    dmu = gradient(data.HofM, 
                   data.mu)

    # Remove -Inf
    idx = findall(x -> !isfinite(x), dmu)
    dmu[idx] .= 0.0
    dmu[idx] .= minimum(dmu)

    data.dmu = dmu
end

# Get magnetic entropy change from magnetization data
function deltaS(data::DATA, 
                Mdata::Matrix{Float64}, # emu/g
                T::Float64,  # K
                H0::Float64, # A/m
                Hf::Float64) # A/m
    #=
        . Mdata is the original dataset, M(H,T) [emu/g]
        . data is the processed dataset with all the info such as 
        the temperature span, the magnetic field range, permeability, etc
        
            dS(H, T) = integral_0^H dM/dT (H', T) dH'

        Can calculate delta S for any limit H0, Hf, where Hf > H0
    =#
    
    # For now, just set a sign to the integral
    sign::Float64 = 1.0
    if Hf < H0 # Switch the limits and invert the sign
        sign = -1.0
        aux = H0
        H0 = Hf
        Hf = aux # H0
    end

    mu0::Float64 = pi*4e-7

    # Set the integral range
    Hrange::Vector{Float64} = [H0]
    Hi::Int32 = 0
    for (i, H) in enumerate(data.HofM)
        
        # Set the start of the range
        if Hi == 0 && H > H0
            Hi = i # The first index where H > H0
        end

        if H >= Hf
            Hrange = vcat(Hrange, data.HofM[Hi:i-1])
            push!(Hrange, Hf)
            break
        end
    end

    # Gradient over the temperature for each magnetic field
    _, dM_dT = gradient(Mdata, data.HofM, data.TofM) # emu/g/K

    # Interpolate for the target temperature, over the magnetic field range
    dM_dT_spline = Spline2D(data.HofM, data.TofM, dM_dT)
    
    dM::Vector{Float64} = zeros(length(Hrange)) # emu/g/K
    for (iH, H) in enumerate(Hrange)
        dM[iH] = dM_dT_spline(H, T) # emu/g/K
    end

    # Entropy change
    dS::Float64 = sign * trapz(Hrange.*mu0, dM)

    return dS
end

# Magnetostatic energy
function getEnergy(mesh::MESH, 
                   data::DATA,
                   H::Vector{Float64},
                   B::Vector{Float64},
                   InsideOnly::Bool = false)

    # Assumes all the materials have the same properties
    # Calculates the magnetostatic energy following the same approach
    # as in FEMM 
    
    # Magnetic permeability in a vacuum
    mu0::Float64 = pi*4e-7

    # Magnetostatic Energy | Non Linear materials, following FEMM
    energy::Float64 = 0.0

    # Energy in free space
    if !InsideOnly
        for k in setdiff(1:mesh.nt, mesh.InsideElements)
            energyDensity::Float64 = 0.5*mu0*H[k]^2
            energy += energyDensity*mesh.VE[k]
        end
    end

    # Energy inside magnetic material
    for k in mesh.InsideElements
        
        # Upper limit of the integral
        Bq::Float64 = B[k]

        # Set the value of H at the integral limit
        Hq::Float64 = interp1(data.B, 
                              data.HofM, 
                              Bq)

        # find last index in B before Bq
        idx = 0
        for (i, v) in enumerate(data.B)
            if v > Bq
                idx = i - 1
                break
            end
        end

        # Set the data from B[1] to Bq and H[1] to Hq 
        xin::Vector{Float64} = [data.B[1:idx]; Bq]
        yin::Vector{Float64} = [data.HofM[1:idx]; Hq]

        # Calculate the energy density for this element
        energyDensity::Float64 = trapz(xin, yin)

        # Multiply by the volume to get the energy
        energy += mesh.VE[k]*energyDensity 

    end

    return energy
end

# Plot magnetization Vs magnetic field
function plotData(data::DATA)
    mu0::Float64 = pi*4e-7

    fig, ax, sct = scatter( mu0.*data.HofM,
                          data.M./(data.rho*1e3))
    ax.xlabel = "mu_0 H (T)"
    ax.ylabel = "M (emu/g)"
    # ax.title  = ""

    return fig, ax, sct
end