# Handles magnetic materials and their datasets

mutable struct DATA
    # Magnetization data
    M # ::Matrix{Float64}

    # Magnetic field H
    HofM::Vector{Float64}
    
    # Temperature
    TofM::Vector{Float64}
    
    # Magnetic Flux
    B::Vector{Float64}

    # Density
    rho::Float64

    # Permeability mu = B/H
    mu::Vector{Float64}

    # d/dH mu (derivative of permeability)
    dmu::Vector{Float64}

    # Constructor
    DATA() = new()
end

function loadMaterial( materialProperties,      # Dict
                       folder::String,          # Folder with materials
                       data::String,            # Data folder of target material
                       key::String,             # Material name
                       density::Float64,        # Density (g/cm3)
                       T::Float64)              # Temperature (K)
    
    # vacuum magnetic permeability
    mu0 = pi*4e-7

    # Load material properties
    materialProperties[key].M = readdlm(folder*"/"*data*"/M.dat")                 # emu/g
    materialProperties[key].HofM = vec(readdlm(folder*"/"*data*"/HofM.dat"))      # Oe
    materialProperties[key].TofM = vec(readdlm(folder*"/"*data*"/TofM.dat"))      # K
    materialProperties[key].rho = density # g/cm3

    # Convert data units
    materialProperties[key].M .*= materialProperties[key].rho*1e3 # A/m
    materialProperties[key].HofM .*= 1e-4/mu0                      # A/m

    # Interpolate data over the target temperature
    spl = Spline2D( materialProperties[key].HofM,
                    materialProperties[key].TofM,
                    materialProperties[key].M)

    M::Vector{Float64} = zeros(length(materialProperties[key].HofM))
    for i in 1:length(M)
        M[i] = spl(materialProperties[key].HofM[i], T)
    end
    materialProperties[key].M = M
    
    # Magnetic flux density
    materialProperties[key].B = mu0.*(materialProperties[key].HofM .+
                                       materialProperties[key].M)

    # Get the permeability and its derivative
    materialPermeability(materialProperties, key)

end

function materialPermeability(materialProperties, key::String)

    # Sets the permeability of the material from the dataset
    # And the derivate of the permeability for the Newton Rapshon methods
    # And cleans NaN and Inf

    # Permeability
    materialProperties[key].mu = materialProperties[key].B./materialProperties[key].HofM
    
    # Remove Inf
    idx = findall(x -> !isfinite(x), materialProperties[key].mu)
    materialProperties[key].mu[idx] .= 0.0
    materialProperties[key].mu[idx] .= maximum(materialProperties[key].mu)

    # d/dH mu
    dmu = gradient(materialProperties[key].HofM, 
                   materialProperties[key].B./materialProperties[key].HofM)

    # Remove -Inf
    idx = findall(x -> !isfinite(x), dmu)
    dmu[idx] .= 0.0
    dmu[idx] .= minimum(dmu)

    materialProperties[key].dmu = dmu
end
