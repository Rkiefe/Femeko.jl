# Defines FEM related methdods and then includes the 2D and 3D FEM logic

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

# 2D Gradient function
# Input 2D matrix, 1D vector for the rows, 1D vector for the columns
function gradient(M::Matrix{Float64}, 
                  H::Vector{Float64}, 
                  T::Vector{Float64})
    
    dM_dT::Matrix{Float64} = zeros(size(M))
    dM_dH::Matrix{Float64} = zeros(size(M))

    for i in 1:length(H)
        dM_dT[i, :] = gradient(T, M[i, :]) # For each field H', calculate dM/dT (H', T)
            # H, T
    end

    for i in 1:length(T)
        dM_dH[:, i] = gradient(H, M[:, i]) # For each temperature T', calculate dM/dH (H, T')
            # H, T
    end

    return dM_dH, dM_dT
end

# 1D numerical integration | trapezoidal integration
function trapz(xin, yin)

    # Apply the trapezoidal numerical integral
    result::Float64 = 0.0
    for i in 1:length(xin)-1
        result += 0.5*(xin[i+1]-xin[i])*(yin[i] + yin[i+1])
    end

    return result
end

# Check if a matrix is approximately symmetric
function isapproxsymmetric(A, tol=1e-10)
    return all(isapprox.(A - A', 0.0; rtol=tol))
end

# Include FEM logic
include("FEM2D.jl")
include("FEM3D.jl")


