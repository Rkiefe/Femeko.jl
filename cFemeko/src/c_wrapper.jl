# Wrapper to the C++ code, to be called in Julia



# Local Stiffness matrix
function CstiffnessMatrix(p::Matrix{Float64}, t::Matrix{Int32}, VE::Vector{Float64}, mu::Vector{Float64})
    # Get the matrix dimensions
    nt::Int32 = size(t,2)
    nv::Int32 = size(p,2)

    # Ready the output
    Ak::Matrix{Float64} = zeros(16,nt)
    
    # Call C++ function
    @ccall "FEMc.so".stiffnessMatrix(Ak::Ptr{Float64}, p::Ptr{Float64}, t::Ptr{Int32}, nv::Int32, nt::Int32, VE::Ptr{Float64}, mu::Ptr{Float64})::Cvoid

    return Ak
end # Wrapper for C++ function to get local stiffness matrix 