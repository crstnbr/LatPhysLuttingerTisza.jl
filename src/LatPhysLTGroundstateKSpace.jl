################################################################################
#
#   definition of Luttinger-Tisza ground state in k-space
#
################################################################################

# LA seems crucial, but is already included in base
# using LinearAlgebra 

struct LTGroundstate

    # the minimal energy modes 
    Modes :: Vector{Vector{Float64}}

    # their respective constraint values
    Constraints :: Vector{Float64} 

end

# function to compute energy from LT interaction matrix
function Energy(
    hamiltonian :: LTHamiltonian{L,U,HB},
    k :: Vector{<:Real}
) :: Float64 where {L,NS,U,HB<:AbstractBondHamiltonian{L,NS}}

    # first compute the matrix at the given k 
    M_k = getMatrixAtK(hamiltonian, k)

    # compute the eigenvalues 
    Eigvals = eigvals(M_k)

    # return minimal value 
    return minimum(Eigvals)

end



