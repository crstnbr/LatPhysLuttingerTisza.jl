################################################################################
#
#   definition of Luttinger-Tisza ground state in k-space
#
################################################################################

struct LTGroundstate

    # the minimal energy
    E_min :: Float64

    # the k_vectors with minimal energy, columns correspond to components
    k_vectors :: Matrix{Float64}

    # their respective constraint values
    constraint_values :: Vector{Float64}

end
