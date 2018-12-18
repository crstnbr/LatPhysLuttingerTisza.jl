################################################################################
#
#   definition of Luttinger-Tisza ground state in k-space
#
################################################################################

struct LTGroundstate

    # the minimal energy modes 
    Modes :: Vector{Vector{Float64}}

    # their respective constraint values
    Constraints :: Vector{Float64} 

end