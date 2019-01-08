################################################################################
#
#   DEFINITION OF LUTTINGER TISZA SPIN CONFIGURATIONS
#   contain Array of spins and lattice on which spins are residing
#
#   FILE CONTAINS
#   - abstract type definition
#   - interface of abstract type
#   - concrete type
#   - interface implementation
#
################################################################################




################################################################################
#
#   ABSTRACT TYPE
#
################################################################################

# Definition of abstract type
abstract type AbstractLTSpinConfiguration{N,L} end

# export the abstract type
export AbstractLTSpinConfiguration



# Interface functions

# obtain the list of spins
function spins(
            ltsc :: LTSC
        ) :: Vector{Vector{Float64}} where {N,L,LTSC <: AbstractLTSpinConfiguration{N,L}}

    # throw an error as this has to implemented by concrete type
    error("Interface function 'spins' not implemented yet for concrete type " * string(LTSC))
end
# obtain a single spin
function spin(
            ltsc  :: LTSC,
            index :: Integer
        ) :: Vector{Vector{Float64}} where {N,L,LTSC <: AbstractLTSpinConfiguration{N,L}}

    # refer to the main interface function
    return spins(ltsc)[index]
end

# obtain the lattice
function lattice(
            ltsc :: LTSC
        ) :: L where {N,L,LTSC <: AbstractLTSpinConfiguration{N,L}}

    # throw an error as this has to implemented by concrete type
    error("Interface function 'lattice' not implemented yet for concrete type " * string(LTSC))
end

# export the interface functions
export spins, spin
export lattice




################################################################################
#
#   CONCRETE TYPE
#
################################################################################

# Definition of concrete type
mutable struct LTSpinConfiguration{N,L} <: AbstractLTSpinConfiguration{N,L}

    # List of spins
    spins :: Vector{Vector{Float64}}

    # Lattice
    lattice :: L

end

# export the concrete type
export LTSpinConfiguration


# implementation of interface

# obtain the list of spins
function spins(
            ltsc :: LTSpinConfiguration{N,L}
        ) :: Vector{Vector{Float64}} where {N,L}

    # return the internal parameter
    return ltsc.spins
end

# obtain the lattice
function lattice(
            ltsc :: LTSpinConfiguration{N,L}
        ) :: L where {N,L}

    # return the internal parameter
    return ltsc.lattice
end
