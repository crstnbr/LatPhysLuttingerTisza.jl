################################################################################
#
#   DEFINITION OF LUTTINGER TISZA STATES AND GROUNDSTATES
#
#   FILE CONTAINS
#   1) Luttinger Tisza states
#       - abstract type definition
#       - interface function definition
#       - concrete type imlementation
#       - interface implementation
#   2) Luttinger Tisza ground states and ground state manifolds
#       - abstract type definition
#       - ground states in real space
#       - ground state manifolds in k space
#
################################################################################




abstract type AbstractLTState end

function eigenvectors(lts) end
function alphas(lts) end
function spinconfiguration(lts, lattice, u=[1,0,0], v=[0,1,0]) :: LTSpinConfiguration{N,L} end

mutable struct LTState <: AbstractLTState
    eigenvectors :: Vector{Vector{Float64}}
    alphas :: Vector{Complex}
end








abstract type AbstractLTGroundstate{LTS <: AbstractLTState} end

function spinconfiguration(ltgs, index=1) end



mutable struct LTGroundstateRealSpace{LTS,L<:AbstractLattice{S,B,U} where {S,B,U},H<:AbstractBondHamiltonian{L,N} where {L,N}} <: AbstractLTGroundstate{LTS}

    # lattice
    lattice :: L
    # bond hamiltonian
    h_bond :: H

    # the minimal energy
    E_min :: Float64

    # the degenerate LT States
    lt_states :: Vector{LTS}

    # their respective constraint values
    constraint_values :: Vector{Float64}

end



mutable struct LTGroundstateKSpace{LTS,U<:AbstractUnitcell{S,B} where {S,B},H<:AbstractBondHamiltonian{L,N} where {L,N}} <: AbstractLTGroundstate{LTS}

    # unitcell
    unitcell :: U
    # bond hamiltonian
    h_bond :: H

    # the minimal energy
    E_min :: Float64

    # the k_vectors with minimal energy, columns correspond to components
    k_vectors :: Vector{Vector{Float64}}
    # the corresponding LT States
    lt_states :: Vector{LTS}

    # their respective constraint values
    constraint_values :: Vector{Float64}

end
