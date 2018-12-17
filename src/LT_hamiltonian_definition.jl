struct LTHamiltonian{
        L,
        U<:AbstractUnitcell{S,B} where {S, B<:AbstractBond{L,N} where {N}},
        HB<:AbstractBondHamiltonian{L,NS} where {NS}
    }

    # the unitcell
    unitcell :: U

    # the bond Hamiltonian
    h_bond   :: HB

    # ...

end


# function to create such an object
function newLTHamiltonian(
            unitcell :: U,
            bond_hamiltonian :: HB
        ) :: LTHamiltonian{L,U,HB} where {L,U<:AbstractUnitcell{S,B} where {S, B<:AbstractBond{L,N} where {N}},HB<:AbstractBondHamiltonian{L,NS} where {NS}}

    # return such an object
    return LTHamiltonian{L,U,HB}(unitcell, bond_hamiltonian)
end

# function to obtain the Hamiltonian matrix at some k
function getHamiltonianAtK(
        hamiltonian :: LTHamiltonian{L,U,HB},
        k :: Vector{<:Real}
    ) :: Matrix{Complex} where {L,U,HB}

end
