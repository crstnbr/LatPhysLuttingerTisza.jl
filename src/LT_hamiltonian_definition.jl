struct LTHamiltonian{
        L,
        U<:AbstractUnitcell{S,B} where {S, B<:AbstractBond{L,N} where {N}},
        HB<:AbstractBondHamiltonian{L,NS} where {NS}
    }

    # the bond Hamiltonian
    h_bond :: HB

    # ...

end


# function to create such an object



# function to obtain the Hamiltonian matrix at some k
function getHamiltonianAtK(
        hamiltonian :: LTHamiltonian{L,U,HB},
        k :: Vector{<:Real}
    ) :: Matrix{Complex} where {L,U,HB}

end
