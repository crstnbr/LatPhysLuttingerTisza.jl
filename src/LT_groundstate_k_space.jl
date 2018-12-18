################################################################################
#
#   definition of Luttinger-Tisza ground state in k-space
#
################################################################################

struct LTGroundstate

    # the minimal energy modes 
    modes :: Vector{Vector{Float64}}

    # their respective constraint values
    constraints :: Vector{Float64} 

end

# function to compute energy from LT interaction matrix
function energy(
    hamiltonian :: LTHamiltonian{L,U,HB},
    k :: Vector{Float64}
) :: Float64 where {L,NS,U,HB<:AbstractBondHamiltonian{L,NS}}

    # first compute the matrix at the given k 
    M_k = getMatrixAtK(hamiltonian, k)

    # return minimal eigenvalue 
    return eigmin(M_k)

end 

# function to compute deviation
function deviation(
    k :: Vector{Float64},
    eigenvectors :: Vector{Vector{Complex}},
    sites :: Vector{S},
    d_spin :: Int64,
    alpha :: Vector{Float64}
) :: Float64 where {L, D, S <: AbstractSite{L, D}}

    d_eig   = length(eigenvectors)
    N_sites = length(sites)

    sum = 0

    for r in 1 : N_sites

        site  = sites[r]
        label = site.label 
        pos   = site.point

        spin_vector = zeros(Float64, d_eig)
        product     = dot(pos, k)

        for s in 1 : d_eig
            spin_vector .+= (alpha[s] * cos(product) + alpha[s + d_eig] * sin(product)) * eigenvectors[s][d_spin * (label - 1) + 1 : d_spin * label]
        end

        sum += abs(norm(spin_vector) - 1)

    end

    return sum / N_sites

end

# function to minimize the constraint
function getLTConstraint(
    k :: Vector{Float64},
    eigenvectors :: Vector{Vector{Complex}},
    sites :: Vector{S}
) :: Float64 where {L, D, S <: AbstractSite{L, D}}

    alpha = zeros(2 * length(eigenvectors))
    constraint = Optim.minimum(Optim.optimize(x->deviation(k, eigenvectors, sites, d_spin, x), alpha))

    return constraint

end

function getLTConstraintDetailed(spin_eigenvectors::Array{Array{Complex{Float64},1},1}, spin_dimension::Int64, lattice_vectors::Array{Array{Float64,1},1}, momentum_vector::Array{Float64, 1})

    full_size = length(spin_eigenvectors[1])
    eigenspace_dimension = length(spin_eigenvectors)
    alpha = zeros(2 * eigenspace_dimension)
    res = Optim.optimize(x->spatial_deviation(spin_eigenvectors, full_size, spin_dimension, eigenspace_dimension, lattice_vectors, momentum_vector, x), alpha, NelderMead(),
    Optim.Options(iterations = 10000))

    constraint_value = Optim.minimum(res)
    parameters = Optim.minimizer(res)

    return constraint_value, parameters

end

function ComputeConstraint(
    momentum_vector::Array{Float64, 1},
    unitcell::Unitcell,
    lattice_vectors::Array{Array{Float64,1},1},
    spin_dimension::Int64,
    epsilon_degenerate::Float64,
    bondInteractionMatrix::Function=getBondInteractionMatrixHeisenberg)

    LT_constraint = 0

    matrix = getSpinInteractionMatrixKSpace(unitcell, momentum_vector, bondInteractionMatrix)

    if size(matrix) == (1,1)

        return LT_constraint

    else

        eigenfactorization = eigfact(matrix)
        eigenvalues  = eigenfactorization[:values]
        eigenvectors = eigenfactorization[:vectors]

        degenerate = zeros(Int64,   length(eigenvalues)) .- 1
        treat      =  ones(Int64,   length(eigenvalues))

        for b in 2:length(eigenvalues)

            if eigenvalues[b] - epsilon_degenerate <= eigenvalues[b-1]

                degenerate[b-1] = b
                treat[b] = 0

            end

        end

        treated_first = false
        for b in 1:length(eigenvalues)

            if treat[b] == 0

                continue

            elseif treated_first

                break

            end

            degenerate_bands = Int64[b]
            while degenerate[degenerate_bands[end]] != -1

                push!(degenerate_bands, degenerate[degenerate_bands[end]])

            end

            LT_constraint = getLTConstraint([eigenvectors[:,j] for j in degenerate_bands], spin_dimension, lattice_vectors, momentum_vector)
            treated_first = true

        end

        return LT_constraint

    end

end



