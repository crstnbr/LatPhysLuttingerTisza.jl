################################################################################
#
#   CONSTRAINT RELATED
#
################################################################################

# function to compute the hard spin constraint
function constraintFunction(
    k :: Vector{Float64},
    eigenvectors :: Vector{Vector{Float64}},
    alpha :: Vector{Float64},
    sites :: Vector{S},
    d_spin :: Int64
    ) :: Float64 where {L<:Integer,D,S<:AbstractSite{L, D}}

    d_eig   = length(eigenvectors)
    N_sites = length(sites)

    sum = 0
    spin_vectors = Vector{Float64}[]

    # iterate over all sites in the lattice
    for r in 1 : N_sites

        site  = sites[r]
        index = label(site)
        pos   = point(site)

        spin_vector = zeros(Float64, d_spin)
        product     = dot(pos, k)

        # iterate over the (degenerate) subset of eigenvectors with minimal energy for vector k
        for s in 1 : d_eig
            spin_vector .+= (alpha[s] * cos(product) + alpha[s + d_eig] * sin(product)) * eigenvectors[s][d_spin * (index - 1) + 1 : d_spin * index]
        end

        # compute deviation from unit length
        sum += abs(norm(spin_vector) - 1.0)

    end

    # average over all sites
    return sum / N_sites

end

# function to compute the constraint for a given k on a finite lattice
function computeConstraint(
    k :: Vector{Float64},
    hamiltonian :: LTHamiltonian{L,U,HB},
    sites :: Vector{S},
    d_spin :: Int64
    ) :: Float64 where {L<:Integer,NS,U,HB<:AbstractBondHamiltonian{L,NS},D,S<:AbstractSite{L, D}}

    # compute the interaction matrix
    M_k = getMatrixAtK(hamiltonian, k)

    # diagonalize it
    eigenfactorization = eigen(M_k)
    eigenvalues        = real(eigenfactorization.values)
    eigenvectors       = real(eigenfactorization.vectors)

    # get (degenerate) subset of minimal energy eigenvectors
    E_min      = minimum(eigenvalues)
    deviations = [abs(E_min - E) for E in eigenvalues]
    degenerate = collect(1 : length(eigenvalues))[deviations .<= 1e-10]

    # minimize the constraint on the lattice
    alpha            = ones(2 * length(degenerate))
    constraint_value = Optim.minimum(Optim.optimize(x->constraintFunction(k, [eigenvectors[:, i] for i in degenerate], x, sites, d_spin), alpha, NelderMead(), Optim.Options(iterations = 10000)))

    return constraint_value

end








################################################################################
#
#   GENERAL FUNCTION
#
################################################################################


# function to compute the groundstate given a mesh of random k_vectors and connections between high symmetry points
function getLTGroundstateKSpace(
        N_random :: Int64,
        symmetric_k_vectors :: Vector{Vector{Float64}},
        hamiltonian :: LTHamiltonian{L,U,HB},
        sites :: Vector{S},
        ruc :: RU
        ;
        initial_tries :: Int64 = 10,
        d_spin :: Int64 = 3,
        epsilon :: Float64 = 1e-6,
        epsilon_k :: Float64 = 1e-10,
        groundstate_energy :: Float64 = Inf
    ) :: LTGroundstate where {L<:Integer,NS,U,HB<:AbstractBondHamiltonian{L,NS},D,S<:AbstractSite{L, D}, P,B,RU<:AbstractReciprocalUnitcell{P, B}}

    d_spatial    = length(point(sites[1]))
    bounds_lower = -2 * pi .* ones(Float64, d_spatial)
    bounds_upper = 2 * pi .* ones(Float64, d_spatial)

    # size of the symmetric mesh
    N_symmetric = length(symmetric_k_vectors)

    # check if the GS energy has to be calculated
    if groundstate_energy == Inf
        # search for 10 random points and compare the found minimum
        for tries in 1 : initial_tries
            # find a suitable starting point for the Newton algorithm
            k_start = zeros(Float64, d_spatial)
            for j in 1 : d_spatial
                comp = rand()
                k_start[j] = comp * bounds_lower[j] + (1 - comp) * bounds_upper[j]
            end
            # optimize the energy
            groundstate_energy = min(Optim.minimum(Optim.optimize(x->energy(hamiltonian, x), k_start)), groundstate_energy)
        end
        # print the groundstate_energy
        println("Groundstate energy is E_min = $(groundstate_energy)")
    end

    # initialize random points, minimize their energies and fold back to first BZ
    function dE(k :: Vector{Float64}) :: Float64
        M_k = getMatrixAtK(hamiltonian, k)
        eigenvalues = eigvals(M_k) .- groundstate_energy
        return minimum(eigenvalues .* eigenvalues)
    end

    random_k_vectors = zeros(Float64, N_random, d_spatial)
    index = 1
    while index <= N_random

        # initialize a random starting point
        k = zeros(Float64, d_spatial)
        for j in 1 : d_spatial
            comp = rand()
            k[j] = comp * bounds_lower[j] + (1 - comp) * bounds_upper[j]
        end

        # start with the initial energy
        e0 = dE(k)
         # iterate i over 100 newton steps (maximum)
        for i in 1 : 100
            # check if the energy is already converged
            if e0 < epsilon
                # fold back to 1.BZ and save the k vector
                random_k_vectors[index, :] = shiftToFirstBZ(ruc, k)
                # increment the index
                index += 1
                # break the newton loop
                break
            end
            # the current energy
            H_0 = e0
            # the gradient of the energy
            H_eps = zeros(Float64, d_spatial)
            for j in 1 : d_spatial
                shift = zeros(Float64, d_spatial)
                shift[j] = epsilon_k
                H_eps[j] = dE(k .+ shift)
            end
            dH = (H_eps .- H_0) ./ epsilon_k
            # absolute value of the gradient
            dHdH = dot(dH, dH)
            # break the Newton loop if the gradient diverges or is flattened too much
            if abs(dHdH) < 1e-20 || abs(dHdH) > 1e20
                break
            end
            # increment k
            dk = dH * (H_0 / dHdH)
            k -= dk
            # calculate a new energy
            e0 = dE(k)
        end
    end

    # check which points in the symmetric mesh also have minimal energy
    symmetric_k_vectors_min = Vector{Float64}[]
    for i in 1 : N_symmetric
        k = symmetric_k_vectors[i]
        if abs(energy(hamiltonian, k) - groundstate_energy) <= epsilon
            push!(symmetric_k_vectors_min, k)
        end
    end

    # iterate over all points and compute the contraint
    N_symmetric        = length(symmetric_k_vectors_min)
    N_total            = N_random + N_symmetric
    k_vectors          = zeros(Float64, N_total, d_spatial)
    constraint_values  = zeros(Float64, N_total)
    for i in 1 : N_total
        if i <= N_random
            k                    = random_k_vectors[i, :]
            constraint_values[i] = computeConstraint(k, hamiltonian, sites, d_spin)
            k_vectors[i, :]      = k
        else
            k                    = symmetric_k_vectors_min[i - N_random]
            constraint_values[i] = computeConstraint(k, hamiltonian, sites, d_spin)
            k_vectors[i, :]      = k
        end
    end

    # return LTGroundstate
    return LTGroundstate(groundstate_energy, k_vectors, constraint_values)

end
