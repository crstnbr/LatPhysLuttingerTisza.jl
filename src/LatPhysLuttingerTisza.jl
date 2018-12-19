################################################################################
#
#   module LatPhysLuttingerTisza
#   -> LatticePhysics_Base
#   -> LatticePhysics_ReciprocalSpace
#   -> LinearAlgebra
#
#   --> TYPE DEFINITIONS
#           - LTBandstructure
#           - LTGroundstateManifold
#
#   --> IO OF TYPES
#
#   --> SPIN INTERACTION MATRICES
#
#   --> CALCULATION OF LT CONSTRAINTS
#
#   --> CALCULATION OF LT BANDSTRUCTURES AND EIGENVALUES
#           - in real space (pure eigenvalues)
#           - in momentum space (along cardinal directions)
#           - in momentum space (along reciprocal lattice directions)
#           - in momentum space (along path)
#
#   --> CALCULATION OF LT GROUNDSTATE MANIFOLDS
#           - calculation for 2D unitcells
#           - calculation for 3D unitcells
#
################################################################################






# start of module
module LatPhysLuttingerTisza

# used libraries
using LatPhysBase
using LatPhysReciprocal
using LinearAlgebra

# Definition of Hamiltonians (only for calculations)
include("LT_hamiltonian_definition.jl")
include("LT_groundstate_k_space.jl")
include("LT_groundstate_plotting.jl")



# end of module
end
