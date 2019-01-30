# LatPhysLuttingerTisza.jl

Luttinger Tisza calculations for [`LatticePhysics.jl`](https://github.com/janattig/LatticePhysics.jl).



## Contents

Calculation of Luttinger Tisza ground states for systems defined by bond Hamiltonians on unitcells (both defined in [`LatPhysBase.jl`](https://github.com/janattig/LatPhysBase.jl)).




## Installation

You can install the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to first install the unregistered
dependencies manually with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysBase.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysLatticeConstruction.jl"
(v1.0) pkg> add "https://github.com/janattig/LatPhysReciprocal.jl"
```
to finally install the main package with
```julia
(v1.0) pkg> add "https://github.com/janattig/LatPhysLuttingerTisza.jl"
```
