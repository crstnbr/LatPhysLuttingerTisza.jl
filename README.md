# LatPhysLuttingerTisza.jl [![pipeline status](http://gitsrv.thp.uni-koeln.de/attig/LatPhysLuttingerTisza.jl/badges/master/pipeline.svg)](http://gitsrv.thp.uni-koeln.de/attig/LatPhysLuttingerTisza.jl/commits/master) [![coverage report](http://gitsrv.thp.uni-koeln.de/attig/LatPhysLuttingerTisza.jl/badges/master/coverage.svg)](http://gitsrv.thp.uni-koeln.de/attig/LatPhysLuttingerTisza.jl/commits/master)

Luttinger Tisza calculations for [`LatticePhysics.jl`](http://gitsrv.thp.uni-koeln.de/attig/LatticePhysics.jl)



## Contents

... TODO




## Installation (usage only):

For usage purposes only, you can install the package via the package mode in Julia (Pkg). However, since the package
is not listed in the Julia package repositories, you have to use
```julia
(v1.0) pkg> add "git@gitsrv.thp.uni-koeln.de:attig/LatPhysLuttingerTisza.jl.git"
```
Note: this can lead to Errors under Windows 10 due to incorrect SSH access. Use the following command instead:
```julia
(v1.0) pkg> add "http://gitsrv.thp.uni-koeln.de/attig/LatPhysLuttingerTisza.jl.git"
```
You will be prompted a username and password validation but it should work the same way.


## Installation (developement):

For developement purposes, it is best to clone the package via git to a developement
git location of your choice and use
```julia
(v1.0) pkg> dev "path/to/the/repository/on/your/machine"
```

Alternatively, you could use
```julia
(v1.0) pkg> dev "git@gitsrv.thp.uni-koeln.de:attig/LatPhysLuttingerTisza.jl.git"
```
or (on Windows)
```julia
(v1.0) pkg> dev "http://gitsrv.thp.uni-koeln.de/attig/LatPhysLuttingerTisza.jl.git"
```
to clone a development version of the package to `~/.julia/dev/`.


Finally, develope the package as you are used to within the editor of your choice.
