## Sources for PO project

The application is composed from the followings classes :

- [SolverNS2](solverns2.hpp) : call the other classes to initialize the objects needed
- [SolverA](solvera.hpp) : solves the problems for the boundary conditions
- [SolverEigenNS2](solvereigenns2.hpp) : solves the eigen problem
- [SolverSpectralProblem](solverspectralproblem.hpp) : solves the spectral problem

### Options

All the options are prefixed by `solverns2.`
Here, a list of the options :
* `verbose` : level of verbosity

#### SolverEigenNS2 options
* `needEigen` : need eigenmodes
* `markerList` : list of markers of the boundary
* `computeEigen` : compute the eigenmodes if true else load them
* `nbMode` : number of modes to load (for computing use `solvereigen.nev`)
* `print` : print matrices to Matlab format

#### SolverA options
* `needA0` : need relief a0
* `computeA0` : compute a0 if true else load it
* `radius` : cylinder's raduis
* `speed` : average speed
* `alpha0` : expression of alpha0, depends on radius and speed

* `needA1` : need relief a1
* `computeA1` : compute a1 if true else load it
* `alpha1` : expression of alpha1

* `needA2` : need relief a2
* `computeA2` : compute a2 if true else load it
* `alpha2` : expression of alpha2

#### SolverSpectralProblem options
* `needSP`: need to run the spectral problem
* `nu` : viscosity
* `computeRijk` : compute Rijk coefficients if true else load them
* `computeRiak` : compute Riak coefficients if true else load them
* `computeRaik` : compute Raik coefficients if true else load them
* `computeRfk` : compute Rfk coefficients if true else load them
* `computeRpk` : compute Rpk coefficients if true else load them
* `f` : right hand side

* `v_ex` : v exact

All the options `compute...` if false need to be use with the exact same mesh than the one use to compute them.

In order to avoid a null pivot during the eigen problem, use the option `-eps_target <shift>`

### Examples

`mpirun -np 2 feelpp_po_app --solverns2.needEigen true --solverns2.computeEigen true --solvereigen.nev 10 --solvereigen.ncv 20 --solverns2.markerList wall1 wall2`

`mpirun -np 3 feelpp_po_app --solverns2.needEigen true --solverns2.computeEigen false --solverns2.nbMode 10 --solverns2.needA0 --solverns2.speed 1 --solverns2.radius 1 --solverns2.alpha0 "2.*speed*(1.-(x*x+y*y)/(radius*radius))"`
