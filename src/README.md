## Sources for PO project

The application is composed from the followings classes :

- [SolverNS2](solverns2.hpp) : call the other classes to initialize the objects needed
- [SolverA](solvera.hpp) : solves the problems for the boundary conditions
- [SolverEigenNS2](solvereigenns2.hpp) : solves the eigen problem
- [SolverSpectralProblem](solverspectralproblem.hpp) : solves the spectral problem

### Options

All the options are prefixed by `solverns2.`
Here, a list of the options :
* `verbose` : level of verbosity (int default 0)

#### SolverEigenNS2 options
* `needEigen` : need eigenmodes (bool default true)
* `markerList` : list of markers of the boundary (list of string)
* `computeEigen` : compute the eigenmodes if true else load them (bool default true)
* `nbMode` : number of modes to load (for computing use `solvereigen.nev`) (int default 1)
* `print` : print matrices to Matlab format (debug) (bool default false)

#### SolverA options
* `needA0` : need relief a0 (bool default false)
* `computeA0` : compute a0 if true else load it (bool default true)
* `radius` : cylinder's raduis (double default 0.5)
* `speed` : average speed (double default 1)
* `alpha0` : expression of alpha0, depends on radius and speed (ginac expression default "2. * speed * (1. - (x*x + y*y) / (radius * radius))")

* `needA1` : need relief a1 (bool default false)
* `computeA1` : compute a1 if true else load it (bool default true)
* `alpha1` : expression of alpha1 (ginac expression default "0.")

* `needA2` : need relief a2 (bool default "false")
* `computeA2` : compute a2 if true else load it (bool default "true")
* `alpha2` : expression of alpha2 (ginac expression default "4.*speed/(radius*radius)")

#### SolverSpectralProblem options
* `needSP`: need to run the spectral problem (bool default true)
* `nu` : viscosity (double default 1)
* `computeRijk` : compute Rijk coefficients if true else load them (bool default false)
* `computeRiak` : compute Riak coefficients if true else load them (bool default false)
* `computeRaik` : compute Raik coefficients if true else load them (bool default true)
* `computeRfk` : compute Rfk coefficients if true else load them (bool default true)
* `computeRpk` : compute Rpk coefficients if true else load them (bool default true)
* `f` : right hand side (ginac expression default "{0,0,1}")

* `v_ex` : v exact (ginac expression default "{0,0,2*(1-4*(x*x + y*y))}:x:y")

All the options `compute...` if false need to be use with the exact same mesh than the one use to compute them.
`nbMode` must be lesser than the number of modes saved on disk.

All these options can be put in a config file. By default, [po_app.cfg](po_app.cfg) will be loaded. You can specify custom config file with `--config-file` option.

For the eigen problem, the number of eigen modes computed is specified by the `--solvereigen.nev` and `--solvereigen.ncv` options. `ncv` must be greater than `nev`.
In order to avoid a null pivot during the eigen problem, use the option `-eps_target <shift>`. This option can not be put in the config file for the moment.

### Examples

To compute 10 eigen modes on a geo file where there is 2 connex boundary marked wall1 and wall2, on 2 proc, use :
`mpirun -np 2 feelpp_po_app --solverns2.needEigen true --solverns2.computeEigen true --solvereigen.nev 10 --solvereigen.ncv 20 -eps_target 35 --solverns2.markerList wall1 wall2`

To load 10 eigen modes from the disk, and compute **a** with alpha0 as a ginac expression, use :
`mpirun -np 3 feelpp_po_app --solverns2.needEigen true --solverns2.computeEigen false --solverns2.nbMode 10 --solverns2.needA0 true --solverns2.computeA0 true --solverns2.speed 1 --solverns2.radius 1 --solverns2.alpha0 "2.*speed*(1.-(x*x+y*y)/(radius*radius))"`
