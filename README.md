PO Project
============

This is the PO project

## Build Status

here is the Travis build status of po:
[![Build Status](https://magnum.travis-ci.com/feelpp/po.png?token=Bxps8gX6edMDEv345qns&branch=master)](https://magnum.travis-ci.com/feelpp/po)



## Branches

Po follows the git flow framework. The default branch is now develop on github and 
the master branch is used for stable releases.

* master : use the method describe in the report
* compZ : use only one component
* strongBC : apply strongly the boundary condition in the cylinder
* g0 : find a basis of H¹_0 and find the corresponding basis P_sigma(g_0)
* darcy : use the Raviart-Thomas elements to compute psi0

## Mailing lists

 - po@cemosis.fr  (discussions)
 - po-commits@cemosis.fr (git commits)

## Programming environment

It is assumed that you have cloned Feel++ (``` git clone http://github.com/feelpp/feelpp.git```).

### Cloning Po

Then issue the following commands in order to clone the po repository in feelpp/ research
```
cd research
git clone http://github.com/feelpp/po.git
```

### Compiling Po

Just do the ```cmake``` Feel++ procedure with the option ```-DFEELPP_ENABLE_RESEARCH_PO=ON```and it will configure Po and generate the necessary Makefiles.
Afterwards just compile the applications you are interested in.

### Running Po

Please see the wiki page.
