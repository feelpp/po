[![Join the chat at https://gitter.im/feelpp/po](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/feelpp/po?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

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

### Running DMD on Atlas with Slurm

```#!/bin/bash

# DMD on 48 CPU

#SBATCH -p public
# number of cores
#SBATCH -n 48
# Hyperthreading is enabled on atlas, if you do not want to use it
# You must specify the following option
#SBATCH --ntasks-per-core 1
# min-max number of nodes
#SBATCH -N 2
# max time of exec (will be killed afterwards)
##SBATCH -t 12:00:00
# number of tasks per node
##SBATCH --tasks-per-node 1
# specify execution constraitns
##SBATCH --constraint \"intel\"
# min mem size
##SBATCH --mem=16684
# display info about cpu binding
##SBATCH --cpu_bind=verbose
# send a mail at the end of the exec
#SBATCH --mail-type=END
#SBATCH --mail-user=login@server.com

export JOBS=48

# data repository
export SCRATCH_DIR=/data/scratch/$USER/cyl.$JOBS
mkdir -p $SCRATCH_DIR

source /etc/profile.d/modules.sh
module load feelpp.profile

# ACUSIM parent directory
export FROM_DIR=/data/plasticomnium/Cylindre_long/20171205_cyl

# move to installation directory
cd ~/opt/po/src

# run the conversion steps
# --idir: input directory (ACUSIM data)
# --runId: ACUSIM run Id (defaults to 1)
# --odir: feelpp data directory
# --iname: simulation data name
# --fields: fields to export
# --DoConversion: ??? but required for some reason
mpirun --bind-to core -x LD_LIBRARY_PATH ./feelpp_po_converter --idir $FROM_DIR --runId=1 --odir $SCRATCH_DIR/db --iname cyl --fields pressure --DoConversion=1 > $SCRATCH_DIR/stdout_conversion_1.log
mpirun --bind-to core -x LD_LIBRARY_PATH ./feelpp_po_converter --idir $FROM_DIR --runId=1 --odir $SCRATCH_DIR/db --iname cyl --fields pressure --DoConversion=0 > $SCRATCH_DIR/stdout_conversion_0.log
# run the DMD step
# --ifile: feelpp db file
# --directory: ensight export directory
# --field: field to process
# --scalar-product: product to use, defaults to L2
# --time-initial: must be large enough for the behavior to be linear
# --export-modes: DMD modes to export
# --v: feelpp verbose option
mpirun --bind-to core -x LD_LIBRARY_PATH ./feelpp_po_dmd --ifile $SCRATCH_DIR/db/mydb.dbfeelpp.json --directory $SCRATCH_DIR --field pressure --scalar-product L2 --time-initial=0.01 --export-modes=32 --v 1 > $SCRATCH_DIR/dmd_$JOBS.log```
