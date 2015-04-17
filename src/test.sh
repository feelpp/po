#!/bin/bash

APP=feelpp_po_basischange

SOLVER=("power" "lapack" "subspace" "arnoldi" "krylovschur")
PROBLEM=("ghep" "gnhep")
TRANSFORM=("shift" "shift_invert" "fold" "cayley")
SPECTRUM=("largest_magnitude" "smallest_magnitude" "largest_real" "smallest_real")

# for np in ${NP[@]};
# do
#     for nev in ${NEV[@]};
#     do
#         ncv=$(( $(($nev + 15)) > $(($nev*2)) ? $(($nev + 15)) : $(($nev*2)) ))
#         mpirun -np $np $APP --solvereigen.nev=$nev --solvereigen.ncv=$ncv -st_redundant_pc_factor_zeropivot 1e-40
#     done
# done

for solver in ${SOLVER[@]};
do
    for problem in ${PROBLEM[@]};
    do
	for transform in ${TRANSFORM[@]};
	do
	    for spectrum in ${SPECTRUM[@]};
	    do
		echo "$APP --gmsh.hsize 0.15 --solvereigen.solver $solver --solvereigen.problem $problem --solvereigen.transform $transform --solvereigen.spectrum $spectrum -st_shift 16"
		./$APP --gmsh.hsize 0.15 --solvereigen.solver $solver --solvereigen.problem $problem --solvereigen.transform $transform --solvereigen.spectrum $spectrum -st_shift 16
		echo "$APP --gmsh.hsize 0.15 --solvereigen.solver $solver --solvereigen.problem $problem --solvereigen.transform $transform --solvereigen.spectrum $spectrum"
		./$APP --gmsh.hsize 0.15 --solvereigen.solver $solver --solvereigen.problem $problem --solvereigen.transform $transform --solvereigen.spectrum $spectrum
	    done
	done
    done
done
