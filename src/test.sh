#!/bin/bash

APP=feelpp_po_app
GEO=cylinderAdim.geo

NP=(1 2 5 10 20 30)
NEV=(10 25 50 100 200)
H=(0.1 0.075 0.05)

# for np in ${NP[@]};
# do
#     for nev in ${NEV[@]};
#     do
#         ncv=$(( $(($nev + 15)) > $(($nev*2)) ? $(($nev + 15)) : $(($nev*2)) ))
#         mpirun -np $np $APP --solvereigen.nev=$nev --solvereigen.ncv=$ncv -st_redundant_pc_factor_zeropivot 1e-40
#     done
# done

for h in ${H[@]};
do
    ./$APP --gmsh.filename $GEO --gmsh.hsize $h
done
