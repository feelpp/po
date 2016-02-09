

## Etape 1 : conversion du maillage AcusimRawMesh et partitionement
Compiler ce programme (se trouvant dans feelpp/applications/mesh/ )
```
make feelpp_mesh_partitioner
```
Lancer l'executable en séquentielle soit directement en ligne de commande ou soit 
avec un fichier de config (en utilisant ```--config-file=myfile.cfg```).

Voici un exemple qui va générer 4 paritionements (on peut en mettre un nombre arbitraire) à partir du fichier .arm
```
dim=3
ifile=/mypath/test.1.arm
odir=/myoutputpath/
part=1
part=4
part=8
part=16
```
En ligne de commande :
```
feelpp_mesh_partitioner --dim 3 --ifile /mypath/test.1.arm --odir /myoutputpath/ --part 1 4 8 16
```

Ceci devrait générer 4 maillages partitionés au format HDF5 dans le repertoire ```myoutputpath``` suffixer par _p1, _p4, _p8, _p16.
Chaque maillage est constititué de 2 fichier (.json, .h5). Le premier est une description des marqueurs, l'autre contient toue la structure de donnée du maillage.

Pour feel++, il faudra ensuite donner le nom du fichier .json (via ```--gmsh.filename=mymesh_p4.json``` par exemple). 


## Etape 2 : conversion des résultats d'Acusim (vitesse P1 et pression P1)
Compiler ce programme (se trouvant dans feelpp/research/po/data_conversion/ )
```
make feelpp_conversion_acusim_fields 
```
Ensuite l'utilsateur devra fournir le .json crée dans l'étape précédente, le fichier de noeud .crb et
une liste de fichier de vitesse et pression (.out en ASCII).
Ce code permettra de sauvegarder au format HDF5 le champ de vitesse et pression sur le partitionnement qui a été choisi (via la .json).

Voici un example de fichier de config :
```
input.mesh.filename=/myoutputpath/test.1_p4.json

input.acusim.nodes=/mypath/MESH.DIR/test.crd
input.acusim.pressure=/mypath/DATA_TXT/data_pression/test_step3712.out
input.acusim.pressure=/mypath/DATA_TXT/data_pression/test_step4000.out
input.acusim.velocity=/mypath/DATA_TXT/data_vitesse/test_step3712.out
input.acusim.velocity=/mypath/DATA_TXT/data_vitesse/test_step4000.out

output.directory=/myoutputpath/fields/

do-export.feel-format=1
do-export.visu-format=0
```

Attention, il faut lancer ce code avec un nombre de processus égale au partitionnement. Par exemple, si on choisit le fichier .json avec 8 paritions,
il faut lancer ce code sur 8 proc ( via par exemple ```mpirun -np 8 ./feelpp_conversion_acusim_fields --config-file myfile.cfg``` )

L'option ```do-export.feel-format``` permet de dire qu'il faut sauvegarder les champs au format HDF5 (et relisable par feel++ par la suite).
L'option ```do-export.visu-format``` permet de dire qu'il faut sauvegarder les champs au format ensightgold (et relisable par Paraview ou HyperView ( mais pas par feel++ ).





| Tables   |      Are      |  Cool |
| -------- | ------------- | ----- |
| col 1 is |  left-aligned | $1600 |
| col 2 is |    centered   |   $12 |
| col 3 is | right-aligned |    $1 |