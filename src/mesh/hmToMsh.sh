#!/bin/bash

dossier='HYPERMESH.DIR'
file='../out.msh'
verbose='false'

while getopts 'hd:o:v' flag; do
    case "${flag}" in
	h) echo "hmToMsh transcode a HyperMesh mesh to a Gmsh mesh"
	    echo "The files .crd, .ebc and .cnn have to be in a HYPERMESH.DIR directory and the output is in out.msh"
	    echo "options:"
	    echo "-h      show this help"
	    echo "-d      change the directory of the HyperMesh files"
	    echo "-o      change the output file"
	    echo "-v      be verbose"
	    exit 0
	    ;;
	d) dossier="${OPTARG}";;
	o) file="../${OPTARG}";;
	v) verbose='true';;
    esac
done

if [ -d $dossier ]
then

cd $dossier
if [ $verbose = 'true' ]
then echo 'debut du transcodage'
fi

resDebut=$(date +%s)

geoEntity=0
phEntity=1

# Format
echo '$MeshFormat' > $file
echo '2.2 0 8' >> $file
echo '$EndMeshFormat' >> $file

# Nodes
echo '$Nodes' >> $file
nbNodes=`cat *.crd | wc -l`
echo $nbNodes >> $file
cat *.crd >> $file
echo '$EndNodes' >> $file

resNodes=$(date +%s)
dt=$(echo "$resNodes-$resDebut" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
if [ $verbose = 'true' ]
then echo "$nbNodes Noeuds : $dm min $ds sec"
fi

# Elements
echo '$Elements' >> $file

# Nombre d'elements
nbTria=`cat *.ebc | wc -l`
nbTetra=`cat *.cnn | wc -l`
let nbElement="$nbTria+$nbTetra" 
echo $nbElement >> $file

# Triangles
arrayEbc=(`ls *.ebc | sort -t "." -k 4`)
for i in ${arrayEbc[*]}; do
    sed -r 's/[ \t]*[0-9]+[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)/\1 2 2 '"$phEntity $geoEntity"' \2 \3 \4/g' $i >> $file
    let geoEntity++;
    let phEntity++;
done

resTria=$(date +%s)
dt=$(echo "$resTria-$resNodes" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
if [ $verbose = 'true' ]
then echo "$nbTria Triangles : $dm min $ds sec"
fi

# Tetra
arrayCnn=(`ls *.cnn`)
for i in ${arrayCnn[*]}; do
    sed -r 's/[ \t]*([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)/\1 4 2 '"$phEntity $geoEntity"' \2 \3 \4 \5/g' $i >> $file
    let geoEntity++;
    let phEntity++;
done

echo '$EndElements' >> $file

resTetra=$(date +%s)
dt=$(echo "$resTetra-$resTria" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
if [ $verbose = 'true' ]
then echo "$nbTetra Tetras : $dm min $ds sec"
fi

# Physical Entity
echo '$PhysicalNames' >> $file

# Nombre d'entite physique
let phEntity--;
echo $phEntity >> $file

# Surface
phEntity=1
arrayEbcName=(`ls *.ebc | sort -t "." -k 4 | sed -r 's/.*\.(fluid_[^.]*)\.tet4\.([^.]+)\.tria3\.ebc/\1_\2/g'`)
for i in ${arrayEbcName[*]}; do
    echo "2 $phEntity \"$i\"" >> $file
    let phEntity++;
done

# Volume
arrayCnnName=(`ls *.cnn | sed -r 's/.*\.([^.]+)\.tet4\.cnn/\1/g'`)
for i in ${arrayCnnName[*]}; do
    echo "3 $phEntity \"$i\"" >> $file
    let phEntity++;
done

echo '$EndPhysicalNames' >> $file

resFin=$(date +%s)
dt=$(echo "$resFin-$resDebut" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
if [ $verbose = 'true' ]
then
echo "$nbElement Elements : $dm min $ds sec"
echo "fin du transcodage"
fi

else
echo "pas de dossier $dossier"
fi
