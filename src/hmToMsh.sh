#!/bin/bash

echo 'debut du transcodage'
resDebut=$(date +%s)

file="test.msh"
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
dt=$(echo "$resDebut - $resNodes" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
echo "$nbNodes Noeuds : $dm min $ds sec"

# Elements
echo '$Elements' >> $file

# Nombre d'elements
nbTria=`cat *.ebc | wc -l`
nbTetra=`cat *.cnn | wc -l`
let nbElement="$nbTria+$nbTetra" 
echo $nbElement >> $file

# Triangles
arrayEbc=(`ls *.ebc`)
for i in ${arrayEbc[*]}; do
    sed -r 's/[ \t]*[0-9]+[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)/\1 2 2 '"$phEntity $geoEntity"' \2 \3 \4/g' $i >> $file
    let geoEntity++;
    let phEntity++;
done

resTria=$(date +%s)
dt=$(echo "$resTria-$resNodes" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
echo "$nbTria Triangles : $dm min $ds sec"

# Tetra
sed -r 's/[ \t]*([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9]+)/\1 4 2 '"$phEntity $geoEntity"' \2 \3 \4 \5/g' *.cnn >> $file

echo '$EndElements' >> $file

resTetra=$(date +%s)
dt=$(echo "$resTetra-$resTria" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
echo "$nbTetra Tetras : $dm min $ds sec"


# Physical Entity
echo '$PhysicalNames' >> $file

# Nombre d'entite physique
echo $phEntity >> $file

# Surface
phEntity=1
arrayEbcName=(`ls *.ebc | sed -r 's/.*\.([^.]+)\.tria3\.ebc/\1/g'`)
for i in ${arrayEbcName[*]}; do
    echo "2 $phEntity \"$i\"" >> $file
    let phEntity++;
done

# Volume
nomVolume=(`ls *.cnn | sed -r 's/.*\.([^.]+)\.tet4\.cnn/\1/g'`)
echo "3 $phEntity \"$nomVolume\"" >> $file

echo '$EndPhysicalNames' >> $file

resFin=$(date +%s)
dt=$(echo "$resFin-$resDebut" | bc)
dm=$(echo "$dt/60" | bc)
ds=$(echo "$dt-$dm*60" | bc)
echo "$nbElement Elements : $dm min $ds sec"
echo "fin du transcodage"
