#!/bin/bash

###Programme qui permet de faire de la conversion des fichiers pdb en format sdf


if  [ ! $1 ]
then
  echo "manque la liste des pdb"
  exit
fi


#file=/home/ampere/regad/Pocket_project/Dataset-Stephanie/liste.txt
file=$1

dir_data=/home/ampere/regad/Pocket_project/Dataset-Stephanie/Complexe_data/
#dir_data=/home/ampere/regad/Pocket_project/Validation_Classif/FromLigand/

dir_work=/home/ampere/regad/Pocket_project/Dataset-Stephanie/
#dir_work=/home/ampere/regad/Pocket_project/Validation_Classif/FromLigand/



rm -rf $dir_work/All_ligand.sdf



for pdb in $(\cat $file)
  do
     echo "   >>>"$pdb
     cd $dir_data/$pdb
     babel -ipdb ligand.pdb -osdf $pdb"_ligand".sdf 
     echo $pdb"_ligand".sdf >> $dir_work/All_ligand.sdf
     cat $pdb"_ligand".sdf  | sed '1d' >> $dir_work/All_ligand.sdf
     cd $dir_work
  done
  




