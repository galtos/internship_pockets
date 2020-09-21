#!/bin/bash




if  [ ! $1 ]
then
  echo "manque la liste des pdb (fichier avec le path) et le path"
  exit
else
  echo "Step 2 : protonation des protéines"
fi


file=$1
#file=liste_prot_with_lignm.csv
#path_work=/home/ampere/regad/Pocket_project/Dataset-Stephanie/Complexe_data   #Pour les data de stephanie
path_work=$2

li1="http://www.pdb.org/pdb/download/downloadLigandFiles.do?ligandIdList="
li2="&structuIsList="
li3="&instanceType=all&excludeUnobserved=false&includeHydrogens=false"



###creation de la liste des ligands

for ligne in $(\cat $file)
  do
  cd $path_work
  pdb=` echo $ligne |  awk 'BEGIN{FS=":"}{print $1}' | tr [a-z] [A-Z] ` 
  pdb2=` echo $ligne |  awk 'BEGIN{FS=":"}{print $1}' | tr [A-Z] [a-z] ` 
  cd $pdb2
  nm_lig=` echo $ligne | awk 'BEGIN{FS=":"}{print $2}' `
  echo $pdb " "$nm_lig
  cmd=` echo $li1$nm_lig$li2$pdb$li3 `
  wget $cmd
  mv downloadLigandFiles* ligand.sdf
  fichier=ligand.sdf 
  if ! [ -f $fichier ]; 
    then
        echo $pdb >> $path_work/liste_pdb_sans_sdfFile.dat
  fi
  done

