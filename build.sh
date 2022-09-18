#!/bin/bash
set -o errexit


NORM="\\033[0;39m" 
BLUE="\\033[1;34m" 

echo -e  "$BLUE" " Copying scripts gene annotation scripts to $genomePath/bin/scripts ..." "$NORM"
echo -e  "$BLUE" "--------------------------" "$NORM"
cp ${genomePath}/src/gene_annotation_scripts/*.sh $genomePath/bin/scripts/
cp ${genomePath}/src/gene_annotation_scripts/*.py $genomePath/bin/scripts/
cp ${genomePath}/src/gene_annotation_scripts/*.p*l $genomePath/bin/scripts/

echo -e  "$BLUE" " Changing group to hillerlab and add write rights ..." "$NORM"
echo -e  "$BLUE" "--------------------------" "$NORM"

#chgrp hillerlab ${genomePath}/bin/scripts/* -R
#chmod g+w ${genomePath}/bin/scripts/* -R

echo -e  "$BLUE" " DONE" "$NORM"
echo -e  "$BLUE" " --------------------------" "$NORM"
 
