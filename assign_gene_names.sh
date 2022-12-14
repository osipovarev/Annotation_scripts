#!/usr/bin/env bash
#

set -eo pipefail

ANNOBED=$1
NAMEDANNO=$2
NAMESFILE=$3
OUTFILE=$4


if [ $# -lt 4 ];
then
	echo -e "Usage: $0 [ANNOBED] [NAMEDANNO] [NAMESFILE] [OUTFILE]
		ANNOBED : annotation.bed , with transcripts, like trans_1
		NAMEDANNO : annotation_with_names.bed with transcripts, like rna-XM_004941692.3
		NAMESFILE : gene - transcript names correspodance, like SUSD6 \t rna-XM_004941692.3
		OUTFILE : out.bed"
	exit 0
fi


WDIR=$(pwd)
ANNODICT="temp.anno_dict.txt"
NAMESDICT="temp.names_dict.txt"


## Check if all input files exist
if [ ! -f $ANNOBED ];
then
	echo -e "$ANNOBED file not found!"
	exit 1
fi

if [ ! -f $NAMEDANNO ];
then
	echo -e "$NAMEDANNO file not found!"
	exit 1
fi


echo -e "Running overlapSelect....."
echo -e "NB: using 80% CDS overlap requirement"
overlapSelect -statsOutput -strand -inCds -overlapThreshold=0.8 -inFmt=bed -selectFmt=bed $NAMEDANNO $ANNOBED stdout | cut -f1-2 | rev | cut -d"." -f2- | rev | awk '{print $2"\t"$1}' > $ANNODICT

echo -e "Running renamig using replace_names_from_dict.py....."
replace_names_from_dict.py -f 1 -a $ANNODICT -d <(awk '{print $2","$1}' $NAMESFILE) | grep -v  ^reg_ | grep -v ^rna- | awk '{print $2","$1}' | sort -u > $NAMESDICT

echo -e "Giving geneName_ prefix to each transcript....."
replace_names_from_dict.py -a $ANNOBED -d $NAMESDICT -p  > $OUTFILE

echo -e "Cleaning up....."
rm $ANNODICT $NAMESDICT

echo -e "All Done! Check $OUTFILE"
