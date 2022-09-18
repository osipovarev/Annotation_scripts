#!/usr/bin/env bash
#


###
db=$1
REFLIST=$2
# number of references it should be found in
SHARED=$3
ANNODIR=$(pwd)
EVMBED=$ANNODIR/$db/$db.evm.join.bed
OUTBED=$ANNODIR/$db/$db.evm.added_togas.bed

TOGAS=""
TOGABED="query_annotation.bed"
REFS=$(echo $REFLIST | sed 's/,/\./g')

for ref in $(echo $REFLIST | sed 's/,/ /g');
do
	TOGASDIR="/projects/hillerlab/genome/gbdb-HL/$ref/TOGA/vs_${db}"
	TOGAS=$TOGAS" $TOGASDIR/$TOGABED"
done

getOverlappingTranscripts.py -n $SHARED -f $TOGAS > $ANNODIR/$db/uniq.overlap.$REFS.bed


remove_duplicated_names.py -f <(getUniqTranscripts.py -f <(cat $EVMBED $ANNODIR/$db/uniq.overlap.$REFS.bed ) 2> /dev/null) > $OUTBED
rm $ANNODIR/$db/uniq.overlap.$REFS.bed

echo "All Done! Added TOGAs to gene models. See results in $OUTBED"


