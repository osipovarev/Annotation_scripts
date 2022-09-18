#!/usr/bin/env bash
#


###
db=$1
ANNODIR=$(pwd)
EVIDIR=$ANNODIR/$db/Evidences/
EVMDIR=$ANNODIR/$db/EVM/
GTFSEVM=$EVMDIR/evm_out.lst
JOINEVM=$db/$db.evm.join.gtf
OUTBED=$db/$db.evm.join.bed


## convert EVM results to gtf
echo "Running evm.out -> gff3 -> gtf....."
for i in $(find $EVMDIR -name \*evm.out | awk -F "/" 'BEGIN { OFS="/";} {$NF=""; print $0 }');
do
	echo -e "$i/$db.evm.out.gtf" >> $GTFSEVM
	EVM_to_GFF3.pl $i/$db.evm.out $i | gt gff3 -tidy -sort | gt gff3_to_gtf |  formatEvmOutGtf.pl > $i/$db.evm.out.gtf
done

## merge results
echo -e "Merging results....."
joingenes -f $GTFSEVM -o $JOINEVM

## clean ./ and / after joingenes if needed
#sed -i 's/\.\///g' $JOINEVM
sed -i "s#$EVMDIR##g" $JOINEVM
sed -i 's/\///g'  $JOINEVM


## gtf -> gff3 -> gp -> bed
echo -e "Running gtf -> gff3 -> gp -> bed......."
gt gtf_to_gff3 $JOINEVM | gt gff3 -tidy -sort -fixregionboundaries | gff3ToGenePred stdin stdout | genePredToBed stdin $OUTBED

## check if the resulting file is ok
if [ ! -s $OUTBED ];
then
	echo -e "File $OUTBED is empty! CHeck your EVM run in $EVMDIR"
	exit 1
fi


## clean up
echo -e "Cleaning up now....."
rm -r $EVMDIR
rm $JOINEVM

echo -e "All Done! Check results in $OUTBED"
