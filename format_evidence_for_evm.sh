#!/usr/bin/env bash
#


# Asumes TOGA dir structure: /projects/hillerlab/genome/gbdb-HL/$ref/TOGA/vs_$db 
# Expects all inputs in gff3 format

set -o pipefail

### Options parser
while getopts "d:t:p:r:a:c" opt;
do
	case $opt in
		d) db=$OPTARG ;;
		t) TOGAREFS=$OPTARG ;;
		p) PROT=$OPTARG ;;
		r) RNA=$OPTARG ;;
		a) DENOVO=$OPTARG ;;
		c) CONTINUE=$OPTARG ;;
	esac
done

### Check if at least DB arguemnt was provided
if [ -z $db ];
then
	echo "Usage: $0 -d DB [OPTIONS]"
	exit 1
fi


### Inititate Evidence directory structure and variables for predictions
ANNODIR=$(pwd)
EVIDIR=$ANNODIR/$db/Evidences/
EVMDIR=$ANNODIR/$db/EVM/
mkdir -p $EVIDIR
mkdir -p $EVMDIR
GENEPREDS=""
EVMPROT=""
rm -f "$EVMDIR/weights.txt"



## Prepare fasta files for the $db
FASTA="$genomePath/gbdb-HL/$db/$db.fa"
if [ ! -f $FASTA ];
then
        twoBitToFa $genomePath/gbdb-HL/$db/$db.2bit $genomePath/gbdb-HL/$db/$db.fa
fi



### Prepare TOGA projections for all references
if [ ! -z $TOGAREFS ];
then
	for ref in $( echo $TOGAREFS | sed 's/,/ /g');
	do
	        TOGADIR="/projects/hillerlab/genome/gbdb-HL/$ref/TOGA/vs_${db}"
		TOGABED="query_annotation.bed"
		echo -e "Running formating of TOGA $ref projections to $db....."

		if [ -f $TOGADIR/$TOGABED ];
		then
			bedToGenePred $TOGADIR/$TOGABED stdout | genePredToGtf file stdin stdout | gt gtf_to_gff3 -tidy | gt gff3 -tidy -sort -setsource $ref | evmFormatIsoSeq.pl > $EVIDIR/$db.evm.toga.$ref.gff3
			GENEPREDS=$GENEPREDS" $EVIDIR/$db.evm.toga.$ref.gff3"
	            	echo -e "OTHER_PREDICTION $ref 8" >> $EVMDIR/weights.txt
        else
			echo -e "File $TOGADIR/$TOGABED does not exist!"
			exit 1
		fi
	done
fi



### Prepare protein
if [ ! -z $PROT ];
then
	if [ -f $PROT ];
	then
		evmFormatAln.pl $PROT > $EVIDIR/$db.evm.protein.gff3
        	EVMPROT=$EVMPROT" --protein $EVIDIR/$db.evm.protein.gff3"
        echo -e "PROTEIN gth 2" >> $EVMDIR/weights.txt
	else
		echo -e "File $PROT does not exist!"
		exit 1
	fi
fi


### Prepare RNAseq
if [ ! -z $RNA ];
then
	if [ -f $RNA ];
	then
		gt gtf_to_gff3 -tidy $RNA | gt gff3 -tidy -sort -setsource rna | evmFormatIsoSeq.pl > $EVIDIR/$db.evm.rna.gff3
        	GENEPREDS=$GENEPREDS" $EVIDIR/$db.evm.rna.gff3"
        echo -e "OTHER_PREDICTION rna 2" >> $EVMDIR/weights.txt
	else
		echo -e "File $RNA does not exist!"
		exit 1
	fi
fi


### Prepare all denovo gene predictions
if [ ! -z $DENOVO ];
then
	if [ -f $DENOVO ];
	then
		cat $DENOVO | gt gff3 -tidy -sort -setsource denovo | evmFormatBraker.pl | grep -v '# ' > $EVIDIR/$db.evm.denovo.gff3
	        GENEPREDS=$GENEPREDS" $EVIDIR/$db.evm.denovo.gff3"
        echo -e "ABINITIO_PREDICTION denovo 1" >> $EVMDIR/weights.txt
    else
		echo -e "File $DENOVO does not exist!"
		exit 1
	fi
fi


### Combine all genePreds for EVM
echo -e "Merging all evidences now: $GENEPREDS"

if [ -z $GENEPREDS ];
then
	echo -e "No gene predictions to merge! Check your input"
	exit 1
else
	gt merge $(echo "$GENEPREDS") > $EVIDIR/$db.evm.genepreds.gff3
fi

if [ ! -s $EVIDIR/$db.evm.genepreds.gff3 ];
then
	echo -e " File $EVIDIR/$db.evm.genepreds.gff3 is empty! You probably want to check your input."
	exit 1
fi


echo -e "$EVMDIR/weights.txt a written in $EVMDIR/weights.txt file; check it and change the default weights if necessary"




### Prepare EVM runs
echo -e "Preparing EVM partitions now......"
cd $EVMDIR

partition_EVM_inputs.pl --genome $FASTA --gene_predictions $EVIDIR/$db.evm.genepreds.gff3 $EVMPROT --segmentSize 1000000 --overlapSize 150000 --partition_listing evm_partitions_list.out

echo -e "Writing EVM commands in evm.$db.jobList......"
write_EVM_commands.pl --genome $FASTA --weights $EVMDIR/weights.txt --gene_predictions $EVIDIR/$db.evm.genepreds.gff3 $EVMPROT --output_file_name $db.evm.out --partitions evm_partitions_list.out > evm.$db.jobList

cd $ANNODIR

echo -e "All done!\n Next steps:\n 1) Check $EVMDIR/weights.txt file and change default weights if necessary\n 2) Run $EVMDIR/evm.$db.jobList with para and after it's finished run a script for processing EVM results"

