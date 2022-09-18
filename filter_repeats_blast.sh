#!/usr/bin/env bash
#

set -euo pipefail

if [ "$#" -lt 2 ];
then
	echo "Usage: $0 [db] [db.denovo.bed]"
	exit 0
fi

db=$1
denovo_gff3=$2
where=$(pwd)
wdir=$where/$db
denovo=${denovo_gff3%gff3}bed
BLASTDIR="BlastFilter_repeats/"

cd $wdir
mkdir -p $BLASTDIR/
norep=${denovo#*/}
norep=${norep%.*}.norepeats.bed
norep_prot=${norep%bed}aa.fas

# gff3 ->gp -> bed
echo  -e "Running now:\n gff3ToGenePred $denovo_gff3 stdout | genePredToBed stdin $denovo"
gff3ToGenePred $denovo_gff3 stdout | genePredToBed stdin $denovo

echo -e "Running now:\n bedtools intersect -v -f 0.1 -split -a $denovo -b repeats.$db.bed -wa > $BLASTDIR/$norep"
bedtools intersect -v -f 0.1 -split -a $denovo -b repeats.$db.bed -wa > $BLASTDIR/$norep

# bed -> gp -> prot
echo -e "Running now:\n bedToGenePred $BLASTDIR/$norep stdout | genePredToProt -includeStop -starForInframeStops stdin /projects/hillerlab/genome/gbdb-HL/$db/$db.2bit $BLASTDIR/$norep_prot"
bedToGenePred $BLASTDIR/$norep stdout | genePredToProt -includeStop -starForInframeStops stdin /projects/hillerlab/genome/gbdb-HL/$db/$db.2bit $BLASTDIR/$norep_prot


# split prot.fastas
cd $BLASTDIR/
mkdir -p split_fasta/
mkdir -p results_blast/
cd split_fasta/
faSplit sequence ../$norep_prot 300 split_fasta
cd $wdir/$BLASTDIR/

# prepare blastp jobList
SWISSPROT="/projects/hillerlab/genome/data/uniref/swissprot_vertebrates"
for i in $( ls split_fasta/); do echo "blastp -evalue 1e-10 -num_threads 24 -db $SWISSPROT -query split_fasta/$i -outfmt 6 -out results_blast/${i%.fa}.BlastHits_out6.tsv"; done > jobList.$db.blastp

echo -e "All done for $db "
echo -e "jobs are written in $wdir/$BLASTDIR/jobList.$db.blastp"
echo -e "Running it now......."
bash jobList.$db.blastp

echo -e "Finished running blatsp jobs."


# merge balst hits
outblast=merged_${db}.denovo.norepeats.BlastHits_out6.tsv
cat results_blast/*.tsv > $outblast

# filter
cd $wdir
speciesdict="/projects/hillerlab/genome/data/uniref/vertebrates.dictionary.lst"
filtered_fa=filtered.$norep_prot
filtered_bed=filtered.$norep
filter_fasta_with_blast.py -l 200 -b $BLASTDIR/$outblast -f $BLASTDIR/$norep_prot -s $speciesdict > $BLASTDIR/$filtered_fa
filter_bed_with_fasta.py -a $BLASTDIR/$norep -f $BLASTDIR/$filtered_fa > $filtered_bed


# go back
cd $where
echo -e "ALL DONE for $db. Check $filtered_bed"

