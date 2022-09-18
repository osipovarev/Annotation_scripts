#!/bin/bash
#SBATCH -J gth               # Job name
#SBATCH -n 1                 # Run a single task        
#SBATCH -c 1                 # Number of CPU cores per task
#SBATCH --time=8:00:00	     # Set walltime limit
#SBATCH --mem-per-cpu=5000  # Job memory request
#SBATCH -o log.gth	     # Standard output and error log

db=$1
echo "Working with assembly $db"

## report what node we're on
echo "Working on $HOSTNAME" 

## make working dir in /tmp
mkdir -p /tmp/$USER/$SLURM_JOB_ID/
echo "Building in /tmp/$USER/$SLURM_JOB_ID/"

## goto working dir
cd /tmp/$USER/$SLURM_JOB_ID/

## put genome here 
cp /projects/hillerlab/genome/gbdb-HL/$db/$db.fa .

## put proteins here
cp /projects/project-osipova/NectarivoryProject/Genome_annotation_2021/split_all_protein/$2 .

## gth command
gth -skipalignmentout -gff3out -paralogs -prseedlength 20 -prminmatchlen 20 -prhdist 2 -minmatchlen 32 -seedlength 32 -gcmincoverage 70 -species chicken -genomic $db.fa -protein $2 -o $db.$2.gthAln.gff3 -force

## move result file back to project space

mv $db.$2.gthAln.gff3 /projects/project-osipova/NectarivoryProject/Genome_annotation_2021/$db/Protein_evidence/results_gth

## clean up on the node! no matter what..
rm -r /tmp/$USER/$SLURM_JOB_ID/
