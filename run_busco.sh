#!/usr/bin/env bash
#

set -eo pipefail

# trap read debug



USAGE=$(cat <<EOF
        USAGE: run_busco.sh -m [MODE] -q [QUERY] -t [TWOBIT] -d [BUSCO ODB] -c [CORES]
        -m : MODE g/genome OR a/annotation
        -q : QUERY; if genome mode it can be: 
                        1) assembly code (mm10 or HLstrHab1) if 2bit file is in standard place;
                        2) genomic 2bit file (.2bit)
                        3) genomic fasta (.fa /.fasta /.faa /.fas)
                   if annotation mode it can be: 
                        1) bed12 file (.bed) 
                        2) genePred file (.gp)
                        3) protein fasta (.fa /.fasta /.faa /.fas)
        -t : TWOBIT required if MODE=annotation and .bed or .gp input is provided 
        -d : BUSCO ODB e.g mammalia_odb10
        -o : OUTPUT optional; default=busco_<mode>
        -c : CORES number of cores to use; default: 12
        -k : KEEPTEMP specify to keep temp files of the busco run
        -h : Prints this help messege

EOF
)


KEEPTEMP="F"


### IDEA:
# Push BUSCO job to cluster via sbatch or so
# or do not run on delta?

### Options parser
while getopts "m:d:q:t:o:c:knh" opt;
do
        case $opt in
                m) MODE=$OPTARG ;;
                d) ODB=$OPTARG ;;
                q) QUERY=$OPTARG ;;
                t) TWOBIT=$OPTARG ;;
                o) OUTPUT=$OPTARG ;;
                c) CORES=$OPTARG ;;
                k) KEEPTEMP="T" ;;
		h)
                        echo "$USAGE"
                        exit 0 ;;
                \?)
                        echo -e "Invalid option: -$OPTARG\n"
                        echo "$USAGE"
                        exit 1 ;;
                :)
                        echo -e "Option -$OPTARG requires an argument\n"
                        echo $USAGE
                        exit  1 ;;
        esac
done

shift $((OPTIND -1))


WDIR=$(pwd)
BUSCO="/beegfs/software/busco/datasets/odb10/"
HLPATH="/projects/hillerlab/genome/gbdb-HL/"


### Check if at least MODE, QUERY and ODB arguemnts were provided
if [ -z $MODE ] || [ -z $ODB ] || [ -z $QUERY ];
then
        echo "$USAGE"
        exit 1
fi


### Check provided MODE argument
if [ $MODE != "a" ] && [ $MODE != "annotation" ] && [ $MODE != "genome" ] && [ $MODE != "g" ];
then
	echo "Use correct MODE flag: a/annotation OR g/genome]"
	exit 1
fi


### Check the requested busco_odb
if [[ "$ODB" == *"/"* ]];
then
        ## it is a full path
        BUSCOODB=$ODB
else
        ## find odb in the default place
        BUSCOODB=$BUSCO$ODB
        echo -e "The default BUSCO ODB path will be used: $BUSCOODB"
fi

if [ ! -d $BUSCOODB ];
then
        echo "ERROR:Dir $BUSCOODB not found!"
        exit 1
fi


### Set the number of cores if not provided
if [ -z $CORES ];
then
        CORES=12
fi



check_input_path() {
        ## Checks input file; make path absolute
        QUERY=$1

        # relative path?
        if [ ! -f "$WDIR/$QUERY" ];
        then
                # absolute path?
                if [ ! -f $QUERY ];
                then
                        echo -e "ERROR:input file $QUERY OR $WDIR/$QUERY not found!" >&2
                        exit 1
                else
                        INPUT=$QUERY
                fi
        else
                INPUT="$WDIR/$QUERY"
        fi   
        echo "$INPUT"
}



######      MODE: genome      ######
if [ $MODE == "g" ] || [ $MODE == "genome" ];
then
        
        ### Check the input
        EXT="${QUERY##*.}"
        if [[ EXT == "2bit" ]];
        then
                ## Working with 2bit
                echo -e "2bit file provided . . .";
                INPUT=$(check_input_path $QUERY)
                FAINPUT="genomic.fasta"
                FA_make_cmd="twoBitToFa -noMask $INPUT $FAINPUT"

        elif [[ $EXT == "fa" ]] || [[ $EXT == "faa" ]] || [[ $EXT == "fasta" ]] || [[ $EXT == "fas" ]];
        then
                ## Working with genomic fasta
                echo -e "genomic fasta is provided . . ."
                INPUT=$(check_input_path $QUERY)
                FAINPUT=$INPUT
                FA_make_cmd=""

        else
                ## Working with assembly 2bit in the standard place
                echo -e "assuming that assembly code is provided . . ."
                INPUT="$HLPATH/$QUERY/$QUERY.2bit"
                if [ ! -f $INPUT ];
                then
                        echo -e "If you meant an assembly code, $INPUT does not exist."
                        exit 1
                fi

                ## Prepare unmasked genomic fasta
                FAINPUT="genomic.fasta"
                FA_make_cmd="twoBitToFa -noMask $INPUT $FAINPUT"
        fi

        ### Run BUSCO -genome
        TMPDIR=$(mktemp -d -t busco_XXXXXXXX)
        cd $TMPDIR
        echo -e "Preparing genomic fasta now . . ."
        eval $FA_make_cmd
        echo -e "Running genome busco now in $TMPRAND . . . "
        busco -i $FAINPUT -o busco_genome -l $BUSCOODB -m genome -c $CORES --offline -f 
fi




######      MODE: annotation      ######
if [ $MODE == "a" ] || [ $MODE == "annotation" ];
then
        ### Check input file
        INPUT=$(check_input_path $QUERY)
        echo -e "THIS IS YOUR INPUT: $INPUT"
        EXT="${QUERY##*.}"

        if [ -z $INPUT ];
        then
                echo "You are running annotation mode - provide annotation.bed/.gp or protein.fa file"
                exit 1
        fi


        if [ $EXT == "bed" ];
        then
                ## Working with annotation bed file
                echo -e "annotation bed file is provided . . ."
                ## check 2bit file
                if [[ -z $TWOBIT ]] || [[ ! -f $TWOBIT ]];
                then
                        echo "ERROR: 2bit file not found! provide a full path to 2bit file with -t flag."
                        exit 1
                fi

                ## assign protin file name and the command to prepare it
                PROT=${INPUT%$EXT}aa.fas 
                PROT=${PROT##*/}      
                # TODO: check how genePredToProt handles stops and FS?
                PROT_make_cmd="bedToGenePred $INPUT stdout | genePredToProt -includeStop -starForInframeStops stdin $TWOBIT $PROT"

        elif [ $EXT == "gp" ];
        then
                ## Working with annotation gp file
                echo -e "annotation gp file is provided . . ."
                ## check 2bit file
                if [[ -z $TWOBIT ]] || [[ ! -f $TWOBIT ]];
                then
                        echo "ERROR: 2bit file not found! provide a full path to 2bit file with -t flag."
                        exit 1
                fi
                ## assign protin file name and the command to prepare it
                PROT=${INPUT%$EXT}aa.fas
                PROT=${PROT##*/}
                # TODO: check how genePredToProt handles stops and FS?
                PROT_make_cmd="genePredToProt -includeStop -starForInframeStops $INPUT $TWOBIT $PROT"


        elif [[ $EXT == "fa" ]] || [[ $EXT == "faa" ]] || [[ $EXT == "fasta" ]] || [[ $EXT == "fas" ]];
        then
                ## Working with protein fasta file
                echo -e "protein fasta file is provided . . ."
                PROT=$INPUT
                PROT_make_cmd=""
        else
               echo -e "Provide correct query file: .bed, .gp, .fa, .fasta, .faa, .fas"
               exit 1 
        fi    

        ### Check the prepared protein input file
        if [ -z $PROT ];
        then
                echo -e "Can't make/find protein file: $PROT !"
                exit 1
        fi



        ### Run BUSCO -annotation
        TMPDIR=$(mktemp -d -t busco_XXXXXXXX)
        cd $TMPDIR

        echo "Preparing protein . . ."
        eval $PROT_make_cmd

        echo -e "this is your protein: $PROT"
        echo -e "Running annotation busco now in $TMPDIR . . . "
        echo "Running busco now . . . "
        busco -i $PROT -o busco_anno -l $BUSCOODB -m protein -c $CORES --offline -f
fi




### Create output directory
# Assign output dir if wasn't provided; check is it was provided
if [ -z $OUTPUT ];
then
        OUTPUT="busco_${MODE}"
fi

if [[ ${OUTPUT:0:1} != '/' ]];
then
        # it is a relative path
        FULL_OUTPUT="$WDIR/$OUTPUT"
else
        # it is an absolute path
        FULL_OUTPUT=$OUTPUT
fi

mkdir -p $FULL_OUTPUT
# check id created successfully
if [ ! -d $FULL_OUTPUT ];
then
        "Can not create output directory: $FULL_OUTPUT !"
        exit 1
fi


### Clean up: keeps short_summary.txt and full_table.tsv unless --keeptemp specified
FILE_SUM=$(find -type f -name "short_summary.txt")
FILE_TAB=$(find -type f -name "full_table.tsv")
mv $FILE_SUM $FULL_OUTPUT/
mv $FILE_TAB $FULL_OUTPUT/
cd ..

if [ $KEEPTEMP == "T" ];
then
        echo -e "All intermediate files are in : $TMPDIR"
else
        echo "Cleaning up now . . ."
        rm -rf $TMPDIR
fi

# add param keep temp files
echo -e "All Done!\nResults written in $FULL_OUTPUT"

### Go back
cd $WDIR
