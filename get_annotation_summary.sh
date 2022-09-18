#!/usr/bin/env bash
#

db=$1
ANNODIR=$2
ANNO=$ANNODIR/$db.final_anno.gp

TRANSCRIPTS=$(wc -l $ANNO | awk '{print $1}')
GENES=$(genePredSingleCover $ANNO stdout | wc -l)

echo -e "$db\t$TRANSCRIPTS\ttranscripts\tstats"
echo -e "$db\t$GENES\tgenes\tstats"

BUSCO=$(grep -P "^\t[0-9]" $ANNODIR/busco_${db}/short_summary.txt | awk '{print $1}')
dblist="$db\n$db\n$db\n$db\n$db\n$db"
blist="busco\nbusco\nbusco\nbusco\nbusco\nbusco"
classes="complete\nsingle\nduplicated\nfragmented\nmissing\ntotal"

paste <(echo -e $dblist) <(echo -e $BUSCO | tr ' ' '\n') <(echo -e $classes) <(echo -e $blist)


