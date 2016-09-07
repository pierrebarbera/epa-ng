#!/bin/bash

ABSPATH=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
OUT=$ABSPATH/runtime_heur
LOG=$OUT/log
TREE=$ABSPATH/data/lucas/20k.newick
MSA=$ABSPATH/data/lucas/1k_100.fasta
REF_MSA=$ABSPATH/data/lucas/1k_reference.fasta
QRY_MSA=$ABSPATH/data/lucas/1k_query_100.fasta
mkdir -p $OUT
rm -rf $OUT/*
touch $LOG

: <<'END'

echo "EPA" >> $LOG
for i in {1..5}
do
	{ time ../bin/epamk -t $TREE -s $REF_MSA -q $QRY_MSA -G 0.02 -w $OUT > /dev/null 2>&1 ; } 2>> $LOG
done

exit

echo "RAXML" >> $LOG
for i in 1
do
	rm -rf $OUT/RAxML*
        { time ../../standard-RAxML/raxmlHPC-AVX -f v -H -s $MSA -G 0.02 -t $TREE -n 1k -m GTRGAMMA -w $OUT > /dev/null 2>&1 ; } 2>> $LOG
done

END

echo "PPLACER" >> $LOG
for i in {1..5}
do
        { time ./pplacer -t $TREE -s 1k/RAxML_info.1k --strike-box 0 --max-strikes 20 --no-pre-mask -j 0 --out-dir $OUT $MSA > /dev/null 2>&1 ; } 2>> $LOG
done

