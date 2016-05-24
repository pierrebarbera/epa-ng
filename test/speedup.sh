#!/bin/bash

ABSPATH=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
OUT=$ABSPATH/speedup
LOG=$OUT/log
TREE=$ABSPATH/data/ref.tre
MSA=$ABSPATH/data/combined.fasta
REF_MSA=$ABSPATH/data/aln.fasta
QRY_MSA=$ABSPATH/data/query.fasta
mkdir -p $OUT
rm -rf $OUT/*
touch $LOG

for i in 1 2 4 8 16 32
do
	export OMP_NUM_THREADS=$i
	echo "$i THREADS:" >> $LOG
	{ time ../bin/epamk -t $TREE -s $REF_MSA -q $QRY_MSA -w $OUT > /dev/null 2>&1 ; } 2>> $LOG
done
