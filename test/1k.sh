#!/bin/bash

ABSPATH=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd)
OUT=$ABSPATH/1k
TREE=$ABSPATH/data/lucas/20k.newick
MSA=$ABSPATH/data/lucas/1k.fasta
mkdir -p $OUT
rm -rf $OUT/*

echo "RUNNING RAXML"
time ../../standard-RAxML/raxmlHPC-SSE3 -f v -G 0.01 -s $MSA -t $TREE -n 1k -m GTRGAMMA -w $OUT

echo "RUNNING EPA"
time ../bin/epa $TREE $MSA -oO -w $OUT

echo "RUNNING JPLACE_COMPARE"
./jplace_compare.py -v $OUT/RAxML_portableTree.1k.jplace $OUT/epa_result.jplace
