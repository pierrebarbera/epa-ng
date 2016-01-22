#!/bin/bash

OUT="20k"
TREE="data/lucas/20k.newick"
MSA="data/lucas/20k.fasta"
mkdir -p $OUT
rm -rf $OUT/*

echo "RUNNING RAXML"
../../standard-RAxML/raxmlHPC -f v -s $MSA -t $TREE -n 20k -m GTRGAMMA -w $OUT

echo "RUNNING EPA"
../bin/epa $TREE $MSA -O -w $OUT

echo "RUNNING JPLACE_COMPARE"
./jplace_compare.py -v $OUT/RAxML_portableTree.20k.jplace $OUT/epa_result.jplace
