#!/bin/bash

OUT="lou"

mkdir -p $OUT
rm -rf $OUT/*
./leave_one_out.py ../../standard-RAxML/raxmlHPC ../bin/epa data/ref.tre data/combined.fasta $OUT $1
