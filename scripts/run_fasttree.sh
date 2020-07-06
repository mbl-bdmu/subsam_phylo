#!/usr/bin/env bash
# syntax: ./run_fasttree.sh <alignment.FASTA>

ALNFILE=$1
TREEFILE=${ALNFILE%.fa}.fasttre.nwk

fasttree -nt -gtr -gamma < $ALNFILE > $TREEFILE

