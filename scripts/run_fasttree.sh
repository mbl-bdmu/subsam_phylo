#!/usr/bin/env bash
# syntax: ./run_fasttree.sh <alignment.FASTA>

ALNFILE=$1
TREEFILE=${ALNFILE%.fa}.fasttre.nwk

/mnt/c/Users/Guest1/fgmp_working_dir/apps/FastTreeMP -nt -gtr -gamma < $ALNFILE > $TREEFILE

