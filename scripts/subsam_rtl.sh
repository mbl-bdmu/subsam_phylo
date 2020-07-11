#!/usr/bin/env bash
NWKFILE=$1
OUTNAME=$2
RTLNUM=$3
CPUNUM=$4

python3 /usr/bin/Treemmer_v0.3.py $NWKFILE -RTL $RTLNUM --cpu $CPUNUM --verbose 1 --plot_always &&
mv ${NWKFILE}_trimmed_list_RTL_${RTLNUM} $OUTNAME

## ADD OPTION TO PROTECT SPECIFIC HEADERS FROM TRIM