#!/usr/bin/env bash

HEADS=$1
NWKFILE=$2
ALN=$3
OUTDIR=$(dirname $1)

DATES=$4
REG=$(basename $5)

# Re-extract dates from tips.
export headers=$HEADS
export dates=$DATES
python - <<EOF
import datetime as dt
import os
import re
import sys

import pandas as pd
from PyAstronomy import pyasl

headers = os.getenv('headers') #snakemake.input[0]
dates = os.getenv('dates') #snakemake.output[0]

# Extract calendar dates, convert to decimal year.
heads_dates = []
pattern = "([0-9]{4}-[0-9]{2}-[0-9]{2})"
with open(headers, "r") as heads:
	for h in heads:
		h = h.strip("\n")
		d = re.search(pattern, h).group().split("-")
		d = [int(d) for d in d]
		d = dt.datetime(d[0],d[1],d[2])
		d = pyasl.decimalYear(d)
		heads_dates.append([h,d])

# Set up dataframe of headers and decimal dates
deciyr_colnames = ["name", "date"]
tab = pd.DataFrame(heads_dates, index=None, columns=deciyr_colnames)
tab.to_csv(dates, index=False, header=True, sep=",")
EOF

# Calculate root-to-tip divergence, plot RTT and tree.
CMD="treetime clock --tree $NWKFILE --dates $DATES --aln $ALN --outdir $OUTDIR --plot-rtt $REG --covariation --clock-filter"
echo "Running $CMD"; eval $CMD

# TO DOs: 
# Use snakemake's logging feature to properly log estimated rate of change, 
# as this file will be used in plotting treetime's rerooted phylo in a time scale. 
