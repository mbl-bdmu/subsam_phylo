#!/usr/bin/env python

import random
import sys

import pandas as pd

# Read input headers, demes, and max size of random sample per deme.
headers = snakemake.input[0]
locations = snakemake.input[1]
size = snakemake.params["size"]
outfile = snakemake.output[0]

#headers = "GEO_1042.headers.txt"
#locations= "GEO_1042.loc.tsv"
#size = 34
#outfile = "test.headers.txt"

# Read all headers.
with open(headers, "r") as heads:
	headers = heads.read().splitlines()

# Filter out headers with no location.
loc = pd.read_csv(locations, sep="\t", header=None)
headers_with_deme = [h for h in headers if h in loc[0].tolist()]
headers_no_deme = [h for h in headers if h not in loc[0].tolist()]
for h in headers_no_deme:
	print(f'Excluding {h}, no deme.')

# Get unique set of locations, make into dictionary.
demes_headers_dict = {d:[] for d in list(set(loc[1].tolist()))}

# Loop filtered headers, look up its deme in loc_dict.
loc_dict = loc.set_index(0, drop=True, inplace=False).to_dict()
for h in headers_with_deme:
	h_deme = loc_dict[1][h]
	# Append header to dictionary by deme.
	for d in demes_headers_dict:
		if h_deme == d:
			demes_headers_dict[d].append(h)

# For sequences under each deme, randomly sample to max size from parameters. 
samples = []
for d in demes_headers_dict:
	sample = random.sample(demes_headers_dict[d], int(size))
	samples.extend(sample)
	print(d, len(sample), sample, "\n")
print("Total", len(samples), samples, "\n")
samples = "".join([h+"\n" for h in samples])

# Write to file. 
with open(outfile, "w") as out:
	out.write(samples)
