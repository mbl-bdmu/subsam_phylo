#!/usr/bin/env python

import random

import pandas as pd


# Read headers names and size of uniform subsample.
headers = snakemake.input[0]
demes = snakemake.input[1]
fixed_percent = snakemake.params[0]
outfile = snakemake.output[0]

# Exclude headers without deme or location metadata. Of headers with 
# deme, assign the deme to the header. Create dataframe from head, deme.
heads_demes = []
deme_colnames = ["header", "deme"]
demes = pd.read_csv(demes, sep="\t", header=None, names=deme_colnames)
with open(headers, "r") as heads:
	for h in heads:
		h = h.strip("\n")
		if h in demes["header"].tolist():
			d = demes[demes.header.str.contains(h)]["deme"].tolist()[0]
			heads_demes.append([h, d])
tab = pd.DataFrame(heads_demes[:], index=None, columns=deme_colnames)

# Randomly select fixed_percent% of headers from each deme stratum.
tab_strats = tab.groupby(["deme"]).agg({"header": lambda x: ','.join(x)}).reset_index()
post_subsampled_heads = []
for i, h in tab_strats.iterrows():
	heads_in_deme = h[1].split(',')
	k = round(len(heads_in_deme)*fixed_percent)
	post_subsampled_heads.extend(random.sample(heads_in_deme, k))
sampled_headers_out = pd.DataFrame(post_subsampled_heads, 
								   index=None, columns=None)

# Save headers to file
sampled_headers_out.to_csv(outfile, index=False, header=False)

