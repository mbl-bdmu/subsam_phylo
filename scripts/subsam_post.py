#!/usr/bin/env python

import datetime as dt
import random
import re

import pandas as pd
from PyAstronomy import pyasl


# Read headers names and size of uniform subsample.
headers = snakemake.input[0]
demes = snakemake.input[1]
fixed_percent = snakemake.params[0]
outfile = snakemake.output[0]

# Exclude headers without deme or location metadata. Of headers with 
# deme, assign the deme to the header.
heads_with_demes = []
deme_colnames = ["header", "deme"]
demes = pd.read_csv(demes, sep="\t", header=None, names=deme_colnames)
with open(headers, "r") as heads:
	for h in heads:
		h = h.strip("\n")
		if h in demes["header"].tolist():
			d = demes[demes.header.str.contains(h)]["deme"].tolist()[0]
			heads_with_demes.append([h, d])

# Among headers with demes, extract calendar dates & convert to decimal
# year; exclude headers without any dates or whose dates do not follow
# the format 'YYYY-MM-DD'. Downstream phylodynamic analysis with BEAST2
# would require all sequences to have dates.
heads_with_demes_and_dates = []
pattern = "([0-9]{4}-[0-9]{2}-[0-9]{2})"
for hwd in heads_with_demes:
	h = hwd[0]
	D = hwd[1]
	d = re.search(pattern, h)
	if d is None:
		print(f'Excluding {h}. Date does not follow YYYY-MM-DD format')
		continue
	else:
		heads_with_demes_and_dates.append([h,D])
tab = pd.DataFrame(heads_with_demes_and_dates[:], index=None, columns=deme_colnames)

# Randomly select set % of headers from each deme stratum.
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

