#!/usr/bin/env python

import random
import sys

# Read input headers and size of random sample.
headers = snakemake.input[0]
size = snakemake.params["size"]
outfile = snakemake.output[0]

# Randomly sample specified number of headers. 
with open(headers, "rU") as heads:
	headers = heads.read().splitlines()
rand_headers = random.sample(headers, int(size))
rand_headers = [h+"\n" for h in rand_headers]
rand_headers = "".join(rand_headers)

# Write to file. 
with open(outfile, "w") as out:
	out.write(rand_headers)
