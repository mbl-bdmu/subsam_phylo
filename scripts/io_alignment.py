#!/usr/bin/env python

import sys

from pyfaidx import Fasta

# Read input alignment and inclusion list.
alignment_file = snakemake.input[0]
inclusion_file = snakemake.input[1]
outfile_name = snakemake.output[0]
orig_aln = Fasta(
	alignment_file, sequence_always_upper=True,
    read_long_names=False
    )

# Loop sequences. Skip those not in inclusion list. 
with open(inclusion_file, "r") as inclu:
	inclusion_file = inclu.read().splitlines()
lengths = []
out_multifasta = []
for record in orig_aln:
	if record.name in inclusion_file:
		record_len = len(record)
		lengths.append(record_len)
		filt_record = ">"+str(record.name)+"\n"+str(record)+"\n"
		out_multifasta.append(filt_record)
out_multifasta = "".join(out_multifasta)

# Write to file. Ensure all are equal length.
if all(l == lengths[0] for l in lengths):
	with open(outfile_name, "w") as out:
		out.write(out_multifasta)
else:
	sys.exit("Error: not all sequences in alignment are of same length.")
 
