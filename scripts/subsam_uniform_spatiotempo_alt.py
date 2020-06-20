#!/usr/bin/env python

import datetime as dt
import random
import re

from PyAstronomy import pyasl
import numpy as np
import pandas as pd


def assign_interval(date, last, intervals):
	"""Assign a date to a time interval.
	"""
	if (date == last):
		value = "-".join(str(_) for _ in intervals[-1])
	else:
		for i in intervals:
			if (date >= i[0]) & (date < i[1]):
				value = "-".join(str(_) for _ in i)
	return value


# Read headers names and size of uniform subsample.
headers = snakemake.input[0]
demes = snakemake.input[1]
outfile = snakemake.output[0]

# Exclude headers without deme or location metadata. Of the included 
# headers, count the number of demes that have an included header.
filt_heads = []
filt_demes = []
deme_colnames = ["header", "deme"]
demes = pd.read_csv(demes, sep="\t", header=None, names=deme_colnames)
with open(headers, "r") as heads:
	for h in heads:
		h = h.strip("\n")
		if h in demes["header"].tolist():
			filt_heads.append(h)
			d = demes[demes.header.str.contains(h)]["deme"]
			d = d.tolist()[0]
			filt_demes.append(d)
heads = filt_heads

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
deciyr_colnames = ["header", "decimal_year"]
tab = pd.DataFrame(heads_dates, index=None, columns=deciyr_colnames)

# Get range of sample dates, divide into intervals, 
# following Hidano & Gates method for getting intervals.
deci_yrs = sorted(tab.decimal_year.tolist(), key=None)
first, last = deci_yrs[0], deci_yrs[-1]
date_range_length = last - first
Tsamp = date_range_length
nr = len(set(filt_demes))
nt = len(filt_heads)
ni = nt/nr+1
interval_length = Tsamp/ni
intervals = np.arange(first, last, interval_length)
intervals = [[i, i+interval_length] for i in intervals]

# Assign each sequence to a time interval and to a deme.
tab["time_interval"] = tab["decimal_year"].apply(
	assign_interval,
	args=(last, intervals)
	)
tab = tab.merge(demes, how='inner', 
				left_on='header', right_on='header')

# Randomly select 1 header from each interval-deme stratum.
tab_strats = tab.groupby(["time_interval", "deme"]).agg({"header": lambda x: ','.join(x)}).reset_index()
uniform_sampled_heads = []
tab_strats["header"].apply(lambda x: uniform_sampled_heads.append(random.sample(x.split(','), 1)))
sampled_headers_out = pd.DataFrame(uniform_sampled_heads, 
								   index=None, columns=None)

# Save headers to file
sampled_headers_out.to_csv(outfile, index=False, header=False)
