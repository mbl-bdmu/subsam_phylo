#!/usr/bin/env python

from collections import Counter
import datetime as dt
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
size = snakemake.params["size"]
num_intervals = snakemake.params["intervals"]
outfile = snakemake.output[0]

# Exclude headers without deme or location metadata.
filt_heads = []
deme_colnames = ["header", "deme"]
demes = pd.read_csv(demes, sep="\t", header=None, names=deme_colnames)
with open(headers, "r") as heads:
	for h in heads:
		h = h.strip("\n")
		if h in demes["header"].tolist():
			filt_heads.append(h)
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

# Get range of sample dates, divide into intervals.
deci_yrs = sorted(tab.decimal_year.tolist(), key=None)
first, last = deci_yrs[0], deci_yrs[-1]
date_range_length = last - first
interval_length = date_range_length/num_intervals
intervals = np.arange(first, last, interval_length)
intervals = [[i, i+interval_length] for i in intervals]

# Assign each sequence to a time interval and to a deme.
tab["time_interval"] = tab["decimal_year"].apply(
	assign_interval,
	args=(last, intervals)
	)
tab = tab.merge(demes, how='inner', 
				left_on='header', right_on='header')

# Count number of sequences present per deme per interval.
tally = Counter(tab.time_interval.tolist())
uniform_prob_interval = 1/len(tally)
tally_interval = pd.DataFrame.from_dict(tally, orient='index', 
							   	  columns=["interval_tally"])
tab = tab.merge(tally_interval, how='left', 
				left_on='time_interval', right_index=True)
uniq_demes = set(demes["deme"].tolist())
interval_deme = [series.time_interval+"_"+series.deme for index,series in tab.iterrows()]
tally_interval_deme = Counter(interval_deme)
tally_interval_deme = [i+"_"+str(d) for i,d in tally_interval_deme.items()]
tally_interval_deme  = [s.split("_") for s in tally_interval_deme]
interval_deme_colnames = ["time_interval", "deme", "interval_deme_tally"]
tally_interval_deme = pd.DataFrame(tally_interval_deme, index=None, columns=interval_deme_colnames)
tab = tab.merge(tally_interval_deme, how='left', on=["time_interval", "deme"])

# Count number of demes present per interval.
tally_deme = tally_interval_deme.groupby(["time_interval"])["deme"].apply(lambda x: ','.join(x))
tally_deme = tally_deme.reset_index()
tally_deme["deme_tally"] = tally_deme.deme.apply(lambda x: len(x.split(",")))
tally_deme = tally_deme[["time_interval", "deme_tally"]]
tab = tab.merge(tally_deme, how='left', on="time_interval")

# Assign a uniform probability for each interval group, a uniform probability 
# for each deme group per interval, and a uniform probability for each sequence
# per deme group per interval. Get the joint probability for each sequence 
# from these three.
tab["uniform_prob_interval_deme"] = tab.interval_deme_tally.apply(lambda x: 1/int(x))
tab["uniform_prob_deme"] = tab.deme_tally.apply(lambda x: 1/int(x))
tab["joint_prob"] = tab["uniform_prob_interval_deme"]*tab["uniform_prob_deme"]
tab["joint_prob"] = tab.joint_prob.apply(lambda x: x*uniform_prob_interval)

# Subsample N sequences headers according to assigned joint probability.
weights = tab.joint_prob.tolist()
sequence_headers = tab.header.tolist()
sampled_headers = np.random.choice(sequence_headers, size, 
								   replace=False, p=weights)
sampled_headers = sampled_headers.tolist()
sampled_headers_out = [h for h in sequence_headers if h in sampled_headers]
sampled_headers_out = pd.DataFrame(sampled_headers_out, 
								   index=None, columns=None)

# Save headers to file
sampled_headers_out.to_csv(outfile, index=False, header=False)

