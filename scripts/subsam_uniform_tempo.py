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
size = snakemake.params["size"]
num_intervals = snakemake.params["intervals"]
outfile = snakemake.output[0]

# Extract calendar dates, convert to decimal year.
heads_dates = []
pattern = "([0-9]{4}-[0-9]{2}-[0-9]{2})"
with open(headers, "r") as heads:
	for h in heads:
		h = h.strip("\n")
		d = re.search(pattern, h)
		if d is None:
			print(d, f'Excluding {h}. Date does not follow YYYY-MM-DD format')
			continue
		else:
			d = d.group().split("-")
			d = [int(d) for d in d]
			d = dt.datetime(d[0],d[1],d[2])
			d = pyasl.decimalYear(d)
			heads_dates.append([h,d])

# Set up dataframe of headers and decimal dates
colnames = ["header", "decimal_year"]
tab = pd.DataFrame(heads_dates, index=None, columns=colnames)

# Get range of sample dates, divide into intervals.
deci_yrs = sorted(tab.decimal_year.tolist(), key=None)
first, last = deci_yrs[0], deci_yrs[-1]
date_range_length = last - first
interval_length = date_range_length/num_intervals
intervals = np.arange(first, last, interval_length)
intervals = [[i, i+interval_length] for i in intervals]

# Assign each sequence to time interval.
tab["time_interval"] = tab["decimal_year"].apply(
	assign_interval,
	args=(last, intervals)
	)

# Count number of sequences per interval.
tally = Counter(tab.time_interval.tolist())
uniform_prob_inter_interval = 1/len(tally)
tally_df = pd.DataFrame.from_dict(tally, orient='index', 
							   	  columns=["interval_tally"])
tab = tab.merge(tally_df, how='left', 
				left_on='time_interval', right_index=True)

# Assign a uniform probability for each interval and a uniform 
# probability for each sequence belonging to the same interval. 
# Get the joint probability for each sequence from these two.
tab["uniform_prob_intra_interval"] = tab.interval_tally.apply(lambda x: 1/x)
tab["joint_prob"] = tab.uniform_prob_intra_interval.apply(
	lambda x: x*uniform_prob_inter_interval
	)

# Subsample N sequences according to assigned joint probability.
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

