#!/usr/bin/env python

import sys, getopt, errno
import re
from collections import defaultdict
from pathlib import Path
import argparse
import gzip
import os
from random import choices
from datetime import timedelta

# non-standard libraries
from dateutil.parser import parse as dparse
from tqdm import tqdm
from Bio import SeqIO

def get_options():
	description = "Splits reads in fastq files based on run time."
	parser = argparse.ArgumentParser(description=description,
									 prog='python analyse_RU.py')
	IO = parser.add_argument_group('Input/options.out')
	IO.add_argument('--dir',
					type=str,
					required=True,
					help='Input directory containing unzipped fastq files.')
	IO.add_argument('--start',
					type=str,
					required=True,
					help='Start time of sequencing run, found in final_summary.txt file (e.g. 2023-04-20T14:29:36.956327+01:00).')
	IO.add_argument('--end',
					type=str,
					required=True,
					help='End time of sequencing run, found in final_summary.txt file (e.g. 2023-04-21T14:30:40.938064+01:00).')
	IO.add_argument('--output',
					type=str,
					default="parsed_output",
					help='Output prefix. '
						 'Default="RU_output"')
	IO.add_argument('--breaks',
					type=int,
					default=6,
					help='Number of time breaks for parsing. '
						'Default=6')
	return parser.parse_args()


def get_fq(directory):
    types = ([".fastq"], [".fastq", ".gz"], [".fq"], [".fq", ".gz"])
    files = (
        str(p.resolve()) for p in Path(directory).glob("**/*") if p.suffixes in types
    )
    yield from files


def main():
	options = get_options()

	fastqDir = options.dir
	init_time = dparse(options.start)
	end_time = dparse(options.end)
	output = options.output
	num_breaks = options.breaks
    
	time_diff = (end_time - init_time) / num_breaks
	
	# generate time list
	time_list = []
	min_time = init_time

	for i in range(1, num_breaks + 1):
		max_time = init_time + (time_diff * i)
		time_list.append((min_time, max_time))
		min_time = max_time

	for f in get_fq(fastqDir):
		
		with open(f, "rb") as f1:
			num_entries = int(sum(1 for _ in f1) / 4)

		# generate seq list for all sequences
		seq_list = [[] for i in range(num_breaks)]
		
		# get filename and extension
		base = os.path.splitext(os.path.basename(f))[0].split("_")

		print("Analysing file: {}".format(os.path.basename(f)))

		if "barcode" in base[1]:
			barcode = base[1]
		else:
			barcode = "NA"

		channel_type = "NA"

		if "adaptive" in f:
			channel_type = "adaptive"
		elif "control" in f:
			channel_type = "control"
		
		for entry in tqdm(SeqIO.parse(open(f),"fastq"), total=num_entries):
			fields = entry.description
			read_time = dparse([i for i in fields.split() if i.startswith('start_time')][0].split('=')[1])
			channel = int([i for i in fields.split() if i.startswith('ch')][0].split('=')[1])

			split_fields = fields.split()
			readName = split_fields[0][1:]
			assert (channel >= 0 and channel <= 512)

			# place in time bucket
			for index, time_range in enumerate(time_list):
				if time_range[0] <= read_time < time_range[1]:
					seq_list[index].append(entry)

		for index, seq_break_list in enumerate(seq_list):
			SeqIO.write(seq_break_list, os.path.join(output, "bc_" + str(barcode) + "_ch_" + channel_type + "_time_" + str(index + 1) + ".fastq"), "fastq")

if __name__ == "__main__":
    main()