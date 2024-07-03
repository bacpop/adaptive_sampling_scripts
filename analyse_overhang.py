#!/usr/bin/env python

import sys, getopt, errno
import re
import mappy as mp
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path
import argparse
import gzip
import os
from random import choices

def get_options():
	description = "Aligns and determines enrichment factor"
	parser = argparse.ArgumentParser(description=description,
									 prog='python analyse_RU.py')
	IO = parser.add_argument_group('Input/options.out')
	IO.add_argument('-f',
					required=True,
					help='Input directory containing fastq files.')
	IO.add_argument('-i',
					required=True,
					help='Reference for minimap2 alignment (can be fasta or mmi).')
	IO.add_argument('-c',
					default="1-256",
					help='channels to separate, in the form \"a-b\". '
						 '[Default=1-256]')
	IO.add_argument('-p',
					default=0.8,
					type=float,
					help='Minimum proportion of bases matching reference in alignment. '
						 'Default=0.8')
	IO.add_argument('-o',
					default="RU_output",
					help='Output prefix. '
						 'Default="RU_output"')
	IO.add_argument('-r',
					default=False,
					action="store_true",
					help='Remove multi-mapping reads.'
						 'Default=False')
	IO.add_argument('-b',
					default=False,
					action="store_true",
					help='Align only pass reads.'
						 'Default=False')
	IO.add_argument('-t',
					default=None,
					help='Specify targets within minimap2 index with comma separated list. Default = None ')
	IO.add_argument('-v',
					default=False,
					action="store_true",
					help='Verbose output.'
						 'Default=False')
	return parser.parse_args()

class ReferenceStats:
	def __init__(self, reference):
		self.totalReads = 0
		self.totalLength = 0
		self.reference = reference

def readfq(fp):  # this is a generator function
    """Read FASTA/Q records from file handle
    https://github.com/lh3/readfq/blob/091bc699beee3013491268890cc3a7cbf995435b/readfq.py
    """
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, "".join(seqs)
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

def get_fq(directory):
    types = ([".fastq"], [".fastq", ".gz"], [".fq"], [".fq", ".gz"])
    files = (
        str(p.resolve()) for p in Path(directory).glob("**/*") if p.suffixes in types
    )
    yield from files

def default_val():
    return (None, None, 0, 0, 0)

# determine best loci alignment in each reference
def get_best_map(index, fasta, cutoff=0.7):
    a = mp.Aligner(index, preset="asm10")

    ref_dict = defaultdict(default_val)

    fasta_sequences = SeqIO.parse(open(fasta), 'fasta')
    for fasta in fasta_sequences:
        id, sequence = fasta.id, str(fasta.seq)

        for hit in a.map(sequence):
            query_hit = hit.blen

            # set cutoff for minimum alignment length
            if query_hit < cutoff * len(sequence):
                continue

            if query_hit > ref_dict[hit.ctg][-1]:
                ref_dict[hit.ctg] = (id, hit.r_st, hit.r_en, query_hit)

    return ref_dict

def main():
	options = get_options()

	fastqDir = options.f
	reference = options.i
	channels = options.c
	fields = channels.split("-")
	min_channel = int(fields[0])
	max_channel = int(fields[1])
	matching_prop = options.p
	output = options.o
	remove_multi = options.r
	only_pass = options.b
	verbose = options.v
	mm_target = options.t

	# split list of targets
	if mm_target != None:
		mm_target = set(mm_target.split(","))

	# initialise results dictionaries
	results_dict = {}
	read_seqs = {"target" : {},
				 "non-target" : {}}
	read_lens = {}

	mapper = mp.Aligner(reference, preset="map-ont")

	# keep list of barcodes to ensure printing in order
	barcode_list = []
	barcode_summary = {}

	for f in get_fq(fastqDir):
		if only_pass:
			if "/fastq_pass/" not in f:
				continue

		# only run on gzipped files
		if ".gz" not in f:
			continue

		if verbose:
			print("Aligning: {}".format(str(f)))
		# get filename and extension
		base = os.path.splitext(os.path.basename(f))[0].split("_")
		# print(base)
		if "barcode" in base[2]:
			barcode = base[2]
		else:
			barcode = "NA"
		
		if barcode not in barcode_summary:
			barcode_summary[barcode] = {}
			read_lens[barcode] = {}
			read_lens[barcode]["adaptive"] = []
			read_lens[barcode]["control"] = []

		with gzip.open(f, "rt") as handle:
			input_sequences = SeqIO.parse(handle, 'fastq')
			for entry in input_sequences:
				fields, seq = entry.description, str(entry.seq)
				fields = fields.split()
				readName = fields[0][1:]
				channel = int(fields[3].split("=")[1])
				assert (fields[3].split("=")[0] == "ch")
				assert (channel >= 0 and channel <= 512)


				# determine is read is in target channel
				target = False
				if channel >= min_channel and channel <= max_channel:
					target = True

				# read_lengths[readName] = length
				# read_seqs[readName] = seq

				ref_align = "unaligned"

				# perc id is alignment coverage of reference, just in case clipping present
				perc_id = 0
				match_len = 0
				overhang_len = 0
				read_len = len(seq)

				align_count = 0

				# Map seq, take alignment with largest matching sequence only
				for r in mapper.map(seq):
					align_count += 1

					len_ref_map = abs(r.r_en - r.r_st)

					# take as best alignment if matching length is longer
					if mm_target is None:
						if match_len < r.mlen:
							match_len = r.mlen
							ref_align = r.ctg
							perc_id = match_len / len_ref_map
							overhang_len = read_len - len_ref_map
					# otherwise make sure read aligns to known target
					else:
						if r.ctg in mm_target:
							if match_len < r.mlen:
								match_len = r.mlen
								ref_align = r.ctg
								perc_id = match_len / len_ref_map
								overhang_len = read_len - len_ref_map

				# if below cutoff, or removing multialigning reads, set reference as unaligned
				if perc_id < matching_prop:
					ref_align = "unaligned"
				if remove_multi and align_count > 1:
					ref_align = "unaligned"


				# take length as alignment length
				total_length = len(seq)
				if ref_align == "unaligned":
					match_len = total_length

				# take length as alignment length
				if target:
					read_lens[barcode]["adaptive"].append((ref_align, overhang_len, total_length))
				else:
					read_lens[barcode]["control"].append((ref_align, overhang_len, total_length))

				if ref_align not in barcode_summary[barcode]:
					barcode_summary[barcode][ref_align] = {}
				
				if target not in barcode_summary[barcode][ref_align]:
					barcode_summary[barcode][ref_align][target] = 0

				if overhang_len > barcode_summary[barcode][ref_align][target]:
					barcode_summary[barcode][ref_align][target] = overhang_len

	# write summary file
	with open(output + "_summary.txt", "w") as o_sum:
		o_sum.write("Channel\tBarcode\tAlignment\tOverhang_len\n")
		# create dictionary to determine enrichment
		for barcode, barcode_dict in barcode_summary.items():
			for ref_align, ref_align_dict in barcode_dict.items():
				for target, overhang_len in ref_align_dict.items():
					target_id = "NAS" if target == True else "Control"

					o_sum.write("{}\t{}\t{}\t{}\n".format(target_id, str(barcode), ref_align, str(overhang_len)))

	with open(output + "_reads.txt", "w") as o_reads:
		o_reads.write("Barcode\tChannel\tAlignment\tOverhang_len\tRead_len\n")
		for barcode, barcode_dict in read_lens.items():
			for target, overhang_list in barcode_dict.items():
				target_id = "NAS" if target == "adaptive" else "Control"
				for entry in overhang_list:
					o_reads.write(str(barcode) + "\t" + str(target_id) + "\t" + str(entry[0]) + "\t" + str(entry[1]) + "\t" + str(entry[2]) + "\n")


if __name__ == "__main__":
    main()
