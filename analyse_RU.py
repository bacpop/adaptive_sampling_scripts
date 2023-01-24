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

def get_options():
	description = "Aligns and determines enrichment factor"
	parser = argparse.ArgumentParser(description=description,
									 prog='python analyse_RU.py')
	IO = parser.add_argument_group('Input/options.out')
	IO.add_argument('-f',
					required=True,
					help='Input directory containing fastq files')
	IO.add_argument('-i',
					required=True,
					help='Reference for minimap2 alignment (can be fasta or mmi)')
	IO.add_argument('-c',
					default="1-256",
					help='channels to separate, in the form \"a-b\".'
						 '[Default = 1-256]')
	IO.add_argument('-p',
					default=0.8,
					type=float,
					help='Minimum proportion of bases matching reference in alignment'
						 'Default=0.8')
	IO.add_argument('-o',
					default="RU_output",
					help='Output prefix'
						 'Default="RU_output"')
	IO.add_argument('-r',
					default=False,
					action="store_true",
					help='Remove multi-mapping reads')
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

def AnalyseForChannels(_readsForChannels, _samFilename, multi_align, matching_prop=0):
	ref_dict = dict()
	ref_dict["unaligned"] = []
	ReferenceStatDict = dict()
	ReferenceStatDict["unaligned"] = ReferenceStats("unaligned")
	try:
		with open(_samFilename, 'r') as infile:
			for line in infile:
				fields = line.split()
				if fields[0] in _readsForChannels:
					reference = fields[2]
					sequenceLength = read_lengths[fields[0]]
					identity = 0
					if reference == "*":
						stats = ReferenceStatDict["unaligned"]
						ref_dict["unaligned"].append((fields[0], reference, identity))
						stats.totalReads += 1
						stats.totalLength += sequenceLength
						continue
					CIGAR = fields[5]
					num_matches = 0
					num_mismatches = 0
					matches = re.findall(r'(\d+)([A-Z]{1})', CIGAR)
					# get number of M and I/D in CIGAR
					for m in matches:
						if m[1] == "M":
							num_matches += int(m[0])
						elif m[1] == "I" or m[1] == "D":
							num_mismatches += 1
					# get number of non-identical 'matches' in CIGAR
					num_mismatches = int(fields[11].split(":")[-1]) - num_mismatches
					# get number of identical bases in alignment
					num_matches -= num_mismatches
					# pass alignment if proportion of matching bases in below threshold
					identity = num_matches / sequenceLength
					if identity < matching_prop or fields[0] in multi_align:
						stats = ReferenceStatDict["unaligned"]
						ref_dict["unaligned"].append((fields[0], reference, identity))
						stats.totalReads += 1
						stats.totalLength += sequenceLength
						continue
					if reference not in ReferenceStatDict.keys():
						ReferenceStatDict[reference] = ReferenceStats(reference)
						ref_dict[reference] = []
					ref_dict[reference].append((fields[0], reference, identity))
					stats = ReferenceStatDict[reference]
					stats.totalReads += 1
					stats.totalLength += sequenceLength
	except (OSError, IOError) as e:
		if getattr(e, 'errno', 0) == errno.ENOENT:
			print("Could not find file " + infile)
		print(e)
		sys.exit(2)

	return ReferenceStatDict, ref_dict

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

	results_dict = {}

	# target_channel_reads = list()
	# non_target_channel_reads = list()
	#read_lengths = dict()
	read_seqs = {"target" : {},
				 "non-target" : {}}

	mapper = mp.Aligner(reference, preset="map-ont")

	# store multi-aligned bases
	multi_align = set()

	for f in get_fq(fastqDir):
		if f.endswith(".gz"):
			fopen = gzip.open
		else:
			fopen = open

		# get filename and extension
		base = os.path.splitext(os.path.basename(f))[0].split("_")
		# print(base)
		if "barcode" in base[2]:
			barcode = base[2]
			file_id = "_".join([base[1], base[2]])
		else:
			barcode = "NA"
			file_id = "_".join([base[1], "NA"])

		if barcode not in results_dict:
			results_dict[barcode] = {"target_channel_bases": 0,
									 "non_target_channel_bases" : 0,
									 "target_channel_reads" : 0,
									 "non_target_channel_reads" : 0,
									 "target_channel_bases_mapped": 0,
									 "non_target_channel_bases_mapped": 0,
									 "target_channel_reads_mapped": 0,
									 "non_target_channel_reads_mapped": 0,
									 "total_bases" : 0,
									 "total_reads" : 0,
									 "ref_dict" : {}}
			read_seqs["target"][barcode] = {}
			read_seqs["non-target"][barcode] = {}


		input_sequences = SeqIO.parse(fopen(f), 'fasta')
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

			align_count = 0

			# Map seq, take alignment with largest matching sequence only
			for r in mapper.map(seq):
				align_count += 1

				len_ref_map = r.r_en - r.r_st

				# take as best alignment if matching length is longer
				if match_len < r.mlen:
					ref_align = r.ctg
					perc_id = r.mlen / len_ref_map

			# add to read_seqs
			if ref_align not in read_seqs["target"][barcode]:
				read_seqs["target"][barcode][ref_align] = []
				read_seqs["non-target"][barcode][ref_align] = []

			if target:
				read_seqs["target"][barcode][ref_align].append((readName, perc_id, seq))
			else:
				read_seqs["non-target"][barcode][ref_align].append((readName, perc_id, seq))

			# if below cutoff, or removing multialigning reads, set reference as unaligned
			if perc_id < matching_prop:
				ref_align = "unaligned"
			if remove_multi and align_count > 1:
				ref_align = "unaligned"

			length = len(seq)

			if ref_align not in results_dict[barcode]["ref_dict"]:
				results_dict[barcode]["ref_dict"][ref_align] = {"target_channel_bases": 0,
																"non_target_channel_bases" : 0,
																"target_channel_reads" : 0,
																"non_target_channel_reads" : 0}

			if target:
				results_dict[barcode]["target_channel_bases"] += length
				results_dict[barcode]["target_channel_reads"] += 1
				results_dict[barcode]["ref_dict"][ref_align]["target_channel_bases"] += length
				results_dict[barcode]["ref_dict"][ref_align]["target_channel_reads"] += 1
				if ref_align != "unaligned":
					results_dict[barcode]["target_channel_bases_mapped"] += length
					results_dict[barcode]["target_channel_reads_mapped"] += 1
			else:
				results_dict[barcode]["non_target_channel_bases"] += length
				results_dict[barcode]["non_target_channel_reads"] += 1
				results_dict[barcode]["ref_dict"][ref_align]["non_target_channel_bases"] += length
				results_dict[barcode]["ref_dict"][ref_align]["non_target_channel_reads"] += 1
				if ref_align != "unaligned":
					results_dict[barcode]["non_target_channel_bases_mapped"] += length
					results_dict[barcode]["non_target_channel_reads_mapped"] += 1

			results_dict[barcode]["total_bases"] += length
			results_dict[barcode]["total_reads"] += 1



	# write summary file
	with open(output + "_summary.txt", "w") as o_sum:
		o_sum.write("Statistic\tChannel\tAlignment\tValue\n")
		# create dictionary to determine enrichment
		enrichment_dict = {}

		# target channels
		print("Reference stats for channels " + str(channels) + ": ")
		for barcode in results_dict.keys():
			print("Barcode: " + str(barcode) + ": ")
			for ref, item in results_dict[barcode]["ref_dict"].items():
					print(ref + "\t" + str(item["target_channel_reads"]) + "\t" + str(item["target_channel_bases"]))
					o_sum.write("Reads_mapped\t{}\t{}\t{}\n".format("Target", ref, str(item["target_channel_reads"])))
					o_sum.write("Bases_mapped\t{}\t{}\t{}\n".format("Target", ref, str(item["target_channel_bases"])))
					enrichment_dict[ref] = {}
					enrichment_dict[ref]["Target_bases_mapped"] = item["target_channel_bases"]

			for barcode, entry in read_seqs["target"].items():
				for ref, read_list in entry:
					with open(output + "_" + str(barcode) + "_target_" + str(ref) + ".fasta", "w") as o:
						for read_tup in read_list:
							read_id, identity, seq = read_tup
							o.write(">" + read_id + "\t" + reference + "\t" + str(identity) + "\n" + seq + "\n")

			print("Total number of reads mapped: " + str(results_dict[barcode]["target_channel_reads_mapped"]) + "/" + str(results_dict[barcode]["target_channel_reads"]))
			print("Total read bases: " + str(results_dict[barcode]["target_channel_bases"]))
			print("Total read bases mapped: " + str(results_dict[barcode]["target_channel_bases_mapped"]))

			# iterate over enrichment dict, calculate enrichment for total bases
			for key in enrichment_dict:
				enrichment_dict[key]["Target_prop_bases"] = enrichment_dict[key]["Target_bases_mapped"] / results_dict[barcode]["target_channel_bases"]

			# write to summary file
			o_sum.write("Reads_total\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["target_channel_reads"])))
			o_sum.write("Reads_mapped\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["target_channel_reads_mapped"])))
			o_sum.write("Bases_total\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["target_channel_bases"])))
			o_sum.write("Bases_mapped\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["target_channel_bases_mapped"])))

		# non target channels
		print("\nReference stats for all other channels: ")
		for barcode in results_dict.keys():
			print("Barcode: " + str(barcode) + ": ")
			for ref, item in results_dict[barcode]["ref_dict"].items():
					print(ref + "\t" + str(item["non_target_channel_reads"]) + "\t" + str(item["non_target_channel_bases"]))
					o_sum.write("Reads_mapped\t{}\t{}\t{}\n".format("Target", ref, str(item["non_target_channel_reads"])))
					o_sum.write("Bases_mapped\t{}\t{}\t{}\n".format("Target", ref, str(item["non_target_channel_bases"])))
					if ref.reference not in enrichment_dict:
						enrichment_dict[ref.reference] = {}
					enrichment_dict[ref.reference]["Nontarget_bases_mapped"] = item["non_target_channel_bases"]

			for barcode, entry in read_seqs["non-target"].items():
				for ref, read_list in entry:
					with open(output + "_" + str(barcode) + "_nontarget_" + str(ref) + ".fasta", "w") as o:
						for read_tup in read_list:
							read_id, identity, seq = read_tup
							o.write(">" + read_id + "\t" + reference + "\t" + str(identity) + "\n" + seq + "\n")

			print("Total number of reads mapped: " + str(results_dict[barcode]["non_target_channel_reads_mapped"]) + "/" + str(results_dict[barcode]["non_target_channel_reads"]))
			print("Total read bases: " + str(results_dict[barcode]["non_target_channel_bases"]))
			print("Total read bases mapped: " + str(results_dict[barcode]["non_target_channel_bases_mapped"]))

			# calculate enrichment for all entries with mappings in target and non-target
			for key in enrichment_dict:
				if "Nontarget_bases_mapped" in enrichment_dict[key] and "Target_prop_bases" in enrichment_dict[key]:
					non_target_prop_bases = enrichment_dict[key]["Nontarget_bases_mapped"] / results_dict[barcode]["non_target_channel_bases"]
					enrichment = enrichment_dict[key]["Target_prop_bases"] / non_target_prop_bases

					o_sum.write("Enrichment\t{}\t{}\t{}\n".format("NA", key, str(enrichment)))

			# write to summary file
			o_sum.write("Reads_total\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["non_target_channel_reads"])))
			o_sum.write("Reads_mapped\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["non_target_channel_reads_mapped"])))
			o_sum.write("Bases_total\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["non_target_channel_bases"])))
			o_sum.write("Bases_mapped\t{}\t{}\t{}\n".format("Target", "Total", str(results_dict[barcode]["non_target_channel_bases_mapped"])))