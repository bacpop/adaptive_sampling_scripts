#!/usr/bin/env python

import sys, getopt, errno
import re

target_channel_reads = list()
non_target_channel_reads = list()
read_lengths = dict()
read_seqs = dict()

class ReferenceStats:
	def __init__(self, reference):
		self.totalReads = 0
		self.totalLength = 0
		self.reference = reference


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

remove_multi = False

try:
	opts, args = getopt.getopt(sys.argv[1:],"hf:s:c:p:o:r:")
except getopt.GetoptError:
	print("Option not recognised.")
	print("python analyse_RU.py -f <fastq file> -s <sam file> -c <channels> -p <min. proportion matching bases")
	print("python analyse_RU.py -h for further usage instructions.")
	sys.exit(2)
for opt, arg in opts:
	if opt == "-h":
		print("analyse_RU.py -i <inputfile>")
		print("-f <fastq file>\t\t\t\t Name of fastq file containing all reads.")
		print("-s <sam file>\t\t\t Name of sam file containing all mappings to references.")
		print("-c <channels>\t channels to separate, in the form \"a-b\".")
		print("-p <proportion>\t Minimum proportion of bases matching reference in alignment")
		print("-o <output-pref>\t Output prefix")
		print("-r <remove-multimappers>\t Remove multi-mapping reads")
		sys.exit()
	elif opt in ("-f"):
		fastqFile = arg
	elif opt in("-s"):
		samFile = arg
	elif opt in ("-c"):
		channels=arg
		fields=arg.split("-")
		min_channel = int(fields[0])
		max_channel = int(fields[1])
	elif opt in ("-p"):
		matching_prop = float(arg)
	elif opt in ("-o"):
		output = str(arg)
	elif opt in ("-r"):
		if str(arg) in ("true", "True", "T", "t"):
			remove_multi = True

if fastqFile == '' or samFile == '' or channels == '':
	print("You must specify -f and -s and -c")
	sys.exit(2)

#split the reads by channel
try:
	with open(fastqFile, 'r') as infile:
		target_channel_bases = 0
		non_target_channel_bases = 0
		count = 0
		for line in infile:
			if count % 4 == 0:
				assert(line[0]=='@')
				fields = line.split()
				readName = fields[0][1:]
				channel = int(fields[3].split("=")[1])
				assert(fields[3].split("=")[0] == "ch")
				assert(channel >= 0 and channel <= 512)
				if channel >= min_channel and channel <= max_channel:
					target_channel_reads.append(readName)
				else:
					non_target_channel_reads.append(readName)
			if count % 4 == 1:
				length = len(line.strip())
				read_lengths[readName] = length
				read_seqs[readName] = line.strip()
				if channel >= min_channel and channel <= max_channel:
					target_channel_bases += length
				else:
					non_target_channel_bases += length
			count += 1
except (OSError, IOError) as e:
	if getattr(e, 'errno', 0) == errno.ENOENT:
		print("Could not find file " + fastqFile)
	print(e)
	sys.exit(2)

# identify multi aligning read_seqs
multi_align = set()
if remove_multi == True:
	single_align = set()
	with open(samFile, 'r') as infile:
		for line in infile:
			fields = line.split()
			if fields[0] not in single_align:
				single_align.add(fields[0])
			else:
				multi_align.add(fields[0])

# target channels
dictionary, output_dict = AnalyseForChannels(target_channel_reads, samFile, multi_align, matching_prop)
print("Reference stats for channels " + str(channels) + ": ")
totalMapped = 0
for ref in dictionary.values():
		print( ref.reference + "\t" + str(ref.totalReads) + "\t" + str(ref.totalLength))
total_read_bases_mapped = 0
prev_reads = set()
for ref, read_list in output_dict.items():
	with open(output + "_target_" + str(ref) + ".fasta", "w") as o:
		for entry in read_list:
			read_id, reference, identity = entry
			o.write(read_id + "\t" + reference + "\t" + str(identity) + "\n" + read_seqs[read_id] + "\n")
			if read_id in prev_reads:
				continue
			prev_reads.add(read_id)
			if ref != "unaligned":
				total_read_bases_mapped += read_lengths[read_id]
				totalMapped += 1
print("Total number of reads mapped: " + str(totalMapped) + "/" + str(len(target_channel_reads)))
print("Total read bases: " + str(target_channel_bases))
print("Total read bases mapped: " + str(total_read_bases_mapped))

# non target channels
dictionary, output_dict = AnalyseForChannels(non_target_channel_reads, samFile, multi_align, matching_prop)
print("\nReference stats for all other channels: ")
totalMapped = 0
prev_reads = set()
for ref in dictionary.values():
		print( ref.reference + "\t" + str(ref.totalReads) + "\t" + str(ref.totalLength))
total_read_bases_mapped = 0
for ref, read_list in output_dict.items():
	with open(output + "_nontarget_" + str(ref) + ".fasta", "w") as o:
		for entry in read_list:
			read_id, reference, identity = entry
			o.write(read_id + "\t" + reference + "\t" + str(identity) + "\n" + read_seqs[read_id] + "\n")
			if read_id in prev_reads:
				continue
			prev_reads.add(read_id)
			if ref != "unaligned":
				total_read_bases_mapped += read_lengths[read_id]
				totalMapped += 1
print("Total number of reads mapped: " + str(totalMapped) + "/" + str(len(non_target_channel_reads)))
print("Total read bases: " + str(non_target_channel_bases))
print("Total read bases mapped: " + str(total_read_bases_mapped))
