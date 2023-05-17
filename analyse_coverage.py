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
import numpy as np

def get_options():
    description = "Aligns and determines coverage of reference"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python analyse_coverage.py')
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
                    default="coverage_output",
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
    summary_dict = {}

    mapper = mp.Aligner(reference, preset="map-ont")

    # keep list of barcodes to ensure printing in order
    barcode_list = []

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

        if barcode not in results_dict:
            barcode_list.append(barcode)
            results_dict[barcode] = {}
            summary_dict[barcode] = {}

        with gzip.open(f, "rt") as handle:
            input_sequences = SeqIO.parse(handle, 'fastq')
            for entry in input_sequences:
                fields, seq = entry.description, str(entry.seq)
                fields = fields.split()
                channel = int(fields[3].split("=")[1])
                assert (fields[3].split("=")[0] == "ch")
                assert (channel >= 0 and channel <= 512)


                # determine is read is in target channel
                target = False
                if channel >= min_channel and channel <= max_channel:
                    target = True

                ref_align = "unaligned"

                # perc id is alignment coverage of reference, just in case clipping present
                perc_id = 0
                match_len = 0

                align_count = 0

                ref_start = 0
                ref_end = 0
                ref_len = 0

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
                            ref_start = r.r_st
                            ref_end = r.r_en
                            ref_len = r.ctg_len
                    # otherwise make sure read aligns to known target
                    else:
                        if r.ctg in mm_target:
                            if match_len < r.mlen:
                                match_len = r.mlen
                                ref_align = r.ctg
                                perc_id = match_len / len_ref_map
                                ref_start = r.r_st
                                ref_end = r.r_en
                                ref_len = r.ctg_len

                # if below cutoff, or removing multialigning reads, set reference as unaligned
                if perc_id < matching_prop:
                    ref_align = "unaligned"
                    ref_start = 0
                    ref_end = 0
                    ref_len = 0
                if remove_multi and align_count > 1:
                    ref_align = "unaligned"
                    ref_start = 0
                    ref_end = 0
                    ref_len = 0

                length = len(seq)

                if ref_align not in summary_dict[barcode]:
                    if ref_align != "unaligned":
                        results_dict[barcode][ref_align] = {"adaptive" : np.zeros(ref_len, dtype=int),
                                                            "control" : np.zeros(ref_len, dtype=int)}
                    # 1st is total, second is reference length, 3rd is average coverage
                    summary_dict[barcode][ref_align] = {"adaptive" : [0, ref_len, 0],
                                                        "control" : [0, ref_len, 0]}

                channel_type = "adaptive" if target else "control"
                summary_dict[barcode][ref_align][channel_type][0] += length
                if ref_align != "unaligned":
                    results_dict[barcode][ref_align][channel_type][ref_start - 1 : ref_end - 1] += 1
                    summary_dict[barcode][ref_align][channel_type][2] = summary_dict[barcode][ref_align][channel_type][0] / summary_dict[barcode][ref_align][channel_type][1]


    for barcode, align_dict in results_dict.items():
        for align, channel_dict in align_dict.items():
            for channel, pos_hist in channel_dict.items():
                outfile = output + "_" + barcode + "_" + align + "_" + channel
                np.savetxt(outfile + "_hist.csv", pos_hist, delimiter=",", fmt='%s')

    with open(output + "_summary.txt", "w") as f:
        f.write("Barcode\tChannel\tAlignment\tRef_length\tTotal_bases\tAvg_cov\n")
        for barcode, align_dict in summary_dict.items():
            for align, channel_dict in align_dict.items():
                for channel, sum_list in channel_dict.items():
                    f.write(barcode + "\t" + channel + "\t" + align + "\t" + str(sum_list[1]) + "\t" + str(sum_list[0]) + "\t" + str(sum_list[2]) + "\n")


if __name__ == "__main__":
    main()