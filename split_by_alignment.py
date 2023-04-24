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
    IO.add_argument('--indir',
                    required=True,
                    help='Input directory containing fastq files.')
    IO.add_argument('--ref',
                    required=True,
                    help='Reference for minimap2 alignment (can be fasta or mmi).')
    IO.add_argument('--channels',
                    default="1-256",
                    help='channels to separate, in the form \"a-b\". '
                         '[Default=1-256]')
    IO.add_argument('--pid',
                    default=0.84,
                    type=float,
                    help='Minimum proportion of bases matching reference in alignment. '
                         'Default=0.84')
    IO.add_argument('--outdir',
                    default="results",
                    help='Output directory. '
                         'Default="results"')
    IO.add_argument('--remove',
                    default=False,
                    action="store_true",
                    help='Remove multi-mapping reads.'
                         'Default=False')
    IO.add_argument('--pass-only',
                    default=False,
                    action="store_true",
                    help='Align only pass reads.'
                         'Default=False')
    IO.add_argument('--target',
                    default=None,
                    help='Specify targets within minimap2 index with comma separated list. Default = None ')
    IO.add_argument('-v',
                    default=False,
                    action="store_true",
                    help='Verbose output.'
                         'Default=False')
    IO.add_argument('--manual',
                    default=False,
                    action="store_true",
                    help='Basecalling done manually.'
                         'Default=False')
    return parser.parse_args()

class ReferenceStats:
    def __init__(self, reference):
        self.totalReads = 0
        self.totalLength = 0
        self.reference = reference

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

    fastqDir = options.indir
    reference = options.ref
    channels = options.channels
    fields = channels.split("-")
    min_channel = int(fields[0])
    max_channel = int(fields[1])
    matching_prop = options.pid
    output = options.outdir
    remove_multi = options.remove
    only_pass = options.pass_only
    verbose = options.v
    mm_target = options.target
    manual = options.manual

    # split list of targets
    if mm_target != None:
        mm_target = set(mm_target.split(","))

    # initialise results dictionaries
    read_seqs = {"adaptive" : {},
                 "control" : {}}

    mapper = mp.Aligner(reference, preset="map-ont")

    # keep list of barcodes to ensure printing in order
    barcode_list = []

    for f in get_fq(fastqDir):
        if only_pass:
            if "pass/" not in f:
                continue

        # only run on gzipped files
        if ".gz" not in f:
            continue

        if verbose:
            print("Aligning: {}".format(str(f)))
        # get filename and extension
        base = os.path.splitext(os.path.basename(f))[0].split("_")
        # print(base)
        if not manual:
            if "barcode" in base[2]:
                barcode = base[2]
            else:
                barcode = "NA"

            if barcode not in read_seqs:
                read_seqs["adaptive"][barcode] = {}
                read_seqs["control"][barcode] = {}

                barcode_list.append(barcode)

        with gzip.open(f, "rt") as handle:
            input_sequences = SeqIO.parse(handle, 'fastq')
            for entry in input_sequences:
                fields, seq = entry.description, str(entry.seq)
                fields = fields.split()
                if manual:
                    channel = int(fields[4].split("=")[1])
                    assert (fields[4].split("=")[0] == "ch")
                    barcode = fields[-1].split("=")[1]

                    if barcode not in read_seqs:
                        read_seqs["adaptive"][barcode] = {}
                        read_seqs["control"][barcode] = {}

                        barcode_list.append(barcode)
                else:
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
                    # otherwise make sure read aligns to known target
                    else:
                        if r.ctg in mm_target:
                            if match_len < r.mlen:
                                match_len = r.mlen
                                ref_align = r.ctg
                                perc_id = match_len / len_ref_map

                # if below cutoff, or removing multialigning reads, set reference as unaligned
                if perc_id < matching_prop:
                    ref_align = "unaligned"
                if remove_multi and align_count > 1:
                    ref_align = "unaligned"

                if ref_align not in read_seqs["adaptive"][barcode]:
                    read_seqs["adaptive"][barcode][ref_align] = []
                    read_seqs["control"][barcode][ref_align] = []

                if target:
                    read_seqs["adaptive"][barcode][ref_align].append(entry)
                else:
                    read_seqs["control"][barcode][ref_align].append(entry)

    if not os.path.exists(output):
        os.mkdir(output)

    # write fastq files
    for channel, barcode_dict in read_seqs.items():
        for barcode, ref_align_dict in barcode_dict.items():
            for ref_align, seq_list in ref_align_dict.items():
                SeqIO.write(seq_list, os.path.join(output, "bc_" + str(barcode) + "_ch_" + str(channel) + "_al_" + str(ref_align) + ".fastq"), "fastq")

if __name__ == "__main__":
    main()