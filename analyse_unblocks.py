"""
Adapted from summarise_fq.py (https://github.com/LooseLab/readfish/blob/master/ru/summarise_fq.py)
"""
import gzip
from pathlib import Path
from statistics import mean, median, stdev
from collections import defaultdict
import sys
import os, glob
import mappy as mp
from Bio import SeqIO
import numpy as np

import argparse

def get_options():
    description = "Parses reads based on unblocks"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python analyse_unblocks.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--indir',
                    required=True,
                    help='Path to read directory')
    IO.add_argument('--summary',
                    default=None,
                    help='Path to run summary file. Default is inferred from --indir.')
    IO.add_argument('--mux-period',
                    type=int,
                    default=480,
                    help='Period of mux scan (in seconds) to ignore reads prior. Default = 480 (8 mins)')
    IO.add_argument('--out',
                    default="result.txt",
                    help='Output file.')
    IO.add_argument('--ref',
                    default=None,
                    help='Specify minimap2 index. No alignment done if not specified. ')
    IO.add_argument('--target',
                    default=None,
                    help='Specify target within minimap2 index. Default = None ')
    IO.add_argument('--loci',
                    default=None,
                    help='Loci to find in minimap2 index ')
    return parser.parse_args()

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
    reference = options.ref
    out = options.out
    indir = options.indir
    summary = options.summary
    mux_period = options.mux_period
    loci = options.loci
    target = options.target

    if summary is None:
        sum_list = glob.glob(os.path.join(indir, "sequencing_summary_*.txt"))
        summary = sum_list[0]

    # create unblocks set
    unblock_dict = {}
    no_unblocks = 0

    # create mux-period set
    mux_set = set()
    with open(summary, "r") as f:
        # ignore header
        next(f)
        for line in f:
            entry = line.strip().split("\t")
            read_id = entry[4]
            unblock = entry[23]
            unblock_dict[read_id] = unblock

            start_time = float(entry[9])
            if start_time < mux_period:
                mux_set.add(read_id)

            # determine type of read
            if unblock == "data_service_unblock_mux_change":
                no_unblocks += 1

    print("Total unblocks: {}".format(str(no_unblocks)))

    if reference is not None:
        mapper = mp.Aligner(reference, preset="map-ont")

        print("Using reference: {}".format(reference), file=sys.stderr)

        if loci is not None:
            print("Aligning loci...")
            ref_dict = get_best_map(reference, loci)

    target_reads_dict = defaultdict(list)
    unblocks_reads_dict = defaultdict(list)
    other_reads_dict = defaultdict(list)

    for f in get_fq(indir):
        if f.endswith(".gz"):
            fopen = gzip.open
        else:
            fopen = open

        # get filename and extension
        base = os.path.splitext(os.path.basename(f))[0].split("_")
        #print(base)
        if "barcode" in base[2]:
            file_id = "_".join([base[1], base[2]])
        else:
            file_id = "_".join([base[1], "NA"])

        with fopen(f, "rt") as fh:
            for name, seq, _ in readfq(fh):
                ref_align = None
                ref_len = 0
                coord_start = "NA"
                coord_end = "NA"
                perc_id = 0
                matched_len = 0
                overlap = 0
                if reference is not None:
                    # Map seq, take alignment with largest mlen
                    for r in mapper.map(seq):
                        # if target defined, only look for those matches
                        if target is not None:
                            if r.ctg != target:
                                continue
                        if r.mlen < matched_len:
                            continue

                        # determine perc_id based on matched alignment length / alignment length in reference (maybe soft clipping)
                        coord_start = r.r_st
                        coord_end = r.r_en

                        len_ref_map = abs(coord_end - coord_start)

                        ref_align = r.ctg
                        perc_id = r.mlen / len_ref_map

                        ref_len = r.ctg_len

                        matched_len = r.mlen

                        # determine whether alignment to loci
                        if loci is not None:
                            ref_start = ref_dict[ref_align][1]
                            ref_end = ref_dict[ref_align][2]

                            overlap = len(range(max(ref_start, coord_start), min(ref_end, coord_end)))


                # check if in mux-period
                mux = 0
                if name in mux_set:
                    mux = 1

                if name in unblock_dict:
                    if unblock_dict[name] == "signal_positive":
                        target_reads_dict[file_id].append((name, len(seq), ref_align, perc_id, overlap, mux, coord_start, coord_end, ref_len))
                    elif unblock_dict[name] == "data_service_unblock_mux_change":
                        unblocks_reads_dict[file_id].append((name, len(seq), ref_align, perc_id, overlap, mux, coord_start, coord_end, ref_len))
                    else:
                        other_reads_dict[file_id].append((name, len(seq), ref_align, perc_id, overlap, mux, coord_start, coord_end, ref_len))
                    #print(name)

    with open(out, "w") as o:
        o.write("Type\tFilter\tBarcode\tRef\tLength\tName\tperc_id\tloci_overlap\tMux\tRef_start\tRef_end\tRef_len\n")
        for file_id, length_list in target_reads_dict.items():
            type = file_id.split("_")
            for len_entry in length_list:
                o.write("Target\t" + type[0] + "\t" + type[1] + "\t" + str(len_entry[2]) + "\t" + str(len_entry[1])
                        + "\t" + len_entry[0] + "\t" + str(len_entry[3]) + "\t" + str(len_entry[4]) + "\t"
                        + str(len_entry[5]) + "\t" + str(len_entry[6]) + "\t" + str(len_entry[7]) + "\t" + str(len_entry[8]) + "\n")
        for file_id, length_list in unblocks_reads_dict.items():
            type = file_id.split("_")
            for len_entry in length_list:
                o.write("Non-target\t" + type[0] + "\t" + type[1] + "\t" + str(len_entry[2]) + "\t" + str(len_entry[1])
                        + "\t" + len_entry[0] + "\t" + str(len_entry[3]) + "\t" + str(len_entry[4]) + "\t"
                        + str(len_entry[5]) + "\t" + str(len_entry[6]) + "\t" + str(len_entry[7]) + "\t" + str(len_entry[8]) + "\n")
        for file_id, length_list in other_reads_dict.items():
            type = file_id.split("_")
            for len_entry in length_list:
                o.write("Other\t" + type[0] + "\t" + type[1] + "\t" + str(len_entry[2]) + "\t" + str(len_entry[1])
                        + "\t" + len_entry[0] + "\t" + str(len_entry[3]) + "\t" + str(len_entry[4]) + "\t"
                        + str(len_entry[5]) + "\t" + str(len_entry[6]) + "\t" + str(len_entry[7]) + "\t" + str(len_entry[8]) + "\n")


if __name__ == "__main__":
    main()
