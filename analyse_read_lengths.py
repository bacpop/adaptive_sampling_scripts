import gzip
from pathlib import Path
import argparse
import os, glob
from collections import defaultdict

def get_options():
    description = "Parses reads based on lengths"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python analyse_read_lengths.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--indir',
                    required=True,
                    help='Path to read directory')
    IO.add_argument('--out',
                    default="result.txt",
                    help='Output file.')
    IO.add_argument('--old',
                    action="store_true",
                    default=False,
                    help='Use is running with old version (v5) of guppy.')
    return parser.parse_args()

def get_fq(directory):
    types = ([".fastq", ".gz"], [".fq", ".gz"])
    files = (
        str(p.resolve()) for p in Path(directory).glob("**/*") if p.suffixes in types
    )
    yield from files

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

def main():
    options = get_options()
    out = options.out
    indir = options.indir
    old = options.old

    target_reads_dict = defaultdict(list)
    unblocks_reads_dict = defaultdict(list)
    other_reads_dict = defaultdict(list)

    sum_list = glob.glob(os.path.join(indir, "sequencing_summary_*.txt"))
    summary = sum_list[0]
    unblock_dict = {}

    with open(summary, "r") as f:
        # ignore header
        next(f)
        for line in f:
            entry = line.strip().split("\t")
            if old:
                read_id = entry[3]
                unblock = entry[22]
                channel = entry[5]
            else:
                read_id = entry[4]
                unblock = entry[23]
                channel = entry[6]
            unblock_dict[read_id] = (unblock, channel)


    for f in get_fq(indir):
        if f.endswith(".gz"):
            fopen = gzip.open
        else:
            fopen = open

        base = os.path.splitext(os.path.basename(f))[0].split("_")
        if "barcode" in base[2]:
            file_id = "_".join([base[1], base[2]])
        else:
            file_id = "_".join([base[1], "NA"])


        with fopen(f, "rt") as fh:
            for name, seq, _ in readfq(fh):

                if name in unblock_dict:
                    unblock_tup = unblock_dict[name]
                    if unblock_tup[0] == "signal_positive":
                        target_reads_dict[file_id].append((name, len(seq), unblock_tup[1]))
                    elif unblock_tup[0] == "data_service_unblock_mux_change":
                        unblocks_reads_dict[file_id].append((name, len(seq), unblock_tup[1]))
                    else:
                        other_reads_dict[file_id].append((name, len(seq), unblock_tup[1]))
                else:
                    other_reads_dict[file_id].append((name, len(seq), "NA"))

    with open(out, "w") as o:
        o.write("Type\tFilter\tBarcode\tName\tLength\tChannel\n")
        for file_id, length_list in target_reads_dict.items():
            type = file_id.split("_")
            for len_entry in length_list:
                o.write("Target\t" + type[0] + "\t" + type[1] + "\t" + str(len_entry[0]) + "\t" + str(len_entry[1]) + "\t" +
                        str(len_entry[2]) + "\n")
        for file_id, length_list in unblocks_reads_dict.items():
            type = file_id.split("_")
            for len_entry in length_list:
                o.write("Non-target\t" + type[0] + "\t" + type[1] + "\t" + str(len_entry[0]) + "\t" + str(len_entry[1]) + "\t" +
                        str(len_entry[2]) + "\n")
        for file_id, length_list in other_reads_dict.items():
            type = file_id.split("_")
            for len_entry in length_list:
                o.write("Other\t" + type[0] + "\t" + type[1] + "\t" + str(len_entry[0]) + "\t" + str(len_entry[1]) + "\t" +
                        str(len_entry[2]) + "\n")


if __name__ == "__main__":
    main()