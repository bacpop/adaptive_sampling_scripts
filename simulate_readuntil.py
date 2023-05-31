import argparse
from ru.basecall import Mapper as MappyMapper
from ru.graphalign.bifrost_mapper import Mapper as GraphMapper
from Bio import SeqIO
from timeit import default_timer as timer
import numpy as np


def get_options():
    description = "Simulates readuntil alignment"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python simulate_readuntil.py')

    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infile',
                    required=True,
                    help='Input fasta file')
    IO.add_argument('--index',
                    required=True,
                    help='Alignment index')
    IO.add_argument('--id',
                    type=float,
                    default=0.75,
                    help='Minimum identity for bifrost alignment. Default = 0.75')
    IO.add_argument('--min-len',
                    default=50,
                    type=int,
                    help='Minimum length for bifrost alignment. Default = 50')
    IO.add_argument('--avg-poi',
                    default=180,
                    type=int,
                    help='Average for random poisson sampling of reads for alignment. '
                         'Default = 180')
    IO.add_argument('--aligner',
                    required=True,
                    choices=["mappy", "graph"],
                    help='Choice of aligner. Either mappy or graph.')
    IO.add_argument('--out',
                    default="parsed",
                    help='Output prefix')
    return parser.parse_args()

# determine if iterable is empy
def peek(iterable):
    try:
        first = next(iterable)
    except StopIteration:
        return False
    return True

def time_alignment(mapper, infile, aligner, id, min_len, avg_poi):
    time_list = []
    rejections = 0

    file_entries = SeqIO.parse(open(infile), 'fasta')
    for entry in file_entries:
        raw_sequence = str(entry.seq)
        reject = 0

        # sample from poisson, add one to avoid zero values
        read_end = np.random.poisson(avg_poi) + 1
        if read_end > len(raw_sequence):
            sequence = raw_sequence
        else:
            sequence = raw_sequence[0:read_end]

        # time alignment only
        t0 = timer()
        result = mapper.map_read(sequence)
        t1 = timer()

        if aligner == "graph":
            if len(sequence) < min_len or result < id:
                reject = 1
                rejections += 1
        else:
            if not peek(result):
                reject = 1
                rejections += 1

        time_list.append(((t1 - t0), len(sequence), reject))

    return time_list, rejections

def main():
    options = get_options()
    infile = options.infile
    out = options.out
    index = options.index
    aligner = options.aligner
    id = options.id
    min_len = options.min_len
    avg_poi = options.avg_poi

    if aligner == "mappy":
        mapper = MappyMapper(index)
    else:
        mapper = GraphMapper(index)

    if mapper.initialised:
        print("{} initialised".format(aligner))

    time_list, rejections = time_alignment(mapper, infile, aligner, id, min_len, avg_poi)

    print("Output file: {}\nRejections: {}".format(out, str(rejections)))

    with open(out, "w") as f:
        f.write("Time\tSeq_len\tRejection\n")
        for entry in time_list:
            time, length, reject = entry
            f.write(str(time) + "\t" + str(length) + "\t" + str(reject) + "\n")

    return 0

if __name__ == '__main__':
    main()