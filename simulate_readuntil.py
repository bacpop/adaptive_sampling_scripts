import argparse
from ru.basecall import Mapper as MappyMapper
from ru.graphalign.bifrost_mapper import Mapper as GraphMapper
from Bio import SeqIO
from timeit import default_timer as timer


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
    IO.add_argument('--aligner',
                    required=True,
                    choices=["mappy", "graph"],
                    help='Choice of aligner. Either mappy or graph.')
    IO.add_argument('--out',
                    default="parsed",
                    help='Output prefix')
    return parser.parse_args()


def time_alignment(mapper, infile, aligner, id, min_len):
    time_list = []
    rejections = 0

    file_entries = SeqIO.parse(open(infile), 'fasta')
    for entry in file_entries:
        sequence = str(entry.seq)

        if aligner == "graph":
            t0 = timer()
            if len(sequence) < min_len:
                rejections += 1
            else:
                result = mapper.map_read(sequence)
                if result < id:
                    rejections += 1
            t1 = timer()
            time_list.append(t1 - t0)

        else:
            t0 = timer()
            results = mapper.map_read(sequence)
            if not results:
                rejections += 1
            t1 = timer()
            time_list.append(t1 - t0)

    return time_list, rejections

def main():
    options = get_options()
    infile = options.infile
    out = options.out
    index = options.index
    aligner = options.aligner
    id = options.id
    min_len = options.min_len

    if aligner == "mappy":
        mapper = MappyMapper(index)
    else:
        mapper = GraphMapper(index)

    if mapper.initialised:
        print("{} initialised".format(aligner))

    time_list, rejections = time_alignment(mapper, infile, aligner, id, min_len)

    print("Output file: {}\nRejections: {}".format(out, str(rejections)))

    with open(out, "w") as f:
        f.writelines(time_list)

    return 0

if __name__ == '__main__':
    main()