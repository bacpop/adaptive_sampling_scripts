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
    IO.add_argument('--graph-index',
                    required=True,
                    help='Graph .gfa index')
    IO.add_argument('--mappy-index',
                    required=True,
                    help='Mappy .mmi index')
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

def time_alignment(graph_index, mappy_index, infile, aligner, id, min_len, avg_poi):
    time_list_graph = []
    rejections_graph = 0
    time_list_mappy = []
    rejections_mappy = 0

    mappy_mapper = MappyMapper(mappy_index)
    if mappy_mapper.initialised:
        print("mappy initialised")

    graph_mapper = GraphMapper(graph_index)
    if graph_mapper.initialised:
        print("graph initialised")

    file_entries = SeqIO.parse(open(infile), 'fasta')
    for entry in file_entries:
        raw_sequence = str(entry.seq)
        mappy_reject = 0
        graph_reject = 0

        # sample from poisson, add one to avoid zero values
        read_end = np.random.poisson(avg_poi) + 1
        if read_end > len(raw_sequence):
            sequence = raw_sequence
        else:
            sequence = raw_sequence[0:read_end]

        # time alignment only for mappy
        t0 = timer()
        mappy_result = mappy_mapper.map_read(sequence)
        t1 = timer()

        if not peek(mappy_result):
            mappy_reject = 1
            rejections_mappy += 1

        time_list_mappy.append(((t1 - t0), len(sequence), mappy_reject))

        # repeat for graph
        t0 = timer()
        graph_result = graph_mapper.map_read(sequence)
        t1 = timer()

        if len(sequence) < min_len or graph_result < id:
            graph_reject = 1
            rejections_graph += 1

        time_list_graph.append(((t1 - t0), len(sequence), graph_reject))


    return time_list_graph, rejections_graph, time_list_mappy, rejections_mappy

def main():
    options = get_options()
    infile = options.infile
    out = options.out
    graph_index = options.graph_index
    mappy_index = options.mappy_index
    id = options.id
    min_len = options.min_len
    avg_poi = options.avg_poi

    time_list_graph, rejections_graph, time_list_mappy, rejections_mappy = time_alignment(graph_index, mappy_index, infile, id, min_len, avg_poi)

    print("Output file: {}\n Mappy Rejections: {}".format(out, str(rejections_mappy)))
    print("Output file: {}\n Graph Rejections: {}".format(out, str(rejections_graph)))

    with open(out, "w") as f:
        f.write("Tool\tTime\tSeq_len\tRejection\n")
        for entry in time_list_graph:
            time, length, reject = entry
            f.write("Graph" + str(time) + "\t" + str(length) + "\t" + str(reject) + "\n")
        for entry in time_list_mappy:
            time, length, reject = entry
            f.write("Mappy" + str(time) + "\t" + str(length) + "\t" + str(reject) + "\n")

    return 0

if __name__ == '__main__':
    main()