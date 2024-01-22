from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import re

def get_options():
    description = "Splits simulated fastq files from NanoSim-H by reference origin"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python split_simulated.py')

    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infile',
                    required=True,
                    help='Input fastq')
    IO.add_argument('--names',
                    required=True,
                    type=str,
                    help='Names of fasta headers to split, separate by comma. Must match order in "--pos"')
    IO.add_argument('--pos',
                    required=True,
                    type=str,
                    help='Positions in form a-b, separate by comma. Must match order in "--names"')
    IO.add_argument('--min-overlap',
                    type=int,
                    default=50,
                    help='Minimum overlap with target locus to be determined as true positive read.'
                         'Default=50')
    IO.add_argument('--out',
                    default="parsed",
                    help='Output prefix')
    return parser.parse_args()

def split_simulated(input_file, out_prefix, name_dict, min_overlap):
    target_file = []
    nontarget_file = []

    file_entries = SeqIO.parse(open(input_file), 'fasta')
    for entry in file_entries:
        target = False
        id, sequence, description = entry.id, entry.seq, entry.description

        split_description = description.split("_")
        ref_name = split_description[0]

        if ref_name in name_dict:
            loci_start, loci_end = name_dict[ref_name]

            read_length = len(sequence)
            read_forward = True if split_description[4] == "F" else False

            # correct position based on read orientation
            if read_forward:
                read_start = int(split_description[1])
                read_end = read_start + read_length
            else:
                read_end = int(split_description[1])
                read_start = read_end - read_length

            overlap = len(range(max(read_start, loci_start), min(read_end, loci_end) + 1))

            if overlap >= min_overlap:
                target = True

        if target:
            target_file.append(entry)
        else:
            nontarget_file.append(entry)

    SeqIO.write(target_file, out_prefix + "_target.fasta", "fasta")
    SeqIO.write(nontarget_file, out_prefix + "_nontarget.fasta", "fasta")

def main():
    options = get_options()
    names = options.names.split(",")
    pos = [tuple([int(b) for b in a.split("-")]) for a in options.pos.split(",")]
    min_overlap = options.min_overlap
    infile = options.infile
    out = options.out

    # replace spaces and underscores
    names = [re.sub(" ", "-", a) for a in names]
    names = [re.sub("_", "-", a) for a in names]

    name_dict = {}
    for i in range(len(names)):
        name_dict[names[i]] = pos[i]

    split_simulated(infile, out, name_dict, min_overlap)

    return 0

if __name__ == '__main__':
    main()