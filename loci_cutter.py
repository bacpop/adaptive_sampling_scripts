import mappy as mp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import argparse
from collections import defaultdict

def get_options():
    description = "Cuts out loci based on alignment of reference sequences"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python loci_cutter.py')
    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infiles',
                    required=True,
                    help='List of file paths to cut, one per line')
    IO.add_argument('--query',
                    help='Fasta file of sequences to align and cut')
    IO.add_argument('--separate',
                    action="store_true",
                    default=False,
                    help='Generate separate fasta files of cut and remaining sequences. ')
    IO.add_argument('--cutoff',
                    type=float,
                    default=0.7,
                    help='Cutoff of alignment length to confirm match. '
                         'Default = 0.7')
    IO.add_argument('--count',
                    action="store_true",
                    default=False,
                    help='Go through cut input files, count different aligning sequences. ')
    IO.add_argument('--outpref',
                    default="result",
                    help='Output prefix ')
    return parser.parse_args()

def get_best_map(index, fasta, cutoff):
    a = mp.Aligner(index, preset="asm10")

    best_hit = (None, None, 0, 0, 0)

    fasta_sequences = SeqIO.parse(open(fasta), 'fasta')
    for fasta in fasta_sequences:
        id, sequence = fasta.id, str(fasta.seq)

        for hit in a.map(sequence):
            if not hit.is_primary:
                continue
            query_hit = hit.blen

            # set cutoff for minimum alignment length
            if query_hit < cutoff * len(sequence):
                continue

            if query_hit > best_hit[-1]:
                best_hit = (id, hit.ctg, hit.r_st, hit.r_en, query_hit)

    return best_hit

def cut_loci(infiles, in_fasta, separate, cutoff, outpref):
    file_list = []

    count_dict = {}

    with open(infiles, "r") as f:
        for line in f:
            file_list.append(line.strip())

    for file in file_list:
        best_map = get_best_map(file, in_fasta, cutoff)

        if best_map[0] not in count_dict:
            count_dict[best_map[0]] = []
        count_dict[best_map[0]].append(best_map[0])

        if best_map[0] != None:
            out_pref = os.path.splitext(file)[0]

            rem_records = []
            cut_records = []

            fasta_sequences = SeqIO.parse(open(file), 'fasta')
            for fasta in fasta_sequences:
                id, sequence = fasta.id, str(fasta.seq)

                # if match to contig, cut out loci
                if id == best_map[1]:
                    cut = sequence[best_map[2]:best_map[3]]
                    pre_remainder = sequence[:best_map[2]]
                    post_remainder = sequence[best_map[3]:]

                    rem_records.append(SeqRecord(Seq(pre_remainder), id=id + "_precut",
                                             description=fasta.description))
                    # print to separate files
                    if separate:
                        cut_records.append(SeqRecord(Seq(cut), id=id + "_" + best_map[0],
                                                 description=fasta.description))
                    else:
                        rem_records.append(SeqRecord(Seq(cut), id=id + "_cut_" + best_map[0],
                                                 description=fasta.description))

                    rem_records.append(SeqRecord(Seq(post_remainder), id=id + "_postcut",
                                             description=fasta.description))
                else:
                    rem_records.append(SeqRecord(fasta.seq, id=fasta.id,
                                             description=fasta.description))

            SeqIO.write(rem_records, out_pref + "_rem.fa", "fasta")
            if separate:
                SeqIO.write(cut_records, out_pref + "_cut.fa", "fasta")
        # if no match, copy and save as rem
        else:
            rem_records = []

            fasta_sequences = SeqIO.parse(open(file), 'fasta')
            for fasta in fasta_sequences:
                rem_records.append(SeqRecord(fasta.seq, id=fasta.id,
                                             description=fasta.description))
            SeqIO.write(rem_records, out_pref + "_nocut.fa", "fasta")


    # print output files
    with open(outpref + "_summary.txt", "w") as o1, open(outpref + ".tsv", "w") as o2:
        o1.write("Assignment\tCount")
        o2.write("File\tAssignment")
        for assignment, file_list in count_dict.items():
            o1.write(str(assignment) + "\t" + str(len(file_list)))
            for file in file_list:
                o2.write(str(file) + "\t" + str(assignment))

def count_loci(infiles, outpref):
    file_list = []
    count_dict = {}

    with open(infiles, "r") as f:
        for line in f:
            # ensure only cut files are analysed
            if "_cut" in line:
                file_list.append(line.strip())

    for file in file_list:
        base = os.path.splitext(os.path.basename(file))[0].split("_cut")[0]
        fasta_sequences = SeqIO.parse(open(file), 'fasta')
        for fasta in fasta_sequences:
            id = fasta.id
            assignment = id.split("_")[-1]
            if assignment not in count_dict:
                count_dict[assignment] = []
            count_dict[assignment].append(base)

    # print output files
    with open(outpref + "_summary.txt", "w") as o1, open(outpref + ".tsv", "w") as o2:
        o1.write("Assignment\tCount")
        o2.write("File\tAssignment")
        for assignment, file_list in count_dict.items():
            o1.write(str(assignment) + "\t" + str(len(file_list)))
            for file in file_list:
                o2.write(str(file) + "\t" + str(assignment))

def main():
    options = get_options()
    infiles = options.infiles
    query = options.query
    separate = options.separate
    cutoff = options.cutoff
    count = options.count
    outpref = options.outpref

    if not count:
        cut_loci(infiles, query, separate, cutoff, outpref)
    else:
        count_loci(infiles, outpref)

    return 0

if __name__ == "__main__":
    main()