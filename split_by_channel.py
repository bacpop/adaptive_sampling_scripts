from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import statistics

def get_options():
    description = "Splits fastq files by channel"
    parser = argparse.ArgumentParser(description=description,
                                     prog='python split_by_channel.py')

    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--infile',
                    help='Input fastq')
    IO.add_argument('--out',
                    help='Output prefix')
    IO.add_argument('--channels',
                    type=str,
                    help='Channels to choose, in form a-b, inclusive. Remainder put in separate file. ')
    return parser.parse_args()

def split_by_channel(input_file, out_prefix, start_channel, end_channel):
    target_channels = []
    nontarget_channels = []

    file_entries = SeqIO.parse(open(input_file), 'fastq')
    for entry in file_entries:
        id, sequence, description = entry.id, entry.seq, entry.description

        channel_index = description.find("ch=")
        channel = int(description[channel_index:].split("=", 1)[-1].split(" ")[0])

        if channel >= start_channel and channel <= end_channel:
            target_channels.append(entry)
        else:
            nontarget_channels.append(entry)

    SeqIO.write(target_channels, out_prefix + "_target.fastq", "fastq")
    SeqIO.write(nontarget_channels, out_prefix + "_nontarget.fastq", "fastq")

def main():
    options = get_options()
    channels = options.channels.split("-")
    start_channel = int(channels[0])
    end_channel = int(channels[-1])

    split_by_channel(options.infile, options.out, start_channel, end_channel)

    return 0

if __name__ == '__main__':
    main()