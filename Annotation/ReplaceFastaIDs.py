from Bio import SeqIO
import argparse


def rename_fasta(fasta_file):
    """
    Function used to rename the entries in a fasta file with just numbers
    :param fasta_file: input fasta file to use
    :return: list of records with the new name
    """

    # Id counter
    count = 1
    sequence_records = []

    for record in SeqIO.parse(open(fasta_file, 'r'), "fasta"):
        new_name = str(count)
        record.id = new_name
        record.name = ''
        record.description = ''
        sequence_records.append(record)
        count += 1

    return sequence_records

program_description = "Quick script to rename the entry names on a list of fasta files. It will replace the" \
                      "current names with sequential numbers"\

parser = argparse.ArgumentParser(description=program_description)
parser.add_argument("-i", "--input_list", type=str, required=True,
                    help="List with the fasta files to process")

parser.add_argument("-o", "--output_prefix", type=str, required=True, help="prefix to use for the output files")
args = parser.parse_args()

processed_files = 1

for line in open(args.input_list, 'r'):

    line = line.rstrip()

    renamed_seqs = rename_fasta(line)
    output_file = args.output_prefix + "_" + line
    SeqIO.write(renamed_seqs, open(output_file, 'w'), "fasta")

    processed_files += 1

print "Total of processed files: " + str(processed_files)
