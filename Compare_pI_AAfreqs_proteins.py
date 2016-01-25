#!/Users/juan/anaconda/bin/python
"""
Script that takes a set of different protein fasta files (one per genome,
or one per metagenome, for example), and calculates isoelectric point and
amino acid frequencies.

It will generate a table summarizing the results, as well as plots of the data.
The idea is based on Fernandez et al, 2014

Usage:

"""
from collections import defaultdict
import os
import argparse
import sys

try:
    import numpy as np
    from Bio import SeqIO
    from pyteomics import electrochem, auxiliary

except ImportError, e:
    print sys.stderr, " does not exist"
    print sys.stderr, " Exception: %s" % str(e)
    sys.exit(1)


program_description = "Script that takes a set of different protein fasta files (one per genome, or one per " \
                      "metagenome, for example), and calculates isoelectric point and amino acid frequencies"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-i", "--input_list", type=str, help="Input file with the name and location of the "
                                                         "protein fasta files", required=True)
parser.add_argument("-o", "--output_folder", type=str, help="Output folder for tables and plots", required=True)

args = parser.parse_args()

# Create output folder
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

# Read input list
fasta_to_process = list()
for line in open(args.input_list, 'r'):
    if line.strip():
        line = line.rstrip()
        name, path = line.split("\t")
        fasta_to_process.append((name, path))

# Iterate on each file and store the results
logfile = open(args.output_folder + "/logfile.txt", 'w')

freq_results = defaultdict()
pI_results = defaultdict()
total_residues = defaultdict(int)

for entry in fasta_to_process:
    aa_freqs = defaultdict(int)
    pI_count = list()
    name, fasta_file = entry
    print "Processing: " + name
    count = 0
    failed_count = 0
    residue_count = 0

    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        aa_sequence = seq_record.seq

        aa_sequence = str(aa_sequence.rstrip())

        count += 1

        if count % 1000 == 0:
            print count

        # Stop at 2,000 for testing purposes:
        if count == 200:
            break

        try:
            pI_count.append(electrochem.pI(aa_sequence))

            for aa in aa_sequence:
                aa_freqs[aa] += 1
                residue_count += 1

        except auxiliary.PyteomicsError:
            failed_count += 1
            continue

    logfile.write("For dataset " + name + "\n")
    logfile.write("A total of %d entries were found\n" % count)
    logfile.write("A total of %d entries had errors\n" % failed_count)
    logfile.write("A total of %d entries were used in the analysis\n\n" % (count-failed_count))

    freq_results[name] = aa_freqs
    pI_results[name] = pI_count
    total_residues[name] = residue_count

# ####Process the frequency results
aa_list = "AVKIGPFHWYCMKRDENQTS"

total_raw_results = []
total_freq_results = []

for aa in aa_list:
    freq_line = [aa]
    raw_line = [aa]
    for entry in freq_results:
        if freq_results[entry][aa]:
            aa_freq = freq_results[entry][aa] / float(total_residues[entry]) * 100
            freq_line.append(str.format('{0:.3f}', aa_freq))
            raw_line.append(freq_results[entry][aa])

        else:
            freq_line.append(str(0))
            raw_line.append(str(0))

    total_raw_results.append(raw_line)
    total_freq_results.append(freq_line)

# Print output tables
raw_output = open(args.output_folder + "/aa_raw_count.txt", 'w')
freq_output = open(args.output_folder + "/aa_freq.txt", 'w')

raw_output.write("Residue")
freq_output.write("Residue")

for entry in freq_results:
    raw_output.write("\t" + entry)
    freq_output.write("\t" + entry)

raw_output.write("\n")
freq_output.write("\n")

for line in total_raw_results:
    raw_output.write("\t".join(map(str, line)) + "\n")

for line in total_freq_results:
    freq_output.write("\t".join(map(str, line)) + "\n")

raw_output.close()
freq_output.close()

# To add, radar plots

# ###Process the pI results

total_pI_results = defaultdict(list)

for entry in pI_results:
    total_pI_results["header"].append(entry)

    hist_sum = np.histogram(np.array(pI_results[entry]), bins=24, range=(2, 14), density=True)

    for hist_bin, value in zip(hist_sum[1], hist_sum[0]):
        total_pI_results[hist_bin].append(value)

pI_output = open(args.output_folder + "/pI_results.txt", 'w')

pI_output.write("pI\t" + "\t".join(total_pI_results["header"]) + "\n")

for data_point in sorted(total_pI_results):

    if data_point == "header":
        continue

    else:
        pI_output.write(str(data_point) + "\t")
        pI_output.write("\t".join(map(str, total_pI_results[data_point])) + "\n")


pI_output.close()

logfile.close()
