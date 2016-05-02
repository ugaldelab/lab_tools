import os
import argparse
import subprocess
import sys
from collections import defaultdict

program_description = "This script compares two fasta files and correlates their ids. Useful to associate" \
                      "two different genomes"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-1", "--fasta_file_1", type=str, required=True, help="Fasta file genome 1")
parser.add_argument("-2", "--fasta_file_2", type=str, required=True, help="Fasta file genome 2")
parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")

args = parser.parse_args()
blastout_format = "6 qseqid qlen sseqid slen evalue length pident"

# Blast 1 vs 2

output_file_1 = "1vs2.blastp"
genome_db_2 = "db_genome1"
blast1_dict = defaultdict()

# Make blastDB
try:
    subprocess.call(["makeblastdb", "-dbtype", "prot", "-in", args.fasta_file_2, "-out", genome_db_2])
except OSError as e:
    if e.errno == os.errno.ENOENT:
        sys.exit("Makeblastdb not found")
    else:
        sys.exit("Error with input files")

# Run Blastp
if not os.path.exists(output_file_1):
    try:
        subprocess.call(["blastp", "-db", genome_db_2, "-query", args.fasta_file_1, "-evalue", "1e-10",
                         "-outfmt", blastout_format, "-out", output_file_1, "-num_threads", "2"])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            sys.exit("Blastp not found")
        else:
            sys.exit("Error with input files")

# Parse hits into dictionary
for line in open(output_file_1, 'r'):
    line = line.rstrip()
    id1, id1_length, id2, id2_length, evalue, aln_length, identity = line.split("\t")

    aln_identity = float(identity)
    coverage = float(aln_length) / float(id1_length)

    # For each hit, store the identity and the alignment length
    if aln_identity >= 90 and coverage >= 0.7:
        if id1 in blast1_dict:
            dict_hit = blast1_dict[id1][0]
            dict_identity = blast1_dict[id1][1]
            dict_coverage = blast1_dict[id1][2]

            if aln_identity > dict_identity and coverage > dict_coverage:
                blast1_dict[id1] = [id2, aln_identity, coverage]
            else:
                continue
        else:
            blast1_dict[id1] = [id2, aln_identity, coverage]

# Blast 2 vs 1
output_file_2 = "2vs1.blastp"
genome_db_1 = "db_genome2"
blast2_dict = defaultdict()

# Make blastDB
try:
    subprocess.call(["makeblastdb", "-dbtype", "prot", "-in", args.fasta_file_1, "-out", genome_db_1])
except OSError as e:
    if e.errno == os.errno.ENOENT:
        sys.exit("Makeblastdb not found")
    else:
        sys.exit("Error with input files")

# Run Blastp
if not os.path.exists(output_file_2):
    try:
        subprocess.call(["blastp", "-db", genome_db_1, "-query", args.fasta_file_2, "-evalue", "1e-10",
                         "-outfmt", blastout_format, "-out", output_file_2, "-num_threads", "2"])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            sys.exit("Blastp not found")
        else:
            sys.exit("Error with input files")

# Parse hits into dictionary
for line in open(output_file_2, 'r'):
    line = line.rstrip()
    id1, id1_length, id2, id2_length, evalue, aln_length, identity = line.split("\t")

    aln_identity = float(identity)
    coverage = float(aln_length) / float(id1_length)

    # For each hit, store the identity and the alignment length
    if aln_identity >= 90 and coverage >= 0.7:
        if id1 in blast1_dict:
            dict_hit = blast2_dict[id1][0]
            dict_identity = blast2_dict[id1][1]
            dict_coverage = blast2_dict[id1][2]

            if aln_identity > dict_identity and coverage > dict_coverage:
                blast2_dict[id1] = [id2, aln_identity, coverage]
            else:
                continue
        else:
            blast2_dict[id1] = [id2, aln_identity, coverage]

# ------------------ #
# Compare the results of the annotation, based on genome 1
output_line = []
for gene_id in blast1_dict:
    hit_id = blast1_dict[gene_id][0]

    try:
        reciprocal_hit = blast2_dict[hit_id][0]
    except KeyError:
        output_line.append([gene_id, ""])
        continue

    if gene_id == reciprocal_hit:
        output_line.append([gene_id, hit_id])
    else:
        output_line.append([gene_id, ""])

output_table = open(args.output_prefix + ".txt", 'w')

for result in output_line:
    output_table.write("\t".join(result) + "\n")

output_table.close()