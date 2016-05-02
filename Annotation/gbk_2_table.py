from Bio import SeqIO
import argparse

program_description = "Parse GBK and creates table"
parser = argparse.ArgumentParser(description=program_description)
parser.add_argument("-g", "--gbk_file", type=str, required=True, help="GBK file")
parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix for the files")

args = parser.parse_args()

# Types CDS, gene, misc_binding, ncRNA, rRNA, repeat_region, source, tRNA, tmRNA

output_table = open(args.output_prefix + ".txt", 'w')
output_fasta = open(args.output_prefix + ".faa", 'w')

table_line = []

for index, record in enumerate(SeqIO.parse(args.gbk_file, "genbank")):
    for feature in record.features:
        note = ""

        if feature.type == "source" or feature.type == "misc_binding" or feature.type== "repeat_region":
            continue

        if feature.type == "CDS":
            contig = record.id
            gene_id = feature.qualifiers["locus_tag"][0]
            start = feature.location.start
            stop = feature.location.end
            strand = feature.location.strand

            try:
                product = feature.qualifiers["product"][0]
            except KeyError:
                product = ""

            try:
                protein_seq = feature.qualifiers["translation"][0]
                note = ""
            except KeyError:
                protein_seq = ""
                note = "pseudogene"

            table_line.append([contig, gene_id, feature.type, str(start), str(stop), str(strand),
                               product, note])

            # Create fasta protein files
            if protein_seq == "":
                continue
            else:
                output_fasta.write(">" + gene_id + "\n")
                output_fasta.write(protein_seq + "\n")

        if feature.type == "ncRNA":

            contig = record.id
            gene_id = feature.qualifiers["locus_tag"][0]
            start = feature.location.start
            stop = feature.location.end
            strand = feature.location.strand
            try:
                product = feature.qualifiers["product"][0]
            except KeyError:
                product = ""

            table_line.append([contig, gene_id, feature.type, str(start), str(stop), str(strand),
                               product, note])

        if feature.type == "rRNA":
            contig = record.id
            gene_id = feature.qualifiers["locus_tag"][0]
            start = feature.location.start
            stop = feature.location.end
            strand = feature.location.strand
            try:
                product = feature.qualifiers["product"][0]
            except KeyError:
                product = ""

            table_line.append([contig, gene_id, feature.type, str(start), str(stop), str(strand),
                               product, note])

        if feature.type == "tRNA":
            contig = record.id
            gene_id = feature.qualifiers["locus_tag"][0]
            start = feature.location.start
            stop = feature.location.end
            strand = feature.location.strand
            try:
                product = feature.qualifiers["product"][0]
            except KeyError:
                product = ""

            table_line.append([contig, gene_id, feature.type, str(start), str(stop), str(strand),
                               product, note])

        if feature.type == "tmRNA":
            contig = record.id
            gene_id = feature.qualifiers["locus_tag"][0]
            start = feature.location.start
            stop = feature.location.end
            strand = feature.location.strand
            try:
                product = feature.qualifiers["product"][0]
            except KeyError:
                product = ""

            table_line.append([contig, gene_id, feature.type, str(start), str(stop), str(strand),
                               product, note])


header=["contig","geneID", "feature Type", "start", "end", "strand", "product", "note"]
output_table.write("\t".join(header) + "\n")
for line in table_line:
    output_table.write("\t".join(line) + "\n")

output_table.close()
output_fasta.close()