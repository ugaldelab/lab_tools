#!/usr/bin/python
"""
This script annotate a protein fasta file using the EggNOG database. It requires the NOG.hmm (~5gb download)
file that can be downloaded from:
http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.hmm.tar.gz

And also requires HMMER3 to do the search against the database
"""


def parse_nog_annotations(input_nog):
    """
    This function reads the NOG annotation file (usually NOG.annotations.tsv),
    and generates a dictionary where the key is the OG ID, and the value is a
    list with the information:
    [Functional category, description]
    :param input_nog: Input file, usually NOG.annotations.tsv (obtained from the Eggnog website)
    :return: nog_annotation: A dictionary. Key is the NOG ID, and the values are a list with the one
    letter functional category and the description.
    """

    from collections import defaultdict

    nog_annotation = defaultdict(list)

    for line in open(input_nog, 'r'):
        line = line.rstrip()
        nog, og_id, members1, members2, func_cat, description = line.split("\t")

        nog_annotation[og_id] = [func_cat, description]

    return nog_annotation


def get_func_cat_description(input_category):
    """
    Get the functional description of a letter category.
    Data obtained from:
    http://eggnogdb.embl.de/download/latest/eggnog4.functional_categories.txt
    The input could be more than one letter
    :param input_category: One letter code for Eggnog
    :return:
    """

    func_cat_description = {
        "J": "Translation, ribosomal structure and biogenesis",
        "A": "RNA processing and modification",
        "K": "Transcription",
        "L": "Replication, recombination and repair",
        "B": "Chromatin structure and dynamics",
        "D": "Cell cycle control, cell division, chromosome partitioning",
        "Y": "Nuclear structure",
        "V": "Defense mechanisms",
        "T": "Signal transduction mechanisms",
        "M": "Cell wall/membrane/envelope biogenesis",
        "N": "Cell motility",
        "Z": "Cytoskeleton",
        "W": "Extracellular structures",
        "U": "Intracellular trafficking, secretion, and vesicular transport",
        "O": "Posttranslational modification, protein turnover, chaperones",
        "C": "Energy production and conversion",
        "G": "Carbohydrate transport and metabolism",
        "E": "Amino acid transport and metabolism",
        "F": "Nucleotide transport and metabolism",
        "H": "Coenzyme transport and metabolism",
        "I": "Lipid transport and metabolism",
        "P": "Inorganic ion transport and metabolism",
        "Q": "Secondary metabolites biosynthesis, transport and catabolism",
        "R": "General function prediction only",
        "S": "Function unknown"
        }

    description_list = []

    for letter in input_category:
        description_list.append(func_cat_description[letter])

    return description_list


def run_hmmsearch(input_fasta, db, num_cpus, tab_output, run_output):
    """

    :param input_fasta:
    :param db:
    :param num_cpus:
    :param tab_output:
    :param run_output:
    :return:
    """
    import subprocess
    import os
    import sys

    try:
        subprocess.call(["hmmsearch", "-E", "0.001", "--cpu", str(num_cpus), "--tblout", tab_output, "-o", run_output,
                        db, input_fasta])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            sys.exit("hmmsearch not found in the PATH")
        else:
            sys.exit("Some error with the input files")


def parse_hmm_results(input_file):
    """

    :param input_file:
    :return: Dictionary with the results from the hmm search
    """
    from collections import defaultdict

    raw_hmm_result = defaultdict(list)

    for line in open(input_file, 'r'):
        line = line.rstrip()
        if line.startswith("#"):
            continue

        query_name = line.split()[0]
        target_name = line.split()[2]
        evalue = float(line.split()[4])

        raw_hmm_result[query_name].append([target_name, evalue])

    # Get the lower evalue

    hmm_result = defaultdict(tuple)
    for query in raw_hmm_result:
        hmm_result[query] = sorted(raw_hmm_result[query], key=lambda x: x[1])[0]

    return hmm_result


def main():
    """
    Starting from a protein fasta file, annotate using hmmsearch against the eggnog database
    """
    import argparse
    from collections import defaultdict

    parser = argparse.ArgumentParser(prog="AnnotateGenomeEggNOG.py", description="Annotate a genome using EggNOG")
    parser.add_argument("-i", "--input_fasta", type=str, required=True, help="Input Fasta protein file")
    parser.add_argument("-a", "--eggnog_description", type=str, required=True, help="File with the Eggnog descriptions")
    parser.add_argument("-e", "--eggnog_db", type=str, required=True, help="Input hmm file to run hmmsearch")
    parser.add_argument("-p", "--cpus", type=int, required=True, help="Number of cpus to use")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix for the output files")

    args = parser.parse_args()

    # Read eggnog information
    eggnog = parse_nog_annotations(args.eggnog_description)

    # Run the hmmsearch
    tab_out = args.output_prefix + "_hmmsearch.txt"
    run_out = args.output_prefix + 'run_out.txt'
    run_hmmsearch(args.input_fasta, args.eggnog_db, args.cpus, tab_out, run_out)

    # Parse results
    hmm_results = parse_hmm_results(tab_out)

    # Create output files
    output_annotation = open(args.output_prefix + ".annotation.txt", 'w')
    output_summary = open(args.output_prefix + ".summary.txt", 'w')

    summary_func_categories = defaultdict(int)

    for entry in hmm_results:
        # get the product description
        target, evalue = hmm_results[entry]

        eggnog_id = target.split(".")[1]

        func_categories, product_description = eggnog[eggnog_id]

        # get the description for the one letter code
        category_descriptions = get_func_cat_description(func_categories)

        output_line = [entry, eggnog_id, str(evalue), product_description,
                       func_categories, ",".join(category_descriptions)]

        output_annotation.write("\t".join(output_line) + "\n")

        for letter in func_categories:
            summary_func_categories[letter] += 1

    for entry in summary_func_categories:
        entry_description = get_func_cat_description(entry)[0]
        count = str(summary_func_categories[entry])

        output_line = [entry, entry_description, count]
        output_summary.write("\t".join(output_line) + "\n")

    output_summary.close()
    output_annotation.close()


if __name__ == "__main__":
    main()
