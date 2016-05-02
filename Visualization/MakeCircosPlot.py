"""
This script will take a GFF file, (optional) an eggnog annotation file, and will
make the basic files for a circos plot
"""


def read_gff(input):
    """
    Read GFF file and store information. The key for the dictionary is the locus tag
    :param input:
    :return:
    """
    from collections import defaultdict

    scaf_dictionary = defaultdict(list)
    scaffold_list = defaultdict(list)

    for line in open(input, 'r'):
        gene_dictionary = defaultdict(list)
        line = line.rstrip()

        if line.startswith("##sequence-region"):
            text, scaffold, scaf_start, scaf_end = line.split()
            scaffold_list[scaffold] = [scaf_start, scaf_end]
            continue

        elif line.startswith("#"):
            continue

        if line.startswith(">"):
            break

        scaf_id = line.split("\t")[0]
        start = line.split("\t")[3]
        end = line.split("\t")[4]
        type = line.split("\t")[2]
        strand = line.split("\t")[6]

        locus_tag = line.split("\t")[8].split(";")[0].split("=")[1]

        gene_dictionary[locus_tag] = [start, end, type, strand]
        #scaf_dictionary[scaf_id] = gene_dictionary
        scaf_dictionary[scaf_id].append(gene_dictionary)

    return scaf_dictionary, scaffold_list


def read_eggnog_annotation(input_file):
    """
    Read the eggnog annotation
    """
    from collections import defaultdict

    out_dict = defaultdict()

    for line in open(input_file, 'r'):
        line = line.rstrip()
        gene = line.split("\t")[0]
        func_letter = line.split("\t")[4]

        out_dict[gene] = func_letter

    return out_dict


def get_colors_category(input_letter):
    """
    Get the colors for each functional category
    Based on this:
    http://www.pnas.org/content/suppl/2009/02/20/0813403106.DCSupplemental/SD4_PDF.pdf
    Colors here are in RGB values
    """
    colors = {
    "A":"205,133,0",
    "B":"black",
    "C":"39,64,139",
    "D":"255,239,219",
    "E":"30,144,255",
    "F":"108,166,205",
    "G":"0,0,255",
    "H":"191,239,255",
    "I":"0,205,205",
    "J":"255,215,0",
    "K":"255,127,0",
    "L":"205,102,0",
    "M":"205,175,149",
    "N":"171,130,255",
    "O":"154,255,154",
    "P":"93,71,139",
    "Q":"69,139,116",
    "R":"229,229,229",
    "S":"179,179,179",
    "T":"255,99,71",
    "U":"255,20,147",
    "V":"255,192,203",
    "W":"0,255,0",
    "Y":"255,255,0",
    "Z":"255,0,0",
    "None":"127,127,127",
    "tRNA":"51,0,51",
    "rRNA":"51,102,51"
    }

    try:
        return colors[input_letter]
    except KeyError:
        return "127,127,127"



def main():
    """
    Generate the highlights and karyotypes files
    :return:
    """
    import argparse
    program_description = " "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-g", "--gff_file", type=str, help="Input GFF file", required=True)
    parser.add_argument("-e", "--eggnog", type=str, help="Eggnog file", required=True)
    parser.add_argument("-o", "--output_prefix", type=str, help="Output prefix for files", required=True)

    args = parser.parse_args()

    gene_info, scaffolds = read_gff(args.gff_file)

    #print the karyotypes
    karyotype_file = open(args.output_prefix + "_karyotype.txt", 'w')

    for scaf in scaffolds:
        start, stop = scaffolds[scaf]
        karyotype_file.write("chr - " + scaf + " " + scaf + " " + "0 " + str(stop) + " black" "\n")

    karyotype_file.close()

    #print the eggnog higlights

    #eggnog_colors = open(args.output_prefix + "_eggnog_colors.txt", 'w')

    eggnog_plus = open(args.output_prefix + "_eggnog_plus.txt", 'w')
    eggnog_negative = open(args.output_prefix + "_eggnog_negative.txt", 'w')

    eggnog_annotation = read_eggnog_annotation(args.eggnog)

    for scaffold in gene_info:
        for locus_information in gene_info[scaffold]:

            for locus in locus_information:

                start, end, type, strand = locus_information[locus]

                if type == "CDS":
                    if locus in eggnog_annotation:
                        func_category = eggnog_annotation[locus][0]
                    else:
                        func_category = None

                elif type == "tRNA":
                    func_category = "tRNA"

                elif type == "rRNA":
                    func_category = "rRNA"
                else:
                    func_category = None

                color = get_colors_category(func_category)

                output_line = [scaffold, str(start), str(end), "fill_color=" + color]

                if strand == "+":
                    eggnog_plus.write("\t".join(output_line) + "\n")

                if strand == "-":
                    eggnog_negative.write("\t".join(output_line) + "\n")

    eggnog_plus.close()
    eggnog_negative.close()


if __name__ == '__main__':
    main()