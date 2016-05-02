"""
This script downloads data from EBI using the ERR record ID.
The input is a list of records, while the output is the FASTQ or SFF associated
with such record
"""

from collections import defaultdict
import argparse
import os
import subprocess
from Bio import SeqIO

def read_xml(ena_url):
    """
    Recieves a plain URL for a nucleotide present in ENA. For example,
    http://www.ebi.ac.uk/ena/data/view/ERR562552
    It will query the ENA database, obtain the result in xml, and parse
    the sequence file, the type and the checksum value

    :param ena_url:
    :return:
    """
    import urllib2
    import xml.etree.ElementTree as ET

    input_file = urllib2.urlopen(ena_url + "&display=xml")
    xml_content = input_file.read()
    input_file.close()

    root = ET.fromstring(xml_content)

    xml_files = list()

    for entry in root.iter("FILE"):
        checksum = entry.attrib["checksum"]  # MD5 checksum to verify
        fastq_name = entry.attrib["filename"]  # Name of the file
        file_type = entry.attrib["filetype"]

        xml_files.append((file_type, checksum, fastq_name))

    return xml_files

# ---------------------- #


program_description = "Script to download data from EBI based on the ERR record. The input is a list with record IDS," \
                      "and the output is the associated sequence file for each record"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-i", "--input_list", type=str, help="Input list with the ID records to get from EBI",
                    required=True)
parser.add_argument("-o", "--output_folder", type=str, help="Output folder to put sequence files", required=True)

args = parser.parse_args()

# Create the output directory
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

# Read the input list
entry_count = 0

for line in open(args.input_list, 'r'):
    if not line.strip():
        continue

    ena_url = "http://www.ebi.ac.uk/ena/data/view/" + line

### Not finished yet...



