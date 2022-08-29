#!/usr/bin/env python3
"""
Title: fix_local_ids.py

Author: Mara Vizitiu

Description:
    This script takes a fasta file as input (interleaved or not). The input fasta file is intended for conversion into a BLAST database using the -parse_seqids option, but contains local IDs longer than the 50 characters accepted by makeblastdb.
    This script converts those headers containing improper local IDs (the string before the first space) into 10-character long alphanumeric strings and outputs a "fixed" fasta file with the converted headers, as well as a tab-separated file of the original headers and corresponding alphanumeric strings.

Usage: python3 fix_local_ids.py input_fasta_file output_fasta_file
"""

#%%Import modules
import random
import string
import sys

#%%Open input and output files

#Open the fasta file that will be turned into the BLAST database, the indexing file and the output file
input_file = open(sys.argv[1], "r")
index_file = open("renamed_headers.tsv", "w")
output_file = open(sys.argv[2], "w")

#%%Parse input fasta file

#List that will store all the new headers to make sure they don't repeat themselves
headers = []

#Parse the fasta file looking for the header lines
for line in input_file:
    if line.startswith(">"):
#Separate the local id from the header (the part of the header until the first space) and remove the > symbol
        local_id = line[1:].split(" ")[0]
#If the local id is longer than 50 characters
        if len(local_id) > 50:
#Generate a random alphanumeric string of length 10 to that header
            new_header = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
#Keep generating new alphanumeric strings if the newly generated one has already been used
            while new_header in headers:
                new_header = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
#Write the new local ID and the corresponding original header to the indexing file
            index_file.write(">{}\t{}\n".format(new_header,line))
#Write the new, randomly generated, header to the output fasta file
            output_file.write(">{}\n".format(new_header))
#If the local ID shorter than 50 characters, write the header as is to the output fasta file
        else:
            output_file.write(line)
#If the line currently being read is not a header, just paste it to the output fasta file
    else:
        output_file.write(line)

#%%Close all files
input_file.close()
index_file.close()
output_file.close()