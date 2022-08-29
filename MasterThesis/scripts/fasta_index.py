#!/usr/bin/env python3

'''
Title: fasta_index.py

Author: Mara Vizitiu

Description: This script can accept one fasta file as input and it will rename the sequence headers contained inside the file with information from the name of the file, while also creating an index file where the new headers and the old ones are saved.

Usage: python3 fasta_index.py -f input_fasta_file -o                                                                               output_fasta_file
'''
import argparse
import os

#Add arguments
parser = argparse.ArgumentParser(description="This script can index and deindex fasta files")
parser.add_argument('-f','--fasta', metavar='fasta input', required = True, help='Input fasta file')
parser.add_argument('-o','--output', metavar='fasta output', help='Output fasta file', required = True)

args = parser.parse_args()

#The name of the file, that will be used for the indexed headers
basename = os.path.basename(args.fasta).split(".")[0]

input_file = open(args.fasta, "r") #Input fasta
output_file = open(args.output, "w") #Output fasta
index_file = open("Index_{}.tsv".format(basename), "w") #Output index file

#Variable for adding a numeric tag to each header
i = 1

for line in input_file:
    if line.startswith(">"): #When finding a header line
        old_header = line.rstrip()[1:] #Save the old header
        new_header = ">{}_{}".format(basename, i) #Create a new header based on the name of the input file and the numeric tag
        index_file.write("{}\t{}\n".format(old_header, new_header[1:])) #Write the old and the corresponding new header to the index file
        output_file.write(new_header+"\n") #Add the new header to the output file
        i = i+1 #Increase the numeric tag by 1
    else: #When finding a sequence line
        line = line.replace(".", "") #Remove dots from then sequences if necessary
        output_file.write(line) #Write the sequence line to the output file

input_file.close()
output_file.close()
index_file.close()