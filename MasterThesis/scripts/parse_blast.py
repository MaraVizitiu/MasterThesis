#!/usr/bin/env python3
"""
Title: parse_blast.py

Author: Mara Vizitiu

Description: This script is used for parsing output files of BLAST analyses. The input file must be formatted in BLAST outfmt 6. The database used to run the analyses must be EukProt.
    A supplementary file is required that contains the taxonomic lineages of the species in EukProt. The taxonomic information for the hits will be added to the output file of this script.
    The information contained in the output file of this script is tab delimited and as follows:
        - name of query sequence
        - name of hit
        - percentage identity
        - length of alignment
        - e-value
        - alignmnet score
        - number of the hit for the current query (e.g. 1st hit, 2nd best hit etc.)
        - lineage of the species the matching sequence belongs to
        - a * is added if a particular hit belongs to a user-specified taxonomic clade and is among the top "i" hits (where the user can also specify the value of "i")

Usage: python3 parse_blast.py [-h] -i blast result file -o parsed output file -e EukProt taxonomic information file -c clade -n number of hits
"""
#%%Importing modules

import argparse

#%%Setting up argparse

parser = argparse.ArgumentParser(description = "This script will parse a BLAST output files formatted in outfmt 6 and add taxonomic information from the EukProt_included_data_sets file.")

parser.add_argument("-i", "--input", metavar = "blast result file", required = True, help = "Input file that is a BLAST result file formatted in outfmt 6")
parser.add_argument("-o", "--output", metavar = "parsed output file", required = True, help = "Output file that is a parsed version of the input file, with added taxonomic information")
parser.add_argument("-e", "--eukprot", metavar = "EukProt taxonomic information file", required = True, help = "File that includes the taxonomic information for the target sequences. File must be tab-separated, formatted as EukProt ID on column 1, semi-colon-separated lineage on column 2.")
parser.add_argument("-c", "--clade", metavar = "clade", required = True, help = "Name of the clade to mark the hits by")
parser.add_argument("-n", "--number", metavar = "number of hits", required = True, help = "Number of top hits to mark the hits by when the hit belongs to the clade of interest")

args = parser.parse_args()

#%%Input and output files

input_file = open(args.input, "r")     #blast result file
output_file = open(args.output, "w", encoding="utf-8") #output file
eukprot_info = open(args.eukprot, "r", encoding="utf-8") #file that includes taxonomic information for the species on EukProt

#%%Variables and asignments

i = 0                                   #counter variable that indicates the number of the hit to a particular query (i.e., no. 1 hit, no. 2 hit, etc.)
current_query = ""                      #variable that stores the current query that the hits have been matched to
headers_of_interest = [0, 1, 2, 3, 10, 11] #columns of the blast result file that we are interested in keeping in the final output file
eukprot_taxonomy = {}                   #dictionary to store the taxonomic information from the input file

#%%EukProt taxonomy file

#Parsing the taxonomy file and saving the lineage associated with each species into the designated dictionary

for line1 in eukprot_info:
    EP = line1.split()[0]               #EukProt accession number
    taxonomy = line1.split("\t")[1]     #lineage
    eukprot_taxonomy[EP] = taxonomy.rstrip() #save the accession numbers as keys and the lineages as corresponding values

#%%Blast result file

#Parsing the blast result file (in outfmt 6), saving the information of interest, associating the taxonomic information, and writing to the output file

for line2 in input_file:
    output_list = []                    #list of information to write on each row of the output file
    #The counter variable is set to 1 whenever a new query is found in the blast results file. The counter variable is increased by 1 for each hit of the same query.
    if line2.split()[0] == current_query: #found another hit for the same query
        i += 1                          #increase counter variable
    else:                               #a new query is found
        current_query = line2.split()[0] #assign query to the corresponding variable
        i = 1                           #restart the counter
    for column in headers_of_interest:  #find the information on each column of interest from the list
        output_list.append(line2.split()[column]) #add that information to the output_list for each row of the output file
    output_list.append(str(i))          #add the counter to the list
    target_id = line2.split()[1][0:7]   #grab the EukProt accession from the hit name to match it with its taxonomic information
    if target_id in eukprot_taxonomy.keys():
        output_list.append(eukprot_taxonomy[target_id]) #add taxonomic information to the list
    #If a "top x" hit to a particular query sequence belongs to a member of the clade decided by the user, mark that row in the output file with a *
        if args.clade in eukprot_taxonomy[target_id] and i <= int(args.number):
            output_list.append("*")
    else:
        output_list.append(" Eukaryota;Opisthokonta;Metazoa;Porifera;Demospongiae;Heteroscleromorpha;Spongillida;Spongillidae") #if a hit belongs to one of the manually freshwater sponge datasets, add this taxonomic information
    output_file.write("\t".join(output_list)+"\n") #write all items stored in the list to the output file

#%%

#Close all input and output files
input_file.close()
output_file.close()
eukprot_info.close()