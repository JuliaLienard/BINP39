#!/usr/bin/ python3

'''
Author Julia Lienard
Date : 2024-10-29

Description: The script uses the table with identified methylation sites (obtained from MicrobeMod pipeline) to parse 
the GFF file and identify whether the methylated site is in an intergenic region, by comparing its genomic position to
starting and ending positions of genes (i.e, position before an end position of gene 2 but before start position of gene 1,
if on strand + or if on strand -, methylation site after end position of gene 1 and before start position of gene 2, in reference
to positions in the gff file). Only the 200 nucleotides before a starting position of genes are considered.

For methylated sites identified in non coding sequences and having the same DNA strand as the CDS, 
the methylated site position in the genome, type of methylation, percent of methylation,  the annotation ID (pgap) of the downstream gene
and the distance between the methylated site and the downstream gene.
Methyl site upstream of gene 1 are omitted here. Manually, can check if some methyl site positions are within (size genome - 200bp)

Usage: python GFFparser_noCDS.py file.gff (with no FASTA sequence) methylated_sites.tsv basenameOutput

The gff file will be parsed only from line n°7 (after line with # and header), if different, change accordingly by adding
or removing "next(GFF)"

'''

import sys
import os
import pandas as pd

GFF_file = sys.argv[1]
Methylation_table = sys.argv[2]
basename = sys.argv[3]

List_methyl_strand_pos = []
List_methyl_strand_neg = []
methylation_strand_pos = {}
methylation_strand_neg = {}
with open(Methylation_table, "r") as MethylTable:
	next(MethylTable)
	for line1 in MethylTable:
		info_data = line1.strip().split("\t")
		index_methyl_site = info_data[1]
		methyl_genom_pos = int(info_data[3]) + 1 # as the MethylTable is 0-based for genomic position and the gff file is 1-based, we correct here
		strand_site = info_data[5]
		Methyl_type = info_data[4]
		Percent_methylated = info_data[12]
		if strand_site == "+":
			methylation_strand_pos = {"Methylation_pos" : methyl_genom_pos, "Methylation_type" : Methyl_type, "Percent_modified" : Percent_methylated, "strand" : strand_site}
			List_methyl_strand_pos.append(methylation_strand_pos)
		if strand_site == "-":
			methylation_strand_neg = {"Methylation_pos" : methyl_genom_pos, "Methylation_type" : Methyl_type, "Percent_modified" : Percent_methylated, "strand" : strand_site}
			List_methyl_strand_neg.append(methylation_strand_neg)


with open(GFF_file, "r") as GFF:
	next(GFF)
	next(GFF)
	next(GFF)
	next(GFF)
	next(GFF)
	next(GFF)
	Gene1_start_pos_strand_pos = [] # list of starting position of all genes on strand + for following gene
	Gene0_end_pos_strand_pos = [] # list of ending position of all genes on strand + for prior gene
	Gene_start_pos_strand_neg = [] # list of starting position (correspond to ending position on gff) of all genes on strand - 
	Gene_end_pos_strand_neg = [] # list of ending position (correspond to start position on gff) of all genes on strand -
	PGAP_ID_strand_pos = []
	PGAP_ID_strand_neg = []
	for line2 in GFF:
		infoline = line2.strip().split("\t")
		type_sequence = infoline[2]
		if type_sequence == "gene":
			start_pos = infoline[3]
			end_pos = infoline[4]
			strand_cds = infoline[6]
			PGAP_line = infoline[8]
			PGAP_IDlong1 = PGAP_line.strip().split(";")
			PGAP_IDlong2 = PGAP_IDlong1[0]
			PGAP_IDlong3 = PGAP_IDlong2.strip().split("=gene-")
			PGAP_ID = PGAP_IDlong3[1] # collected PGAP ID

			if start_pos == "1": # the first gene is always on strand +
				Gene0_end_pos_strand_pos.append(end_pos) # we omit in this script to look if a methyl site in upstream of gene 1. Gff ends with the end position of last gene anyway
			else:
				if strand_cds == "+":
					Gene1_start_pos_strand_pos.append(start_pos)
					Gene0_end_pos_strand_pos.append(end_pos)
					PGAP_ID_strand_pos.append(PGAP_ID)
				if strand_cds == "-":
					Gene_end_pos_strand_neg.append(start_pos) # the start position in gff file for genes on strand - correspond to their end position
					Gene_start_pos_strand_neg.append(end_pos) # the end position in gff file for genes on strand - correspond to their start position
					PGAP_ID_strand_neg.append(PGAP_ID)
print(len(PGAP_ID_strand_pos),len(PGAP_ID_strand_neg))

# Collecting lines for POSITIVE STRAND
with open("Methylation_sites_noCDS_StrandPos.txt", "w") as Methyl_location1: # this file is temporary and will be removed at the end of the script
	Methyl_location1.write(f'Methylation_pos\tMethylation_type\tPercent_modified\tDownstream_GeneID\tDistanceFromGene\tstrand\n')
	for entry in List_methyl_strand_pos:
		for i in range(0,(len(Gene1_start_pos_strand_pos)-1)): # loop over the list of gene end positions, from index 0 to length of the list Gene1_start_pos_strand_pos minus 1
			# here the start position of gene 1 of the gff file has not been collected, so we have the same index for start and end positions to be compared with, to the methylated site:
			if entry["Methylation_pos"] < int(Gene1_start_pos_strand_pos[i]) and entry["Methylation_pos"] > int(Gene0_end_pos_strand_pos[i]):
				distance = int(Gene1_start_pos_strand_pos[i]) - entry["Methylation_pos"]
				if distance <= 200:
					Methyl_location1.write(f'{entry["Methylation_pos"]}\t{entry["Methylation_type"]}\t{entry["Percent_modified"]}\t{PGAP_ID_strand_pos[i]}\t{distance}\t+\n')

# Collecting lines for NEGATIVE STRAND
with open("Methylation_sites_noCDS_StrandNeg.txt", "w") as Methyl_location2: # this file is temporary and will be removed at the end of the script
	Methyl_location2.write(f'Methylation_pos\tMethylation_type\tPercent_modified\tDownstream_GeneID\tDistanceFromGene\tstrand\n')
	for entry1 in List_methyl_strand_neg:
		for i in range(0,(len(Gene_end_pos_strand_neg)-1)): # loop over the list of gene end positions, from index 0 to length of the list Gene_end_pos_strand_neg minus 1
			# here the methylated site is compared to end position of the next gene (i+1) and to the start of the gene before :
			if entry1["Methylation_pos"] < int(Gene_end_pos_strand_neg[i+1]) and entry1["Methylation_pos"] > int(Gene_start_pos_strand_neg[i]):
				distance = entry1["Methylation_pos"] - int(Gene_start_pos_strand_neg[i])
				if distance <= 200:
					Methyl_location2.write(f'{entry1["Methylation_pos"]}\t{entry1["Methylation_type"]}\t{entry1["Percent_modified"]}\t{PGAP_ID_strand_neg[i]}\t{distance}\t-\n')
			# here for the gene, we check if methylated site is not upstream of the last gene on strand - (at max 200nucleotides), meaning placed after the end position on gff file of the last gene		
			if entry1["Methylation_pos"] > int(Gene_start_pos_strand_neg[len(Gene_start_pos_strand_neg)-1]) and entry1["Methylation_pos"] < (int((Gene_start_pos_strand_neg[len(Gene_start_pos_strand_neg)-1])) + 200):
				distance = (Gene_start_pos_strand_neg[len(Gene_start_pos_strand_neg)-1] + 200) - entry1["Methylation_pos"]
				Methyl_location2.write(f'{entry1["Methylation_pos"]}\t{entry1["Methylation_type"]}\t{entry1["Percent_modified"]}\t{PGAP_ID_strand_neg[i]}\t{distance}\t-\n')


# print(len(Gene0_end_pos_strand_pos)) # to check how many genes on strand (+) collected
# print(len(Gene_start_pos_strand_neg)) # to check how many genes on strand (-) collected
# print(int(Gene1_start_pos_strand_pos[len(Gene1_start_pos_strand_pos)-1]))
# print(int(Gene_start_pos_strand_neg[len(Gene_start_pos_strand_neg)-1])) # give the position of the last gene on strand -
# print((int((Gene_start_pos_strand_neg[len(Gene_start_pos_strand_neg)-1])) + 200)) # calculates the last genomic position considered to have a methylated site upstream of a starting position of a gene in strand -

# Concatenate the DataFrames
with open("Methylation_sites_noCDS_StrandPos.txt", "r") as table1_file, open("Methylation_sites_noCDS_StrandNeg.txt", "r") as table2_file:
	# Read the tables into pandas DataFrames (headers are included automatically)
	table1 = pd.read_csv(table1_file, sep='\t')
	table2 = pd.read_csv(table2_file, sep='\t')
	concatenated_table = pd.concat([table1,table2], ignore_index=True)

# Save the concatenated table to a new file, keeping the header from one table
	concatenated_table.to_csv(basename + "_Methylation_sites_noCDS.txt", sep='\t', index=False, header=True) # this table is the concatenation of the 2 temporary table files Methyl_location1 and Methyl_location2

print(f"Results saved to Methylation_sites_noCDS.txt")

os.remove("Methylation_sites_noCDS_StrandPos.txt") # this temporary file is removed 
os.remove("Methylation_sites_noCDS_StrandNeg.txt") # this temporary file is removed 



