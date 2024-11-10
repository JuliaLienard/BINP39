#!/usr/bin/ python3

'''
Author Julia Lienard
Date : 2024-10-10

Description: The script compares two tables from two different samples containing identified methylated 
CDS (called by their pgap annotation ID) and outputs one tables with common methylated CDS between samples.

Usage: python FindUniqueMethylGenes.py Table1_methylated_sites.tsv Table2_methylated_sites.tsv

basename output will be Table1Table2_COMMONmethyl_cds.txt

'''

import sys
import os

Methyl_file1 = sys.argv[1]
Methyl_file2 = sys.argv[2]
File1Name = os.path.basename(Methyl_file1).split("_")
basename1 = File1Name[0]

File2Name = os.path.basename(Methyl_file2).split("_")
basename2 = File2Name[0]


# check correct number of arguments provided according to usage
if basename1 == basename2:
    print("Names of the provided tables are identical, please correct")
    quit()


with open(Methyl_file1, "r") as Methyl1, open(basename1+basename2+"_COMMONmethyl_cds.txt", "w") as output:
	next(Methyl1)
	output.write(f'PGAP_ID\tGene_strand\tGene_start_pos\tGene_end_pos\tIndex_methylated_site({basename1[0:4]})\tMethylation_type\tMethylation_pos\tMethyl_strand\tcoverage\tPercent_modified\tOriginSample\n')
	for line1 in Methyl1:
		methyl_call1 = line1.strip().split("\t")
		pgap_id_1 = methyl_call1[0]
		Methylation1_type = methyl_call1[5]
		methylated_site1_position = methyl_call1[6]
		Gene_strand1 = methyl_call1[1]
		start_pos = methyl_call1[2]
		end_pos = methyl_call1[3]
		index_methyl_site = methyl_call1[4]
		Methyl_strand = methyl_call1[7]
		coverage = methyl_call1[8]
		Percent_methylated = methyl_call1[9]

		with open(Methyl_file2, "r") as Methyl2:
			next(Methyl2)
			for line2 in Methyl2:
				methyl_call2 = line2.strip().split("\t")
				methylated_site2_position = methyl_call2[6]
				pgap_id_2 = methyl_call2[0]
				Methylation2_type = methyl_call2[5]
				Gene_strand2 = methyl_call2[1]
				if pgap_id_1 == pgap_id_2 and methylated_site1_position == methylated_site2_position and Gene_strand1 == Gene_strand2 and Methylation1_type == Methylation2_type:
					#output.write(f'{line1.strip()}\t{basename1+basename2}\n')
					output.write(f'{pgap_id_1}\t{Gene_strand1}\t{start_pos}\t{end_pos}\t{index_methyl_site}\t{Methylation1_type}\t{methylated_site1_position}\t{Methyl_strand}\t{coverage}\t{Percent_methylated}\t{basename1+basename2}\n')

















