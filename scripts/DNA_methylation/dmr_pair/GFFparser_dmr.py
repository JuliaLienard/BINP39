#!/usr/bin/ python3
'''
Author : Julia Lienard 
Date : 2024/11/08

The script parses the output bed file from the modkit dmr pair command (better if already filtered to keep lines for 
highest dmr scores), look whether the methylated sites differentially methylated is inside a coding sequence, based
on the provided gff3 annotation file (used on PGAP gff3 file) and output a txt file, indicating the
GeneID (PGAP), the methylation site (genomic position), methylation type, sample with higher methylation than the other,
dmr score and the gene product according to the gff file.

Usage : GFFparser_dmr.py GFF_file dmr.bed NameofOuputFile


'''

import sys


GFF_file = sys.argv[1]
dmr_table = sys.argv[2]
basename = sys.argv[3]



with open(basename+"_dmr_sites_CDS.txt", "w") as Methyl, open(dmr_table, "r") as dmr:
	Methyl.write(f'Methyl_genome_pos\tMethyl_type\tSample_dmr\tscore_dmr\tPGAP_ID\tstrand\tproduct\n')
	for line in dmr:
		info_line = line.strip().split("\t")
		methyl_genom_pos = int(info_line[1]) + 1 # as the dmr_table is 0-based for genomic position and the gff file is 1-based, we correct here
		sampleA = info_line[9]
		sampleB = info_line[10]
		dmr_score = info_line[4]
		if sampleA[0] == "a":
			type_final = "6mA"
			sampleA_stat = float(sampleA[2:])
			sampleB_stat = float(sampleB[2:])
		if sampleA[0] == "m":
			statA = sampleA.split(",") # m:0.00,21839:18.18
			m5_statA = statA[0] # m:0.00
			m4_statA = statA[1] # 21839:18.18
			m5_statA_block = m5_statA.split(":")
			m5_statA_final = float(m5_statA_block[1]) # 0.00
			m4_statA_block = m4_statA.split(":")
			m4_statA_final = float(m4_statA_block[1]) # 18.18

			statB = sampleB.split(",") # m:0.00,21839:18.18
			m5_statB = statB[0] # m:0.00
			m4_statB = statB[1] # 21839:18.18
			m5_statB_block = m5_statB.split(":")
			m5_statB_final = float(m5_statB_block[1]) # 0.00
			m4_statB_block = m4_statB.split(":")
			m4_statB_final = float(m4_statB_block[1]) # 18.18
			type_final = "" # methylation type
			sample_final = "" # which sample has the highest dmr score
			max_value = 0 # value of percent modified
			if m5_statA_final > max_value:
				type_final = "5mC"
				sample_final = "sampleA"
				max_value = m5_statA_final
			if m4_statA_final > max_value:
				type_final = "4mC"
				sample_final = "sampleA"
				max_value = m4_statA_final
			if m5_statB_final > max_value:
				type_final = "5mC"
				sample_final = "sampleB"
				max_value = m5_statB_final
			if m4_statB_final > max_value:
				type_final = "4mC"
				sample_final = "sampleB"
				max_value = m4_statB_final


		with open(GFF_file, "r") as GFF:
			next(GFF)
			next(GFF)
			next(GFF)
			next(GFF)
			next(GFF)
			next(GFF)
			gene_found = False
			for indexline, line in enumerate(GFF):
				infoline = line.strip().split("\t")
				type_sequence = infoline[2]
				if type_sequence == "CDS":
					start_pos = infoline[3]
					end_pos = infoline[4]
					info_CDS = infoline[8]
					strand = infoline[6]
					# Split the attributes by semicolon to find individual key-value pairs
					for element in info_CDS.split(";"):
           			 # Check if the attribute contains the 'product' key
						if element.startswith("product="):
							product_name = element.split("=")[1] # Extract the product name by removing 'product='
					if (int(methyl_genom_pos) >= int(start_pos) and int(methyl_genom_pos) <= int(end_pos)):
						gene_found = True
						pgap_id_complete = infoline[8].strip().split(";")
						pgap_id_long = pgap_id_complete[0].strip().split("-")
						pgap_id = pgap_id_long[1]
						if sampleA[0] == "a":
							if float(sampleA_stat) > float(sampleB_stat):
								sample_final = "sampleA"
								Methyl.write(f'{methyl_genom_pos}\t{type_final}\t{sample_final}\t{dmr_score}\t{pgap_id}\t{strand}\t{product_name}\n')
							elif float(sampleA_stat) < float(sampleB_stat):
								sample_final = "sampleB"
								Methyl.write(f'{methyl_genom_pos}\t{type_final}\t{sample_final}\t{dmr_score}\t{pgap_id}\t{strand}\t{product_name}\n')	
							else:
								Methyl.write(f'{methyl_genom_pos}\t{type_final}\tunknown\t{dmr_score}\t{pgap_id}\t{strand}\n')	# if the score is the same between the sample
						if sampleA[0] == "m":
							Methyl.write(f'{methyl_genom_pos}\t{type_final}\t{sample_final}\t{dmr_score}\t{pgap_id}\t{strand}\t{product_name}\n')
			if gene_found == False: # if no sites found within a CDS, we print a line with the rest of the info
				Methyl.write(f'{methyl_genom_pos}\t{type_final}\t{sample_final}\t{dmr_score}\tNotInCDS\tNA\tNA\n')













