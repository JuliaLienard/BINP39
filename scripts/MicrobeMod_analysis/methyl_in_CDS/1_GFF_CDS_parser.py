#!/usr/bin/ python3

'''
Author Julia Lienard
Date : 2024-10-10

Description: The script uses the table with identified methylation sites to parse the GFF file
and identify whether the methylated site is within a coding sequence, based on its position and the
positions of genes. For methylated sites identified to be within CDS, are printed the type of CDS (ex: gene, pseudogene..), 
the annotation ID (pgap) of the CDS, the gene positions, strand, index of the methylated site

Usage: python GFFparser.py file.gff SampleA_methylated_sites.tsv basenameForOutput

basename will be used to print the origin of the methylation site index value in the output and for the name of the output file

'''

import sys
import os

GFF_file = sys.argv[1]
Methylation_table = sys.argv[2]
basename = sys.argv[3]


with open(basename+"_Methylation_sites_CDS.txt", "w") as MethylGenes, open(Methylation_table, "r") as MethylTable:
	MethylGenes.write(f'PGAP_ID\tGene_strand\tGene_start_pos\tGene_end_pos\tIndex_methylated_site_{basename}\tMethyl_type\tMethylation_pos\tMethyl_strand\tcoverage\tPercent_modified\n')
	for line in MethylTable:
		info_data = line.strip().split("\t")
		if info_data[0] != "index":
			methylsite = line.strip().split("\t")
			index_methyl_site = methylsite[1]
			methyl_genom_pos = int(methylsite[3]) + 1 # as the MethylTable is 0-based for genomic position and the gff file is 1-based, we correct here
			strand_site = methylsite[5]
			Methyl_type = methylsite[4]
			coverage = methylsite[11]
			Percent_methylated = methylsite[12]

			with open(GFF_file, "r") as GFF:
				next(GFF)
				next(GFF)
				next(GFF)
				next(GFF)
				next(GFF)
				next(GFF)
				for indexline, line in enumerate(GFF):
					infoline = line.strip().split("\t")
					type_sequence = infoline[2]
					if type_sequence == "CDS":
						start_pos = infoline[3]
						end_pos = infoline[4]
						strand_cds = infoline[6]
						if (int(methyl_genom_pos) >= int(start_pos) and int(methyl_genom_pos) <= int(end_pos)):
							pgap_id_complete = infoline[8].strip().split(";")
							pgap_id_long = pgap_id_complete[0].strip().split("-")
							pgap_id = pgap_id_long[1]
							MethylGenes.write(f'{pgap_id}\t{strand_cds}\t{start_pos}\t{end_pos}\t{index_methyl_site}\t{Methyl_type}\t{methyl_genom_pos}\t{strand_site}\t{coverage}\t{Percent_methylated}\n')
						










