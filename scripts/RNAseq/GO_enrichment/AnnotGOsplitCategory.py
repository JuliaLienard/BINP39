#!/usr/bin/ python3

'''
Author Julia Lienard
Date : 2024-05-27
Description: the program takes the annotation file "LU439_SmT_CDS_FULL_GOannot_Filt.txt" which is derived from the original genome annotation file
obtained by Giulia Ribeiro (original file called "blast2go_table_LU439_SmT_CDS_FULL_GOannot_v2_May27.csv") and processed first by the script FormatBlast2go.R, to
reformat the table and filter out all rows with no GO annotation.

The file "LU439_SmT_CDS_FULL_GOannot_Filt.txt" is here parsed to split for each gene, their associated GO term and ID, as they are all in one cell. O
ne GO term/ID per row, so several rows for the same gene if several associated GO term/ID.

Usage: python AnnotGOsplitCategory.py Annot_file.txt

'''


#import sys

#Annot_filt = sys.argv[1]


NumberLineOutput = 0
with open("LU439_SmT_CDS_FULL_GOannot_Filt.txt", "r") as file, open("LU439_SmT_CDS_FULL_GOannot_Filt_GOsplit.txt", "w") as Genes2GO:
	Genes2GO.write(f'pgap_ID\tproduct_PGAP\tGO_ID\tGO_name\tGO_category\n')
	next(file) # skip the header line next(file) # skip the header line
	for line in file:
		row = line.strip().replace("\"", "").split("\t")
		Gene_ID = row[3]
		GO_ID_row = row[12].split(";")
		GO_names_row = row[13].split(";")
		Description_row = row[4].split(";")
		Description = Description_row[0]
		if len(GO_ID_row) != len(GO_names_row):
			print("error") # used to check if there is a discrepency between the number of GO_ID and corresponding GO_names
		NumberLineOutput += int(len(GO_ID_row))
		if len(GO_ID_row) != 0:
			index_element = 0
			while index_element < len(GO_ID_row):
				GO_ID_indiv = GO_ID_row[index_element].strip()
				GO_name_indiv = GO_names_row[index_element].strip()
				if GO_ID_indiv[:2] == "P:":
					category = "Biological Process"
				elif GO_ID_indiv[:2] =="F:":
					category = "Molecular Function"
				elif GO_ID_indiv[:2] =="C:":
					category = "Cellular component"	
				Genes2GO.write(f'{Gene_ID}\t{Description}\t{GO_ID_indiv[2:]}\t{GO_name_indiv[2:]}\t{category}\n')
				index_element += 1

print(f'{NumberLineOutput} lines without header expected in the output table') # check the expected output
