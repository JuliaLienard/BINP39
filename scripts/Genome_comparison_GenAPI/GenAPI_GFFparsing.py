# -*- coding: utf-8 -*-
"""
Author Julia Lienard
Date : 2024-10-26

Description: The script uses the table with identified missing genes ID of genome A compared to genome B, 
using GenAPi. From this table of geneID, the gff annotation file from genome B is parsed to retrieve the
corresponding gene product and output that, with the geneID in a txt file.

Usage: python GenAPI_GFFparsing.py MissingGeneGenomeA.txt GenomeB.gff

"""

import sys

Missing_Gen2 = sys.argv[1]
GFF_Gen1 = sys.argv[2]

missing_set = set()
with open(Missing_Gen2, "r") as Missing:
    for line1 in Missing:
        geneID = line1.strip()
        missing_set.add(geneID)
        
with open(GFF_Gen1, "r") as GFF, open("Missing_Genes_GFFannot.txt", "w") as output:
    output.write(f'#GeneID\tType_sequence\tstart_pos\tend_pos\tProduct\n')
    for line2 in GFF:
        if not line2.startswith("#"):
            annot = line2.strip().split("\t")
            GeneIDline = annot[8] # extract annotation info
            GeneID_full = GeneIDline.strip().split(";")
            GeneID_v1 = GeneID_full[0]
            GeneID_v2 = GeneID_v1.strip().split("=")
            GeneID_v3 = GeneID_v2[1]
            if GeneID_v3 in missing_set:
                next(GFF)
                annot = line2.strip().split("\t")
                type_seq = annot[2]
                start_pos = annot[3]
                end_pos = annot[4]
                Product_v1 = GeneIDline.strip().split("product=")
                Product_v2 = Product_v1[1]
                output.write(f'{GeneID_v3}\t{type_seq}\t{start_pos}\t{end_pos}\t{Product_v2}\n')
        if line2.startswith("##FASTA"):
            sys.exit("Encountered ##FASTA, stopping the script.")