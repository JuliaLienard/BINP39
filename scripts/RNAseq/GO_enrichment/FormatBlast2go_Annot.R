# This script is formatting the genome annotation file of M. avium strain LU439 SmT1, exporting it as a new table called
# "LU439_SmT_CDS_FULL_GOannot_Filt.txt".



# Setting working directory where the genome annotation file is:
setwd("~/Desktop/Courses/BINP39/RNAseq_visualization/LU439/1_data/1_genomeLU439T1_Blast2Goannotation/Pgap_B2GO_InterPro_EggNOG/Filt_by_GOonly_split/")

# Opening the annotation file, into a dataframe, that contains 4704 genes
LU439_SmT_CDS_FULL_GOannot <- read_delim("~/Desktop/Courses/BINP39/RNAseq_visualization/LU439/1_data/1_genomeLU439T1_Blast2Goannotation/Pgap_B2GO_InterPro_EggNOG/blast2go_table_LU439_SmT_CDS_FULL_GOannot_v3_May27.csv", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE) # 4704 genes

#Removing unwanted characters from columns ID (pgap annotation)
LU439_SmT_CDS_FULL_GOannot$ID <- gsub("cds-", "", as.character(LU439_SmT_CDS_FULL_GOannot$ID))

# Renaming column names
colnames(LU439_SmT_CDS_FULL_GOannot)[4] <- "pgap_annot"
colnames(LU439_SmT_CDS_FULL_GOannot)[13] <- "GO_ID"
colnames(LU439_SmT_CDS_FULL_GOannot)[14] <- "GO_name"

# Removing rows with no Gene Ontology terms
library(tidyr)
LU439_SmT_CDS_FULL_GOannot <- LU439_SmT_CDS_FULL_GOannot |> drop_na(GO_ID) # 4054 genes left

# Exporting as a new file
write.table(LU439_SmT_CDS_FULL_GOannot, 
            file = "~/Desktop/Courses/BINP39/RNAseq_visualization/LU439/1_data/1_genomeLU439T1_Blast2Goannotation/Pgap_B2GO_InterPro_EggNOG/Filt_by_GOonly_split/LU439_SmT_CDS_FULL_GOannot_Filt.txt", sep ="\t", row.names = FALSE, quote = FALSE)
