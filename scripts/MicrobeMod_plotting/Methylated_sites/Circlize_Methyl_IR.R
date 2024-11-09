

library(circlize)
library(readr)
library(dplyr)



# Loading data frame with methylated sites within 200bp upstream codon start position

O242_Methyl_IR <- read_delim("O242_Methylation_sites_noCDS.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
O251_Methyl_IR <- read_delim("O251_Methylation_sites_noCDS.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
T257_Methyl_IR <- read_delim("T257_Methylation_sites_noCDS.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
T278_Methyl_IR <- read_delim("T278_Methylation_sites_noCDS.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

# filtering data for 6mA sites
O242_Methyl_IR_6mA <- O242_Methyl_IR %>% filter(Methylation_type == "6mA")
O251_Methyl_IR_6mA <- O251_Methyl_IR %>% filter(Methylation_type == "6mA")
T257_Methyl_IR_6mA <- T257_Methyl_IR %>% filter(Methylation_type == "6mA")
T278_Methyl_IR_6mA <- T278_Methyl_IR %>% filter(Methylation_type == "6mA")
## for promoter region
O242_Methyl_IR_6mA_promot <- O242_Methyl_IR_6mA%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)
O251_Methyl_IR_6mA_promot <- O251_Methyl_IR_6mA%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)
T257_Methyl_IR_6mA_promot <- T257_Methyl_IR_6mA%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)
T278_Methyl_IR_6mA_promot <- T278_Methyl_IR_6mA%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)

# Plotting for all 6mA Methyl IR in the 4 samples

png("circlize_MAC101_methyl_IR_6mA.png", width = 850, height = 800, res = 150)  # Set resolution and dimensions

# Adjust circular plot parameters
circos.par("track.height" = 0.1)
circos.initialize(factors = O242_Methyl_IR_6mA$Methylation_type, x = O242_Methyl_IR_6mA$Methylation_pos)

# track 1 - 6mA IR for MAC101 O242 (higher number of positions among df)
circos.track(O242_Methyl_IR_6mA$Methylation_type, 
             y = O242_Methyl_IR_6mA$Percent_modified,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(6), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.3) # size label genomic positions
             })

# track 2 - 6mA IR for MAC101 O251
circos.track(O251_Methyl_IR_6mA$Methylation_type, 
             y = O251_Methyl_IR_6mA$Percent_modified)

# track 3 - 6mA IR for MAC101 T257
circos.track(T257_Methyl_IR_6mA$Methylation_type, 
             y = T257_Methyl_IR_6mA$Percent_modified)

# track 4 - 6mA IR for MAC101 T278
circos.track(T278_Methyl_IR_6mA$Methylation_type, 
             y = T278_Methyl_IR_6mA$Percent_modified)


circos.lines(O242_Methyl_IR_6mA$Methylation_pos, O242_Methyl_IR_6mA$Percent_modified, type = "h", straight = FALSE,
             col="magenta4", sector.index = get.current.sector.index(),
             track.index = 1)

circos.lines(O251_Methyl_IR_6mA$Methylation_pos, O251_Methyl_IR_6mA$Percent_modified, type = "h", straight = FALSE,
             col="purple", sector.index = get.current.sector.index(),
             track.index = 2)
circos.lines(T257_Methyl_IR_6mA$Methylation_pos, T257_Methyl_IR_6mA$Percent_modified, type = "h", straight = FALSE,
             col="palegreen4", sector.index = get.current.sector.index(),
             track.index = 3)
circos.lines(T278_Methyl_IR_6mA$Methylation_pos, T278_Methyl_IR_6mA$Percent_modified, type = "h", straight = FALSE,
             col="springgreen3", sector.index = get.current.sector.index(),
             track.index = 4)
# Add a title
title("MAC101 \n >= 200bp upstream \n start codon", line = -14, cex.main = 1.5)  # Adjust line for vertical position and size as needed

# Add legend outside the circular plot
legend("topleft", legend = c("O242", "O251", "T257", "T278"),
       col = c("magenta4", "purple","palegreen4","springgreen3"), lty = 1, lwd = 2, cex = 0.8, box.lty = 0)  # Customize position and style

# Clear and close
circos.clear()
dev.off()


# Plotting for all 6mA Methyl promoter region in the 4 samples

png("circlize_MAC101_methyl_promoter_6mA.png", width = 850, height = 800, res = 150)  # Set resolution and dimensions

# Adjust circular plot parameters
circos.par("track.height" = 0.1)
circos.initialize(factors = O242_Methyl_IR_6mA_promot$Methylation_type, x = O242_Methyl_IR_6mA_promot$Methylation_pos)

# track 1 - 6mA IR for MAC101 O242 (higher number of positions among df)
circos.track(O242_Methyl_IR_6mA_promot$Methylation_type, 
             y = O242_Methyl_IR_6mA_promot$Percent_modified,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(6), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.3) # size label genomic positions
             })

# track 2 - 6mA IR for MAC101 O251
circos.track(O251_Methyl_IR_6mA_promot$Methylation_type, 
             y = O251_Methyl_IR_6mA_promot$Percent_modified)

# track 3 - 6mA IR for MAC101 T257
circos.track(T257_Methyl_IR_6mA_promot$Methylation_type, 
             y = T257_Methyl_IR_6mA_promot$Percent_modified)

# track 4 - 6mA IR for MAC101 T278
circos.track(T278_Methyl_IR_6mA_promot$Methylation_type, 
             y = T278_Methyl_IR_6mA_promot$Percent_modified)


circos.lines(O242_Methyl_IR_6mA_promot$Methylation_pos, O242_Methyl_IR_6mA_promot$Percent_modified, type = "h", straight = FALSE,
             col="turquoise4", sector.index = get.current.sector.index(),
             track.index = 1)

circos.lines(O251_Methyl_IR_6mA_promot$Methylation_pos, O251_Methyl_IR_6mA_promot$Percent_modified, type = "h", straight = FALSE,
             col="turquoise3", sector.index = get.current.sector.index(),
             track.index = 2)
circos.lines(T257_Methyl_IR_6mA_promot$Methylation_pos, T257_Methyl_IR_6mA_promot$Percent_modified, type = "h", straight = FALSE,
             col="orangered", sector.index = get.current.sector.index(),
             track.index = 3)
circos.lines(T278_Methyl_IR_6mA_promot$Methylation_pos, T278_Methyl_IR_6mA_promot$Percent_modified, type = "h", straight = FALSE,
             col="orange", sector.index = get.current.sector.index(),
             track.index = 4)
# Add a title
title("MAC101 - potential \n promoter regions", line = -14, cex.main = 1.5)  # Adjust line for vertical position and size as needed

# Add legend outside the circular plot
legend("topleft", legend = c("O242", "O251", "T257", "T278"),
       col = c("turquoise4", "turquoise3","orangered","orange"), lty = 1, lwd = 2, cex = 0.8, box.lty = 0)  # Customize position and style

# Clear and close
circos.clear()
dev.off()


# filtering data for 4mC sites
O242_Methyl_IR_4mC <- O242_Methyl_IR %>% filter(Methylation_type == "4mC")
O251_Methyl_IR_4mC <- O251_Methyl_IR %>% filter(Methylation_type == "4mC")
T257_Methyl_IR_4mC <- T257_Methyl_IR %>% filter(Methylation_type == "4mC")
T278_Methyl_IR_4mC <- T278_Methyl_IR %>% filter(Methylation_type == "4mC")
## for promoter region
O242_Methyl_IR_4mC_promot <- O242_Methyl_IR_4mC%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)
O251_Methyl_IR_4mC_promot <- O251_Methyl_IR_4mC%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)
T257_Methyl_IR_4mC_promot <- T257_Methyl_IR_4mC%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)
T278_Methyl_IR_4mC_promot <- T278_Methyl_IR_4mC%>% filter(DistanceFromGene >=10 & DistanceFromGene <= 35)


# Plotting for all 4mC Methyl IR in the 4 samples

png("circlize_MAC101_methyl_IR_4mC.png", width = 850, height = 800, res = 150)  # Set resolution and dimensions

# Adjust circular plot parameters
circos.par("track.height" = 0.1)
circos.initialize(factors = O242_Methyl_IR_4mC$Methylation_type, x = O242_Methyl_IR_4mC$Methylation_pos)

# track 1 - 4mC IR for MAC101 O242 (higher number of positions among df)
circos.track(O242_Methyl_IR_4mC$Methylation_type, 
             y = O242_Methyl_IR_4mC$Percent_modified,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(6), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.3) # size label genomic positions
             })

# track 2 - 4mC IR for MAC101 O251
circos.track(O251_Methyl_IR_4mC$Methylation_type, 
             y = O251_Methyl_IR_4mC$Percent_modified)

# track 3 - 4mC IR for MAC101 T257
circos.track(T257_Methyl_IR_4mC$Methylation_type, 
             y = T257_Methyl_IR_4mC$Percent_modified)

# track 4 - 4mC IR for MAC101 T278
circos.track(T278_Methyl_IR_4mC$Methylation_type, 
             y = T278_Methyl_IR_4mC$Percent_modified)


circos.lines(O242_Methyl_IR_4mC$Methylation_pos, O242_Methyl_IR_4mC$Percent_modified, type = "h", straight = FALSE,
             col="blue4", sector.index = get.current.sector.index(),
             track.index = 1)

circos.lines(O251_Methyl_IR_4mC$Methylation_pos, O251_Methyl_IR_4mC$Percent_modified, type = "h", straight = FALSE,
             col="blue", sector.index = get.current.sector.index(),
             track.index = 2)
circos.lines(T257_Methyl_IR_4mC$Methylation_pos, T257_Methyl_IR_4mC$Percent_modified, type = "h", straight = FALSE,
             col="black", sector.index = get.current.sector.index(),
             track.index = 3)
circos.lines(T278_Methyl_IR_4mC$Methylation_pos, T278_Methyl_IR_4mC$Percent_modified, type = "h", straight = FALSE,
             col="grey55", sector.index = get.current.sector.index(),
             track.index = 4)
# Add a title
title("MAC101 \n >= 200bp upstream \n start codon", line = -14, cex.main = 1.5)  # Adjust line for vertical position and size as needed

# Add legend outside the circular plot
legend("topleft", legend = c("O242", "O251", "T257", "T278"),
       col = c("blue4", "blue","black","grey55"), lty = 1, lwd = 2, cex = 0.8, box.lty = 0)  # Customize position and style

# Clear and close
circos.clear()
dev.off()

# Plotting for all 4mC Methyl promoter region in the 4 samples

png("circlize_MAC101_methyl_promoter_4mC.png", width = 850, height = 800, res = 150)  # Set resolution and dimensions

# Adjust circular plot parameters
circos.par("track.height" = 0.1, cex = 1)
circos.initialize(factors = O242_Methyl_IR_4mC_promot$Methylation_type, x = O242_Methyl_IR_4mC_promot$Methylation_pos)

# track 1 - 4mC IR for MAC101 O242 (higher number of positions among df)
circos.track(O242_Methyl_IR_4mC_promot$Methylation_type, 
             y = O242_Methyl_IR_4mC_promot$Percent_modified,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(6), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.3) # size label genomic positions
             })

# track 2 - 4mC IR for MAC101 O251
circos.track(O251_Methyl_IR_4mC_promot$Methylation_type, 
             y = O251_Methyl_IR_4mC_promot$Percent_modified)

# track 3 - 4mC IR for MAC101 T257
circos.track(T257_Methyl_IR_4mC_promot$Methylation_type, 
             y = T257_Methyl_IR_4mC_promot$Percent_modified)

# track 4 - 4mC IR for MAC101 T278
circos.track(T278_Methyl_IR_4mC_promot$Methylation_type, 
             y = T278_Methyl_IR_4mC_promot$Percent_modified)


circos.lines(O242_Methyl_IR_4mC_promot$Methylation_pos, O242_Methyl_IR_4mC_promot$Percent_modified, type = "h", straight = FALSE,
             col="slateblue", sector.index = get.current.sector.index(),
             track.index = 1)

circos.lines(O251_Methyl_IR_4mC_promot$Methylation_pos, O251_Methyl_IR_4mC_promot$Percent_modified, type = "h", straight = FALSE,
             col="slateblue4", sector.index = get.current.sector.index(),
             track.index = 2)
circos.lines(T257_Methyl_IR_4mC_promot$Methylation_pos, T257_Methyl_IR_4mC_promot$Percent_modified, type = "h", straight = FALSE,
             col="firebrick", sector.index = get.current.sector.index(),
             track.index = 3)
circos.lines(T278_Methyl_IR_4mC_promot$Methylation_pos, T278_Methyl_IR_4mC_promot$Percent_modified, type = "h", straight = FALSE,
             col="firebrick1", sector.index = get.current.sector.index(),
             track.index = 4)
# Add a title
title("MAC101 - potential \n promoter regions", line = -14, cex.main = 1.5)  # Adjust line for vertical position and size as needed

# Add legend outside the circular plot
legend("topleft", legend = c("O242", "O251", "T257", "T278"),
       col = c("slateblue", "slateblue4","firebrick","firebrick1"), lty = 1, lwd = 2, cex = 0.8, box.lty = 0)  # Customize position and style

# Clear and close
circos.clear()
dev.off()

