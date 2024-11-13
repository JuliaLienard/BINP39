
# Making a barplot with summary statistics of MicrobeMod using cutoff 0.7 (for STREME)

# Loading data 
library(readxl)
library(dplyr)
MicrobeMod_cutoff0_7_summaryMAC101 <- read_excel("MicrobeMod_cutoff0.7_summaryMAC101.xlsx")
MicrobeMod_cutoff0_7_summaryLU439 <- read_excel("MicrobeMod_cutoff0.7_summaryLU439.xlsx")
MicrobeMod_cutoff0_7_summary_MAC101_LU439 <- bind_rows(MicrobeMod_cutoff0_7_summaryMAC101,MicrobeMod_cutoff0_7_summaryLU439)

# 1 ) Methylation detection quality (% map reads)
# making a plot
library(ggplot2)
MethylationQual <- ggplot(MicrobeMod_cutoff0_7_summary_MAC101_LU439, aes(x = Sample, y = Nb_sites, fill = threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Methyl_type) +
  labs(title = "Number of Methylated Sites detected",
       x = "Sample",
       y = "Number of Sites",
       fill = "Methylation \n (%reads)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16) # Increase font size of facet labels
  )
MethylationQual
ggsave("MAC101_LU439_Methylation_qual_stats.png", MethylationQual, width = 7, height = 5, dpi = 300)


# 2) STREME 6mA - Grouping samples by motifs
# STREME MOTIFS 6mA - statistics -Grouping the motifs
MicrobeMod_cutoff0_7_summarySTREMEMAC101 <- read_excel("MicrobeMod_cutoff0.7_summarySTREMEMAC101.xlsx")
MicrobeMod_cutoff0_7_summarySTREMELU439 <- read_excel("MicrobeMod_cutoff0.7_summarySTREMELU439.xlsx")
MicrobeMod_cutoff0_7_summarySTREMEMAC101_LU439 <- bind_rows(MicrobeMod_cutoff0_7_summarySTREMEMAC101,MicrobeMod_cutoff0_7_summarySTREMELU439)
MicrobeMod_cutoff0_7_summarySTREMEMAC101_LU439_filt_6mA <- MicrobeMod_cutoff0_7_summarySTREMEMAC101_LU439 %>% filter((Motif != "nonassigned") & (Methyl_type == "6mA"))

#Removing Motif 5 that has a low E value STREME and counts for 0.3% of 6mA motifs and probable motifs
MicrobeMod_cutoff0_7_summarySTREMEMAC101_LU439_filt_6mA <- filter(MicrobeMod_cutoff0_7_summarySTREMEMAC101_LU439_filt_6mA, Motif_group == "Motif1" | Motif_group == "Motif2" | Motif_group =="Motif3" | Motif_group =="Motif4")


methylMOTIF6mA <- ggplot(MicrobeMod_cutoff0_7_summarySTREMEMAC101_LU439_filt_6mA, aes(x = Sample, y = as.numeric(Nb_sites), fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Motif_group) +
  labs(title = "Methylated motifs identified (6mA)",
       x = "Sample",
       y = "Number of Sites",
       fill = "") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 0, size = 11),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 11, face = "bold") # Increase font size of facet labels
  )

methylMOTIF6mA
ggsave("MAC101_LU439_Methylation_STREME6mA_stats.png", methylMOTIF6mA, width = 7, height = 6, dpi = 300)


# 3) Plot the 6mA motifs sequences
filtered_Motif6mA <- MicrobeMod_cutoff0_7_summarySTREMEMAC101_LU439_filt_6mA %>%
  group_by(Motif_group, Motif) %>%
  slice(1) %>%  # Select only one occurrence of each unique motif per group
  ungroup()

## Create the plot
library(ggrepel)
library(stringr)
# Wrap motif names to avoid overlap
filtered_Motif6mA$Motif_wrapped <- str_wrap(filtered_Motif6mA$Motif, width = 10)

# Create the horizontal bar plot with wrapped labels
Motif_6mA_seq <- ggplot(filtered_Motif6mA, aes(x = Motif_group, y = as.numeric(Nb_sites))) +
  geom_jitter(aes(color = Motif_group), size = 3, width = 0.2) +
  geom_text_repel(aes(label = Motif), size = 5, max.overlaps = 20, angle = 90) +  # Rotate labels to vertical
  labs(
       x = "Motif Group (6mA)",
       y = "Nb_sites") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14), 
        axis.text.y = element_text(angle = 90, hjust = 1, size = 14),
        axis.title.x = element_text(angle = 180, hjust = 0.5, size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none")
Motif_6mA_seq
ggsave("MAC101_LU439_Motif_group6mA.png", Motif_6mA_seq, width = 6, height = 8, dpi = 300)


# 4) STREME MOTIFS 4mC - grouping by samples, and by motig group => just MAC101

MicrobeMod_cutoff0_7_summarySTREMEMAC101_filt_4mC <- MicrobeMod_cutoff0_7_summarySTREMEMAC101 %>% filter((Motif != "nonassigned") & (Methyl_type == "4mC"))

methylMOTIF4mC <-
  ggplot(MicrobeMod_cutoff0_7_summarySTREMEMAC101_filt_4mC, aes(x = Sample, y = as.numeric(Nb_sites), fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  #facet_wrap(~Sample) +
  labs(title = "Methylated motifs identified (4mC)",
       x = "Sample",
       y = "Number of Sites",
       fill = "") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold") # Increase font size of facet labels
  )

methylMOTIF4mC
ggsave("MAC101_Methylation_STREME4mC_stats.png", methylMOTIF4mC, width = 7, height = 6, dpi = 300)


# 5) Plot the 4mC motifs sequences
filtered_Motif4mC <- MicrobeMod_cutoff0_7_summarySTREMEMAC101_filt_4mC %>%
  group_by(Motif_group, Motif) %>%
  slice(1) %>%  # Select only one occurrence of each unique motif per group
  ungroup()

## Create the plot
library(ggrepel)
library(stringr)
# Wrap motif names to avoid overlap
filtered_Motif4mC$Motif_wrapped <- str_wrap(filtered_Motif4mC$Motif, width = 10)

# Create the horizontal bar plot with wrapped labels
Motif_4mC_seq <- ggplot(filtered_Motif4mC, aes(x = Motif_group, y = as.numeric(Nb_sites))) +
  geom_jitter(aes(color = Motif_group), size = 3, width = 0.2) +
  geom_text_repel(aes(label = Motif), size = 5, max.overlaps = 20, angle = 90) +  # Rotate labels to vertical
  labs(
    x = "Motif Group (4mC)",
    y = "Nb_sites") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14), 
        axis.text.y = element_text(angle = 90, hjust = 1, size = 14),
        axis.title.x = element_text(angle = 180, hjust = 0.5, size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none")
Motif_4mC_seq
ggsave("MAC101_Motif_group4mC.png", Motif_4mC_seq, width = 6, height = 8, dpi = 300)
