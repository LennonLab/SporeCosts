
# Load necessary libraries
library(dplyr)

# Load the data files
presence_absence <- read.csv("/Users/canankarakoccanan/genomes_v2/nucleotide_presence_absence.csv")
genome_sizes <- read.csv("/Users/canankarakoccanan/genomes_v2/genome_sizes_and_noncoding_nucleotides.csv")
phylo_names <- read.csv("/Users/canankarakoccanan/genomes_v2/phylo_names.csv")

# Merge the phylo_names with presence_absence to get the Taxid
names(presence_absence)[1] <- "Name"
presence_absence <- presence_absence %>%
  left_join(phylo_names, by = "Name")%>%
  mutate(across(Spo0A:CgeB.YkvP, as.numeric, .names = "numeric_{col}")) %>%
  rowwise() %>%
  mutate(numGenes = sum(c_across(numeric_Spo0A:numeric_CgeB.YkvP), na.rm = TRUE)) %>%
  ungroup() %>%
  select(-starts_with("numeric_"))

# Merge the presence_absence with genome_sizes
merged_data <- presence_absence %>%
  left_join(genome_sizes, by = "Taxid")

# Plotting example
library(ggplot2)

# Ensure noncoding_nucleotides is numeric
merged_data$Noncoding_Nucleotides <- as.numeric(merged_data$Noncoding_Nucleotides)

distance <- read.table("/Users/canankarakoccanan/genomes_v2/distances_to_closest_sporeformer_bootstrap_WRS.csv", sep = ",", dec = ".", header = T, stringsAsFactors = F)
distance$Taxid <- gsub("genome_", "", as.character(distance$Lost_Taxon))

merged_data$Taxid <- as.character(merged_data$Taxid)

merged_data_distance <- distance %>%
  left_join(merged_data, by = "Taxid")


ggplot(merged_data_distance, aes(x = as.numeric(Distance), y = (as.numeric(lost_nucleotides))/numGenes)) +
  geom_point(shape = 21, size = 3)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  theme_bw()+
  labs(x = "distance", y = "number spore genes")
