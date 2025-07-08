library(phyloseq)
library(decontam)
library(tidyverse)
library(Matrix)
library(readxl) 
library(patchwork)
library(vegan)
library(ggforce)
library(lme4) 
library(lmerTest)
library(performance)
library(ggforce)
#remotes::install_github("david-barnett/microViz")
library(microViz)

# #replicating Alice's code

#importing cleaned bacteria data
VagData <- readRDS("/Users/caoyang/Desktop/Tetel Lab/datasets/vaginal_bacteria_cleanedv3.rds")
GutData <- readRDS("/Users/caoyang/Desktop/Tetel Lab/datasets/fecal_bacteria_cleanedv3.rds")
########################################################################
#cleaning vaginal data
otu_mat <- VagData %>%
  otu_table() %>%
  as("matrix")


# Transpose if taxa are rows
if (taxa_are_rows(VagData)) {
  otu_mat <- t(otu_mat)
}

most_abundant_taxa <- apply(otu_mat, 1, function(x) names(x)[which.max(x)])
tax_mat <- tax_table(VagData) %>% as("matrix")
most_abundant_species <- paste(
  tax_mat[most_abundant_taxa, "Genus"],
  tax_mat[most_abundant_taxa, "Species_exact"]
)
sample_data(VagData)$most_abundant_species <- most_abundant_species

species_to_cst <- data.frame(
  Species = c("crispatus", "gasseri", "iners", "jensenii"),
  CST = c("I", "II", "III", "V") # IV is anaerobic/diverse cluster
)

VagSample.df <- data.frame(sample_data(VagData))

VagSample.df <- VagSample.df %>%
  mutate(
    CST = case_when(
      str_detect(most_abundant_species, "gasseri") ~ "II",
      str_detect(most_abundant_species, "crispatus") ~ "I",
      str_detect(most_abundant_species, "iners") ~ "III",
      str_detect(most_abundant_species, "jensenii") ~ "V",
      TRUE ~ "IV" # Assign "IV" for diverse/anaerobic or unclassified species
    )
  )


VagSample.df.cl <- VagSample.df  %>%
  select(biome_id, logDate, CST, qr, SampleID, sampleType, most_abundant_species) %>%
  arrange(biome_id)

#VagSample.df.cl is organized by biome_id, but due to phyloseq constraint the sample_data is still organized by sample ID
sample_data(VagData) <- VagSample.df.cl #there are 1571samples

CST_NA <- VagSample.df.cl %>% 
  filter(is.na(CST))


########################################################################
#cleaning gut data
otu_matG <- GutData %>%
  otu_table() %>%
  as("matrix") %>%
  { if (taxa_are_rows(GutData)) t(.) else . }

most_abundant_taxaG <- apply(otu_matG, 1, function(x) names(x)[which.max(x)])
tax_matG <- tax_table(GutData) %>% as("matrix") 
most_abundant_speciesG <- sapply(most_abundant_taxaG, function(taxa) {
  if (taxa %in% rownames(tax_matG)) {
    paste(tax_matG[taxa, "Genus"], tax_matG[taxa, "Species"])
  } else {
    NA
  }
})
sample_data(GutData)$most_abundant_species <- most_abundant_speciesG

#adding b.vulg abundance column to Gut
tax_df_G <- data.frame(tax_table(GutData))
bvulg_asvs <- rownames(tax_df_G)[tax_df_G$Species == "vulgatus"]
bvulg_abund <- colSums(otu_table(GutData)[bvulg_asvs, , drop = FALSE])
total_reads <- colSums(otu_table(GutData))
bvulg_rel_abund <- bvulg_abund / total_reads

#add total_reads to sample_data
sample_data(GutData)$total_reads <- total_reads

# Add relative abundance to sample_data
sample_data(GutData)$gut_bvulg_rel_abundance <- bvulg_rel_abund

#adding lacto abundance column to Gut
lacto_asvs_G <- rownames(tax_df_G)[tax_df_G$Genus == "Lactobacillus"]
lacto_abund_G <- colSums(otu_table(GutData)[lacto_asvs_G, , drop = FALSE])
total_reads_G <- colSums(otu_table(GutData))
lacto_rel_abund_G <- lacto_abund_G / total_reads_G

# Add relative abundance to sample_data
sample_data(GutData)$lacto_rel_abundance_gut <- lacto_rel_abund_G

#adding Gardnerella
gard_asvs_G <- rownames(tax_df_G)[tax_df_G$Genus == "Gardnerella"]
gard_abund_G <- colSums(otu_table(GutData)[gard_asvs_G, , drop = FALSE])
gard_rel_abund_G <- gard_abund_G / total_reads_G
sample_data(GutData)$gard_rel_abundance_gut <- gard_rel_abund_G

#adding Mobiluncus
mobi_asvs_G <- rownames(tax_df_G)[tax_df_G$Genus == "Mobiluncus"]
mobi_abund_G <- colSums(otu_table(GutData)[mobi_asvs_G, , drop = FALSE])
mobi_rel_abund_G <- mobi_abund_G / total_reads_G
sample_data(GutData)$mobi_rel_abundance_gut <- mobi_rel_abund_G

#adding Megasphaera
mega_asvs_G <- rownames(tax_df_G)[tax_df_G$Genus == "Megasphaera"]
mega_abund_G <- colSums(otu_table(GutData)[mega_asvs_G, , drop = FALSE])
mega_rel_abund_G <- mega_abund_G / total_reads_G
sample_data(GutData)$mega_rel_abundance_gut <- mega_rel_abund_G

#adding Sneathia
snea_asvs_G <- rownames(tax_df_G)[tax_df_G$Genus == "Sneathia"]
snea_abund_G <- colSums(otu_table(GutData)[snea_asvs_G, , drop = FALSE])
snea_rel_abund_G <- snea_abund_G / total_reads_G
sample_data(GutData)$snea_rel_abundance_gut <- snea_rel_abund_G

#adding Prevotella abundance column to Gut
prevo_asvs_G <- rownames(tax_df_G)[tax_df_G$Genus %in% c("Prevotella", "Prevotella_7", "Prevotella_9")] 
prevo_abund_G <- colSums(otu_table(GutData)[prevo_asvs_G, , drop = FALSE])
prevo_rel_abund_G <- prevo_abund_G / total_reads_G

# Add relative abundance to sample_data
sample_data(GutData)$prevo_rel_abundance_gut <- prevo_rel_abund_G


#Add Alistipes putredinis (a species) relative abundance to Gut
putredinis_asvs_G <- rownames(tax_df_G)[tax_df_G$Species == "putredinis"]
putredinis_abund_G <- colSums(otu_table(GutData)[putredinis_asvs_G, , drop = FALSE])
putredinis_rel_abund_G <- putredinis_abund_G / total_reads_G
sample_data(GutData)$putredinis_rel_abundance_gut <- putredinis_rel_abund_G


#add Staphylococcus epidermidis to check whether a sample has been mislabeled
staphy_asvs_G <- rownames(tax_df_G)[tax_df_G$Species == "epidermidis"]
staphy_abund_G <- colSums(otu_table(GutData)[staphy_asvs_G, , drop = FALSE])
staphy_rel_abund_G <- staphy_abund_G / total_reads_G
sample_data(GutData)$staphy_rel_abundance_gut <- staphy_rel_abund_G

staphyG <- tax_df_G %>% filter(Genus == "Staphylococcus")

#add Lactobacillus iners
Liners_asvs_G <- rownames(tax_df_G)[tax_df_G$Species == "iners"]
Liners_abund_G <- colSums(otu_table(GutData)[Liners_asvs_G, , drop = FALSE])
Liners_rel_abund_G <- Liners_abund_G / total_reads_G
sample_data(GutData)$Liners_rel_abundance_gut <- Liners_rel_abund_G


####################################################################################
#check for vaginal samples in gut
Gutsample <- data.frame(sample_data(GutData))
suspicious <- Gutsample %>%
    filter(str_detect(most_abundant_species, "Lactobacillus"))

colnames(suspicious)

dim(suspicious)

#running PCA
SuspiciousGut <- subset_samples(GutData, sample_names(GutData) %in% suspicious$SampleID)
CombinedData <- merge_phyloseq(SuspiciousGut, VagData)
sample_data(CombinedData)$Source <- ifelse(
  sample_names(CombinedData) %in% sample_names(SuspiciousGut),
  "SuspiciousGut", "Vagina"
)
CombinedData <- CombinedData %>% tax_fix()
CombinedData <- CombinedData %>%
  tax_fix(unknowns = c(
    "acidifaciens", "anthropi", "asaccharolytica", "avium", "bacterium",
    "beijerinckii", "buccalis", "caccae", "caecimuris", "denticola",
    "epidermidis", "faecalis", "faecis", "faecium", "finegoldii",
    "gilardii", "gordonii", "haemolyticus", "hominis", "ignava",
    "intestinalis", "johnsonii", "koreensis", "lactis", "massiliense",
    "massiliensis", "obesi", "oralis", "oris", "pneumoniae", "rhizophila",
    "rhizosphaerae", "sanguinis", "shahii", "simulans", "soli", "sputigena",
    "stercoris", "taiwanensis", "terrae", "timonensis", "urealyticum",
    "vaginalis", "veronii", "vulgaris"
  ))
CombinedData %>%
  tax_transform("clr") %>%  # No rank = ASV-level
  ord_calc(method = "PCA") %>%
  ord_plot(color = "Source", size = 2)


##############################################################
# otu_df <- as.data.frame(otu_table(GutData))
# tax_df <- as.data.frame(tax_table(GutData))
# 
# # If taxa are columns (not rows), transpose
# if(taxa_are_rows(GutData)) {
#   otu_df <- t(otu_df)
# }
# 
# # Combine counts with taxonomy
# otu_tax_df <- otu_df %>%
#   as_tibble(rownames = "ASV") %>%
#   left_join(tax_df %>% as_tibble(rownames = "ASV"), by = "ASV")
# 
# # Filter for Staphylococcus genus
# staph_df <- otu_tax_df %>%
#   filter(Genus == "Staphylococcus") %>%
#   rowwise() %>%
#   mutate(Total_Count = sum(c_across(where(is.numeric)))) %>%
#   ungroup()
# 
# # Summarize by species
# staph_species_counts <- staph_df %>%
#   group_by(Species) %>%
#   summarise(Total_Count = sum(Total_Count)) %>%
#   arrange(desc(Total_Count))
# 
# # View result
# print(staph_species_counts)





######################
#vaginal data

#adding lactobacillus rel abundance column to Vag
tax_df <- data.frame(tax_table(VagData))
Lacto_asvs <- rownames(tax_df)[tax_df$Genus == "Lactobacillus"] 
Lacto_abund <- colSums(otu_table(VagData)[Lacto_asvs, , drop = FALSE])
total_reads_V <- colSums(otu_table(VagData))
Lacto_rel_abund <- Lacto_abund / total_reads_V

#add total reads to sample data
sample_data(VagData)$total_reads_V <- total_reads_V

# Add relative abundance to sample_data
sample_data(VagData)$lacto_rel_abundance_vag <- Lacto_rel_abund

#adding Gardnerella
gard_asvs_V <- rownames(tax_df)[tax_df$Genus == "Gardnerella"] 
gard_abund_V <- colSums(otu_table(VagData)[gard_asvs_V, , drop = FALSE])
gard_rel_abund_V <- gard_abund_V / total_reads_V
sample_data(VagData)$gard_rel_abundance_vag <- gard_rel_abund_V

#adding Mobiluncus
mobi_asvs_V <- rownames(tax_df)[tax_df$Genus == "Mobiluncus"] 
mobi_abund_V <- colSums(otu_table(VagData)[mobi_asvs_V, , drop = FALSE])
mobi_rel_abund_V <- mobi_abund_V / total_reads_V
sample_data(VagData)$mobi_rel_abundance_vag <- mobi_rel_abund_V

#adding Megasphaera
mega_asvs_V <- rownames(tax_df)[tax_df$Genus == "Megasphaera"] 
mega_abund_V <- colSums(otu_table(VagData)[mega_asvs_V, , drop = FALSE])
mega_rel_abund_V <- mega_abund_V / total_reads_V
sample_data(VagData)$mega_rel_abundance_vag <- mega_rel_abund_V

#adding Sneathia
snea_asvs_V <- rownames(tax_df)[tax_df$Genus == "Sneathia"] 
snea_abund_V <- colSums(otu_table(VagData)[snea_asvs_V, , drop = FALSE])
snea_rel_abund_V <- snea_abund_V / total_reads_V
sample_data(VagData)$snea_rel_abundance_vag <- snea_rel_abund_V

#add Prevotella to sample data
prev_asvs <- rownames(tax_df)[tax_df$Genus %in% c("Prevotella", "Prevotella_7", "Prevotella_9")] 
prev_abund <- colSums(otu_table(VagData)[prev_asvs, , drop = FALSE])
prev_rel_abund <- prev_abund / total_reads_V
# Add relative abundance to sample_data
sample_data(VagData)$prevo_rel_abundance_vag <- prev_rel_abund


#add Alistipes putredinis 
putredinis_asvs_V <- rownames(tax_df)[tax_df$Species == "putredinis"] 
putredinis_abund_V <- colSums(otu_table(VagData)[putredinis_asvs_V, , drop = FALSE])
putredinis_rel_abund_V <- putredinis_abund_V / total_reads_V
sample_data(VagData)$putredinis_rel_abundance_vag <- putredinis_rel_abund_V

#add Staphylococcus epidermidis 
staphy_asvs_V <- rownames(tax_df)[tax_df$Species == "epidermidis"] 
staphy_abund_V <- colSums(otu_table(VagData)[staphy_asvs_V, , drop = FALSE])
staphy_rel_abund_V <- staphy_abund_V / total_reads_V
sample_data(VagData)$staphy_rel_abundance_vag <- staphy_rel_abund_V

staphy <- tax_df %>% filter(Genus == "Staphylococcus")

table(staphyG$Species)
table(staphy$Species)

#add L.iners
Liners_asvs_V <- rownames(tax_df)[tax_df$Species == "iners"] 
Liners_abund_V <- colSums(otu_table(VagData)[Liners_asvs_V, , drop = FALSE])
Liners_rel_abund_V <- Liners_abund_V / total_reads_V
sample_data(VagData)$Liners_rel_abundance_vag <- Liners_rel_abund_V




########################################################################
#matching these sample data of gut and vaginal data

GutforMatch <- data.frame(sample_data(GutData))
VagforMatch <- data.frame(sample_data(VagData))

GutforMatch_trimmed <- GutforMatch %>%
  group_by(biome_id, logDate) %>%
  slice(1) %>%  # keep only the first sample in each group
  ungroup()

VagforMatch_trimmed <- VagforMatch %>%
  group_by(biome_id, logDate) %>%
  slice(1) %>%
  ungroup()

cross.df <- left_join(GutforMatch_trimmed, VagforMatch_trimmed, by = c("biome_id", "logDate"))

cross.CSTNA <- cross.df %>% filter(is.na(CST)) #NA comes from the gut data that did not have matched vaginal data

#rename and remove meaningless columns
cross.df <- cross.df %>%
  rename(
    SampleID_gut = SampleID.x,
    qr_gut = qr.x,
    SampleID_vaginal = SampleID.y,
    qr_vaginal = qr.y,
    most_abundant_species_gut = most_abundant_species.x,
    most_abundant_species_vaginal = most_abundant_species.y
  ) %>%
  select(
    -is_blank,
    -status,
    -sampleType.x,
    -sampleType.y
  )
######################################################################
#mislabel check
summary(cross.df$staphy_rel_abundance_gut)
summary(cross.df$staphy_rel_abundance_vag)

########################################################################
cross.df_CST <- cross.df %>% filter(CST %in% c("I", "II", "III", "IV", "V")) #don't omit NA
#659 has CSTI, 70 has CSTII, 191 has CSTIII, 54 has CSTIV, 30 has CST V
dim(cross.df)
########################################################################
#finding correlation of gut and vaginal microbiome 

#Lactobacillus
#CST coloring
ggplot(cross.df_CST, aes(x = sqrt(lacto_rel_abundance_gut), 
                     y = sqrt(lacto_rel_abundance_vag), 
                     color = CST)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +  # non-linear curves
  labs(
    title = "Relationship Between Gut and Vaginal Lactobacillus Abundance",
    x = "Gut Lactobacillus Relative Abundance (sqrt)",
    y = "Vaginal Lactobacillus Relative Abundance (sqrt)",
    color = "CST"
  ) +
  theme_minimal()

#no CST coloring
ggplot(cross.df, aes(x = sqrt(lacto_rel_abundance_gut), 
                     y = sqrt(lacto_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal Lactobacillus Abundance",
    x = "Gut Lactobacillus Relative Abundance (sqrt)",
    y = "Vaginal Lactobacillus Relative Abundance (sqrt)"
  ) +
  theme_minimal()

cor.test(cross.df$lacto_rel_abundance_gut, 
         cross.df$lacto_rel_abundance_vag, 
         method = "spearman")


#Prevotella
ggplot(cross.df_CST, aes(x = sqrt(prevo_rel_abundance_gut), 
                         y = sqrt(prevo_rel_abundance_vag), 
                         color = CST)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +  # non-linear curves
  labs(
    title = "Relationship Between Gut and Vaginal Prevotella Abundance",
    x = "Gut Prevotella Relative Abundance (sqrt)",
    y = "Vaginal Pre Relative Abundance (sqrt)",
    color = "CST"
  ) +
  theme_minimal()

ggplot(cross.df, aes(x = sqrt(prevo_rel_abundance_gut), 
                     y = sqrt(prevo_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal Prevotella Abundance",
    x = "Gut Prevotella Relative Abundance (sqrt)",
    y = "Vaginal Prevotella Relative Abundance (sqrt)"
  ) +
  theme_minimal()

cor.test(cross.df$prevo_rel_abundance_gut, 
         cross.df$prevo_rel_abundance_vag, 
         method = "spearman")

#Gardnerella
ggplot(cross.df, aes(x = log10(gard_rel_abundance_gut), 
                     y = log10(gard_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal Gardnerella Abundance",
    x = "Gut Gardnerella Relative Abundance (log10)",
    y = "Vaginal Gardnerella Relative Abundance (log10)"
  ) +
  theme_minimal()

cor.test(cross.df$gard_rel_abundance_gut, 
         cross.df$gard_rel_abundance_vag, 
         method = "spearman")

#Mobiluncus
ggplot(cross.df, aes(x = log10(mobi_rel_abundance_gut), 
                     y = log10(mobi_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal Mobiluncus Abundance",
    x = "Gut Mobiluncus Relative Abundance (log10)",
    y = "Vaginal Mobiluncus Relative Abundance (log10)"
  ) +
  theme_minimal()

cor.test(cross.df$mobi_rel_abundance_gut, 
         cross.df$mobi_rel_abundance_vag, 
         method = "spearman")

#Megasphaera #too many zeros, no meaning
ggplot(cross.df, aes(x = mega_rel_abundance_gut, 
                     y = mega_rel_abundance_vag)) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal Megasphaera Abundance",
    x = "Gut Megasphaera Relative Abundance",
    y = "Vaginal Megasphaera Relative Abundance"
  ) +
  theme_minimal()

#Sneathia #too many zeros
ggplot(cross.df, aes(x = snea_rel_abundance_gut, 
                     y = snea_rel_abundance_vag)) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal Sneathia Abundance",
    x = "Gut Sneathia Relative Abundance",
    y = "Vaginal Sneathia Relative Abundance"
  ) +
  theme_minimal()

cor.test(cross.df$snea_rel_abundance_gut, 
         cross.df$snea_rel_abundance_vag, 
         method = "spearman")
########################################################################
#trying to see lactobacillus relative abundanec fluctuation over time
cross.df_clean <- cross.df %>%
  filter(!is.na(logDate), !is.na(lacto_rel_abundance_vag))

ggplot(cross.df_clean, aes(x = logDate, y = lacto_rel_abundance_vag)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  labs(title = "Lactobacillus Relative Abundance Over Time",
       x = "Log Date",
       y = "Relative Abundance") +
  theme_minimal()
########################################################################
#trying to remake the barplot of CST and gut species
species_freq <- cross.df_CST %>%
  group_by(CST, most_abundant_species_gut) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(CST) %>%
  mutate(percent = count / sum(count) * 100)

top_species <- species_freq %>%
  group_by(CST) %>%
  slice_max(order_by = percent, n = 1)

cst_counts <- cross.df_CST %>%
  group_by(CST) %>%
  summarise(n = n(), .groups = "drop")

cst_labels <- setNames(
  paste0(cst_counts$CST, " (n = ", cst_counts$n, ")"),
  cst_counts$CST
)

variance_df <- species_freq %>%
  group_by(CST) %>%
  summarise(variance = var(percent), .groups = "drop")


ggplot(top_species, aes(x = reorder(most_abundant_species_gut, -percent), y = percent, fill = CST)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ CST, scales = "free_x", labeller = labeller(CST = cst_labels)) +
  labs(
    x = "Most Abundant Gut Species",
    y = "Percentage of Samples (%)",
    title = "Most Common Gut Species by Vaginal CST"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(variance_df, aes(x = CST, y = variance, fill = CST)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Vaginal CST",
    y = "Variance of % Gut Species",
    title = "Variance in Gut Species Composition by Vaginal CST"
  ) +
  theme_minimal()


ggplot(species_freq, aes(x = CST, y = percent, fill = CST)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black", size = 1) +
  labs(
    title = "Distribution of Gut Species Percentages by Vaginal CST",
    x = "Vaginal CST",
    y = "Percentage of Dominant Gut Species"
  ) +
  theme_minimal()


bvulg_data <- cross.df_CST %>%
  filter(most_abundant_species_gut == "Bacteroides vulgatus")

# Plot abundance of B. vulgatus across CSTs
ggplot(bvulg_data, aes(x = CST, y = gut_bvulg_rel_abundance, fill = CST)) +
  geom_boxplot() +
  labs(
    title = "Relative Abundance of B. vulgatus by Vaginal CST (samples)",
    y = "Relative Abundance",
    x = "Vaginal CST"
  ) +
  theme_minimal()

###############################################################
#taking medical history into account
medhis <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/cleaned_Report 9-Volunteer Medical History - cleaned_Report 9-Volunteer Medical History.csv")

cross.df <- cross.df %>%
  left_join(medhis, by = "biome_id") %>%
  select(-logDate.y) %>%  
  rename_with(~ str_replace_all(., "\\.x$|\\.y$", ""))

colnames(cross.df)
######################################################################
#making stacked bar


# target_species <- c(
#   "Bacteroides vulgatus", "Finegoldia magna", "Collinsella aerofaciens",
#   "Fenollaria", "Lactobacillus helveticus", "Bacteroides dorei"
# )

#creating a new dataframe that categorize the dominant species into either one of the above of "other"
# plot_df <- cross.df %>%
#   mutate(species_group = ifelse(most_abundant_species_gut %in% target_species,
#                                 most_abundant_species_gut,
#                                 "Other"))


plot_df <- cross.df %>%
  mutate(
    species_group = case_when(
      str_detect(most_abundant_species_gut, "Bacteroides vulgatus") ~ "Bacteroides vulgatus",
      str_detect(most_abundant_species_gut, "Finegoldia magna") ~ "Finegoldia magna",
      str_detect(most_abundant_species_gut, "Collinsella aerofaciens") ~ "Collinsella aerofaciens",
      str_detect(most_abundant_species_gut, "Fenollaria") ~ "Fenollaria",
      str_detect(most_abundant_species_gut, "Lactobacillus helveticus") ~ "Lactobacillus helveticus",
      str_detect(most_abundant_species_gut, "Bacteroides dorei") ~ "Bacteroides dorei",
      TRUE ~ "Other"
    )
  )



#compute the frequency (percentage) of each most dominant species 
species_freq <- plot_df %>%
  filter(!is.na(CST)) %>%
  group_by(CST, species_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(CST) %>%
  mutate(percent = count / sum(count) * 100)

ggplot(species_freq, aes(x = CST, y = percent, fill = species_group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Distribution of Dominant Gut Microbiome Species by CST (samples)",
    x = "Vaginal CST",
    y = "Percentage of Samples",
    fill = "Gut Species"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Bacteroides vulgatus" = "#1f77b4",
      "Finegoldia magna" = "#ff7f0e",
      "Collinsella aerofaciens" = "#2ca02c",
      "Fenollaria" = "#d62728",
      "Lactobacillus helveticus" = "#9467bd",
      "Bacteroides dorei" = "#8c564b",
      "Other" = "gray80"
    )
  )


cst_labels <- data.frame(
  CST = c("I", "II", "III", "IV", "V"),
  label = c("n = 659", "n = 70", "n = 191", "n = 54", "n = 30"),
  y = 105  # y-position above the 100% bar
)

# Plot with manual sample size labels
ggplot(species_freq, aes(x = CST, y = percent, fill = species_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(data = cst_labels,
            aes(x = CST, y = y, label = label),
            inherit.aes = FALSE,
            size = 3.5) +
  labs(
    title = "Distribution of Dominant Gut Microbiome Species by CST (samples)",
    x = "Vaginal CST",
    y = "Percentage of Samples",
    fill = "Gut Species"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Bacteroides vulgatus" = "#1f77b4",
      "Finegoldia magna" = "#ff7f0e",
      "Collinsella aerofaciens" = "#2ca02c",
      "Fenollaria" = "#d62728",
      "Lactobacillus helveticus" = "#9467bd",
      "Bacteroides dorei" = "#8c564b",
      "Other" = "gray80"
    )
  )

contingency_tbl <- table(plot_df$CST, plot_df$species_group)

# Run Chi-squared test
chisq.test(contingency_tbl)

#is CST IV significantly different from otehrs?
plot_dfIV <- plot_df %>%
  mutate(CST_group = ifelse(CST == "IV", "CST_IV", "Other_CST"))

# Contingency table: CST IV vs. others across species_group
cst_iv_tbl <- table(plot_dfIV$CST_group, plot_dfIV$species_group)

# View the table (optional)
print(cst_iv_tbl)

# Run Chi-squared test
chisq_result <- chisq.test(cst_iv_tbl)
print(chisq_result)
#yes CST IV is indeed different from others

#what about CST V?
plot_dfV <- plot_df %>%
  mutate(CST_group = ifelse(CST == "V", "CST_V", "Other_CST"))
cst_v_tbl <- table(plot_dfV$CST_group, plot_dfV$species_group)
print(cst_v_tbl)
chisq_result <- chisq.test(cst_v_tbl)
print(chisq_result)
#CST V difference significance

chi <- chisq.test(table(plot_df$CST, plot_df$species_group))

# Extract standardized residuals
std_res <- chi$stdres
# 
# # View as a heatmap
# pheatmap(std_res,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          display_numbers = TRUE,
#          main = "Standardized Residuals (Chi-squared)")
# 



######################################################################
genera <- c("lacto", "gard", "mobi", "mega", "snea", "prevo")
pairs_df <- expand.grid(vag = genera, gut = genera, stringsAsFactors = FALSE)
get_colname <- function(genus, site) {
  paste0(genus, "_rel_abundance_", site)
}
corr_results <- pairs_df %>%
  rowwise() %>%
  mutate(
    vag_col = get_colname(vag, "vag"),
    gut_col = get_colname(gut, "gut"),
    # Compute Spearman correlation and p-value
    cor_test = list(cor.test(
      sqrt(cross.df[[gut_col]]),  # you can transform if you want sqrt
      sqrt(cross.df[[vag_col]]),
      method = "spearman",
      exact = FALSE, # avoids ties warning
      use = "complete.obs"
    )),
    rho = cor_test$estimate,
    p_value = cor_test$p.value
  ) %>%
  ungroup() %>%
  select(vag, gut, rho, p_value)

ggplot(corr_results, aes(x = gut, y = vag, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = sprintf("%.2f", rho)), color = "black", size = 4) +
  labs(
    title = "Spearman Correlations Between Vaginal and Gut Species Abundances",
    x = "Gut Species",
    y = "Vaginal Species",
    fill = "Spearman rho"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

############################################################
#trying to create a bunch of scatterplot
# Generate scatter data
scatter_data <- pairs_df %>%
  rowwise() %>%
  mutate(
    gut_val = list(sqrt(cross.df[[paste0(gut, "_rel_abundance_gut")]])),
    vag_val = list(sqrt(cross.df[[paste0(vag, "_rel_abundance_vag")]])),
    id = paste(gut, vag, sep = "_")
  ) %>%
  unnest(cols = c(gut_val, vag_val)) %>%
  ungroup()

# Rename columns for clarity
scatter_data <- scatter_data %>%
  rename(gut_abund = gut_val, vag_abund = vag_val)

# Plot
ggplot(scatter_data, aes(x = gut_abund, y = vag_abund)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "red", size = 0.5) +
  facet_grid(rows = vars(vag), cols = vars(gut)) +
  theme_minimal() +
  labs(
    x = "Gut genus abundance (sqrt)",
    y = "Vaginal genus abundance (sqrt)",
    title = "Scatterplots with Spline: Gut vs. Vaginal Genus Abundances"
  ) +
  theme(
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 5)
  )

##################################################
#creating a scatterplot of shannon diversity correlation between vaginal and gut per person
# plot(x = cross.df$vaginal_shannon, y = cross.df$gut_shannon,
#      main = "Average Shannnon Diversity Correlation",
#      xlab = "Vaginal Microbiome",
#      ylab = "Gut Microbiome",
#      pch = 16)           
# 
# spline_fit <- smooth.spline(cross.df$vaginal_shannon, cross.df$gut_shannon)
# lines(spline_fit, col = "blue", lwd = 2)

################################################################
#picking one person with the most samples
sample_counts <- cross.df %>%
  count(biome_id, name = "num_samples") %>%
  arrange(desc(num_samples))
sample_counts #id 65 had 51 samples; 66 had 48 samples; 11 had 47 samples

#looking at how participant 66 Lacto fluctuate over time

df_66 <- cross.df %>%
  filter(biome_id == 66)

df_66_long <- df_66 %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, CST) %>%
  pivot_longer(cols = starts_with("lacto"), names_to = "Site", values_to = "Abundance")

df_66_clean <- df_66_long %>%
  filter(!is.na(logDate) & !is.na(Abundance) & is.finite(Abundance))

ggplot(df_66_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1, linetype = "solid") +
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 66)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site"
  ) +
  theme_minimal()

ggplot(df_66, aes(x = lacto_rel_abundance_gut, y = lacto_rel_abundance_vag)) +
  geom_point(size = 3, color = "steelblue") +
  geom_path(arrow = arrow(length = unit(0.2, "cm")), color = "gray40") +
  labs(
    title = "Gut vs Vaginal Lactobacillus Over Time (Participant 66)",
    x = "Gut Lactobacillus Relative Abundance",
    y = "Vaginal Lactobacillus Relative Abundance"
  ) +
  theme_minimal()

#65
df_65 <- cross.df %>%
  filter(biome_id == 65)

df_65_long <- df_65 %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, CST) %>%
  pivot_longer(cols = starts_with("lacto"), names_to = "Site", values_to = "Abundance")

df_65_clean <- df_65_long %>%
  filter(!is.na(logDate) & !is.na(Abundance) & is.finite(Abundance))

ggplot(df_65_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1, linetype = "solid") +
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 65)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site"
  ) +
  theme_minimal()

#11
df_11 <- cross.df %>%
  filter(biome_id == 11)

df_11_long <- df_11 %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, CST) %>%
  pivot_longer(cols = starts_with("lacto"), names_to = "Site", values_to = "Abundance")

df_11_clean <- df_11_long %>%
  filter(!is.na(logDate) & !is.na(Abundance) & is.finite(Abundance))

ggplot(df_11_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1, linetype = "solid") +
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 11)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site"
  ) +
  theme_minimal()

df_11_clean <- df_11_clean %>% arrange(logDate)

# Identify stretches where CST is constant
cst_spans <- df_11_clean %>%
  select(logDate, CST) %>%
  distinct() %>%
  mutate(
    change = CST != lag(CST, default = first(CST)),
    phase_id = cumsum(change)
  ) %>%
  group_by(phase_id, CST) %>%
  summarise(start = min(logDate), end = max(logDate), .groups = "drop")

ggplot(df_11_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_rect(data = cst_spans, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = CST),
            inherit.aes = FALSE, alpha = 0.2) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 11)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST Phase"
  ) +
  theme_minimal()


cst_spans <- df_66_clean %>%
  select(logDate, CST) %>%
  distinct() %>%
  mutate(
    change = CST != lag(CST, default = first(CST)),
    phase_id = cumsum(change)
  ) %>%
  group_by(phase_id, CST) %>%
  summarise(start = min(logDate), end = max(logDate), .groups = "drop")

ggplot(df_66_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_rect(data = cst_spans, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = CST),
            inherit.aes = FALSE, alpha = 0.2) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 66)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST Phase"
  ) +
  theme_minimal()


cst_spans <- df_65_clean %>%
  select(logDate, CST) %>%
  distinct() %>%
  mutate(
    change = CST != lag(CST, default = first(CST)),
    phase_id = cumsum(change)
  ) %>%
  group_by(phase_id, CST) %>%
  summarise(start = min(logDate), end = max(logDate), .groups = "drop")

ggplot(df_65_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_rect(data = cst_spans, 
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = CST),
            inherit.aes = FALSE, alpha = 0.2) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 65)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST Phase"
  ) +
  theme_minimal()


cst_rects <- df_66_clean %>%
  distinct(logDate, CST)

# 2. Plot with CST-colored backgrounds, including NA
ggplot(df_66_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  # CST background per day
  geom_tile(data = cst_rects,
            aes(x = logDate, y = 0, fill = CST),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  
  # Foreground data: points, lines, smooth
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  
  # Titles and themes
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 66)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal() +
  
  # 3. Assign white color to NA CST
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"  # This goes outside the values list
  
    )

################################################################
#locating the participants with the most gut dominant species shift
shifts_per_participant <- cross.df %>%
  arrange(biome_id, logDate) %>%  # sort by participant and time
  group_by(biome_id) %>%
  mutate(
    prev_species = lag(most_abundant_species_gut),                  # previous time point
    species_changed = most_abundant_species_gut != prev_species     # TRUE when species changes
  ) %>%
  summarise(num_shifts = sum(species_changed, na.rm = TRUE)) %>%  # count changes
  arrange(desc(num_shifts))

#who shifted the most
head(shifts_per_participant, 10) #11 #66 #64 #41 #12 #332

#trying participant 64
df_64 <- cross.df %>%
  filter(biome_id == 64)

df_64_long <- df_64 %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, CST) %>%
  pivot_longer(cols = starts_with("lacto"), names_to = "Site", values_to = "Abundance")

df_64_clean <- df_64_long %>%
  filter(!is.na(logDate) & !is.na(Abundance) & is.finite(Abundance))

cst_rects64 <- df_64_clean %>%
  distinct(logDate, CST)

ggplot(df_65_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(data = cst_rects64,
            aes(x = logDate, y = 0, fill = CST),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 64)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"  # This goes outside the values list
    
  )

#trying participant 41
df_41 <- cross.df %>%
  filter(biome_id == 41)

df_41_long <- df_41 %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, CST) %>%
  pivot_longer(cols = starts_with("lacto"), names_to = "Site", values_to = "Abundance")

df_41_clean <- df_41_long %>%
  filter(!is.na(logDate) & !is.na(Abundance) & is.finite(Abundance))

cst_rects41 <- df_41_clean %>%
  distinct(logDate, CST)

ggplot(df_41_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(data = cst_rects41,
            aes(x = logDate, y = 0, fill = CST),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 41)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"  # This goes outside the values list
    
  )
##########################################################
#viewing the cross site correlation of different species in details
ggplot(cross.df, aes(x = sqrt(snea_rel_abundance_gut), 
                     y = sqrt(lacto_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut Sneathia and Vaginal Lactobacillus Abundance",
    x = "Gut Sneathia Relative Abundance (sqrt)",
    y = "Vaginal Lactobacillus Relative Abundance (sqrt)"
  ) +
  theme_minimal()

ggplot(cross.df, aes(x = sqrt(snea_rel_abundance_gut), 
                     y = sqrt(prevo_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut Sneathia and Vaginal Prevotella Abundance",
    x = "Gut Sneathia Relative Abundance (sqrt)",
    y = "Vaginal Prevotella Relative Abundance (sqrt)"
  ) +
  theme_minimal()

ggplot(cross.df, aes(x = sqrt(mobi_rel_abundance_gut), 
                     y = sqrt(lacto_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  #geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut Mobiluncus and Vaginal Lactobacillus Abundance",
    x = "Gut Mobiluncus Relative Abundance (sqrt)",
    y = "Vaginal Lactobacillus Relative Abundance (sqrt)"
  ) +
  theme_minimal()

##################################################################
#writing the merged dataset into csv
write.csv(cross.df, file = "/Users/caoyang/Desktop/Tetel Lab/datasets/microbiome_crosstalk_merged_abund.csv", row.names = TRUE)

##################################################################
#taking antibiotics into account to lactobacillus correlation
cross_lac <- cross.df %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, taken_antibiotics) %>%
  pivot_longer(cols = starts_with("lacto_rel_abundance"),
               names_to = "Site",
               values_to = "Abundance") %>%
  mutate(Site = recode(Site,
                       "lacto_rel_abundance_gut" = "Gut",
                       "lacto_rel_abundance_vag" = "Vagina"))
# no Antibiotic group
crossAnti1<- cross_lac %>% filter(taken_antibiotics == 1)
# antibiotic group
crossAnti0 <- cross_lac %>% filter(taken_antibiotics == 0)

# Plot for antibiotic group
p1 <- ggplot(crossAnti1, aes(x = logDate, y = Abundance, color = Site)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5) +
  labs(title = "Without Antibiotic Use",
       x = "Time", y = "Lactobacillus Relative Abundance") +
  theme_minimal()

# Plot for no antibiotic group
p2 <- ggplot(crossAnti0, aes(x = logDate, y = Abundance, color = Site)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5) +
  labs(title = "With Antibiotic Use",
       x = "Time", y = "Lactobacillus Relative Abundance") +
  theme_minimal()

# Show them side-by-side
p1 + p2
############################################################################
#temporal lacto correlation by participant without antibiotics
df_43 <- cross.df %>%
  filter(biome_id == 43)
df_43_long <- df_43 %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, CST) %>%
  pivot_longer(cols = starts_with("lacto"), names_to = "Site", values_to = "Abundance")
df_43_clean <- df_43_long %>%
  filter(!is.na(logDate) & !is.na(Abundance) & is.finite(Abundance))
cst_rects43 <- df_43_clean %>%
  distinct(logDate, CST)
ggplot(df_43_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(data = cst_rects41,
            aes(x = logDate, y = 0, fill = CST),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 43, antibiotics)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"  # This goes outside the values list
    
  )


###################################################################
#creating a new column that calculates the difference between gut and vaginal lactobacillus abundance
cross.df$lacto_dif <- cross.df$lacto_rel_abundance_vag - cross.df$lacto_rel_abundance_gut

#visualize the difference fluctuation over time
ggplot(cross.df, aes(x = logDate, y = lacto_dif)) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_point(color = "steelblue", alpha = 0.6) +
  labs(
    title = "Difference in Lactobacillus Relative Abundance (Gut - Vaginal) Over Time",
    x = "Date",
    y = "Lactobacillus Abundance Difference"
  ) +
  theme_minimal()
###################################################################
#PCOA
# abundance_matrix <- cross.df %>%
#   select(-biome_id, -site)
# 
# # Bray-Curtis distance
# dist_matrix <- vegdist(abundance_matrix, method = "bray")

#################################################################
#the participant who's been consistently taking SSRI
df_60 <- cross.df %>%
  filter(biome_id == 60)
df_60_long <- df_60 %>%
  select(logDate, lacto_rel_abundance_gut, lacto_rel_abundance_vag, CST) %>%
  pivot_longer(cols = starts_with("lacto"), names_to = "Site", values_to = "Abundance")
df_60_clean <- df_60_long %>%
  filter(!is.na(logDate) & !is.na(Abundance) & is.finite(Abundance))
cst_rects60 <- df_60_clean %>%
  distinct(logDate, CST)
ggplot(df_60_clean, aes(x = logDate, y = Abundance, color = Site, group = Site)) +
  geom_tile(data = cst_rects41,
            aes(x = logDate, y = 0, fill = CST),
            width = 1, height = Inf, alpha = 0.2, inherit.aes = FALSE) +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
  labs(
    title = "Lactobacillus Relative Abundance Over Time (Participant 60, SSRI)",
    x = "Date",
    y = "Relative Abundance",
    color = "Site",
    fill = "CST"
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "I" = "#1b9e77",
      "II" = "#d95f02",
      "III" = "#7570b3",
      "IV" = "#e7298a",
      "V" = "#66a61e"
    ),
    na.value = "white"  # This goes outside the values list
    
  )

################################################################
#ssri and lacto
ssri_users <- c(11,14,16,24,27,28,30,33,51,60,73)
ssridf <- cross.df %>%
  mutate(
    SSRI_status = if_else(biome_id %in% ssri_users, "SSRI User", "Non-User")
  )

all_days <- seq.Date(as.Date("2022-10-13"), as.Date("2022-12-16"), by = "day")
ssridf$study_day <- match(as.Date(ssridf$logDate), all_days) -1

ssridf <- ssridf %>%
  mutate(day_c = scale(study_day, center = TRUE, scale = FALSE))

ggplot(ssridf, aes(x = logDate, y = lacto_rel_abundance_vag, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Vaginal Lactobacillus Abundance Over Time",
    x = "Date",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()

ggplot(ssridf, aes(x = logDate, y = lacto_rel_abundance_gut, color = SSRI_status)) +
  geom_jitter(width = 0.5, height = 0, alpha = 0.6, size = 2) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(
    title = "Gut Lactobacillus Abundance Over Time",
    x = "Date",
    y = "Relative Abundance",
    color = "SSRI Use"
  ) +
  theme_minimal()

################################################################
ggplot(ssridf, aes(x = SSRI_status, y = Liners_rel_abundance_vag, fill = SSRI_status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # hide default outliers
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    title = "L.iners Relative Abundance by SSRI Use",
    x = "SSRI Status",
    y = "Abundance"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(ssridf, aes(x = SSRI_status, y = Liners_rel_abundance_vag, color = SSRI_status)) +
  geom_sina(alpha = 0.7, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +  # show median line
  labs(
    title = "L.iners Relative Abundance by SSRI Use",
    x = "SSRI Status",
    y = "Abundance"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ssri_liners <- lmer(Liners_rel_abundance_vag ~ SSRI_status + day_c + I(day_c^2) + (1| biome_id), data = ssridf)
summary(ssri_liners)














