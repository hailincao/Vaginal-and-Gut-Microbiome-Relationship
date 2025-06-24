library(phyloseq)
library(decontam)
library(tidyverse)
library(Matrix)
library(readxl) 

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
######################
#vaginal data

#adding lactobacillus rel abundance column to Vag
tax_df <- data.frame(tax_table(VagData))
Lacto_asvs <- rownames(tax_df)[tax_df$Genus == "Lactobacillus"] 
Lacto_abund <- colSums(otu_table(VagData)[Lacto_asvs, , drop = FALSE])
total_reads_V <- colSums(otu_table(VagData))
Lacto_rel_abund <- Lacto_abund / total_reads_V

# Add relative abundance to sample_data
sample_data(VagData)$Lacto_rel_abundance_vag <- Lacto_rel_abund

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


########################################################################
cross.df_CST <- cross.df %>% filter(CST %in% c("I", "II", "III", "IV", "V"))
#659 has CSTI, 70 has CSTII, 191 has CSTIII, 54 has CSTIV, 30 has CST V

########################################################################
#finding correlation of gut and vaginal microbiome 

#Lactobacillus
#CST coloring
ggplot(cross.df_CST, aes(x = sqrt(lacto_rel_abundance_gut), 
                     y = sqrt(Lacto_rel_abundance_vag), 
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
                     y = sqrt(Lacto_rel_abundance_vag))) +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 1.2) +
  labs(
    title = "Overall Relationship Between Gut and Vaginal Lactobacillus Abundance",
    x = "Gut Lactobacillus Relative Abundance (sqrt)",
    y = "Vaginal Lactobacillus Relative Abundance (sqrt)"
  ) +
  theme_minimal()

cor.test(cross.df$lacto_rel_abundance_gut, 
         cross.df$Lacto_rel_abundance_vag, 
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
    title = "Relative Abundance of B. vulgatus by Vaginal CST",
    y = "Relative Abundance",
    x = "Vaginal CST"
  ) +
  theme_minimal()

# #read in the ones that Alice merged with other lifestyle factor
# gut.data <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/gut.lifestyle.merged.csv")
# vaginal.data <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/vaginal.lifestyle.csv")
# 
# #indicating which are vaginal which are gut
# vaginal.data <- vaginal.data %>% 
#   select(-X) %>% 
#   rename(vaginal_shannon = shannon,
#          # vaginal_max_taxa = max_taxa,
#          vaginal_OTU = OTU,
#          vaginal_sampleID = SampleID
#   )
# nrow(vaginal.data) #1415 samples
# 
# gut.data <- gut.data %>% 
#   select(-c(X)) %>% 
#   rename(gut_shannon=shannon,
#          # gut_max_taxa = max_taxa,
#          gut_OTU = OTU,
#          gut_sampleID = SampleID)
# nrow(sample_data(gut.data)) #1249 samples
# sum(!is.na(gut.data$gut_shannon)) #none of the rows in gut_shannon are NA
# 
# cross.df <- left_join(gut.data, vaginal.data, by = c("biome_id", "logDate"))

# #manually checking if there are any duplicating columns
# names(cross.df)
# all(cross.df$cholesterol_prop.x == cross.df$cholesterol_prop.y, na.rm = TRUE)
# cross.df <- cross.df %>% select(-cholesterol_prop.y)
# all(cross.df$caloriesall_avg.x == cross.df$caloriesall_avg.y, na.rm = TRUE)
# which(cross.df$caloriesall_avg.x != cross.df$caloriesall_avg.y)
# cross.df %>%
#   filter(caloriesall_avg.x != caloriesall_avg.y)
# 
# #View(cross.df)

# samp <- data.frame(sample_data(VagData))
# sum(samp$SampleID %in% cross.df$vaginal_sampleID)
# # merge(samp, cross.df, by.x="SampleID", by.y="vaginal_sampleID")
# merged_df <- samp %>% 
#   left_join(as.data.frame(cross.df), by = c("SampleID" = "vaginal_sampleID"))
# colnames(merged_df)
# 
# #merged cross.df(what Alice have cleaned) back to the VagData file
# merged_df <- as.data.frame(merged_df)
# rownames(merged_df) <- merged_df$SampleID
# new_sample_data <- phyloseq::sample_data(merged_df)
# sample_data(VagData) <- new_sample_data
# head(sample_data(VagData))

#how lacto fluctuate in vagina

# samp <- data.frame(sample_data(VagData))
# matched_samples <- samp$SampleID[samp$SampleID %in% cross.df$SampleID_vaginal]
# length(matched_samples)
# VagData <- prune_samples(matched_samples, VagData)
# merged_df <- data.frame(sample_data(VagData)) %>%
#   left_join(data.frame(cross.df), 
#             by = c("SampleID" = "SampleID_vaginal"))
# rownames(merged_df) <- merged_df$SampleID
# sample_data(VagData) <- sample_data(merged_df)
# head(sample_data(VagData))
# nsamples(VagData)

#removing duplicating columns
# samp_df <- as.data.frame(sample_data(VagData))
# duplicate_cols <- sapply(samp_df, function(x) any(duplicated(as.list(samp_df))))
# samp_df[, duplicate_cols]
# 
# dup_content <- duplicated(as.list(samp_df))
# samp_df_clean <- samp_df[, !dup_content]
# 
# #putting the cleaned columns back to the phyloseq
# rownames(samp_df_clean) <- samp_df_clean$SampleID
# sample_data(VagData) <- sample_data(samp_df_clean)
# 
# #verifying
# colnames(sample_data(VagData))
# 
# #cleaning the name
# samp_df <- as.data.frame(sample_data(VagData))
# colnames(samp_df) <- gsub("\\.x$", "", colnames(samp_df)) 
# colnames(samp_df)
# rownames(samp_df) <- samp_df$SampleID
# sample_data(VagData) <- sample_data(samp_df)
# colnames(sample_data(VagData))

#relative abundance
# vaginal_phyloseq_rel <- transform_sample_counts(VagData, function(x) x / sum(x))


#add Lactobacillus rel. abundance to sample data
# Lacto_phy <- subset_taxa(vaginal_phyloseq_rel, Genus == "Lactobacillus")
# Lacto_abund <- rowSums( otu_table(Lacto_phy)[ , , drop = FALSE ] )
# sample_data(vaginal_phyloseq_rel)$Lacto_abundance_Vag <- Lacto_abund[ sample_names(vaginal_phyloseq_rel) ]
# 
# colnames(sample_data(vaginal_phyloseq_rel))



# #how Prevotella abundance fluctuate over time in each CST
# meta_df <- data.frame(sample_data(vaginal_phyloseq_rel)) %>%
#   select(SampleID, timestamp, CST, Prev_abundance_Vag)
# 
# meta_df$timestamp <- as.Date(meta_df$timestamp)
# 
# meta_df <- meta_df %>%
#   filter(!is.na(CST), !is.na(Prev_abundance_Vag), !is.na(timestamp))
# 
# ggplot(meta_df, aes(x = timestamp, y = Prev_abundance_Vag, color = CST)) +
#   geom_point(alpha = 0.6) +  # Scatter plot
#   geom_smooth(method = "loess", se = FALSE) +  # Trend line
#   labs(
#     title = "Prevotella Abundance (Vaginal) Over Time by CST",
#     x = "Time",
#     y = "Prevotella Relative Abundance"
#   ) +
#   theme_minimal() +
#   facet_wrap(~CST, scales = "free_y")  # Separate plots per CST


# #exploring gut microbiome in each vaginal CST
# CST_Gut <- cross.df %>%
#   group_by(CST, gut_OTU) %>%
#   summarise(count = n(), .groups = "drop") %>%
#   group_by(CST) %>%
#   mutate(total_samples = sum(count),
#          percentage = count / total_samples * 100) %>%
#   arrange(CST, desc(count))
# 
# topspecies <- CST_Gut %>%
#   group_by(CST) %>%
#   slice_max(count, n = 5)
# 
# print(topspecies) #the top 5 abundant species
# 
# #plotting them on a bar plot
# ggplot(topspecies, 
#        aes(x = reorder(gut_OTU, -percentage), y = percentage, fill = CST)) +
#   geom_col() +
#   facet_wrap(~ CST, scales = "free_x") +  
#   labs(x = "Top Gut Species", y = "Percentage per CST", 
#        title = "Top 5 Gut Species by Vaginal CST") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none")
# 
# top1species <- CST_Gut %>%
#   group_by(CST) %>%
#   slice_max(count, n = 1)
# 
# print(top1species)
# 
# ggplot(top1species, aes(x = reorder(CST, -count), y = count, fill = gut_OTU)) +
#   geom_col() +
#   geom_text(aes(label = gut_OTU), 
#             position = position_stack(vjust = 0.5), 
#             size = 3, color = "white") +
#   labs(x = "Vaginal CST", y = "Number of Samples", 
#        title = "Most Abundant Gut Species by Vaginal CST") +
#   theme_minimal() +
#   theme(legend.position = "none") 
# 
# table(cross.df$CST)
# cst_totals <- c(I = 220, II = 27, III = 61, IV = 23, V = 13)
# 
# top1species <- top1species %>%
#   mutate(total_samples = cst_totals[CST])
# 
# ggplot(top1species, aes(x = reorder(CST, -count), y = count, fill = gut_OTU)) +
#   geom_col() +
#   geom_text(aes(y = 0, label = paste0("n=", total_samples)),
#             vjust = 1.2, color = "black", size = 3.5) +
#   geom_text(aes(label = gut_OTU),
#             position = position_stack(vjust = 0.5),
#             color = "white", size = 3.5) +
#   geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
#   labs(x = "Vaginal CST", y = "Count of Top Gut Species",
#        title = "Most Abundant Gut Species by Vaginal CST") +
#   theme_minimal() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

# #how lacto fluctuate in gut
# samp2 <- data.frame(sample_data(GutData))
# matched_samples <- samp2$SampleID[samp2$SampleID %in% cross.df$gut_sampleID]
# length(matched_samples)
# GutData <- prune_samples(matched_samples, GutData)
# otu_table(GutData) <- otu_table(t(otu_table(GutData)))
# merged_df <- data.frame(sample_data(GutData)) %>%
#   left_join(data.frame(cross.df), 
#             by = c("SampleID" = "gut_sampleID"))
# rownames(merged_df) <- merged_df$SampleID
# sample_data(GutData) <- sample_data(merged_df)
# colnames(sample_data(GutData))
# nsamples(GutData)
# #View(sample_data(GutData))
# 
# #removing duplicating columns
# samp_df <- data.frame(sample_data(GutData))
# colnames(samp_df)
# samp_df_clean <- samp_df %>% 
#   select(SampleID, is_blank, qr, biome_id.x, logDate.x, timestamp, sampleType, CST, vaginal_shannon)
# rownames(samp_df_clean) <- samp_df_clean$SampleID
# sample_data(GutData) <- sample_data(samp_df_clean)
# 
# #View(sample_data(GutData))
# colnames(sample_data(GutData)) 
# 
# samp_df <- data.frame(sample_data(GutData))
# colnames(samp_df) <- gsub("\\.x$|\\.y$", "", colnames(samp_df))
# #samp_df <- samp_df[, !duplicated(colnames(samp_df))]
# rownames(samp_df) <- samp_df$SampleID
# sample_data(GutData) <- sample_data(samp_df)
# colnames(sample_data(GutData))
# 
# 
# #View(otu_table(GutData))
# 
# #relative abundance
# gut_phyloseq_rel <- transform_sample_counts(GutData, function(x) x / sum(x))
# 
# #add Lactobacillus rel. abundance to sample data
# Lacto_phy <- subset_taxa(gut_phyloseq_rel, Genus == "Lactobacillus")
# Lacto_abund <- rowSums( otu_table(Lacto_phy)[ , , drop = FALSE ] )
# summary(Lacto_abund)
# sample_data(gut_phyloseq_rel)$Lacto_abundance_Gut <- Lacto_abund[ sample_names(gut_phyloseq_rel) ]
# 
# #add Prevotella rel.abundance to sample data
# Prev_phy <- subset_taxa(gut_phyloseq_rel,
#                         Genus %in% c("Prevotella", "Prevotella_7", "Prevotella_9"))
# Prev_abund <- rowSums( otu_table(Prev_phy)[ , , drop = FALSE ] )
# summary(Prev_abund)
# sample_data(gut_phyloseq_rel)$Prev_abundance_Gut <- Prev_abund[ sample_names(gut_phyloseq_rel) ]
# 
# 
# #View((sample_data(gut_phyloseq_rel)))
# colnames(sample_data(gut_phyloseq_rel))
# 
# #how Lacto abundance fluctuate over time in each CST
# meta_df <- data.frame(sample_data(gut_phyloseq_rel)) %>%
#   select(SampleID, timestamp, CST, Lacto_abundance_Gut)
# 
# meta_df$timestamp <- as.Date(meta_df$timestamp)
# 
# meta_df <- meta_df %>%
#   filter(!is.na(CST), !is.na(Lacto_abundance_Gut), !is.na(timestamp))


# ggplot(meta_df, aes(x = timestamp, y = Lacto_abundance_Gut, color = CST)) +
#   geom_point(alpha = 0.6) +  # Scatter plot
#   geom_smooth(method = "loess", se = FALSE) +  # Trend line
#   labs(
#     title = "Lactobacillus Abundance (Gut) Over Time by CST",
#     x = "Time",
#     y = "Lacto_abundance_Gut"
#   ) +
#   theme_minimal() +
#   facet_wrap(~CST, scales = "free_y")  # Separate plots per CST
# 
# #how prevotella abundance fluctuate over time in each CST
# meta_df <- data.frame(sample_data(gut_phyloseq_rel)) %>%
#   select(SampleID, timestamp, CST, Prev_abundance_Gut)
# 
# meta_df$timestamp <- as.Date(meta_df$timestamp)
# 
# meta_df <- meta_df %>%
#   filter(!is.na(CST), !is.na(Prev_abundance_Gut), !is.na(timestamp))
# 
# 
# ggplot(meta_df, aes(x = timestamp, y = Prev_abundance_Gut, color = CST)) +
#   geom_point(alpha = 0.6) +  # Scatter plot
#   geom_smooth(method = "loess", se = FALSE) +  # Trend line
#   labs(
#     title = "Prevotella Abundance (Gut) Over Time by CST",
#     x = "Time",
#     y = "Prevotella Relative Abundance"
#   ) +
#   theme_minimal() +
#   facet_wrap(~CST, scales = "free_y")  # Separate plots per CST
# 
# 
# #how gut and vaginal lacto correlate with each other in CST 
# vag_meta <- vaginal_phyloseq_rel %>%
#   sample_data() %>%       
#   data.frame() %>%      
#   select(
#     SampleID, 
#     qr,
#     biome_id,
#     CST_vag = CST, 
#     Lacto_abundance_Vag, 
#     logDate
#   )
# 
# 
# gut_meta <- gut_phyloseq_rel %>%
#   sample_data() %>%       
#   data.frame() %>%      
#   select(
#     SampleID, 
#     qr,
#     biome_id,
#     CST_gut = CST, 
#     Lacto_abundance_Gut, 
#     logDate
#   )
# 
# colnames(gut_meta)
# colnames(vag_meta)
# 
# merged_df <- left_join(gut_meta, 
#                        vag_meta %>% select(biome_id, qr, logDate, CST_vag, Lacto_abundance_Vag), 
#                        by = c("biome_id","logDate"))
# 
# ggplot(merged_df, aes(x = Lacto_abundance_Vag, y = Lacto_abundance_Gut, color = CST_vag)) +
#   geom_point(alpha = 0.7) +
#   geom_smooth(method = "loess", se = FALSE, span = 0.75, linewidth = 0.8) +
#   facet_wrap(~CST_vag) +
#   labs(
#     x = "Vaginal Lactobacillus Relative Abundance (%)",
#     y = "Gut Lactobacillus Relative Abundance (%)",
#     title = "Vaginal-Gut Lactobacillus Correlation by CST"
#   ) +
#   theme_bw()
# 
# #how Prevotella correlate with each other
# vag_meta2 <- vaginal_phyloseq_rel %>%
#   sample_data() %>%       
#   data.frame() %>%      
#   select(
#     SampleID, 
#     qr,
#     biome_id,
#     CST_vag = CST, 
#     Prev_abundance_Vag, 
#     logDate
#   )
# 
# 
# gut_meta2 <- gut_phyloseq_rel %>%
#   sample_data() %>%       
#   data.frame() %>%      
#   select(
#     SampleID, 
#     qr,
#     biome_id,
#     CST_gut = CST, 
#     Prev_abundance_Gut, 
#     logDate
#   )
# 
# colnames(gut_meta2)
# colnames(vag_meta2)
# 
# merged_df2 <- left_join(gut_meta2, 
#                        vag_meta2 %>% select(biome_id, qr, logDate, CST_vag, Prev_abundance_Vag), 
#                        by = c("biome_id","logDate"))
# 
# ggplot(merged_df2, aes(x = Prev_abundance_Vag, y = Prev_abundance_Gut, color = CST_vag)) +
#   geom_point(alpha = 0.7) +
#   geom_smooth(method = "loess", se = FALSE, span = 0.75, linewidth = 0.8)  +
#   facet_wrap(~CST_vag) +
#   labs(
#     x = "Vaginal Prevotella Relative Abundance (%)",
#     y = "Gut Prevotella Relative Abundance (%)",
#     title = "Vaginal-Gut Prevotella Correlation by CST"
#   ) +
#   theme_bw()








#View(cross.df)

##################################################
#creating a scatterplot of shannon diversity correlation between vaginal and gut per person
plot(x = cross.df$vaginal_shannon, y = cross.df$gut_shannon,
     main = "Average Shannnon Diversity Correlation",
     xlab = "Vaginal Microbiome",
     ylab = "Gut Microbiome",
     pch = 16)           

spline_fit <- smooth.spline(cross.df$vaginal_shannon, cross.df$gut_shannon)
lines(spline_fit, col = "blue", lwd = 2)

###################################################
#replicating Alice's analysis
vag.bacterial.data <- readRDS("/Users/caoyang/Desktop/Tetel Lab/datasets/vaginal_cleaned_max_taxa.rds")
gut.bacterial.data <- readRDS("/Users/caoyang/Desktop/Tetel Lab/datasets/fecal_cleaned_max_taxa.rds")

vag_ids_to_keep <- cross.df$vaginal_sampleID
gut_ids_to_keep <- cross.df$gut_sampleID

# Subset the phyloseq dfs
vag.bacterial.subset <- prune_samples(vag_ids_to_keep, vag.bacterial.data)
gut.bacterial.subset <- prune_samples(gut_ids_to_keep, gut.bacterial.data)

# metadata df
vag.meta <- sample_data(vag.bacterial.subset)
gut.meta <- sample_data(gut.bacterial.subset)

# merge bacterial data of sites
merge.site.bacterial <- merge_phyloseq(vag.bacterial.subset, gut.bacterial.subset)

#View(merge.site.bacterial)
head(tax_table(merge.site.bacterial))
colnames(tax_table(merge.site.bacterial))
colnames(sample_data(merge.site.bacterial))

# bray curtis distance of pairs (beta diversity)
bray_dist <- phyloseq::distance(merge.site.bacterial, method = "bray")

class(bray_dist)

bray_df <- as.matrix(bray_dist) %>%
  as.data.frame() %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "bray") %>%
  filter(Sample1 < Sample2) 

# join df
meta <- as(sample_data(merge.site.bacterial), "data.frame") %>% 
  select(biome_id, SampleID, logDate, sampleType)
bray_meta <- bray_df %>%
  left_join(meta, by = c("Sample1" = "SampleID")) %>%
  rename(biome_id_1 = biome_id, site1 = sampleType, logDate_1 = logDate) %>%
  left_join(meta, by = c("Sample2" = "SampleID")) %>%
  rename(biome_id_2 = biome_id, site2 = sampleType, logDate_2 = logDate)

# Filter to within-person, cross-site comparisons
bray_cross_site <- bray_meta %>%
  filter(biome_id_1 == biome_id_2,
         site1 != site2,
         logDate_1 == logDate_2)
# Join meta data
bray_cross_site.full <- bray_cross_site %>% 
  left_join(cross.df, by = c("logDate_1" = "logDate", "biome_id_1" = "biome_id"))
names(bray_cross_site.full)

colnames(bray_cross_site.full)

# gut vs vaginal Shannon
ggplot(bray_cross_site.full, aes(x = vaginal_shannon, y = gut_shannon, col=as.factor(biome_id_1))) +
  geom_point(alpha = 0.6) +
  # geom_smooth(method = "loess", se = FALSE) +
  geom_smooth(method = "loess", se = FALSE, color="blue") +
  labs(
    x = "Vaginal Shannon Diversity",
    y = "Gut Shannon Diversity",
    title = " "
  ) +
  theme_minimal() +
  ylim(0,5)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0),
        text=element_text(size=16),
        legend.position="none")

cor.test(bray_cross_site.full$vaginal_shannon, 
         bray_cross_site.full$gut_shannon, use = "complete.obs")

## Bray

# Beta Diversity Between Gut and Vaginal Microbiomes pairs
ggplot(bray_cross_site.full, aes(x = as.factor(biome_id_1), y = bray)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "steelblue") +
  labs(
    x = "Participant ID",
    y = "Bray-Curtis Dissimilarity (Gut vs Vaginal)",
    title = " "
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 14))

#View(bray_cross_site.full)
dim(bray_cross_site.full)

# Cross-site Bray-Curtis by CST
ggplot(bray_cross_site.full, aes(x = CST, y = bray, fill = CST)) +
  geom_boxplot() +
  geom_jitter(color="steelblue", width = 0.2, alpha = 0.5) +
  labs(
    x = "CST",
    y = "Bray-Curtis Dissimilarity",
    title = ""
  ) +
  theme_minimal() +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 14)) +
  guides(fill = guide_legend(title = "CST"))

## Edit df
bray_cross_site.full.filtered <- bray_cross_site.full %>% 
  select(-c(Sample1, Sample2, logDate_2, biome_id_2, site1, site2, vaginal_sampleID, gut_sampleID)) %>% 
  rename(logDate = logDate_1,
         biome_id = biome_id_1)
colSums(is.na(bray_cross_site.full.filtered))

bray_cross_site.full.filtered <- bray_cross_site.full.filtered %>% 
  mutate(stress_severity = case_when(
    stress_score <= 14 ~ "Normal",
    stress_score >= 15 & stress_score <= 18 ~ "Mild",
    stress_score >= 19 & stress_score <= 25 ~ "Moderate",
    stress_score >= 26 & stress_score <= 33 ~ "Severe",
    stress_score >= 34 ~ "Extremely Severe",
    TRUE ~ NA_character_
  )) %>% 
  select(-stress_score)
colSums(is.na(bray_cross_site.full.filtered))

# filter all NA cols and all NA rows
bray_cross_site.full.filtered <- bray_cross_site.full.filtered %>%
  select(where(~ !all(is.na(.)))) %>%
  filter(if_any(everything(), ~ !is.na(.))) 
dim(bray_cross_site.full.filtered)
colSums(is.na(bray_cross_site.full.filtered))

################################################################
#this part I did not understand what gut vs. vaginal bray between-person is

# bray_combined <- bind_rows(
#   bray_cross_site %>% mutate(comparison = "within-person"),
#   bray_meta %>%
#     filter(biome_id_1 != biome_id_2, site1 != site2, logDate_1 == logDate_2) %>%
#     mutate(comparison = "between-person")
# )
# 
# # Plot with a side-by-side boxplot
# ggplot(bray_combined, aes(x = comparison, y = bray, fill = comparison)) +
#   geom_boxplot(alpha = 0.7) +
#   # geom_jitter(color="orchid", alpha=0.3) +
#   labs(x = "", y = "Bray-Curtis Dissimilarity", title = "Cross-site Dissimilarity (Gut vs Vaginal)") +
#   scale_fill_manual(values = c("within-person" = "#1f77b4", "between-person" = "#ff7f0e")) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
#         text=element_text(size=16))
# 
# #violin plot
# ggplot(bray_combined, aes(x = comparison, y = bray, fill = comparison)) +
#   geom_violin(alpha = 0.7) +
#   geom_boxplot(width = 0.1, alpha = 0.3) + 
#   labs(x = "", y = "Bray-Curtis Dissimilarity", title = "Cross-site Dissimilarity (Gut vs Vaginal)") +
#   scale_fill_manual(values = c("within-person" = "#1f77b4", "between-person" = "#ff7f0e")) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
#         text=element_text(size=16))

################################################################






