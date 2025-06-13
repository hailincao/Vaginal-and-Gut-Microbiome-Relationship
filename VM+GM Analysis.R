library(phyloseq)
library(decontam)
library(tidyverse)
library(Matrix)
library(readxl) 

#replicating Alice's code

bacterial.data <- readRDS("/Users/caoyang/Desktop/Tetel Lab/Walther-Antonio_Project_022_ITS2.rds")
samples.data <- read_excel("/Users/caoyang/Desktop/Tetel Lab/cleaned_samplesv2.xlsx") #changing original code's sample data into version 2
uminn_data <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/cleaned_uminn_data.csv", header=TRUE)                   

bacterial_otu_table <- otu_table(bacterial.data)
bacterial_tax_table <- tax_table(bacterial.data)

dim(samples.data)
length(unique(samples.data$biome_id))
dim(uminn_data)

## Metadata
uminn_data <- uminn_data %>% 
  select(Sample.ID, Special.Notes) %>% 
  mutate(qr=sub("_.*", "", Sample.ID))

metadata_bacteria <- data.frame(
  SampleID = sample_names(bacterial_otu_table),
  is_blank = grepl("BLANK", colnames(bacterial_otu_table)),
  stringsAsFactors = FALSE) %>%
  mutate(qr = sub("\\_.*", "", SampleID))

metadata_bacteria <- metadata_bacteria %>% 
  left_join(samples.data, by="qr")
rownames(metadata_bacteria) <- metadata_bacteria$SampleID

dupes <- names(table(metadata_bacteria$SampleID))[table(metadata_bacteria$SampleID) > 1]
# [1] "S1519_2_V3V5_S1487"       "S1519_V3V5_S1199"         "S1587_2_V3V5_S1535"      
# [4] "S1587_V3V5_S1247"         "S2558_2_V3V5_S1499"       "S2558_V3V5_S1211"        
# [7] "S3244_2_V3V5_S1475"       "S3244_V3V5_S1187"         "S3791_2_V3V5_S1511"      
# [10] "S3791_V3V5_S1223"         "S3977_2_V3V5_S1523"       "S3977_V3V5_S1235"        
# [13] "S4495_loose_2_V3V5_S1547" "S4495_loose_V3V5_S1259" 

# Set sample data
sample_data_obj <- sample_data(metadata_bacteria)

## Saving as new obj
bacteria_physeq <- phyloseq(bacterial_otu_table, sample_data_obj, bacterial_tax_table)
table(metadata_bacteria$is_blank)

bacteria_physeq_otu <- otu_table(bacteria_physeq)

## Filter out error data - filter samples
uminn_data_keep <- uminn_data %>% 
  filter(str_detect(Special.Notes, "error") | str_detect(Special.Notes, "No swab in tube")) %>% 
  select(Sample.ID, Special.Notes, qr)
otu_table_bacteria <- as.data.frame(t(bacteria_physeq_otu))
dim(otu_table_bacteria) # 3659 122403

sample_ids <- rownames(otu_table_bacteria)
sample_ids <- sub("_.*", "", sample_ids)
length(sample_ids) # 3659

# Prune samples that had errors
samples_to_keep <- !sample_ids %in% uminn_data_keep$qr

t_bacterial_otu_table <- t(bacterial_otu_table)
physeq_no_error <- prune_samples(samples_to_keep, t_bacterial_otu_table)

# Prep for decontam
physeq_no_error_otu <- otu_table(physeq_no_error)
physeq_no_error_otu_df <- as.data.frame(physeq_no_error_otu)

# Save new obj
# save.image("/Volumes/T7/microbiome_data/R_environments/microbiome_cleaningv1.RData")
# load("/Volumes/T7/microbiome_data/R_environments/microbiome_cleaningv1.RData")

# Sample data
metadata_bacteria <- data.frame(
  SampleID = rownames(physeq_no_error_otu_df),
  is_blank = grepl("BLANK", rownames(physeq_no_error_otu_df)),
  stringsAsFactors = FALSE
) %>%
  mutate(qr = sub("\\_.*", "", SampleID))
metadata_bacteria <- metadata_bacteria %>% 
  left_join(samples.data, by="qr")
rownames(metadata_bacteria) <- metadata_bacteria$SampleID

# Convert back to phyloseq obj
otu_table_obj <- otu_table(physeq_no_error_otu_df, taxa_are_rows = FALSE)
sample_data_obj <- sample_data(metadata_bacteria)

## Saving as new obj
bacteria_physeq <- phyloseq(otu_table_obj, sample_data_obj, bacterial_tax_table)

# Save new obj
saveRDS(bacteria_physeq, file = "/Users/caoyang/Desktop/Tetel Lab/datasets/bacteria_intermediary2.rds")

gut.data <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/gut.lifestyle.merged.csv")
vaginal.data <- read.csv("/Users/caoyang/Desktop/Tetel Lab/datasets/vaginal.lifestyle.csv")

vaginal.data <- vaginal.data %>% 
  select(-X) %>% 
  rename(vaginal_shannon = shannon,
         # vaginal_max_taxa = max_taxa,
         vaginal_OTU = OTU,
         vaginal_sampleID = SampleID
  )


gut.data <- gut.data %>% 
  select(-c(X)) %>% 
  rename(gut_shannon=shannon,
         # gut_max_taxa = max_taxa,
         gut_OTU = OTU,
         gut_sampleID = SampleID)

View(gut.data)

cross.df <- vaginal.data %>% 
  left_join(gut.data) %>% 
  filter(!is.na(gut_shannon))

names(cross.df)
dim(cross.df)

View(cross.df)

##################################################
#creating a scatterplot of shannon diversity correlation between vaginal and gut per person
plot(x = cross.df$vaginal_shannon, y = cross.df$gut_shannon,
     main = "Average Shannnon Diversity Correlation",
     xlab = "Vaginal Microbiome",
     ylab = "Gut Microbiome",
     pch = 16)           

?abline

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

View(bray_cross_site.full)
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






