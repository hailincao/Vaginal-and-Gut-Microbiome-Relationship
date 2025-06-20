library(phyloseq)
library(decontam)
library(tidyverse)
library(Matrix)
library(readxl) 

# #replicating Alice's code

#importing cleaned bacteria data
VagData <- readRDS("/Users/caoyang/Desktop/Tetel Lab/datasets/vaginal_bacteria_filter2_cleanedv3.rds")
GutData <- readRDS("/Users/caoyang/Desktop/Tetel Lab/datasets/fecal_bacteria_cleanedv3.rds")



# bacterial.data <- readRDS("/Users/caoyang/Desktop/Tetel Lab/datasets/vaginal_cleaned_max_taxa.rds")
# head(otu_table(bacterial.data))


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

#View(gut.data)

cross.df <- vaginal.data %>% 
  left_join(gut.data) %>% 
  filter(!is.na(gut_shannon))


#View(cross.df)

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

samp <- data.frame(sample_data(VagData))
matched_samples <- samp$SampleID[samp$SampleID %in% cross.df$vaginal_sampleID]
length(matched_samples)
VagData <- prune_samples(matched_samples, VagData)
merged_df <- data.frame(sample_data(VagData)) %>%
  left_join(data.frame(cross.df), 
            by = c("SampleID" = "vaginal_sampleID"))
rownames(merged_df) <- merged_df$SampleID
sample_data(VagData) <- sample_data(merged_df)
head(sample_data(VagData))
nsamples(VagData)

#removing duplicating columns
samp_df <- as.data.frame(sample_data(VagData))
duplicate_cols <- sapply(samp_df, function(x) any(duplicated(as.list(samp_df))))
samp_df[, duplicate_cols]

dup_content <- duplicated(as.list(samp_df))
samp_df_clean <- samp_df[, !dup_content]

#putting the cleaned columns back to the phyloseq
rownames(samp_df_clean) <- samp_df_clean$SampleID
sample_data(VagData) <- sample_data(samp_df_clean)

#verifying
colnames(sample_data(VagData))

#cleaning the name
samp_df <- as.data.frame(sample_data(VagData))
colnames(samp_df) <- gsub("\\.x$", "", colnames(samp_df)) 
colnames(samp_df)
rownames(samp_df) <- samp_df$SampleID
sample_data(VagData) <- sample_data(samp_df)
colnames(sample_data(VagData))

#relative abundance
vaginal_phyloseq_rel <- transform_sample_counts(VagData, function(x) x / sum(x))


#add Lactobacillus rel. abundance to sample data
Lacto_phy <- subset_taxa(vaginal_phyloseq_rel, Genus == "Lactobacillus")
Lacto_abund <- rowSums( otu_table(Lacto_phy)[ , , drop = FALSE ] )
sample_data(vaginal_phyloseq_rel)$Lacto_abundance_Vag <- Lacto_abund[ sample_names(vaginal_phyloseq_rel) ]

colnames(sample_data(vaginal_phyloseq_rel))

#add Prevotella to sample data
Prev_phy <- subset_taxa(vaginal_phyloseq_rel,
                        Genus %in% c("Prevotella", "Prevotella_7", "Prevotella_9"))
Prev_abund <- rowSums( otu_table(Prev_phy)[ , , drop = FALSE ] )
sample_data(vaginal_phyloseq_rel)$Prev_abundance_Vag <- Prev_abund[ sample_names(vaginal_phyloseq_rel) ]


#how Prevotella abundance fluctuate over time in each CST
meta_df <- data.frame(sample_data(vaginal_phyloseq_rel)) %>%
  select(SampleID, timestamp, CST, Prev_abundance_Vag)

meta_df$timestamp <- as.Date(meta_df$timestamp)

meta_df <- meta_df %>%
  filter(!is.na(CST), !is.na(Prev_abundance_Vag), !is.na(timestamp))

ggplot(meta_df, aes(x = timestamp, y = Prev_abundance_Vag, color = CST)) +
  geom_point(alpha = 0.6) +  # Scatter plot
  geom_smooth(method = "loess", se = FALSE) +  # Trend line
  labs(
    title = "Prevotella Abundance (Vaginal) Over Time by CST",
    x = "Time",
    y = "Prevotella Relative Abundance"
  ) +
  theme_minimal() +
  facet_wrap(~CST, scales = "free_y")  # Separate plots per CST


#exploring gut microbiome in each vaginal CST
CST_Gut <- cross.df %>%
  group_by(CST, gut_OTU) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(CST) %>%
  mutate(total_samples = sum(count),
         percentage = count / total_samples * 100) %>%
  arrange(CST, desc(count))

topspecies <- CST_Gut %>%
  group_by(CST) %>%
  slice_max(count, n = 5)

print(topspecies) #the top 5 abundant species

#plotting them on a bar plot
ggplot(topspecies, 
       aes(x = reorder(gut_OTU, -percentage), y = percentage, fill = CST)) +
  geom_col() +
  facet_wrap(~ CST, scales = "free_x") +  
  labs(x = "Top Gut Species", y = "Percentage per CST", 
       title = "Top 5 Gut Species by Vaginal CST") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

top1species <- CST_Gut %>%
  group_by(CST) %>%
  slice_max(count, n = 1)

print(top1species)

ggplot(top1species, aes(x = reorder(CST, -count), y = count, fill = gut_OTU)) +
  geom_col() +
  geom_text(aes(label = gut_OTU), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "white") +
  labs(x = "Vaginal CST", y = "Number of Samples", 
       title = "Most Abundant Gut Species by Vaginal CST") +
  theme_minimal() +
  theme(legend.position = "none") 

table(cross.df$CST)
cst_totals <- c(I = 220, II = 27, III = 61, IV = 23, V = 13)

top1species <- top1species %>%
  mutate(total_samples = cst_totals[CST])

ggplot(top1species, aes(x = reorder(CST, -count), y = count, fill = gut_OTU)) +
  geom_col() +
  geom_text(aes(y = 0, label = paste0("n=", total_samples)),
            vjust = 1.2, color = "black", size = 3.5) +
  geom_text(aes(label = gut_OTU),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3.5) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5) +
  labs(x = "Vaginal CST", y = "Count of Top Gut Species",
       title = "Most Abundant Gut Species by Vaginal CST") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

#how lacto fluctuate in gut
samp2 <- data.frame(sample_data(GutData))
matched_samples <- samp2$SampleID[samp2$SampleID %in% cross.df$gut_sampleID]
length(matched_samples)
GutData <- prune_samples(matched_samples, GutData)
otu_table(GutData) <- otu_table(t(otu_table(GutData)))
merged_df <- data.frame(sample_data(GutData)) %>%
  left_join(data.frame(cross.df), 
            by = c("SampleID" = "gut_sampleID"))
rownames(merged_df) <- merged_df$SampleID
sample_data(GutData) <- sample_data(merged_df)
colnames(sample_data(GutData))
nsamples(GutData)
#View(sample_data(GutData))

#removing duplicating columns
samp_df <- data.frame(sample_data(GutData))
colnames(samp_df)
samp_df_clean <- samp_df %>% 
  select(SampleID, is_blank, qr, biome_id.x, logDate.x, timestamp, sampleType, CST, vaginal_shannon)
rownames(samp_df_clean) <- samp_df_clean$SampleID
sample_data(GutData) <- sample_data(samp_df_clean)

#View(sample_data(GutData))
colnames(sample_data(GutData)) 

samp_df <- data.frame(sample_data(GutData))
colnames(samp_df) <- gsub("\\.x$|\\.y$", "", colnames(samp_df))
#samp_df <- samp_df[, !duplicated(colnames(samp_df))]
rownames(samp_df) <- samp_df$SampleID
sample_data(GutData) <- sample_data(samp_df)
colnames(sample_data(GutData))


#View(otu_table(GutData))

#relative abundance
gut_phyloseq_rel <- transform_sample_counts(GutData, function(x) x / sum(x))

#add Lactobacillus rel. abundance to sample data
Lacto_phy <- subset_taxa(gut_phyloseq_rel, Genus == "Lactobacillus")
Lacto_abund <- rowSums( otu_table(Lacto_phy)[ , , drop = FALSE ] )
summary(Lacto_abund)
sample_data(gut_phyloseq_rel)$Lacto_abundance_Gut <- Lacto_abund[ sample_names(gut_phyloseq_rel) ]

#add Prevotella rel.abundance to sample data
Prev_phy <- subset_taxa(gut_phyloseq_rel,
                        Genus %in% c("Prevotella", "Prevotella_7", "Prevotella_9"))
Prev_abund <- rowSums( otu_table(Prev_phy)[ , , drop = FALSE ] )
summary(Prev_abund)
sample_data(gut_phyloseq_rel)$Prev_abundance_Gut <- Prev_abund[ sample_names(gut_phyloseq_rel) ]


#View((sample_data(gut_phyloseq_rel)))
colnames(sample_data(gut_phyloseq_rel))

#how Lacto abundance fluctuate over time in each CST
meta_df <- data.frame(sample_data(gut_phyloseq_rel)) %>%
  select(SampleID, timestamp, CST, Lacto_abundance_Gut)

meta_df$timestamp <- as.Date(meta_df$timestamp)

meta_df <- meta_df %>%
  filter(!is.na(CST), !is.na(Lacto_abundance_Gut), !is.na(timestamp))


ggplot(meta_df, aes(x = timestamp, y = Lacto_abundance_Gut, color = CST)) +
  geom_point(alpha = 0.6) +  # Scatter plot
  geom_smooth(method = "loess", se = FALSE) +  # Trend line
  labs(
    title = "Lactobacillus Abundance (Gut) Over Time by CST",
    x = "Time",
    y = "Lacto_abundance_Gut"
  ) +
  theme_minimal() +
  facet_wrap(~CST, scales = "free_y")  # Separate plots per CST

#how prevotella abundance fluctuate over time in each CST
meta_df <- data.frame(sample_data(gut_phyloseq_rel)) %>%
  select(SampleID, timestamp, CST, Prev_abundance_Gut)

meta_df$timestamp <- as.Date(meta_df$timestamp)

meta_df <- meta_df %>%
  filter(!is.na(CST), !is.na(Prev_abundance_Gut), !is.na(timestamp))


ggplot(meta_df, aes(x = timestamp, y = Prev_abundance_Gut, color = CST)) +
  geom_point(alpha = 0.6) +  # Scatter plot
  geom_smooth(method = "loess", se = FALSE) +  # Trend line
  labs(
    title = "Prevotella Abundance (Gut) Over Time by CST",
    x = "Time",
    y = "Prevotella Relative Abundance"
  ) +
  theme_minimal() +
  facet_wrap(~CST, scales = "free_y")  # Separate plots per CST


#how gut and vaginal lacto correlate with each other in CST 
vag_meta <- vaginal_phyloseq_rel %>%
  sample_data() %>%       
  data.frame() %>%      
  select(
    SampleID, 
    qr,
    biome_id,
    CST_vag = CST, 
    Lacto_abundance_Vag, 
    logDate
  )


gut_meta <- gut_phyloseq_rel %>%
  sample_data() %>%       
  data.frame() %>%      
  select(
    SampleID, 
    qr,
    biome_id,
    CST_gut = CST, 
    Lacto_abundance_Gut, 
    logDate
  )

colnames(gut_meta)
colnames(vag_meta)

merged_df <- left_join(gut_meta, 
                       vag_meta %>% select(biome_id, qr, logDate, CST_vag, Lacto_abundance_Vag), 
                       by = c("biome_id","logDate"))

ggplot(merged_df, aes(x = Lacto_abundance_Vag, y = Lacto_abundance_Gut, color = CST_vag)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~CST_vag) +
  labs(
    x = "Vaginal Lactobacillus Relative Abundance (%)",
    y = "Gut Lactobacillus Relative Abundance (%)",
    title = "Vaginal-Gut Lactobacillus Correlation by CST"
  ) +
  theme_bw()

#how Prevotella correlate with each other
vag_meta2 <- vaginal_phyloseq_rel %>%
  sample_data() %>%       
  data.frame() %>%      
  select(
    SampleID, 
    qr,
    biome_id,
    CST_vag = CST, 
    Prev_abundance_Vag, 
    logDate
  )


gut_meta2 <- gut_phyloseq_rel %>%
  sample_data() %>%       
  data.frame() %>%      
  select(
    SampleID, 
    qr,
    biome_id,
    CST_gut = CST, 
    Prev_abundance_Gut, 
    logDate
  )

colnames(gut_meta2)
colnames(vag_meta2)

merged_df2 <- left_join(gut_meta2, 
                       vag_meta2 %>% select(biome_id, qr, logDate, CST_vag, Prev_abundance_Vag), 
                       by = c("biome_id","logDate"))

ggplot(merged_df2, aes(x = Prev_abundance_Vag, y = Prev_abundance_Gut, color = CST_vag)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~CST_vag) +
  labs(
    x = "Vaginal Prevotella Relative Abundance (%)",
    y = "Gut Prevotella Relative Abundance (%)",
    title = "Vaginal-Gut Prevotella Correlation by CST"
  ) +
  theme_bw()








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






