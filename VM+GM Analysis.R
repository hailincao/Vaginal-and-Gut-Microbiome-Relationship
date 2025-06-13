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

View(vaginal.data)


gut.data <- gut.data %>% 
  select(-c(X)) %>% 
  rename(gut_shannon=shannon,
         # gut_max_taxa = max_taxa,
         gut_OTU = OTU,
         gut_sampleID = SampleID)

cross.df <- vaginal.data %>% 
  left_join(gut.data) %>% 
  filter(!is.na(gut_shannon))

names(cross.df)
dim(cross.df)














