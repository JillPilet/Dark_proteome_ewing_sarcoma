# The goal of this script is to merge all filtered isoforms into one

### 0. Load libraries-----------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(ggridges)
library(scales)
library(ggrepel)
library(tidyverse)
library(biomaRt)
library(openxlsx)
library(tidyr)

# Set directories  
data_dir = ""

### 1. Load data and annotate --------------------------------------------------
  # 1.1 Load all filtered files ------------------------------------------------
A673 <- read_xlsx(file.path(data_dir, "A673_filtered_all_filters.xlsx"))
TC71 <- read_xlsx(file.path(data_dir, "TC71_filtered_all_filters.xlsx"))
EW7 <- read_xlsx(file.path(data_dir, "EW7_filtered_all_filters.xlsx"))
EW16 <- read_xlsx(file.path(data_dir, "EW16_filtered_all_filters.xlsx"))
Ewima1 <- read_xlsx(file.path(data_dir, "Ewima1_filtered_all_filters.xlsx"))

  # 1.2 Convert appropriate columns to numeric ---------------------------------

numeric_variables <- c("UGS_Median", 
                       "Median_EWING_Cell_line(33)",
                       "Median_Transcript_EWING_Cell_line(33)",
                       "Median_Transcript_Ewima1(3)",
                       "Median_Transcript_MSC(6)",
                       "Median_Ewima1(3)",
                       "Median_MSC(6)",
                       "A673_log2FC_High_Low",
                       "ASP14_log2FC_High_Low",
                       "TC71_72h_log2FC_High_Low",
                       "TC71_96h_log2FC_High_Low",
                       "Nb_reads_Nanopore",
                       "Ewing_sarcoma_UGS_Median", 
                       "TCGA_Primary_Tumor_Median", 
                       "TCGA_Tissue_Normal_Median",
                       "GTEX_Median",
                       "HPA_Median",
                       "UGS_Max",
                       "Ewing_sarcoma_UGS_Max",
                       "TCGA_Primary_Tumor_Max",
                       "TCGA_Tissue_Normal_Max",
                       "GTEX_Max",
                       "HPA_Max",
                       "Other_Median",
                       "Other_Max",
                       "Delta_tss",
                       "Length_Ref_Transcript",
                       "Length_Transcript...7",
                       "Length_Transcript...64",
                       "Nb_Annot_membrane_3sources")

A673[,numeric_variables] = lapply(A673[,numeric_variables], as.numeric)
TC71[,numeric_variables] = lapply(TC71[,numeric_variables], as.numeric)
EW7[,numeric_variables] = lapply(EW7[,numeric_variables], as.numeric)
EW16[,numeric_variables] = lapply(EW16[,numeric_variables], as.numeric)
Ewima1[,numeric_variables] = lapply(Ewima1[,numeric_variables], as.numeric)

A673 = A673 %>%
  rename(closest_FLI1_GGAAm = A673_FLI1_GGAAm)
TC71 = TC71 %>%
  rename(closest_FLI1_GGAAm = TC71_FLI1_GGAAm)
EW7 = EW7 %>%
  rename(closest_FLI1_GGAAm = EW7_FLI1_GGAAm)
EW16 = EW16 %>%
  rename(closest_FLI1_GGAAm = TC71_FLI1_GGAAm)
Ewima1 = Ewima1 %>%
  rename(closest_FLI1_GGAAm = `MSC7-EF1_FLI1_GGAAm`)

  # 1.3 Merge all cell lines----------------------------------------------------
all.models <- full_join(A673, TC71)
all.models <- full_join(all.models, EW7)
all.models <- full_join(all.models, EW16)
all.models <- full_join(all.models, EW7)
all.models <- full_join(all.models, Ewima1)

write.xlsx(all.models, file = ".../all.models_merged.xlsx", sep="\t")

  # 1.4 Automatically annotate Long isoforms and truncated isoforms ------------

# Here I manually changed ENSG_id to remove the transcript version
# Manually change the TSS to remove chromosome

# Load the updated file
all.models <- read.delim(".../all.models_merged_modified.txt", sep="\t")

# Add genomic start and end
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://feb2014.archive.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

all.models = all.models %>% 
  rename(ensembl_gene_id = ENSG_id)

ENSG_ids <- unique(all.models$ensembl_gene_id)

t2g<-getBM(attributes=c('ensembl_gene_id','start_position', 'end_position', 'transcript_start'),  values = ENSG_ids, mart = ensembl)
dim(t2g)

# Remove duplicates of genes in BiomaRt database
t2g = t2g %>% filter(!duplicated(ensembl_gene_id),
                     ensembl_gene_id %in% all.models$ensembl_gene_id)
                                 
dim(t2g)

# Join info from BiomaRt and novel isoforms table
all.models = left_join(all.models, t2g, by="ensembl_gene_id")
dim(all.models)

# Annotate intragenic TSS isoforms
all.models = all.models %>%
         mutate(GGAA_FLI_pos = case_when(Strand=='+' ~ TSS + closest_FLI1_GGAAm,
                                  Strand=='-' ~ TSS - closest_FLI1_GGAAm),
         TSS_in_Genomic_seq = ifelse(TSS < end_position & TSS > start_position, "Yes", "No"), # TSS in the genomic sequence
         GGAA_in_Genomic_seq = ifelse(GGAA_FLI_pos < end_position & GGAA_FLI_pos > start_position, "Yes", "No"), # GGAAm and FLI1 peak in the genomic sequence
         Shorter_transcript = ifelse(Length_Transcript...64 < Length_Ref_Transcript, 'Yes', "No"), # Shorter transcript
         Short.isoform = ifelse(TSS_in_Genomic_seq =="Yes" & 
                                  GGAA_in_Genomic_seq == "Yes" &
                                  Shorter_transcript == "Yes", 
                                "Short_isoform", "No"))

#  Annotate Long isoforms
all.models = all.models %>%
  mutate(TSS_before_gene = case_when(Strand=="+" & TSS < start_position ~ "Yes",
                                     Strand=="+" & TSS > start_position ~ "No",
                                     Strand=="-" & TSS > end_position ~ "Yes",
                                     Strand=="-" & TSS < end_position ~ "No"),
         GGAA_before_gene = case_when(Strand=="+" & GGAA_FLI_pos < start_position ~ "Yes",
                                      Strand=="+" & GGAA_FLI_pos > start_position ~ "No",
                                      Strand=="-" & GGAA_FLI_pos > end_position ~ "Yes",
                                      Strand=="-" & GGAA_FLI_pos < end_position ~ "No"),
         Long.isoform = ifelse(TSS_before_gene == "Yes" & 
                                 GGAA_before_gene == "Yes" & 
                                 Length_Transcript...64 > Length_Ref_Transcript + 100, # Transcript longer than reference
                               "Long.isoform", "No")) 

  # 1.5 Export table to get unique isoforms------------------------------------
write.xlsx(all.models, file = ".../all.models_merged_before_gff_compare.xlsx", sep="\t")

### 2. Last filtering steps ----------------------------------------------------
# Load table with unique isoforms after merging all cell lines
all.models.unique <- read.delim(".../merged_dataset_new_filter_rescue_isoforms.txt", sep="\t")

# Create columns to know in which cell line this isoform was detected
all.models.unique = all.models.unique %>% 
                        mutate(Detected.in.TC71 = ifelse(grepl("TC71", all.models.unique$merged_ids)==TRUE, 1, 0),
                               Detected.in.A673 = ifelse(grepl("A673", all.models.unique$merged_ids)==TRUE, 1, 0),
                               Detected.in.EW7 = ifelse(grepl("EW7", all.models.unique$merged_ids)==TRUE, 1, 0),
                               Detected.in.Ewima1 = ifelse(grepl("MSC-EF1", all.models.unique$merged_ids)==TRUE,1, 0),
                               Detected.in.EW16 = ifelse(grepl("EW16", all.models.unique$merged_ids)==TRUE, 1, 0),
                               Detected.in.PacBio = ifelse(PacBioId!="-", 1, 0),
                               Detected.in.Nanopore = ifelse(NanoporeId!="-", 1, 0))

# Novel isoforms were detected in Nanopore and PacBio ?
all.models.unique = all.models.unique %>% 
  group_by(new_id) %>%
  mutate(PacBio.nb.detections = sum(Detected.in.PacBio)) %>%
  ungroup() %>% 
  group_by(new_id) %>%
  mutate(Nanopore.nb.detections = sum(Detected.in.Nanopore)) %>%
  ungroup()

all.models.unique = all.models.unique %>% 
  mutate(Nb.CL = rowSums(dplyr::select(., Detected.in.TC71, Detected.in.A673, Detected.in.EW7, Detected.in.Ewima1, Detected.in.EW16)),
         Nb.CL.T.F = ifelse(rowSums(dplyr::select(., Detected.in.TC71, Detected.in.A673, Detected.in.EW7, Detected.in.Ewima1, Detected.in.EW16)) >= 2, TRUE, FALSE),
         Detected.in.Nanopore.and.PacBio = ifelse(PacBio.nb.detections>=1 & Nanopore.nb.detections >=1, "yes", "no"))

dim(all.models.unique)

  # 2.2 Filter out lowly expressed transcripts----------------------------------
all.models.unique = all.models.unique %>%
  filter(`Median_Transcript_EWING_Cell_line.33.` > 1.5 |
           Median_Ewima1.3. > 1.2 | 
           Nb_reads_Nanopore > 10)

dim(all.models.unique)

  # 2.3 Create a unique column for NG annotation--------------------------------
all.models.unique  = all.models.unique %>%
  mutate(isoform_type = case_when(Short.isoform=="Short_isoform" ~ "Truncated_isoform",
                                  Long.isoform=="Short.isoform" ~ "Long_isoform",
                                  grepl("Ew_NG", Gene_Annot_if_different) ~ "Neogene"))

# Export files
write.xlsx(all.models.unique, file = ".../all.models_merged_after_gff_compare.xlsx", sep="\t")

### 3. Annotate with Sqanti3----------------------------------------------------
# Load sqanti annotation 
sqanti.annot <- read.delim(".../subset_unique_neoisoforms_classification.txt",
                           sep="\t", header = T)
# Load manually curated isoforms
iso <- read.xlsx(".../all.models_merged_after_gff_compare_curated.xlsx")
iso = iso %>%
  rename(isoform = new_id)

iso <- left_join(iso, sqanti.annot, by="isoform")

table(iso$structural_category)

# Export files
write.xlsx(iso, file = ".../all.models_merged_after_gff_compare_curated_sqanti_annot.xlsx", sep="\t")