################################################################################
# Filter isoforms identified in PacBio, Nanopore :
#     1. Expression level (TPM) in Illumina RNA-seq
#     2. Presence of a GGAA-microsatellite + FLI1 ChIP-seq peak within the TAD
#     3. Transcriptional modulation upon EWS-FLI1 depletion
#     4. Proximity to open chromatin marks (H3K27ac / H3K4me3)
#     5. Low/no expression in mesenchymal stem cells (MSCs)

# The same approach was used for each cell line : A-673, TC-71, EW-7, EW-16, Ewima1

# Load packages ----------------------------------------------------------------
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2) 

# Create data_dir & OUTPUT_DIR -------------------------------------------------
# Directories containing the input annotation file and where outputs are written.

# Choose filters applied -------------------------------------------------------

gene_exp_sup_15TPM <- c("yes", "no")[1]
GGAAm_FLI1_in_TAD <- c("yes", "no")[1]
Modulated_by_EF1 <- c("yes", "no")[1]
Open_chromatin <- c("yes", "no")[1]
Expression_MSC <- c("yes", "no")[1]

# Load data --------------------------------------------------------------------
A673 <- read_xlsx(file.path(data_dir, "A673_NANO_vs_PACBIO_merge_all_annot_20231219_tad_jill.xlsx"))

# (1) Filter most expressed genes in Illumina ----------------------------------

if(gene_exp_sup_15TPM=="yes") {
A673.f = A673 %>% 
  filter(Expression_A673_R4_D066T33 >15)
}

# (2) Filter transcripts with GGAAm and FL1 peak in the TAD --------------------

if(GGAAm_FLI1_in_TAD=="yes") {
A673.f = A673.f %>% 
  filter(FLI1_GGAAm_in_TAD=="Yes")
}

# (3) Filter transcripts downregulated upon EWS-FLI1 depletion -----------------

if(Modulated_by_EF1=="yes") {
  A673.f = A673.f %>% 
    filter(A673_log2FC_High_Low>0.5)
}


# (4) Filter transcripts next to open chromatin state --------------------------
A673.f$A673_H3K27ac = as.numeric(A673.f$A673_H3K27ac)
A673.f$A673_H3K4me3 = as.numeric(A673.f$A673_H3K4me3)

if(Open_chromatin=="yes") { 
A673.f = A673.f %>% 
  filter(A673_H3K27ac==0 | A673_H3K4me3==0)
}

dim(A673.f)

# number of genes
length(unique(A673.f$Gene_Id))

# (5) Filter out transcripts expressed in MSCs ---------------------------------
A673.f$`Median_Transcript_MSC(6)` <- as.numeric(A673.f$`Median_Transcript_MSC(6)`)

if(Expression_MSC=="yes") { 
A673.f = A673.f %>% 
  filter(`Median_Transcript_MSC(6)`<=0.03)
}

dim(A673.f)

# (6) Export new table filtered  -----------------------------------------------

if(gene_exp_sup_15TPM=="yes" & 
   GGAAm_FLI1_in_TAD =="yes"&
   Modulated_by_EF1 =="yes"& 
   Open_chromatin =="yes"&
   Expression_MSC=="yes") { 
write.xlsx(A673.f, file = "OUTPUT_DIR/A673_filtered_all_filters.xlsx", sep="\t")
}

