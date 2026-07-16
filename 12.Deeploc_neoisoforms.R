# Load libraries ---------------------------------------------------------------
library(dplyr)
library(readxl)
library(openxlsx)
library(biomaRt)
library(tidyr)
library(ggplot2)

# Load data --------------------------------------------------------------------
deep_neo <- read.csv(".../neotranscript_deeploc.csv")
annot <- read.csv(".../neoisoform_orf_table_with_dups.csv")

# Create a fasta file for canonical proteins with MANE transcript --------------
# Connect to Ensembl

# Disable curie proxy
Sys.unsetenv("http_proxy")
Sys.unsetenv("https_proxy")
Sys.unsetenv("HTTP_PROXY")
Sys.unsetenv("HTTPS_PROXY")


mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# Your ENSG list
genes <- unique(annot$gene_id)
length(genes)

# Query
res <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "transcript_mane_select",
    "peptide"
  ),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)

# Keep only MANE Select
mane <- res[res$transcript_mane_select != "", ]

write_fasta <- function(df, file="mane_proteins.fasta") {
  con <- file(file, "w")
  
  for (i in 1:nrow(df)) {
    header <- paste0(
      ">",
      df$ensembl_gene_id[i], "|",
      df$ensembl_transcript_id[i], "|",
      df$transcript_mane_select[i]
    )
    
    seq <- df$peptide[i]
    
    writeLines(header, con)
    writeLines(seq, con)
  }
  
  close(con)
}

write_fasta(mane, ".../canonical_proteins.fa")

# Add the annotation on deep_neo file-------------------------------------------
deep_neo = deep_neo %>%
  dplyr::rename(known_orf_id = Protein_ID)
deep_neo_annotated <- left_join(deep_neo, annot, by="known_orf_id")
dim(deep_neo_annotated)

# Now match with the canonical deeploc -----------------------------------------
deep_canonical <- read.csv(".../canonical_deeploc.csv")
deep_canonical <- deep_canonical %>%
  separate(Protein_ID, into = c("gene_id", "Transcript", "RefSeq1", "RefSeq2"), sep = "_")

colnames(deep_canonical) = paste0("canonical_", colnames(deep_canonical))
head(colnames(deep_canonical))

deep_canonical = deep_canonical %>% dplyr::rename(gene_id = canonical_gene_id)

setdiff(deep_neo_annotated$gene_id, deep_canonical$gene_id)

deep_neo_canonical <- left_join(deep_neo_annotated, deep_canonical, by="gene_id")
dim(deep_neo_canonical)

# Now annotate the 2'244 isoforms-----------------------------------------------
iso <- read_xlsx(".../neoisoforms_ORFs_with_dups.xlsx")

iso = iso %>% dplyr::select(Transcript_Id,
                     isoform_version, 
                     isoform_final, 
                     isoform_final_subtype, 
                     Gene_automatic_annot, 
                     Gene_curated, 
                     ensembl_gene_id,
                     neotranscript_orf_type,
                     known_orf_id) %>%
  filter(neotranscript_orf_type %in% c("In-frame ORF", "Out-of-frame ORF")) # keep only isoforms encoding for translons

dim(iso)
length(unique(iso$known_orf_id))


# Finally merge isoforms with deeploc info -------------------------------------
df <- left_join(iso, deep_neo_canonical, by="known_orf_id")
table(df$isoform_final)

write.table(df, ".../neoisoforms_canonical_deeploc.txt", sep="\t")

# Calculate the number of isoforms ---------------------------------------------
# Number of neoexons encoding for EwS-specific translons for figure

neoexon = df %>% filter(isoform_final %in% c("Neoexon", "Neoexon; Intragenic TSS isoform"))
length(unique(neoexon$isoform_version))

neoexon_only = df %>% filter(isoform_final %in% c("Neoexon"))
length(unique(neoexon_only$isoform_version))

intragenic = df %>% filter(isoform_final %in% c("Intragenic TSS isoform", "Neoexon; Intragenic TSS isoform"))
length(unique(intragenic$isoform_version))

intragenic_only = df %>% filter(isoform_final %in% c("Intragenic TSS isoform"))
length(unique(intragenic_only$isoform_version))

`%nin%` = Negate(`%in%`)

# Create a diff localization flag per row
df2 <- df %>%
  mutate(diff_loc = Localizations != canonical_Localizations,
         gained_extra_membrane = case_when(Localizations %in% c("Cell membrane", "Extracellular", "Extracellular|Endoplasmic reticulum")  &  
                                             canonical_Localizations %nin% c("Cell membrane", "Extracellular", "Extracellular|Endoplasmic reticulum") ~ "TRUE",
                                           TRUE ~ "FALSE"))

table(df2$gained_extra_membrane)
table(df2$diff_loc)

# Number of different localizations compared to canonical ----------------------
# Collapse to isoform level
isoform_summary <- df2 %>%
  group_by(isoform_version, isoform_final) %>%
  summarise(
    any_diff = any(diff_loc, na.rm = TRUE),
    .groups = "drop"
  )

# Count per isoform type 
counts_all <- isoform_summary %>%
  group_by(isoform_final) %>%
  summarise(
    total_isoforms = n(),
    diff_isoforms = sum(any_diff),
    .groups = "drop"
  )


p <- ggplot(counts_all, aes(x = isoform_final, y = total_isoforms)) +
  geom_bar(stat = "identity", fill = "grey80", color = "black", linewidth = 0.3) +
  geom_bar(aes(y = diff_isoforms), stat = "identity", fill = "steelblue", color = "black", linewidth = 0.3) +
  labs(
    y = "Number of isoforms"
  ) +
  theme_classic() +
  theme(legend.text = element_text(size=8, color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=8, color="black"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size=8, color="black"),
        axis.line = element_line(size = 0.3),
        axis.text.x = element_text(size=8, angle = 45, hjust =1, color="black"),
        axis.text.y = element_text(size=8, color="black", margin = margin(5, 5, 5, 5))) 

p

ggsave(".../Deeploc_prediction_neoisoforms.pdf", plot = p, 
       device = "pdf", width = 6, height = 9, units = "cm", dpi = 600)

# Statistics to make a sentence in the manuscript ------------------------------
neo_stats <- isoform_summary %>%
  filter(isoform_final %in% c("Intragenic TSS isoform", "Neoexon", "Neoexon; Intragenic TSS isoform")) %>%
  summarise(
    total = n(),
    diff = sum(any_diff)
  )

neo_pct <- 100 * neo_stats$diff / neo_stats$total
neo_pct

# Number of neoproteins with gained domains ------------------------------------
class(df2$gained_extra_membrane)

df3 <- df2 %>%
  group_by(isoform_version, isoform_final) %>%
  summarise(
    known_orf_id = paste(unique(known_orf_id), collapse = ";"),
    any_gained_extra_membrane = any(gained_extra_membrane, na.rm = TRUE),
    .groups = "drop"
  )

# Count per isoform type 
counts_all <- df3 %>%
  group_by(isoform_final) %>%
  summarise(
    total_isoforms = n(),
    diff_isoforms = sum(any_gained_extra_membrane),
    .groups = "drop"
  )

# List of isoforms with gained extracellular or TM localization ----------------
gained_extra_membrane <- df3 %>%
  filter(any_gained_extra_membrane == TRUE)

isoforms_with_extra_TM <- gained_extra_membrane$isoform_version
isoforms_with_extra_TM

# Number of translons with gained TM or extracellular prediction
extra_TM_translons <- unique(gained_extra_membrane$known_orf_id)
extra_TM_translons

extra_TM_translons <- unlist(strsplit(extra_TM_translons, ";"))
extra_TM_translons = unique(extra_TM_translons)
extra_TM_translons

if_oof_table <- annot %>% filter(known_orf_id %in% extra_TM_translons,
                                 !duplicated(known_orf_id))
  
dim(if_oof_table)
table(if_oof_table$neotranscript_orf_type)
