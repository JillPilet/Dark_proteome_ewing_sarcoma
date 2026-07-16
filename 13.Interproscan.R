# Load libraries ---------------------------------------------------------------
library(dplyr)
library(readxl)
library(openxlsx)
library(biomaRt)
library(tidyr)
library(ggplot2)
library(Biostrings)
library(paletteer)

# Load interproscan outputs for neoisoforms and annot---------------------------
inframe_1 <- read.delim(".../Neoisoform_if_interproscan1.tsv", header = F, na.strings = "-")
inframe_2 <- read.delim(".../Neoisoform_if_interproscan3.tsv", header = F, na.strings = "-")
out_of_frame <- read.delim(".../Neoisoform_oof_interproscan.tsv", header = F, na.strings = "-")
annot <- read.csv(".../neoisoform_orf_table_with_dups.csv")
iso <- read_xlsx(".../neoisoforms_ORFs_with_dups.xlsx")

colnames_interpro <- c("orf_id", "md5", "length", "analysis", "sig_accession", 
                                    "sig_description", "start", "stop", "score", "status", 
                                    "date", "ipr_accession", "ipr_description", "go_terms", "pathways")

colnames(inframe_1) <- colnames_interpro
colnames(inframe_2) <- colnames_interpro
colnames(out_of_frame) <- colnames_interpro

# Merge interproscan outputs ---------------------------------------------------
df = rbind(inframe_1, inframe_2)
df = rbind(df, out_of_frame)
dim(df)

# Classify analysis sources ----------------------------------------------------
functional_dbs <- c("Pfam", "SMART", "CDD", "SUPERFAMILY", "Gene3D", 
                    "PANTHER", "ProSiteProfiles", "ProSitePatterns", 
                    "PRINTS", "PIRSF", "Hamap", "TIGRFAM")
structural_dbs <- c("Phobius", "SignalP", "TMHMM", "SignalP_EUK", "SignalP_GRAM_POSITIVE")
disorder_dbs   <- c("MobiDBLite", "Coils")

df <- df %>%
  mutate(hit_category = case_when(
    analysis %in% functional_dbs ~ "Functional domain",
    analysis %in% structural_dbs ~ "Signal peptide / TM",
    analysis %in% disorder_dbs   ~ "Disordered / coiled-coil",
    TRUE                         ~ "Other"
  ))

all_orfs_in_results <- unique(df$orf_id)

orf_summary <- df %>%
  group_by(orf_id) %>%
  summarise(
    has_functional  = any(hit_category == "Functional domain"),
    has_structural  = any(hit_category == "Signal peptide / TM"),
    has_disorder    = any(hit_category == "Disordered / coiled-coil"),
    .groups = "drop"
  ) %>%
  mutate(orf_class = case_when(
    has_functional  ~ "Functional domain",
    has_structural  ~ "Signal peptide / TM only",
    has_disorder    ~ "Disordered only",
    TRUE            ~ "Other"
  )) 


# Add no-hit ORFs (requires all_orfs vector)
no_hit <- tibble(orf_id = setdiff(iso$known_orf_id, orf_summary$orf_id), 
                 has_functional = FALSE, 
                 has_structural = FALSE,
                 has_disorder = FALSE,
                 orf_class = "No domain detected")
no_hit = no_hit %>% filter(!is.na(orf_id))

orf_summary = rbind(orf_summary, no_hit, by="orf_id")
orf_summary = orf_summary %>% filter(orf_id!="orf_id")

# Add the annotation with gene_id file------------------------------------------
orf_summary = orf_summary %>%
  dplyr::rename(known_orf_id = orf_id)
orf_summary <- left_join(orf_summary, annot, by="known_orf_id")

# Load interproscan output for canonical isoforms ------------------------------
canonical_1 <- read.delim(".../interproscan_canonical_1.tsv", header = F, na.strings = "-")
canonical_2 <- read.delim(".../interproscan_canonical_2.tsv", header = F, na.strings = "-")

# Merge canonical 1 and 2 ------------------------------------------------------
colnames(canonical_1) <- colnames_interpro
colnames(canonical_2) <- colnames_interpro
canonical_df = rbind(canonical_1, canonical_2)

dim(canonical_df)

# Annotate the type of domain --------------------------------------------------
canonical_df <- canonical_df %>%
  mutate(hit_category = case_when(
    analysis %in% functional_dbs ~ "Functional domain",
    analysis %in% structural_dbs ~ "Signal peptide / TM",
    analysis %in% disorder_dbs   ~ "Disordered / coiled-coil",
    TRUE                         ~ "Other"
  ))

all_orfs_in_results <- unique(canonical_df$orf_id)

orf_summary_canonical <- canonical_df %>%
  group_by(orf_id) %>%
  summarise(
    has_functional  = any(hit_category == "Functional domain"),
    has_structural  = any(hit_category == "Signal peptide / TM"),
    has_disorder    = any(hit_category == "Disordered / coiled-coil"),
    .groups = "drop"
  ) %>%
  mutate(orf_class = case_when(
    has_functional  ~ "Functional domain",
    has_structural  ~ "Signal peptide / TM only",
    has_disorder    ~ "Disordered only",
    TRUE            ~ "Other"
  ))

orf_summary_canonical <- orf_summary_canonical %>%
  separate(orf_id, into = c("gene_id", "Transcript", "RefSeq1"), sep = "\\|")

# Add no-hit ORFs (requires all_orfs vector)
no_hit <- tibble(gene_id = setdiff(iso$gene_id, orf_summary_canonical$gene_id), 
                 Transcript = NA, 
                 RefSeq1 = NA,
                 has_functional = FALSE, 
                 has_structural = FALSE,
                 has_disorder = FALSE,
                 orf_class = "No domain detected")

orf_summary_canonical = rbind(orf_summary_canonical, no_hit, by="gene_id")
orf_summary_canonical = orf_summary_canonical %>% filter(gene_id!="gene_id")

# Add a prefix to distinct canonical vs neoisforms translons--------------------
colnames(orf_summary_canonical) = paste0("canonical_", colnames(orf_summary_canonical))
head(colnames(orf_summary_canonical))

orf_summary_canonical = orf_summary_canonical %>% dplyr::rename(gene_id = canonical_gene_id)

inter_results <- left_join(orf_summary, orf_summary_canonical, by="gene_id")

# Now annotate the 2'244 isoforms-----------------------------------------------
iso = iso %>% dplyr::select(Transcript_Id,
                            isoform_version, 
                            isoform_final, 
                            isoform_final_subtype, 
                            Gene_automatic_annot, 
                            Gene_curated, 
                            ensembl_gene_id,
                            neotranscript_orf_type,
                            known_orf_id) %>%
  filter(neotranscript_orf_type %in% c("In-frame ORF", "Out-of-frame ORF"))

dim(iso)
length(unique(iso$known_orf_id))


# Finally merge isoforms with deeploc info -------------------------------------
inter_results_annotated <- left_join(inter_results, iso, by="known_orf_id")
table(inter_results_annotated$isoform_final)

write.table(inter_results_annotated, ".../Interproscan_if_oof_results.txt", sep="\t")

# Annotate if there is gain or loss domains ------------------------------------
gain_isoforms <- inter_results_annotated %>%
  mutate(
    has_functional = as.logical(has_functional),
    canonical_has_functional = as.logical(canonical_has_functional),
    
    has_structural = as.logical(has_structural),
    canonical_has_structural = as.logical(canonical_has_structural),
    
    has_disorder = as.logical(has_disorder),
    canonical_has_disorder = as.logical(canonical_has_disorder),
    
    gain_function   = has_functional & !canonical_has_functional,
    gain_structural = has_structural & !canonical_has_structural,
    gain_disorder   = has_disorder & !canonical_has_disorder
  ) %>%
  group_by(isoform_version) %>%
  summarise(
    gain_function   = any(gain_function, na.rm = TRUE),
    gain_structural = any(gain_structural, na.rm = TRUE),
    gain_disorder   = any(gain_disorder, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("gain_"),
    names_to = "category",
    values_to = "gain"
  ) %>%
  filter(gain)

counts <- gain_isoforms %>%
  count(category)

cols <- paletteer_d("ggthemes::wsj_rgby", n=3) 

p <- ggplot(counts, aes(x = category, y = n, fill = category)) +
  scale_fill_manual(values=rep("#D5695DFF",3)) +
  geom_bar(stat = "identity", color = "black", width = 0.9, linewidth = 0.3) +
  labs(
    x = NULL,
    y = "Number of Isoforms"
  ) +
  theme_classic() +
  labs(x="", y="Number of isoforms with at \n
       least one gained domain") +
  theme(legend.text = element_text(size=8, color="black"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.position = "none",
        axis.title = element_text(size=8, color="black"),
        axis.line = element_line(size = 0.3),
        axis.text.x = element_text(size=8, angle = 45, hjust =0.9, color="black"),
        axis.text.y = element_text(size=8, color="black", margin = margin(5, 5, 5, 5)))

p

ggsave(".../gained_domains_Interproscan.pdf", plot = p, 
       device = "pdf", width = 6, height = 8, units = "cm", dpi = 600)

# Annotate loss domains --------------------------------------------------------
loss_isoforms <- inter_results_annotated %>%
  mutate(
    has_functional = as.logical(has_functional),
    canonical_has_functional = as.logical(canonical_has_functional),
    
    has_structural = as.logical(has_structural),
    canonical_has_structural = as.logical(canonical_has_structural),
    
    has_disorder = as.logical(has_disorder),
    canonical_has_disorder = as.logical(canonical_has_disorder),
    
    loss_function   = !has_functional & canonical_has_functional,
    loss_structural = !has_structural & canonical_has_structural,
    loss_disorder   = !has_disorder & canonical_has_disorder
  ) %>%
  group_by(isoform_version) %>%
  summarise(
    loss_function   = any(loss_function, na.rm = TRUE),
    loss_structural = any(loss_structural, na.rm = TRUE),
    loss_disorder   = any(loss_disorder, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("loss_"),
    names_to = "category",
    values_to = "loss"
  ) %>%
  filter(loss)

loss_counts <- loss_isoforms %>%
  count(category)


# Plot lost domains ------------------------------------------------------------

p <- ggplot(loss_counts, aes(x = category, y = n, fill = category)) +
  scale_fill_manual(values=rep("#D3BA68FF", 3)) +
  geom_bar(stat = "identity", color = "black", width = 0.9, linewidth = 0.3) +
  labs(
    x = NULL,
    y = "Number of Isoforms"
  ) +
  theme_classic() +
  labs(x="", y="Number of isoforms with at \n
       least one lost domain") +
  theme(legend.text = element_text(size=8, color="black"),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "none",
      axis.title = element_text(size=8, color="black"),
      axis.line = element_line(size = 0.3),
      axis.text.x = element_text(size=8, angle = 45, hjust =0.9, color="black"),
      axis.text.y = element_text(size=8, color="black", margin = margin(5, 5, 5, 5)))

p

ggsave(".../lost_domains_Interproscan.pdf", plot = p, 
       device = "pdf", width = 6, height = 8, units = "cm", dpi = 600)

# Number of genes with gained domains ------------------------------------------
genes_corres <- iso %>% dplyr::select(isoform_version, Gene_curated) %>%
  filter(!duplicated(isoform_version))
dim(genes_corres)

nb_genes_gain <- left_join(gain_isoforms, genes_corres, by="isoform_version")
dim(nb_genes_gain)
length(unique(nb_genes_gain$Gene_curated))
unique(nb_genes_gain$Gene_curated)
length(unique(nb_genes_gain$isoform_version))

# Number of genes with lost domains --------------------------------------------
nb_genes_lost <- left_join(loss_isoforms, genes_corres, by="isoform_version")
length(unique(nb_genes_lost$Gene_curated))
length(unique(nb_genes_lost$isoform_version))

# What are the domains gained --------------------------------------------------
gain_orf_level <- inter_results_annotated %>%
  mutate(
    has_functional = as.logical(has_functional),
    canonical_has_functional = as.logical(canonical_has_functional),
    
    has_structural = as.logical(has_structural),
    canonical_has_structural = as.logical(canonical_has_structural),
    
    has_disorder = as.logical(has_disorder),
    canonical_has_disorder = as.logical(canonical_has_disorder),
    
    gain_function   = has_functional & !canonical_has_functional,
    gain_structural = has_structural & !canonical_has_structural,
    gain_disorder   = has_disorder & !canonical_has_disorder
  )

gained_orfs <- gain_orf_level %>%
  filter(gain_function | gain_structural | gain_disorder) %>%
  dplyr::select(known_orf_id, isoform_version, gene_id)

gained_domains <- df %>%
  filter(orf_id %in% gained_orfs$known_orf_id) %>%
  left_join(gained_orfs, by = c("orf_id" = "known_orf_id"))

canonical_domains <- canonical_df %>%
  dplyr::select(gene_id = orf_id, ipr_accession)

gained_domains_filtered <- gained_domains %>%
  filter(!ipr_accession %in% canonical_domains$ipr_accession | is.na(ipr_accession))

domain_table <- gained_domains_filtered %>%
  dplyr::select(isoform_version, gene_id, ipr_accession, ipr_description, analysis, hit_category)

write.table(domain_table, ".../Gains_interproscan_details.txt", sep="\t")

