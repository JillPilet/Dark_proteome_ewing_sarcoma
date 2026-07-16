# Load libraries 
library(readxl)     
library(dplyr)      
library(tidyr)      
library(stringr)     
library(ggplot2)    
library(openxlsx)    
library(UpSetR)      
library(extrafont)  
loadfonts()

### 0. Load tables -------------------------------------------------------------
# Table with ORF prediction but when found from multiple isoforms, only one isoform referenced
ORFs <- read.csv(".../neoisoform_orf_table.csv")
#write.xlsx(ORFs, ".../ORFs_Riboseq/neoisoform_orf_table.xlsx")

# Table with ORF prediction but when found from multiple isoforms, all isoforms annotated
ORFs_with_dups <- read.csv(".../neoisoform_orf_table_with_dups.csv")
#write.xlsx(ORFs_with_dups, ".../ORFs_Riboseq/neoisoform_orf_table_with_dups.xlsx")

# Table with neogenes and isoforms identified in Curie
curie_neo <- readxl::read_xlsx(".../neogenes_isoforms_annot_20260320.xlsx")
curie_neo = curie_neo %>% dplyr::rename(Gene_curated = Gene_curated.x)

### 1. Do some basic statistics to know how many neogenes are there ------------
  # 1.1 Truncated isoforms -----------------------------------------------------
table(curie_neo$isoform_final)
truncated = curie_neo %>% filter(grepl("Intragenic", isoform_final))
dim(truncated)
unique(truncated$Gene_curated)

  # 1.2 Non-canonical isoforms -------------------------------------------------
table(curie_neo$isoform_final)
Noncanonical = curie_neo %>% filter(isoform_final =="Non-canonical")
dim(Noncanonical)
unique(Noncanonical$Gene_curated)

  # 1.3 Neogenes isoforms ------------------------------------------------------
table(curie_neo$isoform_final)
Neogene = curie_neo %>% filter(isoform_final =="Neogene")
dim(Neogene)
length(unique(Neogene$Gene_curated))

# Filter only isoforms from known genes because here we focus on isoforms
curie_neo = curie_neo %>% filter(isoform_final %in% c("Intragenic TSS isoform", "Neoexon", "Neoexon; Intragenic TSS isoform", "Non-canonical"))
dim(curie_neo)

table(curie_neo$structural_category)
table(curie_neo$subcategory)


### 2. Upset plots -------------------------------------------------------------
  # 2.1 Genes and isoforms detected in different cell lines --------------------
# Per isoform
curie_neo_upset = curie_neo %>% select(c("Detected.in.A673", "Detected.in.TC71", "Detected.in.EW7", "Detected.in.EW16", "Detected.in.Ewima1"))
curie_neo_upset = as.data.frame(curie_neo_upset)
cols <-c("Detected.in.A673", "Detected.in.TC71", "Detected.in.EW7", "Detected.in.EW16", "Detected.in.Ewima1")

Isoforms_CL <- upset(
  curie_neo_upset,
  sets = cols,
  nsets = length(cols),       # number of sets shown
  nintersects = NA,           # show all intersections
  order.by = "freq",          # order intersections by frequency
  keep.order = TRUE,          # maintain order of sets
  sets.bar.color = "#34bfc7",   # color for set bars
  main.bar.color = "#34bfc7",   # color for intersection bars
  text.scale = c(1, 1, 1, 1, 1, 1), # font sizes
  mainbar.y.label = "Number of Isoforms",      # label for intersection bars
  sets.x.label = "Detected in Sample", # label for set bars
  line.size = 0.5,
  shade.alpha = 0, matrix.dot.alpha = 0, point.size = 3, matrix.color = "gray80"
)

Isoforms_CL

pdf(".../Upset_Nb_isoforms_cell_lines.pdf",
    width = 18/2.54,  # convertir cm en inches
    height = 8/2.54)  # convertir cm en inches

# Dessiner le plot UpSet
Isoforms_CL

# Fermer le device
dev.off()


  # Per Gene 
curie_neo_genes <- curie_neo %>%
  group_by(Gene_curated) %>%
  summarise(
    Gene.Detected.in.A673  = as.integer(any(Detected.in.A673 == 1)),
    Gene.Detected.in.TC71  = as.integer(any(Detected.in.TC71 == 1)),
    Gene.Detected.in.EW7   = as.integer(any(Detected.in.EW7 == 1)),
    Gene.Detected.in.EW16  = as.integer(any(Detected.in.EW16 == 1)),
    Gene.Detected.in.Ewima1 = as.integer(any(Detected.in.Ewima1 == 1))
  ) %>%
  ungroup()

curie_neo_genes = as.data.frame(curie_neo_genes)

rownames(curie_neo_genes) <- curie_neo_genes$Gene_curated

curie_neo_genes$Gene_curated = NULL

cols <- c("Gene.Detected.in.A673", "Gene.Detected.in.TC71", "Gene.Detected.in.EW7", "Gene.Detected.in.EW16", "Gene.Detected.in.Ewima1")
genes_CL <- upset(
  curie_neo_genes,
  sets = cols,
  nsets = length(cols),       # number of sets shown
  nintersects = NA,           # show all intersections
  order.by = "freq",          # order intersections by frequency
  keep.order = TRUE,          # maintain order of sets
  sets.bar.color = "#34bfc7",   # color for set bars
  main.bar.color = "#34bfc7",   # color for intersection bars
  text.scale = c(1, 1, 1, 1, 1, 1), # font sizes
  mainbar.y.label = "Number of Genes",      # label for intersection bars
  sets.x.label = "Detected in Sample", # label for set bars
  line.size = 0.5,
  shade.alpha = 0, matrix.dot.alpha = 0, point.size = 3, matrix.color = "gray80"
)

genes_CL

pdf(".../Upset_Nb_genes_cell_lines.pdf",
    width = 18/2.54,  # convertir cm en inches
    height = 8/2.54)  # convertir cm en inches

# Dessiner le plot UpSet
genes_CL

# Fermer le device
dev.off()

  # 2.2 Compare genes retrieved in Pacbio vs Nanopore---------------------------
# Per isoform
curie_neo_upset = curie_neo %>% select(c("Detected.in.Nanopore", "Detected.in.PacBio"))
curie_neo_upset = as.data.frame(curie_neo_upset)
cols <-c("Detected.in.Nanopore", "Detected.in.PacBio")

Isoforms_technique <- upset(
  curie_neo_upset,
  sets = cols,
  nsets = length(cols),       # number of sets shown
  nintersects = NA,           # show all intersections
  order.by = "freq",          # order intersections by frequency
  keep.order = TRUE,          # maintain order of sets
  sets.bar.color = "#cfc353",   # color for set bars
  main.bar.color = "#cfc353",   # color for intersection bars
  text.scale = c(1, 1, 1, 1, 1, 1), # font sizes
  mainbar.y.label = "Number of Isoforms",      # label for intersection bars
  sets.x.label = "Detected in Sample", # label for set bars
  line.size = 0.5,
  shade.alpha = 0, matrix.dot.alpha = 0, point.size = 3, matrix.color = "gray80"
)

Isoforms_technique

pdf(".../Upset_Nb_isoforms_technique.pdf",
    width = 11/2.54,  # convertir cm en inches
    height = 8/2.54)  # convertir cm en inches

# Dessiner le plot UpSet
Isoforms_technique

# Fermer le device
dev.off()


# Per Gene 
curie_neo_genes <- curie_neo %>%
  group_by(Gene_curated) %>%
  summarise(
    Genes.Detected.in.Nanopore  = as.integer(any(Detected.in.Nanopore == 1)),
    Gene.Detected.in.PacBio  = as.integer(any(Detected.in.PacBio == 1)))%>%
  ungroup()

curie_neo_genes = as.data.frame(curie_neo_genes)

rownames(curie_neo_genes) <- curie_neo_genes$Gene_curated

curie_neo_genes$Gene_curated = NULL

cols <- c("Genes.Detected.in.Nanopore", "Gene.Detected.in.PacBio")
genes_technique <- upset(
  curie_neo_genes,
  sets = cols,
  nsets = length(cols),       # number of sets shown
  nintersects = NA,           # show all intersections
  order.by = "freq",          # order intersections by frequency
  keep.order = TRUE,          # maintain order of sets
  sets.bar.color = "#cfc353",   # color for set bars
  main.bar.color = "#cfc353",   # color for intersection bars
  text.scale = c(1, 1, 1, 1, 1, 1), # font sizes
  mainbar.y.label = "Number of Genes",      # label for intersection bars
  sets.x.label = "Detected in Sample", # label for set bars
  line.size = 0.5,
  shade.alpha = 0, matrix.dot.alpha = 0, point.size = 3, matrix.color = "gray80"
)

genes_technique

pdf(".../Upset_Nb_genes_technique.pdf",
    width = 11/2.54,  # convertir cm en inches
    height = 8/2.54)  # convertir cm en inches

# Dessiner le plot UpSet
genes_technique

# Fermer le device
dev.off()


### 3. Number of coding  isoforms per category ---------------------------------
  # 3.1 Recode annotation ------------------------------------------------------
table(ORFs_with_dups$neotranscript_orf_type)

# Expand 'curie_neo' so each merged_id gets its own row
curie_neo_expanded <- curie_neo %>%
  mutate(merged_ids = str_split(merged_ids, ",")) %>%  # split into list
  unnest(merged_ids) %>%                              # expand into rows
  mutate(merged_ids = str_trim(merged_ids))            # remove extra spaces

# Create identifier columns
curie_neo_expanded <- curie_neo_expanded %>%
  mutate(identifier = merged_ids)

ORFs_with_dups <- ORFs_with_dups %>%
  mutate(identifier = transcript_id.x)

# Perform the join
annot_with_ORFs_with_dups <- left_join(curie_neo_expanded, ORFs_with_dups, by = "identifier")

# Reannotation of ORF categories
annot_with_ORFs_with_dups = annot_with_ORFs_with_dups %>%
  mutate(neotranscript_orf_type = case_when(neotranscript_orf_type == "In-frame non-canonical" ~ "In-frame ORF",
                                            neotranscript_orf_type == "In-frame truncated" ~ "In-frame ORF",
                                            neotranscript_orf_type == "Out-of-frame novel" ~ "Out-of-frame ORF",
                                            is.na(neotranscript_orf_type) == TRUE ~ "No ORF"))

# Collapse back to original structure, if desired
annot_with_ORFs_with_dups_collapsed <- annot_with_ORFs_with_dups %>%
  group_by(Transcript_Id) %>%
  summarise(across(everything(), ~ paste(unique(.x), collapse = "; ")))

# Inspect result
dim(annot_with_ORFs_with_dups)
dim(annot_with_ORFs_with_dups_collapsed)
table(annot_with_ORFs_with_dups_collapsed$neotranscript_orf_type)

# Recode an ORF type per isoform otherwise when we want to plot the number of isoforms (2244) we have multiple ORFs and multiple categories
annot_with_ORFs_with_dups_collapsed = annot_with_ORFs_with_dups_collapsed %>%
  mutate(neotranscript_orf_type_per_isoform = case_when(neotranscript_orf_type == "No ORF; In-frame ORF" ~ "In-frame ORF",
                                                        neotranscript_orf_type == "No ORF" ~ "No ORF",
                                                        neotranscript_orf_type == "In-frame ORF" ~ "In-frame ORF",
                                                        neotranscript_orf_type == "In-frame ORF; No ORF" ~ "In-frame ORF",
                                                        neotranscript_orf_type == "Out-of-frame ORF" ~ "Out-of-frame ORF",
                                                        neotranscript_orf_type == "Out-of-frame ORF; In-frame ORF" ~ "Out-of-frame ORF; In-frame ORF",
                                                        neotranscript_orf_type == "Out-of-frame ORF; In-frame ORF; No ORF" ~ "Out-of-frame ORF; In-frame ORF",
                                                        neotranscript_orf_type == "Out-of-frame ORF; No ORF" ~ "Out-of-frame ORF",
                                                        neotranscript_orf_type == "No ORF; Out-of-frame ORF" ~ "Out-of-frame ORF",
                                                        neotranscript_orf_type == "No ORF; Out-of-frame ORF; In-frame ORF" ~ "Out-of-frame ORF; In-frame ORF",
                                                        neotranscript_orf_type == "No ORF; In-frame ORF; Out-of-frame ORF" ~ "Out-of-frame ORF; In-frame ORF"))     

table(annot_with_ORFs_with_dups_collapsed$neotranscript_orf_type)
table(annot_with_ORFs_with_dups_collapsed$neotranscript_orf_type_per_isoform)

# Create a new annotation column to have splice type, neoexons and truncated
annot_with_ORFs_with_dups_collapsed = annot_with_ORFs_with_dups_collapsed %>% 
  mutate(isoform_final_subtype = case_when(isoform_final=="Intragenic TSS isoform" ~ "Intragenic TSS isoform",
                                           isoform_final=="Neoexon" ~ "Neoexon",
                                           isoform_final=="Neoexon; Intragenic TSS isoform" ~ "Neoexon; Intragenic TSS isoform",
                                           isoform_final=="Non-canonical" ~ subcategory))

table(annot_with_ORFs_with_dups_collapsed$isoform_final)

# Also do the same for the table with duplicates
annot_with_ORFs_with_dups = annot_with_ORFs_with_dups %>% 
  mutate(isoform_final_subtype = case_when(isoform_final=="Intragenic TSS isoform" ~ "Intragenic TSS isoform",
                                           isoform_final=="Neoexon" ~ "Neoexon",
                                           isoform_final=="Neoexon; Intragenic TSS isoform" ~ "Neoexon; Intragenic TSS isoform",
                                           isoform_final=="Non-canonical" ~ subcategory))

# number of oof and i-f translons for Neoexons----------------------------------
test = annot_with_ORFs_with_dups %>%
  filter(isoform_final %in% c("Neoexon", "Neoexon; Intragenic TSS isoform"),
         neotranscript_orf_type=="Out-of-frame ORF")
unique(test$known_orf_id)

test = annot_with_ORFs_with_dups %>%
  filter(isoform_final %in% c("Neoexon", "Neoexon; Intragenic TSS isoform"),
         neotranscript_orf_type=="In-frame ORF")
unique(test$known_orf_id)

test = annot_with_ORFs_with_dups %>%
  filter(isoform_final %in% c("Neoexon"),
         neotranscript_orf_type=="Out-of-frame ORF")
unique(test$known_orf_id)

test = annot_with_ORFs_with_dups %>%
  filter(isoform_final %in% c("Neoexon; Intragenic TSS isoform"),
         neotranscript_orf_type=="Out-of-frame ORF")
unique(test$known_orf_id)

test = annot_with_ORFs_with_dups %>%
  filter(isoform_final %in% c("Intragenic TSS isoform"),
         neotranscript_orf_type=="In-frame ORF")
unique(test$known_orf_id)


# number of oof and i-f translons for non-canonical isoforms -------------------
test = annot_with_ORFs_with_dups %>%
  filter(isoform_final %in% c("Non-canonical"),
         neotranscript_orf_type=="Out-of-frame ORF")
unique(test$known_orf_id)

test = annot_with_ORFs_with_dups %>%
  filter(isoform_final %in% c("Non-canonical"),
         neotranscript_orf_type=="In-frame ORF")
unique(test$known_orf_id)

# Number of genes per category -------------------------------------------------
# In nc-isoforms
annot_with_ORFs_with_dups_collapsed %>%
  filter(isoform_final %in% c("Non-canonical"),
         !duplicated(Gene_curated)) %>% pull(Gene_curated)

# In Neoexon
annot_with_ORFs_with_dups_collapsed %>%
  filter(isoform_final %in% c("Neoexon"),
         !duplicated(Gene_curated)) %>% pull(Gene_curated)

# In Intragenic TSS with neoexon 
annot_with_ORFs_with_dups_collapsed %>%
  filter(isoform_final %in% c("Neoexon; Intragenic TSS isoform"),
         !duplicated(Gene_curated)) %>% pull(Gene_curated)

# In Intragenic TSS with neoexon 
annot_with_ORFs_with_dups_collapsed %>%
  filter(isoform_final %in% c("Intragenic TSS isoform"),
         !duplicated(Gene_curated)) %>% pull(Gene_curated)

# Here we can get multiple ORFs per isoform so we need to use the _with_dup table
write.xlsx(annot_with_ORFs_with_dups_collapsed, ".../neoisoforms_ORFs_with_dups_collasped.xlsx")

# However, for plots in which we count isoforms we need to have no isoforms duplicates
write.xlsx(annot_with_ORFs_with_dups, ".../neoisoforms_ORFs_with_dups.xlsx")

  # 3.2 Number of translons derived from isoforms by subcategory ---------------
table(annot_with_ORFs_with_dups$neotranscript_orf_type)
table(annot_with_ORFs_with_dups$isoform_final)


annot <- annot_with_ORFs_with_dups %>%
  filter(neotranscript_orf_type != "No ORF", # Remove isoforms without translons
         isoform_final == "Non-canonical") %>% # Keep only non-canonical isoforms
  distinct(subcategory, known_orf_id, neotranscript_orf_type) %>%
  count(subcategory, neotranscript_orf_type)


vec <- c("Out-of-frame ORF" = "#D0E8F8",
         "In-frame ORF" = "#6EA9D2",
         "Out-of-frame ORF; In-frame ORF" = "#ff9700")

annot$subcategory = factor(annot$subcategory, levels = c(
                                                   "mono-exon_by_intron_retention",
                                                   "mono-exon", "multi-exon",
                                                   "internal_fragment", "3prime_fragment",
                                                   "combination_of_known_splicesites", "5prime_fragment",
                                                   "combination_of_known_junctions", "intron_retention",
                                                   "at_least_one_novel_splicesite"))


p <- ggplot(annot, aes(fill=neotranscript_orf_type, y=n, x=subcategory)) + 
  geom_bar(stat="identity", color = "black", width = 0.9, linewidth = 0.3) + 
  scale_fill_manual(values = vec) + 
  theme_classic() + 
 # scale_y_continuous(limits=c(0,380))+
  theme(legend.text = element_text(size=8, color="black", family = "Helvetica"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size=8, color="black", family = "Helvetica"),
        axis.line = element_line(size = 0.3),
        axis.text.x = element_text(size=8, angle = 45, hjust =0.9, color="black", family = "Helvetica"),
        axis.text.y = element_text(size=8, color="black", family = "Helvetica", margin = margin(5, 5, 5, 5))) + 
  labs(x="", y="Number of translons") + coord_flip()
  

p

ggsave(".../nb_isoforms_non_canonical_per_structural_category_ORFs.pdf", plot = p, 
       device = "pdf", width = 15, height = 6, units = "cm", dpi = 600)


  # 3.3 Number of translons derived from isoforms by in-house annotation -------

annot <- annot_with_ORFs_with_dups %>%
  filter(neotranscript_orf_type != "No ORF") %>%
  distinct(isoform_final, known_orf_id, neotranscript_orf_type) %>%
  count(isoform_final, neotranscript_orf_type)

colnames(annot) <- c("isoform_final", "ORF_type", "value")

annot$isoform_final = factor(annot$isoform_final, c("Neoexon; Intragenic TSS isoform",
                                                    "Intragenic TSS isoform",
                                                    "Neoexon",
                                                    "Non-canonical"))

vec <- c("Out-of-frame ORF" = "#D0E8F8",
         "In-frame ORF" = "#6EA9D2",
         "Out-of-frame ORF; In-frame ORF" = "#ff9700")


p <- ggplot(annot, aes(fill=ORF_type, y=value, x=isoform_final)) + 
  geom_bar(stat="identity", color = "black", width = 0.9, linewidth = 0.3) + 
  scale_fill_manual(values = vec) + 
  theme_classic() + 
  theme(legend.text = element_text(size=8, color="black", family = "Helvetica"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size=8, color="black", family = "Helvetica"),
        axis.line = element_line(size = 0.3),
        axis.text.x = element_text(size=8, angle = 45, hjust =0.9, color="black", family = "Helvetica"),
        axis.text.y = element_text(size=8, color="black", family = "Helvetica", margin = margin(5, 5, 5, 5))) + 
  labs(x="", y="Number of translons") + coord_flip()


p

ggsave(".../nb_isoforms_isoform_final_ORFs.pdf", plot = p, 
       device = "pdf", width = 12, height = 4, units = "cm", dpi = 600)

  # 3.4 Number of non-canonical isoforms ---------------------------------------
dim(annot_with_ORFs_with_dups_collapsed)
table(annot_with_ORFs_with_dups_collapsed$isoform_final)

subset_non_canonical = annot_with_ORFs_with_dups_collapsed %>%
  filter(isoform_final=="Non-canonical")
  
table(subset_non_canonical$subcategory)

annot <- as.data.frame(table(subset_non_canonical$subcategory))
colnames(annot) <- c("subcategory","value")

annot$subcategory = factor(annot$subcategory, levels = c(
  "mono-exon_by_intron_retention",
  "mono-exon", "multi-exon",
  "internal_fragment", "3prime_fragment",
  "combination_of_known_splicesites", "5prime_fragment",
  "combination_of_known_junctions", "intron_retention",
  "at_least_one_novel_splicesite"))


p <- ggplot(annot, aes(y=value, x=subcategory)) + 
  geom_bar(stat="identity", color = "black", width = 0.9, linewidth = 0.3) + 
  scale_fill_manual(values = "black") + 
  scale_y_continuous(limits=c(0,600))+
  theme_classic() + 
  theme(legend.text = element_text(size=8, color="black", family = "Helvetica"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size=8, color="black", family = "Helvetica"),
        axis.line = element_line(size = 0.3),
        axis.text.x = element_text(size=8, angle = 45, hjust =0.9, color="black", family = "Helvetica"),
        axis.text.y = element_text(size=8, color="black", family = "Helvetica", margin = margin(5, 5, 5, 5))) + 
  labs(x="", y="Number of isoforms") + coord_flip()


p

ggsave(".../nb_non_canonical_isoforms_structural_category.pdf", plot = p, 
       device = "pdf", width = 10, height = 6, units = "cm", dpi = 600)

  # 3.5  Number of neoisoforms by in-house annotation --------------------------

table(annot_with_ORFs_with_dups_collapsed$isoform_final)

annot <- as.data.frame(table(annot_with_ORFs_with_dups_collapsed$isoform_final))
colnames(annot) <- c("isoform_final","value")

annot$isoform_final = factor(annot$isoform_final, c("Neoexon; Intragenic TSS isoform",
                                                    "Intragenic TSS isoform",
                                                    "Neoexon",
                                                    "Non-canonical"))

vec = c("Neoexon" = "#EA7580FF", 
        "Neoexon; Intragenic TSS isoform" = "#088BBEFF",
       "Intragenic TSS isoform"= "#F8CD9CFF",
       "Non-canonical"= "#1BB6AFFF")

p <- ggplot(annot, aes(y=value, x=isoform_final, fill=isoform_final)) + 
  geom_bar(stat="identity", color = "black", width = 0.9, linewidth = 0.3) + 
  scale_fill_manual(values = vec) + 
  theme_classic() + 
  theme(legend.text = element_text(size=8, color="black", family = "Helvetica"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_text(size=8, color="black", family = "Helvetica"),
        axis.line = element_line(size = 0.3),
        axis.text.x = element_text(size=8, angle = 45, hjust =0.9, color="black", family = "Helvetica"),
        axis.text.y = element_text(size=8, color="black", family = "Helvetica", margin = margin(5, 5, 5, 5))) + 
  labs(x="", y="Number of isoforms") + coord_flip()


p

ggsave(".../nb_isoforms_isoform_final.pdf", plot = p, 
       device = "pdf", width = 14, height = 4, units = "cm", dpi = 600)
